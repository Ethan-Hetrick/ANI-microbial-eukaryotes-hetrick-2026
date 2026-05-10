#!/usr/bin/env python3
"""Render and assemble Figure 2 SVG panels.

Panel A is produced by ``bin/figure2_panel_a_rank_probabilities.R``. This script renders the
Python panels B/C and composes A/B/C into a single SVG without requiring
svgutils or other extra composition packages.
"""

from __future__ import annotations

import argparse
import copy
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from publication_figure_style import PUBLICATION_FONT_CSS, configure_matplotlib

configure_matplotlib(style="white")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import seaborn as sns
from scipy.ndimage import gaussian_filter


SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
PANEL_B_GRID_RIGHT = 0.83306
PANEL_B_COLORBAR_LEFT = 0.858
PANEL_B_MARGINAL_KDE_MAX_POINTS = 250_000
PANEL_C_DATA_RIGHT = 0.7430
PANEL_C_HIST_BINS = 420
PANEL_C_KDE_SIGMA = 2.35
PANEL_C_DENSITY_QUANTILES = (0.70, 0.78, 0.84, 0.89, 0.93, 0.965, 0.985, 0.995, 0.998)
PANEL_C_DENSITY_ALPHA = 0.62
PANEL_C_SCATTER_MAX_POINTS_PER_RANK = 9_000
PANEL_C_SCATTER_ALPHA = 0.10
PANEL_C_SCATTER_SIZE = 1.0
PANEL_A_RIGHT_PADDING = 24.0
SHARED_LEGEND_HEIGHT = 66.0
PANEL_AXIS_LABEL_SIZE = 17
PANEL_AXIS_TICK_SIZE = 14
COLORBAR_LABEL_SIZE = 15
COLORBAR_TICK_SIZE = 13
RANK_PALETTE = {"species": "#0072B2", "genus": "#009E73", "family-phylum": "#CC79A7"}
ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for path in [start, *start.parents]:
        if (path / "data" / "all_tables_processed").exists():
            return path
    raise FileNotFoundError("Could not find the publication repository root.")


def parse_args() -> argparse.Namespace:
    repo_default = find_repo_root(Path(__file__).resolve())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=repo_default)
    parser.add_argument("--panel-a", type=Path, default=None)
    parser.add_argument("--panel-b", type=Path, default=None)
    parser.add_argument("--panel-c", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--panel-gap", type=float, default=2.0, help="Vertical gap, in SVG points, between stacked panels.")
    parser.add_argument(
        "--only",
        choices=("all", "panels", "compose"),
        default="all",
        help="Render B/C and compose, render only B/C, or compose existing panel SVGs.",
    )
    return parser.parse_args()


def resolve_output(path: Path | None, repo_root: Path, default_name: str) -> Path:
    if path is None:
        path = repo_root / "assets" / default_name
    elif not path.is_absolute():
        path = repo_root / path
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def simplify_lstr(rank: object) -> str | None:
    if rank == "species":
        return "species"
    if rank == "genus":
        return "genus"
    if pd.isna(rank):
        return None
    return "family-phylum"


def load_plot_data(repo_root: Path) -> pd.DataFrame:
    dataset = ds.dataset(repo_root / "data" / "all_tables_processed", format="parquet")
    table = dataset.to_table(
        columns=["ANI", "AF", "LSTR"],
        filter=ds.field("ANI").is_valid() & ds.field("AF").is_valid() & ds.field("LSTR").is_valid(),
    )
    plot_df = table.to_pandas()
    plot_df["ANI"] = pd.to_numeric(plot_df["ANI"], errors="coerce")
    plot_df["AF"] = pd.to_numeric(plot_df["AF"], errors="coerce")
    plot_df = plot_df[plot_df["ANI"].between(72.5, 100)].copy()
    plot_df["AF_pct"] = plot_df["AF"] * 100
    plot_df["rank_simple"] = plot_df["LSTR"].map(simplify_lstr)
    plot_df = plot_df.dropna(subset=["ANI", "AF_pct", "rank_simple"])
    return plot_df[["ANI", "AF_pct", "rank_simple"]]


def save_svg(fig: matplotlib.figure.Figure, output: Path) -> None:
    fig.savefig(output, format="svg", facecolor="white", bbox_inches=None, metadata={"Date": None})
    plt.close(fig)


def render_panel_b(plot_df: pd.DataFrame, output: Path) -> None:
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(7.5, 1.05),
        height_ratios=(1.15, 5.85),
        left=0.10,
        right=PANEL_B_GRID_RIGHT,
        bottom=0.16,
        top=0.94,
        hspace=0.0,
        wspace=0.0,
    )
    ax_marg_x = fig.add_subplot(gs[0, 0])
    ax_joint = fig.add_subplot(gs[1, 0], sharex=ax_marg_x)
    ax_marg_y = fig.add_subplot(gs[1, 1], sharey=ax_joint)

    hb = ax_joint.hexbin(
        plot_df["ANI"],
        plot_df["AF_pct"],
        gridsize=50,
        cmap="viridis",
        mincnt=1,
        bins="log",
        linewidths=0,
    )

    for group, color in RANK_PALETTE.items():
        subset = plot_df[plot_df["rank_simple"] == group]
        if subset.empty:
            continue
        kde_subset = (
            subset.sample(n=PANEL_B_MARGINAL_KDE_MAX_POINTS, random_state=2026)
            if len(subset) > PANEL_B_MARGINAL_KDE_MAX_POINTS
            else subset
        )
        sns.kdeplot(
            data=kde_subset,
            x="ANI",
            fill=True,
            alpha=0.30,
            color=color,
            linewidth=1.2,
            ax=ax_marg_x,
            warn_singular=False,
        )
        sns.kdeplot(
            data=kde_subset,
            y="AF_pct",
            fill=True,
            alpha=0.30,
            color=color,
            linewidth=1.2,
            ax=ax_marg_y,
            warn_singular=False,
        )

    ax_joint.set_xlim(72.5, 100)
    ax_joint.set_ylim(0, 100)
    ax_joint.set_xlabel("ANI", fontsize=PANEL_AXIS_LABEL_SIZE)
    ax_joint.set_ylabel("Aligned Fraction (%)", fontsize=PANEL_AXIS_LABEL_SIZE)
    ax_joint.tick_params(axis="both", labelsize=PANEL_AXIS_TICK_SIZE, width=1.0)
    ax_joint.set_axisbelow(True)
    ax_joint.grid(True, color="#e6e6e6", linewidth=0.6, alpha=0.9)

    for ax in (ax_marg_x, ax_marg_y):
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.grid(False)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        sns.despine(ax=ax, left=True, bottom=True)

    sns.despine(ax=ax_joint)

    cbar = fig.colorbar(hb, cax=fig.add_axes([PANEL_B_COLORBAR_LEFT, 0.29, 0.026, 0.45]))
    cbar.set_label("log10(count)", fontsize=COLORBAR_LABEL_SIZE)
    cbar.ax.tick_params(labelsize=COLORBAR_TICK_SIZE)

    save_svg(fig, output)


def render_panel_c(plot_df: pd.DataFrame, output: Path) -> None:
    kde_df = plot_df.copy()

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.subplots_adjust(left=0.10, right=PANEL_C_DATA_RIGHT, bottom=0.14, top=0.925)

    for group, color in RANK_PALETTE.items():
        subset = kde_df[kde_df["rank_simple"] == group]
        if subset.empty:
            continue

        scatter_subset = (
            subset.sample(n=PANEL_C_SCATTER_MAX_POINTS_PER_RANK, random_state=2026)
            if len(subset) > PANEL_C_SCATTER_MAX_POINTS_PER_RANK
            else subset
        )
        ax.scatter(
            scatter_subset["ANI"],
            scatter_subset["AF_pct"],
            s=PANEL_C_SCATTER_SIZE,
            color=color,
            alpha=PANEL_C_SCATTER_ALPHA,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

        x_edges = np.linspace(72.5, 100, PANEL_C_HIST_BINS + 1)
        y_edges = np.linspace(0, 100, PANEL_C_HIST_BINS + 1)
        hist, _, _ = np.histogram2d(subset["ANI"], subset["AF_pct"], bins=(x_edges, y_edges), density=False)
        density = gaussian_filter(hist.T, sigma=PANEL_C_KDE_SIGMA)
        positive_density = density[density > 0]
        if positive_density.size < 2:
            continue

        fill_levels = np.quantile(positive_density, PANEL_C_DENSITY_QUANTILES)
        fill_levels = np.unique(fill_levels)
        if fill_levels.size < 2:
            continue

        x_centers = (x_edges[:-1] + x_edges[1:]) / 2
        y_centers = (y_edges[:-1] + y_edges[1:]) / 2
        transparent_color = mcolors.to_rgba(color, 0.00)
        filled_color = mcolors.to_rgba(color, PANEL_C_DENSITY_ALPHA)
        cmap = mcolors.LinearSegmentedColormap.from_list(f"{group}_density", [transparent_color, filled_color])

        ax.contourf(
            x_centers,
            y_centers,
            density,
            levels=fill_levels,
            cmap=cmap,
            antialiased=True,
            zorder=2,
        )
        ax.contour(
            x_centers,
            y_centers,
            density,
            levels=fill_levels,
            colors=color,
            linewidths=1.0,
            alpha=0.85,
            zorder=3,
        )

    ax.set_xlim(72.5, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("ANI", fontsize=PANEL_AXIS_LABEL_SIZE)
    ax.set_ylabel("Aligned Fraction (%)", fontsize=PANEL_AXIS_LABEL_SIZE)
    ax.tick_params(axis="both", labelsize=PANEL_AXIS_TICK_SIZE, width=1.0)
    ax.set_axisbelow(True)
    ax.grid(True, color="#e6e6e6", linewidth=0.6, alpha=0.9)
    sns.despine(ax=ax)
    save_svg(fig, output)


def parse_viewbox(root: ET.Element) -> tuple[float, float, float, float]:
    viewbox = root.attrib.get("viewBox")
    if viewbox:
        values = [float(value) for value in re.split(r"[,\s]+", viewbox.strip()) if value]
        if len(values) == 4:
            return values[0], values[1], values[2], values[3]
    width = float(re.sub(r"[^0-9.]", "", root.attrib.get("width", "0")))
    height = float(re.sub(r"[^0-9.]", "", root.attrib.get("height", "0")))
    return 0.0, 0.0, width, height


def prefix_svg_ids(root: ET.Element, prefix: str) -> None:
    id_map: dict[str, str] = {}
    for element in root.iter():
        old_id = element.attrib.get("id")
        if old_id:
            new_id = f"{prefix}_{old_id}"
            id_map[old_id] = new_id
            element.attrib["id"] = new_id

    if not id_map:
        return

    url_pattern = re.compile(r"url\(#([^)]+)\)")
    hash_pattern = re.compile(r"^#(.+)$")

    def replace_value(value: str) -> str:
        value = url_pattern.sub(lambda match: f"url(#{id_map.get(match.group(1), match.group(1))})", value)
        hash_match = hash_pattern.match(value)
        if hash_match:
            old = hash_match.group(1)
            value = f"#{id_map.get(old, old)}"
        return value

    for element in root.iter():
        for key, value in list(element.attrib.items()):
            if "#" in value:
                element.attrib[key] = replace_value(value)
        if element.text and "#" in element.text:
            element.text = replace_value(element.text)


def make_unfilled_elements_explicit(root: ET.Element) -> None:
    for element in root.iter():
        local_name = element.tag.rsplit("}", 1)[-1]
        if local_name not in {"line", "polyline", "rect"}:
            continue
        style = element.attrib.get("style", "")
        has_inline_fill = "fill:" in style
        has_fill_attr = "fill" in element.attrib
        if not has_inline_fill and not has_fill_attr:
            element.attrib["style"] = f"fill: none; {style}".strip()
            element.attrib["fill"] = "none"
            style = element.attrib["style"]
        if local_name in {"line", "polyline"} and "stroke:" not in style and "stroke" not in element.attrib:
            element.attrib["style"] = f"stroke: #000000; {style}".strip()


def load_svg(path: Path, prefix: str) -> tuple[ET.Element, tuple[float, float, float, float]]:
    root = ET.parse(path).getroot()
    prefix_svg_ids(root, prefix)
    make_unfilled_elements_explicit(root)
    return root, parse_viewbox(root)


def append_nested_svg(
    parent: ET.Element,
    source: ET.Element,
    viewbox: tuple[float, float, float, float],
    x: float,
    y: float,
    width: float,
    height: float,
    preserve_aspect_ratio: str | None = None,
) -> None:
    min_x, min_y, vb_width, vb_height = viewbox
    attrs = {
        "x": f"{x:.2f}",
        "y": f"{y:.2f}",
        "width": f"{width:.2f}",
        "height": f"{height:.2f}",
        "viewBox": f"{min_x:.2f} {min_y:.2f} {vb_width:.2f} {vb_height:.2f}",
    }
    if preserve_aspect_ratio:
        attrs["preserveAspectRatio"] = preserve_aspect_ratio
    nested = ET.SubElement(
        parent,
        f"{{{SVG_NS}}}svg",
        attrs,
    )
    for child in list(source):
        nested.append(copy.deepcopy(child))


def crop_viewbox(
    viewbox: tuple[float, float, float, float],
    crop: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    left, top, right, bottom = crop
    min_x, min_y, width, height = viewbox
    return (
        min_x + left,
        min_y + top,
        width - left - right,
        height - top - bottom,
    )


def extract_x_tick_positions(root: ET.Element, tick_values: tuple[str, ...] = ("75", "100")) -> dict[str, float]:
    positions: dict[str, list[float]] = {tick: [] for tick in tick_values}
    for element in root.iter():
        if not element.tag.endswith("text"):
            continue
        label = "".join(element.itertext()).strip()
        if label not in positions or "x" not in element.attrib:
            continue
        try:
            x_value = float(str(element.attrib["x"]).split()[0])
        except ValueError:
            continue
        positions[label].append(x_value)

    # If a tick value appears in both axes, the x-axis tick is the rightmost one.
    return {tick: max(values) for tick, values in positions.items() if values}


def mean_tick_positions(tick_sets: list[dict[str, float]], tick_values: tuple[str, ...] = ("75", "100")) -> dict[str, float]:
    means: dict[str, float] = {}
    for tick in tick_values:
        values = [ticks[tick] for ticks in tick_sets if tick in ticks]
        if values:
            means[tick] = sum(values) / len(values)
    return means


def align_viewbox_x_to_ticks(
    viewbox: tuple[float, float, float, float],
    source_ticks: dict[str, float],
    target_ticks: dict[str, float],
    display_width: float,
    left_tick: str = "75",
    right_tick: str = "100",
) -> tuple[float, float, float, float]:
    if left_tick not in source_ticks or right_tick not in source_ticks:
        return viewbox
    if left_tick not in target_ticks or right_tick not in target_ticks:
        return viewbox

    source_span = source_ticks[right_tick] - source_ticks[left_tick]
    target_span = target_ticks[right_tick] - target_ticks[left_tick]
    if source_span <= 0 or target_span <= 0:
        return viewbox

    x_scale = target_span / source_span
    min_x = source_ticks[left_tick] - (target_ticks[left_tick] / x_scale)
    _, min_y, _, height = viewbox
    return min_x, min_y, display_width / x_scale, height


def add_panel_label(root: ET.Element, label: str, x: float, y: float) -> None:
    text = ET.SubElement(
        root,
        f"{{{SVG_NS}}}text",
        {
            "x": f"{x:.2f}",
            "y": f"{y:.2f}",
            "font-family": PUBLICATION_FONT_CSS,
            "font-size": "34",
            "font-weight": "700",
        },
    )
    text.text = label


def add_svg_text(
    root: ET.Element,
    text: str,
    x: float,
    y: float,
    size: float = 12.5,
    weight: str = "400",
) -> None:
    element = ET.SubElement(
        root,
        f"{{{SVG_NS}}}text",
        {
            "x": f"{x:.2f}",
            "y": f"{y:.2f}",
            "font-family": PUBLICATION_FONT_CSS,
            "font-size": f"{size:.1f}",
            "font-weight": weight,
            "fill": "#222222",
        },
    )
    element.text = text


def add_shared_legend(root: ET.Element, total_width: float) -> None:
    box_margin = 34.0
    box_x = box_margin
    box_y = 10.0
    box_width = total_width - (box_margin * 2)
    box_height = 43.0
    center_y = box_y + 25.0

    ET.SubElement(
        root,
        f"{{{SVG_NS}}}rect",
        {
            "x": f"{box_x:.2f}",
            "y": f"{box_y:.2f}",
            "width": f"{box_width:.2f}",
            "height": f"{box_height:.2f}",
            "rx": "5.5",
            "ry": "5.5",
            "fill": "#ffffff",
            "stroke": "#d0d0d0",
            "stroke-width": "0.9",
        },
    )

    add_svg_text(root, "Figure 2 legend", box_x + 17.0, center_y + 4.5, size=13.8, weight="700")
    add_svg_text(root, "LSTR", box_x + 157.0, center_y + 4.5, size=12.8, weight="700")

    x = box_x + 207.0
    rank_steps = {"species": 106.0, "genus": 98.0, "family-phylum": 156.0}
    for label, color in RANK_PALETTE.items():
        ET.SubElement(
            root,
            f"{{{SVG_NS}}}line",
            {
                "x1": f"{x:.2f}",
                "x2": f"{x + 36.0:.2f}",
                "y1": f"{center_y:.2f}",
                "y2": f"{center_y:.2f}",
                "stroke": color,
                "stroke-width": "4.2",
                "stroke-linecap": "round",
            },
        )
        add_svg_text(root, label, x + 47.0, center_y + 4.5, size=13.0)
        x += rank_steps[label]

    ET.SubElement(
        root,
        f"{{{SVG_NS}}}rect",
        {
            "x": f"{x + 7.0:.2f}",
            "y": f"{center_y - 8.0:.2f}",
            "width": "17.0",
            "height": "17.0",
            "fill": "#d9d9d9",
            "fill-opacity": "0.35",
            "stroke": "none",
        },
    )
    add_svg_text(root, "ANI discontinuity", x + 36.0, center_y + 4.5, size=13.0)
    x += 178.0

    ET.SubElement(
        root,
        f"{{{SVG_NS}}}line",
        {
            "x1": f"{x:.2f}",
            "x2": f"{x + 36.0:.2f}",
            "y1": f"{center_y:.2f}",
            "y2": f"{center_y:.2f}",
            "stroke": "#222222",
            "stroke-width": "2.2",
            "stroke-linecap": "round",
            "stroke-dasharray": "2.4 5.0",
        },
    )
    add_svg_text(root, "Intersection", x + 47.0, center_y + 4.5, size=13.0)


def extract_font_face_style(root: ET.Element) -> str | None:
    for style in root.iter(f"{{{SVG_NS}}}style"):
        if style.text and "@font-face" in style.text:
            return style.text
    return None


def compose_figure(panel_a: Path, panel_b: Path, panel_c: Path, output: Path, panel_gap: float = 2.0) -> None:
    panels = [
        ("A", panel_a, "panel_a"),
        ("B", panel_b, "panel_b"),
        ("C", panel_c, "panel_c"),
    ]
    loaded = [(label, *load_svg(path, prefix)) for label, path, prefix in panels]
    panel_crops = {
        "A": (0.0, 0.0, 0.0, 8.0),
        "B": (0.0, 32.0, 0.0, 42.0),
        "C": (0.0, 30.0, 0.0, 0.0),
    }

    cropped_loaded = [
        (label, root, crop_viewbox(viewbox, panel_crops.get(label, (0.0, 0.0, 0.0, 0.0))))
        for label, root, viewbox in loaded
    ]

    panel_width = max(viewbox[2] for _, _, viewbox in cropped_loaded)
    viewport_widths = {
        label: panel_width + (PANEL_A_RIGHT_PADDING if label == "A" else 0.0)
        for label, _, _ in cropped_loaded
    }
    x_ticks = {label: extract_x_tick_positions(root) for label, root, _ in cropped_loaded}
    target_ticks = mean_tick_positions([x_ticks.get("B", {}), x_ticks.get("C", {})])
    cropped_loaded = [
        (
            label,
            root,
            align_viewbox_x_to_ticks(viewbox, x_ticks.get(label, {}), target_ticks, viewport_widths[label])
            if label == "A"
            else viewbox,
        )
        for label, root, viewbox in cropped_loaded
    ]

    label_gutter = 54.0
    gap = panel_gap
    y = SHARED_LEGEND_HEIGHT
    placements: list[tuple[str, ET.Element, tuple[float, float, float, float], float, float, float, float, str | None]] = []

    for label, root, viewbox in cropped_loaded:
        width = viewport_widths[label]
        height = viewbox[3] * (width / viewbox[2])
        preserve_aspect_ratio = None
        placements.append((label, root, viewbox, label_gutter, y, width, height, preserve_aspect_ratio))
        y += height + gap

    total_width = max(viewport_widths.values()) + label_gutter
    total_height = y - gap
    out_root = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "width": f"{total_width:.2f}pt",
            "height": f"{total_height:.2f}pt",
            "viewBox": f"0 0 {total_width:.2f} {total_height:.2f}",
            "version": "1.1",
        },
    )
    ET.SubElement(
        out_root,
        f"{{{SVG_NS}}}rect",
        {"x": "0", "y": "0", "width": "100%", "height": "100%", "fill": "white"},
    )
    font_face_style = extract_font_face_style(loaded[0][1])
    if font_face_style:
        defs = ET.SubElement(out_root, f"{{{SVG_NS}}}defs")
        style = ET.SubElement(defs, f"{{{SVG_NS}}}style", {"type": "text/css"})
        style.text = font_face_style

    add_shared_legend(out_root, total_width)

    for label, root, viewbox, x, y0, width, height, preserve_aspect_ratio in placements:
        append_nested_svg(out_root, root, viewbox, x, y0, width, height, preserve_aspect_ratio)
        add_panel_label(out_root, f"({label})", 4.0, y0 + 42.0)

    ET.ElementTree(out_root).write(output, encoding="utf-8", xml_declaration=True)


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    panel_a = resolve_output(args.panel_a, repo_root, "figure2_panel_a.svg")
    panel_b = resolve_output(args.panel_b, repo_root, "figure2_panel_b.svg")
    panel_c = resolve_output(args.panel_c, repo_root, "figure2_panel_c.svg")
    output = resolve_output(args.output, repo_root, "figure2.svg")

    if args.only in {"all", "panels"}:
        plot_df = load_plot_data(repo_root)
        print(f"Loaded {len(plot_df):,} rows for Figure 2B/C", flush=True)
        render_panel_b(plot_df, panel_b)
        print(f"Wrote {panel_b}", flush=True)
        render_panel_c(plot_df, panel_c)
        print(f"Wrote {panel_c}", flush=True)

    if args.only in {"all", "compose"}:
        missing = [path for path in (panel_a, panel_b, panel_c) if not path.exists()]
        if missing:
            raise FileNotFoundError(
                "Missing panel SVG(s): "
                + ", ".join(str(path) for path in missing)
                + ". Render Panel A with bin/figure2_panel_a_rank_probabilities.R first."
            )
        compose_figure(panel_a, panel_b, panel_c, output, panel_gap=args.panel_gap)
        print(f"Wrote {output}", flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
