"""Shared plotting style for publication figures."""

from __future__ import annotations

import matplotlib
import seaborn as sns


PUBLICATION_FONT = "Nimbus Sans"
PUBLICATION_FONT_STACK = [
    PUBLICATION_FONT,
    "Arial",
    "Helvetica",
    "DejaVu Sans",
    "Liberation Sans",
    "sans-serif",
]
PUBLICATION_FONT_CSS = '"Nimbus Sans", Arial, Helvetica, "DejaVu Sans", "Liberation Sans", sans-serif'


def configure_matplotlib(style: str = "whitegrid") -> None:
    """Apply the shared font and baseline plotting style."""
    rc = {
        "font.family": "sans-serif",
        "font.sans-serif": PUBLICATION_FONT_STACK,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    }
    matplotlib.rcParams.update(rc)
    sns.set_theme(style=style, font=PUBLICATION_FONT, rc=rc)
