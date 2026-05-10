#!/usr/bin/env Rscript

# Render Figure 3: ridgeline plot of interspecies and intraspecies ANI
# comparisons by genus, grouped by phylum.

prepend_default_user_library <- function() {
  version <- paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  candidates <- c(
    Sys.getenv("R_LIBS_USER", unset = NA_character_),
    file.path(path.expand("~"), "R", "x86_64-pc-linux-gnu-library", version)
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  for (candidate in candidates) {
    if (dir.exists(candidate) && !candidate %in% .libPaths()) {
      .libPaths(c(candidate, .libPaths()))
    }
  }
}

parse_args <- function(args) {
  opts <- list(
    repo_root = getwd(),
    output_prefix = NULL
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript bin/figure3_ani_ridgeline_by_genus.R [--repo-root PATH] [--output-prefix PATH]\n",
        "\n",
        "Defaults:\n",
        "  --repo-root      current working directory\n",
        "  --output-prefix  <repo-root>/assets/figures/figure3_ridgeline_ani\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg == "--repo-root") {
      i <- i + 1
      opts$repo_root <- args[[i]]
    } else if (arg == "--output-prefix") {
      i <- i + 1
      opts$output_prefix <- args[[i]]
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1
  }

  opts
}

prepend_default_user_library()
required_packages <- c("DBI", "duckdb", "dplyr", "ggplot2", "ggridges", "ggh4x")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(ggh4x)
})

opts <- parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- normalizePath(opts$repo_root, mustWork = TRUE)
parquet_glob <- file.path(repo_root, "data", "all_tables_processed", "*.parquet")

if (is.null(opts$output_prefix)) {
  output_prefix <- file.path(repo_root, "assets", "figures", "figure3_ridgeline_ani")
} else if (grepl("^/", opts$output_prefix)) {
  output_prefix <- opts$output_prefix
} else {
  output_prefix <- file.path(repo_root, opts$output_prefix)
}
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

min_pairs_per_group <- 20
x_min <- 75
x_max <- 100
x_ticks <- seq(75, 100, by = 5)
shade_xmin <- 95
shade_xmax <- 100
base_size <- 11
base_family <- "Helvetica"
dpi <- 300
width_in <- 8
height_in <- 14

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
df <- dbGetQuery(
  con,
  sprintf(
    "
    SELECT ANI, LSTR, GEN1_GENUS, GEN1_SPECIES, GEN2_SPECIES, GEN1_PHYLUM
    FROM read_parquet('%s')
    WHERE ANI IS NOT NULL
      AND LSTR IN ('species', 'genus')
    ",
    parquet_glob
  )
)
dbDisconnect(con, shutdown = TRUE)

df <- df %>%
  mutate(
    GENUS = GEN1_GENUS,
    PHYLUM = ifelse(is.na(GEN1_PHYLUM) | GEN1_PHYLUM == "", "Unknown", GEN1_PHYLUM),
    SP_CANON = pmin(GEN1_SPECIES, GEN2_SPECIES),
    group = case_when(
      LSTR == "species" ~ "intraspecies",
      LSTR == "genus" ~ "interspecies",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group), !is.na(GENUS), GENUS != "")

species_richness <- df %>%
  filter(group == "intraspecies") %>%
  distinct(PHYLUM, GENUS, SP_CANON) %>%
  count(PHYLUM, GENUS, name = "n_species")

df <- df %>%
  left_join(species_richness, by = c("PHYLUM", "GENUS")) %>%
  mutate(n_species = ifelse(is.na(n_species), 0L, n_species)) %>%
  filter(n_species > 1)

keep <- df %>%
  group_by(PHYLUM, GENUS, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= min_pairs_per_group)

df <- df %>%
  semi_join(keep, by = c("PHYLUM", "GENUS", "group"))

if (nrow(df) == 0) {
  stop("No data left after filtering.", call. = FALSE)
}

pair_counts <- df %>%
  count(PHYLUM, GENUS, name = "n_total_pairs")

df <- df %>%
  left_join(pair_counts, by = c("PHYLUM", "GENUS")) %>%
  mutate(
    count_label = formatC(n_total_pairs, format = "e", digits = 1),
    GENUS_LABEL = paste0(GENUS, "  (", count_label, ")")
  )

total_pairs <- nrow(df)
total_label <- formatC(total_pairs, format = "e", digits = 1)

df_all <- df %>%
  mutate(
    PHYLUM = "All phyla",
    GENUS_LABEL = paste0("All comparisons  (", total_label, ")")
  )

df <- bind_rows(df_all, df)

phylum_levels <- c("All phyla", sort(unique(df$PHYLUM[df$PHYLUM != "All phyla"])))
df <- df %>%
  mutate(PHYLUM = factor(PHYLUM, levels = phylum_levels)) %>%
  group_by(PHYLUM) %>%
  mutate(GENUS_LABEL = factor(GENUS_LABEL, levels = rev(sort(unique(GENUS_LABEL))))) %>%
  ungroup()

phylum_sizes <- df %>%
  distinct(PHYLUM, GENUS_LABEL) %>%
  count(PHYLUM, name = "n_genera")

facet_heights <- phylum_sizes$n_genera
names(facet_heights) <- phylum_sizes$PHYLUM
facet_heights["All phyla"] <- facet_heights["All phyla"] * 2
facet_heights <- as.list(facet_heights)

genus_lines <- df %>% distinct(PHYLUM, GENUS_LABEL)

p <- ggplot(df, aes(x = ANI, y = GENUS_LABEL, fill = group)) +
  annotate(
    "rect",
    xmin = shade_xmin, xmax = shade_xmax,
    ymin = -Inf, ymax = Inf,
    fill = "lightsteelblue1",
    alpha = 0.25
  ) +
  geom_hline(
    data = genus_lines,
    aes(yintercept = GENUS_LABEL),
    inherit.aes = FALSE,
    color = "grey88",
    linewidth = 0.25
  ) +
  geom_density_ridges(
    alpha = 0.6,
    scale = 1.0,
    rel_min_height = 0.02,
    color = "black",
    linewidth = 0.2,
    quantile_lines = TRUE,
    quantiles = c(0.25, 0.5, 0.75)
  ) +
  ggh4x::facet_wrap2(
    ~ PHYLUM,
    scales = "free_y",
    ncol = 1,
    strip.position = "top"
  ) +
  scale_fill_manual(values = c("interspecies" = "#F4A09C", "intraspecies" = "#55C9C7")) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = x_ticks,
    expand = c(0, 0)
  ) +
  labs(
    x = "ANI",
    y = "Genus (total pairs)",
    fill = NULL
  ) +
  theme_bw(base_size = base_size, base_family = base_family) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = base_size - 3),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  ggh4x::force_panelsizes(rows = facet_heights)

ggsave(
  filename = paste0(output_prefix, ".png"),
  plot = p,
  width = width_in,
  height = height_in,
  units = "in",
  dpi = dpi,
  bg = "white"
)

ggsave(
  filename = paste0(output_prefix, ".pdf"),
  plot = p,
  width = width_in,
  height = height_in,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(output_prefix, ".svg"),
  plot = p,
  width = width_in,
  height = height_in,
  units = "in",
  bg = "white"
)

message("Saved Figure 3 ridgeline files to ", dirname(output_prefix))
