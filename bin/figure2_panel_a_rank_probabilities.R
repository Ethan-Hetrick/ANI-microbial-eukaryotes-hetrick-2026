#!/usr/bin/env Rscript

# Render Figure 2A: ANI histogram, taxonomic-rank proportions, and
# multinomial-logistic probability curves.

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
    repo_root = NULL,
    output = NULL,
    width = 12,
    height = 8.35
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript bin/figure2_panel_a_rank_probabilities.R [--repo-root PATH] [--output PATH] [--width N] [--height N]\n",
        "\n",
        "Defaults:\n",
        "  --repo-root  auto-detected from the script location or working directory\n",
        "  --output     <repo-root>/assets/figure2_panel_a.svg\n",
        "  --width      12\n",
        "  --height     8.35\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg == "--repo-root") {
      i <- i + 1
      opts$repo_root <- args[[i]]
    } else if (arg == "--output") {
      i <- i + 1
      opts$output <- args[[i]]
    } else if (arg == "--width") {
      i <- i + 1
      opts$width <- as.numeric(args[[i]])
    } else if (arg == "--height") {
      i <- i + 1
      opts$height <- as.numeric(args[[i]])
    } else if (startsWith(arg, "--")) {
      stop("Unknown option: ", arg, call. = FALSE)
    } else if (is.null(opts$output)) {
      opts$output <- arg
    } else {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }
    i <- i + 1
  }

  opts
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(NULL)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

find_repo_root <- function(starts) {
  for (start in starts) {
    if (is.null(start) || is.na(start) || !nzchar(start)) {
      next
    }
    path <- normalizePath(start, mustWork = TRUE)
    if (!dir.exists(path)) {
      path <- dirname(path)
    }

    repeat {
      has_data <- dir.exists(file.path(path, "data", "all_tables_processed"))
      if (has_data) {
        return(path)
      }

      parent <- dirname(path)
      if (identical(parent, path)) {
        break
      }
      path <- parent
    }
  }

  stop("Could not find the repository root.", call. = FALSE)
}

require_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0) {
    stop(
      "Install required R package(s): ",
      paste(missing, collapse = ", "),
      "\nActive R library paths:\n  ",
      paste(.libPaths(), collapse = "\n  "),
      call. = FALSE
    )
  }
}

find_intersection <- function(x, y1, y2) {
  delta <- y1 - y2
  index <- which(diff(sign(delta)) != 0)

  if (length(index) == 0) {
    return(NA_real_)
  }

  i <- index[1]
  x[i] - delta[i] * (x[i + 1] - x[i]) / (delta[i + 1] - delta[i])
}

prepend_default_user_library()

required_packages <- c(
  "DBI",
  "dplyr",
  "duckdb",
  "ggplot2",
  "nnet",
  "patchwork",
  "scales",
  "svglite",
  "systemfonts",
  "tibble"
)
require_packages(required_packages)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(nnet)
  library(patchwork)
  library(scales)
})

opts <- parse_args(commandArgs(trailingOnly = TRUE))
script_file <- script_path()
repo_root <- if (!is.null(opts$repo_root)) {
  normalizePath(opts$repo_root, mustWork = TRUE)
} else {
  find_repo_root(c(dirname(script_file), getwd()))
}

output <- opts$output
if (is.null(output)) {
  output <- file.path(repo_root, "assets", "figure2_panel_a.svg")
}
if (!grepl("^/", output)) {
  output <- file.path(repo_root, output)
}
output <- normalizePath(output, mustWork = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "bin", "publication_figure_style.R"))
all_tables_con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(DBI::dbDisconnect(all_tables_con, shutdown = TRUE), add = TRUE)

all_tables_glob <- normalizePath(
  file.path(repo_root, "data", "all_tables_processed", "*.parquet"),
  mustWork = FALSE
)
DBI::dbExecute(
  all_tables_con,
  sprintf(
    "CREATE VIEW all_tables AS SELECT * FROM read_parquet('%s')",
    gsub("'", "''", all_tables_glob, fixed = TRUE)
  )
)

rank_palette <- c(
  "species" = "#0072B2",
  "genus" = "#009E73",
  "family_phylum" = "#CC79A7"
)

figure_font <- PUBLICATION_FONT
svg_fonts <- publication_svg_fonts()

ani_rank <- DBI::dbGetQuery(
  all_tables_con,
  "
  SELECT ANI, LSTR
  FROM all_tables
  WHERE ANI IS NOT NULL
    AND LSTR IS NOT NULL
  "
) %>%
  mutate(
    ANI = as.numeric(ANI),
    LSTR = as.character(LSTR)
  )

probability_data <- ani_rank %>%
  mutate(
    rank_simple = case_when(
      LSTR == "species" ~ "species",
      LSTR == "genus" ~ "genus",
      TRUE ~ "family_phylum"
    ),
    rank_simple = factor(
      rank_simple,
      levels = c("species", "genus", "family_phylum")
    )
  )

probability_fit <- multinom(rank_simple ~ ANI, data = probability_data, trace = FALSE)

data_x_min <- min(probability_data$ANI, na.rm = TRUE)
data_x_max <- max(probability_data$ANI, na.rm = TRUE)
x_min <- 72.5
x_max <- 100
x_breaks <- seq(75, 100, by = 5)
x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
ani_x_scale <- function() {
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = x_breaks,
    expand = expansion(mult = 0, add = 0)
  )
}

prediction_grid <- tibble::tibble(
  ANI = seq(
    data_x_min,
    data_x_max,
    length.out = 2000
  )
)

predicted_probabilities <- as.data.frame(
  predict(probability_fit, newdata = prediction_grid, type = "probs")
)

for (rank in names(rank_palette)) {
  if (!rank %in% names(predicted_probabilities)) {
    predicted_probabilities[[rank]] <- 0
  }
}

prediction_grid <- bind_cols(prediction_grid, predicted_probabilities)

intersection_ani_left <- find_intersection(
  prediction_grid$ANI,
  prediction_grid$genus,
  prediction_grid$family_phylum
)

intersection_ani_right <- find_intersection(
  prediction_grid$ANI,
  prediction_grid$species,
  prediction_grid$genus
)

discontinuity_band <- tibble::tibble(
  xmin = 80,
  xmax = 96,
  ymin = -Inf,
  ymax = Inf,
  band = "ANI discontinuity"
)

intersection_lines <- tibble::tibble(
  xint = c(intersection_ani_left, intersection_ani_right)
) %>%
  filter(!is.na(xint))

binwidth <- (x_max - x_min) / 100

figure_base_size <- 18
figure_axis_text_size <- 14.6
figure_axis_title_size <- 18
figure_y_axis_title_size <- 16.2
figure_legend_text_size <- 13
figure_legend_title_size <- 14
figure_annotation_size <- 5.2

shared_theme <- theme_minimal(base_size = figure_base_size, base_family = figure_font) +
  theme(
    text = element_text(family = figure_font),
    axis.text = element_text(size = figure_axis_text_size),
    axis.title = element_text(size = figure_axis_title_size),
    axis.title.y = element_text(size = figure_y_axis_title_size),
    legend.text = element_text(size = figure_legend_text_size),
    legend.title = element_text(size = figure_legend_title_size),
    panel.grid = element_blank(),
    plot.margin = margin(2, 10, 2, 10)
  )

ani_histogram <- ggplot(probability_data, aes(x = ANI)) +
  geom_rect(
    data = discontinuity_band,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.25
  ) +
  geom_histogram(
    binwidth = binwidth,
    boundary = x_min,
    color = "black",
    fill = "white"
  ) +
  geom_vline(
    data = intersection_lines,
    aes(xintercept = xint),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  ani_x_scale() +
  scale_y_continuous(labels = label_scientific()) +
  labs(x = NULL, y = "ANI Comparisons") +
  shared_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(12, 10, 7, 10)
  )

rank_proportions <- ggplot(probability_data, aes(x = ANI, fill = rank_simple)) +
  geom_rect(
    data = discontinuity_band,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.25
  ) +
  geom_histogram(
    aes(y = after_stat(count)),
    binwidth = binwidth,
    boundary = x_min,
    position = "fill",
    color = "black",
    linewidth = 0.15
  ) +
  geom_vline(
    data = intersection_lines,
    aes(xintercept = xint),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  ani_x_scale() +
  scale_y_continuous(
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = rank_palette,
    breaks = c("family_phylum", "genus", "species"),
    labels = c("family-phylum", "genus", "species")
  ) +
  labs(x = NULL, y = "Proportion") +
  shared_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(6, 10, 5, 10)
  )

probability_curves <- ggplot(prediction_grid, aes(x = ANI)) +
  geom_rect(
    data = discontinuity_band,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = band),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  geom_line(aes(y = family_phylum, color = "family_phylum"), linewidth = 1) +
  geom_line(aes(y = genus, color = "genus"), linewidth = 1) +
  geom_line(aes(y = species, color = "species"), linewidth = 1) +
  geom_vline(
    data = intersection_lines,
    aes(xintercept = xint, linetype = "Intersection"),
    inherit.aes = FALSE,
    linewidth = 0.6,
    key_glyph = draw_key_path
  ) +
  annotate(
    "text",
    x = intersection_ani_left - 0.4,
    y = 0.55,
    label = round(intersection_ani_left, 2),
    hjust = 1.2,
    vjust = 2,
    size = figure_annotation_size
  ) +
  annotate(
    "text",
    x = intersection_ani_right + 0.4,
    y = 0.55,
    label = round(intersection_ani_right, 2),
    hjust = -0.2,
    vjust = 2,
    size = figure_annotation_size
  ) +
  ani_x_scale() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = rank_palette,
    breaks = c("family_phylum", "genus", "species"),
    labels = c("family-phylum", "genus", "species"),
    name = "Lowest Shared Taxonomic Rank"
  ) +
  scale_fill_manual(
    values = c("ANI discontinuity" = "grey85"),
    name = " "
  ) +
  scale_linetype_manual(
    values = c("Intersection" = "dotted"),
    name = "  "
  ) +
  labs(x = "ANI", y = "Probability") +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(alpha = 0.25, colour = NA, linetype = 0)
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(fill = NA, colour = "black", linewidth = 0.6, alpha = 1)
    ),
    color = guide_legend(
      order = 3,
      override.aes = list(
        colour = unname(rank_palette[c("family_phylum", "genus", "species")]),
        fill = NA,
        linetype = c(1, 1, 1),
        linewidth = c(1.2, 1.2, 1.2),
        alpha = 1
      )
    )
  ) +
  shared_theme +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = NA, colour = NA)
  )

probability_plot <- ani_histogram / rank_proportions / probability_curves +
  plot_layout(heights = c(1.15, 1.15, 3.25))

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    output,
    probability_plot,
    width = opts$width,
    height = opts$height,
    device = svglite::svglite,
    user_fonts = svg_fonts$user_fonts,
    web_fonts = svg_fonts$web_fonts
  )
} else {
  warning("Package 'svglite' is not available; using base R svg() device.")
  ggsave(output, probability_plot, width = opts$width, height = opts$height, device = "svg")
}

message("Wrote ", output)
message("Rows plotted: ", format(nrow(probability_data), big.mark = ","))
message("Intersections: ", paste(round(c(intersection_ani_left, intersection_ani_right), 2), collapse = ", "))
