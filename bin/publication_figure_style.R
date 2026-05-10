# Shared plotting style for publication figures.

PUBLICATION_FONT <- "Nimbus Sans"
PUBLICATION_SVG_FONT_ALIAS <- PUBLICATION_FONT

publication_theme_minimal <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size, base_family = PUBLICATION_FONT) +
    ggplot2::theme(text = ggplot2::element_text(family = PUBLICATION_FONT))
}

publication_theme_bw <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size, base_family = PUBLICATION_FONT) +
    ggplot2::theme(text = ggplot2::element_text(family = PUBLICATION_FONT))
}

publication_font_file <- function() {
  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    stop("Install required R package: systemfonts", call. = FALSE)
  }

  font_file <- systemfonts::match_fonts(PUBLICATION_FONT)$path[[1]]
  if (is.na(font_file) || !file.exists(font_file)) {
    stop("Could not find the ", PUBLICATION_FONT, " font file.", call. = FALSE)
  }

  font_file
}

publication_svg_fonts <- function() {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Install required R package: svglite", call. = FALSE)
  }

  font_file <- publication_font_file()
  extension <- tolower(tools::file_ext(font_file))
  font_args <- list(
    family = PUBLICATION_SVG_FONT_ALIAS,
    local = c(PUBLICATION_FONT, "Helvetica", PUBLICATION_SVG_FONT_ALIAS),
    embed = TRUE
  )

  if (identical(extension, "otf")) {
    font_args$otf <- font_file
  } else {
    font_args$ttf <- font_file
  }

  list(
    user_fonts = setNames(
      list(
        list(
          plain = list(alias = PUBLICATION_SVG_FONT_ALIAS, file = font_file)
        )
      ),
      PUBLICATION_FONT
    ),
    web_fonts = list(do.call(svglite::font_face, font_args))
  )
}
