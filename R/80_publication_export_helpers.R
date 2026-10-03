############################################################
## 80_publication_export_helpers.R
##
## Export-only helpers for publication-ready main figure PDFs.
## These functions do not alter statistical objects or plotting data.
############################################################

mr_mm_to_in <- function(mm) {
  mm / 25.4
}

mr_full_width_mm <- function() {
  170
}

mr_max_height_mm <- function() {
  225
}

mr_capture_output_dirs <- function(cfg) {
  list(
    output_root = cfg$output_root,
    table_dir = cfg$table_dir,
    plot_dir = cfg$plot_dir,
    object_dir = cfg$object_dir,
    pathway_dir = cfg$pathway_dir
  )
}

mr_restore_output_dirs <- function(cfg, dirs) {
  for (nm in names(dirs)) {
    cfg[[nm]] <- dirs[[nm]]
  }

  for (nm in c("output_root", "table_dir", "plot_dir", "object_dir", "pathway_dir")) {
    if (!is.null(cfg[[nm]]) && !is.na(cfg[[nm]])) {
      dir.create(cfg[[nm]], recursive = TRUE, showWarnings = FALSE)
    }
  }

  cfg
}

mr_use_staged_output_dirs <- function(cfg, staging_root) {
  cfg$output_root <- staging_root
  cfg$table_dir <- file.path(staging_root, "tables")
  cfg$plot_dir <- file.path(staging_root, "plots")
  cfg$object_dir <- file.path(staging_root, "objects")
  cfg$pathway_dir <- file.path(staging_root, "pathway_gene_sets_all_clients")

  for (path in c(cfg$output_root, cfg$table_dir, cfg$plot_dir, cfg$object_dir, cfg$pathway_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  cfg
}

mr_new_staging_root <- function() {
  root <- file.path(
    tempdir(),
    paste0("rosmap_manuscript_ready_side_effects_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root
}

mr_manuscript_ready_output_dir <- function(cfg) {
  if (is.null(cfg$plot_dir) || is.na(cfg$plot_dir)) {
    stop("cfg$plot_dir is missing; cannot create manuscript-ready output directory.", call. = FALSE)
  }

  out_dir <- file.path(cfg$plot_dir, "manuscript_ready_pdf")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

mr_child_cleanup_theme <- function() {
  ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank()
  )
}

mr_outer_cleanup_theme <- function(tag_size = 12, figure = NA_integer_) {
  outer_margin <- switch(
    as.character(figure),
    `1` = ggplot2::margin(10, 14, 10, 14),
    `2` = ggplot2::margin(10, 12, 10, 12),
    `3` = ggplot2::margin(9, 11, 9, 11),
    `4` = ggplot2::margin(14, 18, 14, 18),
    `5` = ggplot2::margin(14, 18, 16, 18),
    `6` = ggplot2::margin(14, 18, 16, 18),
    ggplot2::margin(8, 10, 8, 10)
  )

  ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank(),
    plot.margin = outer_margin,
    plot.tag = ggplot2::element_text(face = "bold", size = tag_size, color = "black"),
    plot.tag.position = "topleft"
  )
}

mr_make_manuscript_ready_plot <- function(plot, tag_levels = "A", tag_size = 12, figure = NA_integer_) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for manuscript-ready cleanup.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for manuscript-ready cleanup.", call. = FALSE)
  }

  child_cleanup <- mr_child_cleanup_theme()
  outer_cleanup <- mr_outer_cleanup_theme(tag_size = tag_size, figure = figure)

  plot_clean <- tryCatch(
    plot & ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL),
    error = function(e) {
      plot + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
    }
  )

  plot_clean <- tryCatch(
    plot_clean & child_cleanup,
    error = function(e) {
      plot_clean + child_cleanup
    }
  )

  plot_clean <- mr_apply_figure_export_layout(plot_clean, figure = figure)

  plot_clean +
    patchwork::plot_annotation(
      title = NULL,
      subtitle = NULL,
      caption = NULL,
      tag_levels = tag_levels,
      tag_prefix = "",
      tag_suffix = "",
      theme = outer_cleanup
    )
}

mr_apply_figure_export_layout <- function(plot, figure = NA_integer_) {
  if (is.na(figure)) {
    return(plot)
  }

  extra_theme <- switch(
    as.character(figure),
    `1` = ggplot2::theme(
      strip.text = ggtext::element_markdown(
        size = 5.25,
        lineheight = 0.84,
        margin = ggplot2::margin(6, 2, 6, 2)
      ),
      strip.clip = "off",
      panel.spacing.x = grid::unit(0.45, "lines"),
      panel.spacing.y = grid::unit(0.9, "lines")
    ),
    `2` = ggplot2::theme(
      plot.margin = ggplot2::margin(8, 10, 8, 10)
    ),
    `3` = ggplot2::theme(
      axis.title = ggplot2::element_text(
        face = "bold",
        size = 8.3
      ),
      axis.text = ggplot2::element_text(
        size = 7.4
      ),
      legend.margin = ggplot2::margin(
        1, 1, 1, 1
      ),
      legend.text = ggplot2::element_text(
        size = 7.2
      ),
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 7.4
      ),
      legend.key.size = grid::unit(
        0.17,
        "cm"
      ),
      legend.spacing.x = grid::unit(
        1.5,
        "pt"
      )
    ),
    `4` = ggplot2::theme(
      legend.margin = ggplot2::margin(4, 4, 4, 4),
      legend.text = ggplot2::element_text(size = 5.8, lineheight = 0.86),
      legend.title = ggplot2::element_text(face = "bold", size = 6.8),
      legend.key.size = grid::unit(0.20, "cm"),
      legend.box = "vertical",
      legend.spacing.x = grid::unit(2, "pt"),
      legend.spacing.y = grid::unit(1, "pt")
    ),
    `5` = ggplot2::theme(),
    `6` = ggplot2::theme(
      legend.margin = ggplot2::margin(4, 4, 4, 4),
      legend.text = ggplot2::element_text(size = 5.8, lineheight = 0.86),
      legend.title = ggplot2::element_text(face = "bold", size = 6.8),
      legend.key.size = grid::unit(0.20, "cm"),
      legend.spacing.x = grid::unit(2, "pt"),
      legend.spacing.y = grid::unit(1, "pt")
    ),
    ggplot2::theme(
      plot.margin = ggplot2::margin(7, 9, 7, 9)
    )
  )

  tryCatch(
    {
      plot <- plot & extra_theme
      if (identical(as.integer(figure), 3L)) {
        plot <- plot & ggplot2::guides(color = "none")
      }
      if (identical(as.integer(figure), 4L)) {
        plot <- plot & ggplot2::guides(
          size = "none"
        )
      }
      plot
    },
    error = function(e) {
      plot + extra_theme
    }
  )
}

mr_first_existing_tool <- function(tool_names) {
  hits <- Sys.which(tool_names)
  hits <- unname(hits[nzchar(hits)])
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[[1]]
}

mr_try_embed_fonts <- function(pdf_path) {
  gs <- mr_first_existing_tool(c("gs", "gswin64c", "gswin32c"))

  if (is.na(gs)) {
    return(list(
      attempted = FALSE,
      status = "skipped",
      message = "Ghostscript was not found; skipped optional font embedding."
    ))
  }

  embedded_path <- tempfile(fileext = ".pdf")

  gs_args <- c(
    "-q",
    "-dNOPAUSE",
    "-dBATCH",
    "-dSAFER",
    "-sDEVICE=pdfwrite",
    "-dPDFSETTINGS=/prepress",
    paste0(
      "-sOutputFile=",
      shQuote(embedded_path)
    ),
    "-c",
    shQuote(
      "<</EmbedAllFonts true /SubsetFonts true /NeverEmbed []>> setdistillerparams"
    ),
    "-f",
    shQuote(pdf_path)
  )

  output <- tryCatch(
    system2(
      gs,
      args = gs_args,
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) e
  )

  if (inherits(output, "error")) {
    return(list(
      attempted = TRUE,
      status = "failed",
      message = paste(
        "Ghostscript font embedding failed:",
        conditionMessage(output)
      )
    ))
  }

  status <- attr(output, "status")

  if (is.null(status)) {
    status <- 0L
  }

  if (!identical(as.integer(status), 0L)) {
    return(list(
      attempted = TRUE,
      status = "failed",
      message = paste(
        "Ghostscript font embedding exited with status",
        status,
        paste(output, collapse = " | ")
      )
    ))
  }

  if (
    !file.exists(embedded_path) ||
    is.na(file.info(embedded_path)$size) ||
    file.info(embedded_path)$size <= 0
  ) {
    return(list(
      attempted = TRUE,
      status = "failed",
      message = "Ghostscript did not produce a usable embedded PDF."
    ))
  }

  copied <- file.copy(
    embedded_path,
    pdf_path,
    overwrite = TRUE
  )

  if (!isTRUE(copied)) {
    return(list(
      attempted = TRUE,
      status = "failed",
      message = "Could not replace the source PDF with the embedded PDF."
    ))
  }

  list(
    attempted = TRUE,
    status = "completed",
    message = paste(
      "Ghostscript font embedding completed with",
      "EmbedAllFonts=true, SubsetFonts=true, and NeverEmbed=[]."
    )
  )
}


mr_audit_pdfinfo <- function(pdf_path) {
  pdfinfo <- mr_first_existing_tool("pdfinfo")

  if (is.na(pdfinfo)) {
    return(list(
      available = FALSE,
      status = "skipped",
      message = "pdfinfo was not found; skipped PDF metadata audit.",
      page_size = NA_character_
    ))
  }

  output <- tryCatch(
    system2(
      pdfinfo,
      args = shQuote(pdf_path),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) e
  )

  if (inherits(output, "error")) {
    return(list(
      available = TRUE,
      status = "failed",
      message = conditionMessage(output),
      page_size = NA_character_
    ))
  }

  page_size <- output[grepl("^Page size:", output)]
  page_size <- if (length(page_size) > 0) sub("^Page size:\\s*", "", page_size[[1]]) else NA_character_

  list(
    available = TRUE,
    status = "completed",
    message = "pdfinfo audit completed.",
    page_size = page_size
  )
}

mr_audit_pdffonts <- function(pdf_path) {
  pdffonts <- mr_first_existing_tool("pdffonts")

  if (is.na(pdffonts)) {
    return(list(
      available = FALSE,
      status = "skipped",
      message = "pdffonts was not found; skipped font audit.",
      fonts_not_embedded = NA_integer_
    ))
  }

  output <- tryCatch(
    system2(
      pdffonts,
      args = shQuote(pdf_path),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) e
  )

  if (inherits(output, "error")) {
    return(list(
      available = TRUE,
      status = "failed",
      message = conditionMessage(output),
      fonts_not_embedded = NA_integer_
    ))
  }

  if (length(output) < 3) {
    return(list(
      available = TRUE,
      status = "completed",
      message = "pdffonts audit completed; no font rows were reported.",
      fonts_not_embedded = 0L
    ))
  }

  header <- output[[1]]
  emb_start <- regexpr("\\bemb\\b", header)[[1]]
  sub_start <- regexpr("\\bsub\\b", header)[[1]]

  font_rows <- output[-c(1, 2)]
  font_rows <- font_rows[nzchar(trimws(font_rows))]

  fonts_not_embedded <- NA_integer_
  if (emb_start > 0 && sub_start > emb_start && length(font_rows) > 0) {
    emb_values <- trimws(substr(font_rows, emb_start, sub_start - 1))
    fonts_not_embedded <- sum(emb_values == "no", na.rm = TRUE)
  }

  msg <- if (is.na(fonts_not_embedded)) {
    "pdffonts audit completed, but embedding status could not be parsed."
  } else if (fonts_not_embedded == 0) {
    "pdffonts audit completed; all reported font rows are embedded or subset."
  } else {
    paste("pdffonts audit completed;", fonts_not_embedded, "font row(s) reported as not embedded.")
  }

  list(
    available = TRUE,
    status = "completed",
    message = msg,
    fonts_not_embedded = fonts_not_embedded
  )
}

mr_file_size_mb <- function(path) {
  file.info(path)$size / (1024^2)
}

mr_export_png <- function(plot, png_path, width_in, height_in, dpi = 600) {
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width_in,
    height = height_in,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )

  png_info <- file.info(png_path)
  if (!file.exists(png_path) || is.na(png_info$size) || png_info$size <= 0) {
    stop(
      "PNG export failed or produced an empty file: ",
      png_path,
      call. = FALSE
    )
  }

  invisible(png_path)
}

mr_export_pdf <- function(
    plot,
    filename,
    output_dir,
    figure,
    object_name,
    width_mm = mr_full_width_mm(),
    height_mm,
    embed_fonts = TRUE) {

  if (height_mm > mr_max_height_mm()) {
    stop(
      "Requested height for Figure ", figure, " is ", height_mm,
      " mm, which exceeds the manuscript maximum of ", mr_max_height_mm(), " mm.",
      call. = FALSE
    )
  }

  if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".pdf")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(output_dir, filename)

  width_in <- mr_mm_to_in(width_mm)
  height_in <- mr_mm_to_in(height_mm)

  opened <- FALSE
  grDevices::pdf(
    pdf_path,
    width = width_in,
    height = height_in,
    useDingbats = FALSE,
    onefile = FALSE
  )
  opened <- TRUE
  on.exit({
    if (opened) {
      grDevices::dev.off()
    }
  }, add = TRUE)

  print(plot)
  grDevices::dev.off()
  opened <- FALSE

  png_filename <- sub("\\.pdf$", ".png", filename, ignore.case = TRUE)
  png_path <- file.path(output_dir, png_filename)
  mr_export_png(
    plot = plot,
    png_path = png_path,
    width_in = width_in,
    height_in = height_in
  )

  embed_result <- if (embed_fonts) {
    mr_try_embed_fonts(pdf_path)
  } else {
    list(
      attempted = FALSE,
      status = "skipped",
      message = "Optional font embedding disabled."
    )
  }

  pdfinfo_result <- mr_audit_pdfinfo(pdf_path)
  pdffonts_result <- mr_audit_pdffonts(pdf_path)
  size_mb <- mr_file_size_mb(pdf_path)
  png_size_mb <- mr_file_size_mb(png_path)

  data.frame(
    figure = figure,
    object_name = object_name,
    file = filename,
    path = pdf_path,
    png_file = png_filename,
    png_path = png_path,
    width_mm = width_mm,
    height_mm = height_mm,
    width_in = width_in,
    height_in = height_in,
    file_size_mb = size_mb,
    png_file_size_mb = png_size_mb,
    under_10mb = is.finite(size_mb) && size_mb < 10,
    font_embedding_attempted = embed_result$attempted,
    font_embedding_status = embed_result$status,
    font_embedding_message = embed_result$message,
    pdfinfo_status = pdfinfo_result$status,
    pdfinfo_message = pdfinfo_result$message,
    pdfinfo_page_size = pdfinfo_result$page_size,
    pdffonts_status = pdffonts_result$status,
    pdffonts_message = pdffonts_result$message,
    pdffonts_not_embedded = pdffonts_result$fonts_not_embedded,
    stringsAsFactors = FALSE
  )
}

mr_write_audit_csv <- function(audit_df, output_dir) {
  out <- file.path(output_dir, "manuscript_ready_pdf_audit.csv")
  utils::write.csv(audit_df, out, row.names = FALSE)
  message("Wrote manuscript-ready audit: ", out)
  invisible(out)
}

mr_install_ggsave_skipper <- function(env = .GlobalEnv) {
  old <- list(
    had_global_ggsave = exists("ggsave", envir = env, inherits = FALSE),
    global_ggsave = NULL
  )

  if (old$had_global_ggsave) {
    old$global_ggsave <- get("ggsave", envir = env, inherits = FALSE)
  }

  assign(
    "ggsave",
    function(filename, plot = NULL, ...) {
      target <- tryCatch(as.character(filename)[[1]], error = function(e) "<unknown file>")
      message("Manuscript-ready wrapper skipped original ggsave side effect: ", target)
      invisible(filename)
    },
    envir = env
  )

  old
}

mr_restore_ggsave <- function(state, env = .GlobalEnv) {
  if (isTRUE(state$had_global_ggsave)) {
    assign("ggsave", state$global_ggsave, envir = env)
  } else if (exists("ggsave", envir = env, inherits = FALSE)) {
    rm("ggsave", envir = env)
  }
  invisible(TRUE)
}

mr_source_script <- function(script, env = .GlobalEnv) {
  caller <- parent.frame()
  if (exists("cfg", envir = caller, inherits = FALSE)) {
    assign("cfg", get("cfg", envir = caller, inherits = FALSE), envir = env)
  }

  message("\n------------------------------------------------------------")
  message("Sourcing: ", script)
  message("------------------------------------------------------------")

  tryCatch(
    source(script, local = env),
    error = function(e) {
      stop(
        "\nPipeline failed while sourcing: ", script,
        "\n\nOriginal error:\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

message("Loaded 80_publication_export_helpers.R")
