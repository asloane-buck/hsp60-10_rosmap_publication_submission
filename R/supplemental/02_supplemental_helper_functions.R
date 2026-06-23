# ============================================================
# General helpers
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

pick_col <- function(df, candidates, required = TRUE, label = NULL) {
  if (is.null(label)) {
    label <- "requested column"
  }

  available_cols <- colnames(df)
  selected <- candidates[candidates %in% available_cols][1]

  if (!is.na(selected)) {
    return(selected)
  }

  if (isTRUE(required)) {
    stop(
      paste0(
        "Missing required column for conceptual variable: ", label, "\n",
        "Candidate column names tried: ", paste(candidates, collapse = ", "), "\n",
        "Available column names: ", paste(available_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  NA_character_
}

clean_gene <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_to_upper(x)
  x <- stringr::str_replace(x, "\\.\\d+$", "")
  x[x == ""] <- NA_character_
  x
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

# ============================================================
# Plot saving helpers
# ============================================================

is_final_supplemental_figure <- function(filename_base, output_dir) {
  if (!exists("SUPP_FINAL_FIGURE_HEIGHT_MM", inherits = TRUE)) {
    return(FALSE)
  }

  filename_base %in% names(SUPP_FINAL_FIGURE_HEIGHT_MM) &&
    identical(
      normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      normalizePath(figures_dir, winslash = "/", mustWork = FALSE)
    )
}

supplemental_final_export_size <- function(filename_base) {
  if (!filename_base %in% names(SUPP_FINAL_FIGURE_HEIGHT_MM)) {
    stop("No final supplemental export height configured for: ", filename_base, call. = FALSE)
  }

  height_mm <- unname(SUPP_FINAL_FIGURE_HEIGHT_MM[[filename_base]])
  if (height_mm > SUPP_MAX_HEIGHT_MM) {
    stop(
      "Configured final supplemental export height exceeds journal maximum for ",
      filename_base,
      ": ",
      height_mm,
      " mm.",
      call. = FALSE
    )
  }

  list(
    width = SUPP_FULL_WIDTH_IN,
    height = height_mm * MM_TO_IN,
    width_mm = SUPP_FULL_WIDTH_MM,
    height_mm = height_mm
  )
}

supplemental_cairo_pdf_available <- local({
  available <- NULL

  function() {
    if (!is.null(available)) {
      return(available)
    }

    if (!isTRUE(capabilities("cairo"))) {
      available <<- FALSE
      return(available)
    }

    test_pdf <- tempfile(fileext = ".pdf")
    warned <- FALSE
    available <<- tryCatch(
      {
        withCallingHandlers(
          {
            grDevices::cairo_pdf(test_pdf, width = 1, height = 1)
            grDevices::dev.off()
          },
          warning = function(w) {
            warned <<- TRUE
            invokeRestart("muffleWarning")
          }
        )
        !warned && file.exists(test_pdf)
      },
      error = function(e) FALSE
    )
    if (file.exists(test_pdf)) {
      unlink(test_pdf)
    }

    available
  }
})

supplemental_pdf_device <- function() {
  if (isTRUE(supplemental_cairo_pdf_available())) {
    return(grDevices::cairo_pdf)
  }

  "pdf"
}

save_plot_set <- function(plot, filename_base, width, height, output_dir) {
  ensure_dir(output_dir)

  if (is_final_supplemental_figure(filename_base, output_dir)) {
    final_size <- supplemental_final_export_size(filename_base)
    width <- final_size$width
    height <- final_size$height
    message(
      "Applying journal final export size for ",
      filename_base,
      ": ",
      final_size$width_mm,
      " x ",
      final_size$height_mm,
      " mm."
    )
  }

  pdf_path <- file.path(output_dir, paste0(filename_base, ".pdf"))
  png_path <- file.path(output_dir, paste0(filename_base, ".png"))

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = supplemental_pdf_device(),
    bg = "white"
  )

  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "png",
    dpi = 500,
    bg = "white"
  )

  message("Saved plot set: ", pdf_path, " and ", png_path)
  invisible(c(pdf = pdf_path, png = png_path))
}

find_poppler_tool <- function(tool_name) {
  from_path <- Sys.which(tool_name)
  candidates <- c(
    unname(from_path),
    file.path(
      path.expand("~"),
      ".cache",
      "codex-runtimes",
      "codex-primary-runtime",
      "dependencies",
      "bin",
      tool_name
    )
  )

  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates[file.exists(candidates)][1] %||% NA_character_
}

read_pdf_page_size_mm <- function(pdf_path) {
  pdfinfo_bin <- find_poppler_tool("pdfinfo")
  if (is.na(pdfinfo_bin)) {
    return(list(width_mm = NA_real_, height_mm = NA_real_, method = "pdfinfo unavailable"))
  }

  info <- tryCatch(
    system2(pdfinfo_bin, args = shQuote(pdf_path), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  page_size_line <- info[stringr::str_detect(info, "^Page size:")][1]
  dims <- stringr::str_match(page_size_line, "Page size:\\s+([0-9.]+)\\s+x\\s+([0-9.]+)\\s+pts")

  if (any(is.na(dims[1, 2:3]))) {
    return(list(width_mm = NA_real_, height_mm = NA_real_, method = "pdfinfo parse failed"))
  }

  list(
    width_mm = as.numeric(dims[1, 2]) * 25.4 / 72,
    height_mm = as.numeric(dims[1, 3]) * 25.4 / 72,
    method = "pdfinfo"
  )
}

read_pdf_font_embedding_status <- function(pdf_path) {
  pdffonts_bin <- find_poppler_tool("pdffonts")
  if (is.na(pdffonts_bin)) {
    return("font verification unavailable: pdffonts not found")
  }

  fonts <- tryCatch(
    system2(pdffonts_bin, args = shQuote(pdf_path), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  font_rows <- fonts[!stringr::str_detect(fonts, "^name\\s+type|^-+$")]
  font_rows <- font_rows[nzchar(stringr::str_trim(font_rows))]

  if (length(font_rows) == 0) {
    return("no fonts detected by pdffonts")
  }

  embedded <- vapply(
    font_rows,
    function(row_i) {
      emb_match <- stringr::str_match(
        row_i,
        "\\s+(yes|no)\\s+(yes|no)\\s+(yes|no)\\s+[0-9]+\\s+[0-9]+\\s*$"
      )
      if (is.na(emb_match[1, 2])) {
        return(NA_character_)
      }
      emb_match[1, 2]
    },
    character(1)
  )

  if (all(embedded == "yes", na.rm = TRUE) && !any(is.na(embedded))) {
    "all fonts embedded"
  } else {
    paste0("font embedding warning: emb=", paste(unique(embedded), collapse = "/"))
  }
}

audit_manuscript_ready_supplemental_pdfs <- function(pdf_dir = manuscript_ready_pdf_dir,
                                                     audit_path = file.path(audits_dir, "manuscript_ready_supplemental_pdf_dimension_audit.csv"),
                                                     stop_on_fail = TRUE) {
  final_pdf_names <- paste0(names(SUPP_FINAL_FIGURE_HEIGHT_MM), ".pdf")
  pdf_paths <- file.path(pdf_dir, final_pdf_names)
  missing_paths <- pdf_paths[!file.exists(pdf_paths)]

  if (length(missing_paths) > 0) {
    stop(
      "Cannot audit manuscript-ready supplemental PDFs; missing files: ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
  }

  audit <- purrr::map_dfr(pdf_paths, function(path_i) {
    dims <- read_pdf_page_size_mm(path_i)
    font_status <- read_pdf_font_embedding_status(path_i)
    size_mb <- file.info(path_i)$size / (1024^2)
    width_pass <- !is.na(dims$width_mm) &&
      dims$width_mm >= (SUPP_FULL_WIDTH_MM - 1) &&
      dims$width_mm <= (SUPP_FULL_WIDTH_MM + 1)
    height_pass <- !is.na(dims$height_mm) && dims$height_mm <= SUPP_MAX_HEIGHT_MM
    size_pass <- !is.na(size_mb) && size_mb < 10

    tibble::tibble(
      pdf_name = basename(path_i),
      width_mm = dims$width_mm,
      height_mm = dims$height_mm,
      file_size_MB = size_mb,
      dimension_method = dims$method,
      font_embedding_status = font_status,
      pass_fail = ifelse(width_pass && height_pass && size_pass, "PASS", "FAIL")
    )
  })

  readr::write_csv(audit, audit_path)

  if (isTRUE(stop_on_fail) && any(audit$pass_fail != "PASS")) {
    print(audit, n = Inf)
    stop("One or more manuscript-ready supplemental PDFs failed journal dimension audit.", call. = FALSE)
  }

  invisible(audit)
}

save_panel_set <- function(plot, filename_base, width = 6, height = 4, output_dir = panels_dir) {
  save_plot_set(
    plot = plot,
    filename_base = filename_base,
    width = width,
    height = height,
    output_dir = output_dir
  )
}

theme_supplement <- function(base_size = 8.5) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = 0.3),
      axis.text = ggplot2::element_text(size = base_size * 0.85, color = "grey20"),
      axis.title = ggplot2::element_text(size = base_size * 0.95, color = "grey15"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = base_size * 0.9, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size * 0.85),
      strip.text = ggplot2::element_text(size = base_size * 0.9, face = "bold")
    )
}

panel_label <- function(plot, label) {
  if (inherits(plot, "patchwork") && requireNamespace("cowplot", quietly = TRUE)) {
    return(
      cowplot::ggdraw() +
        cowplot::draw_plot(plot, x = 0.035, y = 0, width = 0.965, height = 0.98) +
        cowplot::draw_label(
          label,
          x = 0,
          y = 1,
          hjust = 0,
          vjust = 1,
          fontface = "bold",
          size = 9.5
        )
    )
  }

  plot +
    ggplot2::labs(tag = label) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 9.5),
      plot.tag.position = c(0.01, 0.99)
    )
}

# ============================================================
# Matrix helpers
# ============================================================

coerce_expression_matrix <- function(object, label = "expression matrix") {
  if (is.null(object)) {
    return(NULL)
  }

  if (inherits(object, "matrix")) {
    mat <- object
  } else if (is.data.frame(object)) {
    df <- as.data.frame(object)
    gene_col <- pick_col(
      df,
      c(
        "gene", "gene_symbol", "symbol", "hgnc_symbol",
        "Gene", "GeneSymbol", "SYMBOL", "GENE",
        "protein", "Protein"
      ),
      required = FALSE,
      label = paste(label, "gene identifier")
    )

    if (!is.na(gene_col)) {
      row_ids <- df[[gene_col]]
      df[[gene_col]] <- NULL
    } else {
      row_ids <- rownames(df)
    }

    numeric_cols <- vapply(df, is.numeric, logical(1))
    if (!any(numeric_cols)) {
      converted_df <- lapply(df, function(col_i) {
        suppressWarnings(as.numeric(as.character(col_i)))
      })
      numeric_cols <- vapply(converted_df, function(col_i) {
        sum(!is.na(col_i)) > 0
      }, logical(1))
      df[numeric_cols] <- converted_df[numeric_cols]
    }

    if (!any(numeric_cols)) {
      stop(
        "Could not identify numeric sample columns in ", label, ".",
        call. = FALSE
      )
    }

    mat <- as.matrix(df[, numeric_cols, drop = FALSE])
    rownames(mat) <- make.unique(as.character(row_ids))
  } else if (
    requireNamespace("SummarizedExperiment", quietly = TRUE) &&
      methods::is(object, "SummarizedExperiment")
  ) {
    mat <- SummarizedExperiment::assay(object)
  } else {
    stop(
      "Unsupported object type for ", label, ": ",
      paste(class(object), collapse = ", "),
      call. = FALSE
    )
  }

  suppressWarnings(storage.mode(mat) <- "numeric")
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("row_", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("sample_", seq_len(ncol(mat)))
  }

  mat
}

collapse_duplicate_gene_rows <- function(mat, label = "matrix") {
  if (is.null(mat) || nrow(mat) == 0) {
    return(mat)
  }

  duplicated_genes <- unique(rownames(mat)[duplicated(rownames(mat))])
  if (length(duplicated_genes) == 0) {
    return(mat)
  }

  message(
    "Collapsing duplicate gene rows in ", label, ": ",
    length(duplicated_genes), " duplicated gene symbols"
  )

  row_groups <- split(seq_len(nrow(mat)), rownames(mat))
  collapsed <- do.call(
    rbind,
    lapply(row_groups, function(idx) {
      values <- colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
      values[is.nan(values)] <- NA_real_
      values
    })
  )

  rownames(collapsed) <- names(row_groups)
  collapsed
}

infer_and_orient_matrix <- function(
    object,
    gene_symbols,
    sample_ids = character(),
    label = "matrix",
    required = FALSE) {
  if (is.null(object)) {
    if (isTRUE(required)) {
      stop(label, " is required but was not provided.", call. = FALSE)
    }
    return(NULL)
  }

  mat <- coerce_expression_matrix(object, label = label)

  clean_interactor_genes <- unique(stats::na.omit(clean_gene(gene_symbols)))
  clean_row_ids <- clean_gene(rownames(mat))
  clean_col_ids <- clean_gene(colnames(mat))
  sample_ids <- unique(stats::na.omit(as.character(sample_ids)))

  row_gene_hits <- sum(clean_row_ids %in% clean_interactor_genes)
  col_gene_hits <- sum(clean_col_ids %in% clean_interactor_genes)
  row_sample_hits <- if (length(sample_ids) > 0) sum(rownames(mat) %in% sample_ids) else 0
  col_sample_hits <- if (length(sample_ids) > 0) sum(colnames(mat) %in% sample_ids) else 0

  genes_as_rows_score <- row_gene_hits + col_sample_hits
  genes_as_cols_score <- col_gene_hits + row_sample_hits

  orientation <- NA_character_
  if (genes_as_rows_score > genes_as_cols_score) {
    orientation <- "genes_rows_samples_columns"
  } else if (genes_as_cols_score > genes_as_rows_score) {
    orientation <- "genes_columns_samples_rows"
  } else if (row_gene_hits > col_gene_hits) {
    orientation <- "genes_rows_samples_columns"
  } else if (col_gene_hits > row_gene_hits) {
    orientation <- "genes_columns_samples_rows"
  }

  diagnostics <- tibble::tibble(
    matrix_label = label,
    n_rows = nrow(mat),
    n_columns = ncol(mat),
    row_gene_hits = row_gene_hits,
    column_gene_hits = col_gene_hits,
    row_sample_hits = row_sample_hits,
    column_sample_hits = col_sample_hits,
    selected_orientation = orientation
  )

  message("Orientation diagnostics for ", label, ":")
  print(diagnostics)

  if (is.na(orientation)) {
    message("Could not infer orientation for ", label, "; matrix will be unavailable.")
    return(list(matrix = NULL, diagnostics = diagnostics))
  }

  if (identical(orientation, "genes_columns_samples_rows")) {
    mat <- t(mat)
  }

  rownames(mat) <- clean_gene(rownames(mat))
  keep_rows <- !is.na(rownames(mat)) & nzchar(rownames(mat))
  mat <- mat[keep_rows, , drop = FALSE]
  mat <- collapse_duplicate_gene_rows(mat, label = label)

  message(
    "Selected orientation for ", label, ": ", orientation,
    "; final dimensions: ", paste(dim(mat), collapse = " x ")
  )

  list(matrix = mat, diagnostics = diagnostics)
}

write_skip_audit <- function(path, panel, reason) {
  readr::write_csv(
    tibble::tibble(
      panel = panel,
      status = "skipped",
      reason = reason
    ),
    path
  )
  message("Wrote skip audit: ", path)
  invisible(path)
}
