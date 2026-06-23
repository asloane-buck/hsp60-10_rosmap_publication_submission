############################################################
## 01_utils.R
## Reusable functions only. No analysis objects are built here.
############################################################

clean_gene_values <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "\\.\\d+$", "")
  x <- stringr::str_replace(x, "\\|.*$", "")
  x <- stringr::str_trim(x)
  x
}

canonical_gene_symbol <- function(x) {
  x <- clean_gene_values(x)
  dplyr::case_when(
    x == "ATP5B" ~ "ATP5F1B",
    x == "ATP5A1" ~ "ATP5F1A",
    x == "ATP5C1" ~ "ATP5F1C",
    TRUE ~ x
  )
}

display_gene_symbol <- function(x) {
  x <- canonical_gene_symbol(x)
  dplyr::case_when(
    x == "ATP5F1B" ~ "ATP5B",
    x == "ATP5F1A" ~ "ATP5A1",
    x == "ATP5F1C" ~ "ATP5C1",
    TRUE ~ x
  )
}

clean_gene_symbols <- function(x) {
  x <- canonical_gene_symbol(x)
  sort(unique(x[!is.na(x) & x != ""]))
}

collapse_duplicate_rows <- function(mat) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  rn <- canonical_gene_symbol(rownames(mat))
  keep <- !is.na(rn) & rn != ""
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]
  rownames(mat) <- rn

  split_idx <- split(seq_along(rn), rn)
  out <- lapply(split_idx, function(idx) {
    if (length(idx) == 1) {
      mat[idx, , drop = FALSE]
    } else {
      matrix(colMeans(mat[idx, , drop = FALSE], na.rm = TRUE), nrow = 1)
    }
  })
  out <- do.call(rbind, out)
  rownames(out) <- names(split_idx)
  as.matrix(out)
}

read_gene_matrix_csv <- function(file, gene_col = 1) {
  x <- readr::read_csv(file, show_col_types = FALSE, name_repair = "minimal")
  x <- as.data.frame(x, check.names = FALSE)

  if (is.character(gene_col)) {
    if (!gene_col %in% colnames(x)) {
      stop(
        "Missing gene column: ", gene_col,
        ". Available columns: ", paste(colnames(x), collapse = ", "),
        call. = FALSE
      )
    }
    gene <- x[[gene_col]]
    x[[gene_col]] <- NULL
  } else {
    gene <- x[[gene_col]]
    x[[gene_col]] <- NULL
  }

  x[] <- lapply(x, function(v) suppressWarnings(as.numeric(v)))
  mat <- as.matrix(x)
  rownames(mat) <- clean_gene_values(gene)
  collapse_duplicate_rows(mat)
}

first_existing <- function(paths, label = "file") {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) {
    stop("Could not find ", label, ". Tried:\n", paste(paths, collapse = "\n"), call. = FALSE)
  }
  hit
}

safe_num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))

available_cols_msg <- function(df) paste(colnames(df), collapse = ", ")

available_objects_msg <- function(env = .GlobalEnv) paste(ls(envir = env), collapse = ", ")

require_objects <- function(objects, context = "object check", env = .GlobalEnv) {
  missing <- objects[!vapply(objects, exists, logical(1), envir = env)]
  if (length(missing) > 0) {
    stop(
      context, " missing required object(s): ", paste(missing, collapse = ", "),
      ". Available objects: ", available_objects_msg(env),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

require_columns <- function(df, cols, context = "column check") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(
      context, " missing required column(s): ", paste(missing, collapse = ", "),
      ". Available columns: ", available_cols_msg(df),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

require_values <- function(df, col, values, context = "value check") {
  require_columns(df, col, context)
  missing <- setdiff(values, unique(as.character(df[[col]])))
  if (length(missing) > 0) {
    stop(
      context, " missing expected ", col, " value(s): ", paste(missing, collapse = ", "),
      ". Available values: ", paste(sort(unique(as.character(df[[col]]))), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

alias_sets <- list(
  sample_id = c("sample_id", "SampleID", "sampleid", "specimen_id"),
  individual_id = c("individual_id", "IndividualID", "individualID", "individualid", "projid"),
  age_death = c("age_death", "age", "age_at_death", "age_meta"),
  pmi = c("pmi", "pmi_meta", "postmortem_interval", "pmi_hours", "pmihours"),
  sex = c("sex", "sex_label", "msex"),
  braak = c("braak_num", "braak", "braaksc", "Braak"),
  cerad = c("cerad_num", "cerad", "ceradsc", "CERAD"),
  diagnosis = c("diagnosis", "EmoryStrictDx.2019", "emory_strict_dx_2019", "cogdx"),
  pathway = c("Pathway", "pathway", "PathwayName", "pathway_name", "gene_set", "set_name")
)

pick_existing <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

require_alias_col <- function(df, candidates, label, context = "alias check") {
  hit <- pick_existing(df, candidates)
  if (is.na(hit)) {
    stop(
      context, " could not find ", label, " column. Tried: ",
      paste(candidates, collapse = ", "),
      ". Available columns: ", available_cols_msg(df),
      call. = FALSE
    )
  }
  hit
}

add_alias_column <- function(df, target, candidates, context = "alias check", required = TRUE) {
  if (target %in% colnames(df)) return(df)
  hit <- pick_existing(df, candidates)
  if (is.na(hit)) {
    if (isTRUE(required)) {
      stop(
        context, " could not create ", target, ". Tried aliases: ",
        paste(candidates, collapse = ", "),
        ". Available columns: ", available_cols_msg(df),
        call. = FALSE
      )
    }
    df[[target]] <- NA
    return(df)
  }
  df[[target]] <- df[[hit]]
  df
}

coalesce_alias_chr <- function(data, candidates) {
  present <- candidates[candidates %in% colnames(data)]
  if (length(present) == 0) return(rep(NA_character_, nrow(data)))
  out <- rep(NA_character_, nrow(data))
  for (cc in present) out <- dplyr::coalesce(out, as.character(data[[cc]]))
  out
}

coalesce_alias_num <- function(data, candidates) {
  safe_num(coalesce_alias_chr(data, candidates))
}

standardize_common_aliases <- function(df, context = "alias standardization",
                                       require_ids = FALSE, require_diagnosis = FALSE) {
  df <- add_alias_column(df, "sample_id", alias_sets$sample_id, context, required = require_ids)
  df <- add_alias_column(df, "individual_id", alias_sets$individual_id, context, required = require_ids)
  df <- add_alias_column(df, "diagnosis", alias_sets$diagnosis, context, required = require_diagnosis)
  df <- add_alias_column(df, "age_death", alias_sets$age_death, context, required = FALSE)
  df <- add_alias_column(df, "pmi", alias_sets$pmi, context, required = FALSE)
  df <- add_alias_column(df, "sex", alias_sets$sex, context, required = FALSE)
  df <- add_alias_column(df, "braak", alias_sets$braak, context, required = FALSE)
  df <- add_alias_column(df, "cerad", alias_sets$cerad, context, required = FALSE)
  df
}

join_key_pairs <- function(by) {
  if (is.null(names(by))) {
    tibble::tibble(left = by, right = by)
  } else {
    right <- unname(by)
    left <- names(by)
    left[left == ""] <- right[left == ""]
    tibble::tibble(left = left, right = right)
  }
}

checked_left_join <- function(x, y, by, label = "left_join", suffix = c(".x", ".y"), ...) {
  keys <- join_key_pairs(by)
  require_columns(x, keys$left, paste0(label, " left input"))
  require_columns(y, keys$right, paste0(label, " right input"))

  message(label, " join keys: ", paste(paste0(keys$left, "=", keys$right), collapse = ", "))
  message(label, " left dimensions before join: ", paste(dim(x), collapse = " x "))
  message(label, " right dimensions before join: ", paste(dim(y), collapse = " x "))

  left_keys <- dplyr::distinct(x, dplyr::across(dplyr::all_of(keys$left)))
  right_keys <- dplyr::distinct(y, dplyr::across(dplyr::all_of(keys$right)))
  names(right_keys) <- keys$left
  names(left_keys) <- keys$left

  unmatched_left <- dplyr::anti_join(left_keys, right_keys, by = keys$left)
  unmatched_right <- dplyr::anti_join(right_keys, left_keys, by = keys$left)
  message(label, " unmatched left keys: ", nrow(unmatched_left))
  message(label, " unmatched right keys: ", nrow(unmatched_right))

  out <- dplyr::left_join(x, y, by = by, suffix = suffix, ...)
  message(label, " dimensions after join: ", paste(dim(out), collapse = " x "))
  out
}

safe_z <- function(x) {
  x <- safe_num(x)
  if (sum(is.finite(x)) < 5 || length(unique(x[is.finite(x)])) < 3) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

zscore_rows <- function(mat) {
  z <- t(scale(t(as.matrix(mat))))
  z[!is.finite(z)] <- NA_real_
  z
}

score_pathway_mean_z <- function(mat, genes, min_genes = 2) {
  genes_use <- intersect(clean_gene_symbols(genes), rownames(mat))
  if (length(genes_use) < min_genes) return(rep(NA_real_, ncol(mat)))
  colMeans(zscore_rows(mat[genes_use, , drop = FALSE]), na.rm = TRUE)
}

compute_pathway_scores <- function(mat, pathway_list, sample_col = "sample_id", method = "mean_z") {
  if (!identical(method, "mean_z")) {
    stop("Unsupported pathway score method: ", method, ". Available method: mean_z", call. = FALSE)
  }
  scores <- purrr::map(pathway_list, ~ score_pathway_mean_z(mat, .x))
  out <- as.data.frame(scores, check.names = FALSE)
  out[[sample_col]] <- colnames(mat)
  out[, c(sample_col, names(pathway_list)), drop = FALSE]
}

summarize_pathway_coverage <- function(mat, pathway_list) {
  tibble(
    Pathway = names(pathway_list),
    pathway = Pathway,
    n_total = purrr::map_int(pathway_list, ~ length(clean_gene_symbols(.x))),
    n_detected = purrr::map_int(pathway_list, ~ length(intersect(clean_gene_symbols(.x), rownames(mat)))),
    pct_detected = round(100 * n_detected / pmax(n_total, 1), 1)
  )
}

percentile01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) <= 1) return(out)
  out[ok] <- 100 * dplyr::percent_rank(x[ok])
  out
}

covar_is_usable_numeric <- function(x, min_n = 20) {
  x <- safe_num(x)
  sum(is.finite(x)) >= min_n && length(unique(x[is.finite(x)])) >= 3
}

covar_is_usable_factor <- function(x, n_total) {
  x <- as.factor(x)
  tab <- table(x, useNA = "no")
  length(tab) >= 2 &&
    length(tab) <= max(2, floor(n_total / 4)) &&
    all(tab >= 3)
}

prep_covariates <- function(meta_df) {
  df <- meta_df

  sex_col <- pick_existing(df, c("sex_label", alias_sets$sex, "gender"))
  age_col <- pick_existing(df, c(alias_sets$age_death, "age_num"))
  pmi_col <- pick_existing(df, c(alias_sets$pmi, "pmi_num"))
  rin_col <- pick_existing(df, c("rin", "rin_score", "rna_integrity_number", "rin_num"))
  batch_col <- pick_existing(df, c("batch", "batch_id", "batch_channel", "batch.channel", "study_batch"))

  df$sex_factor <- if (!is.na(sex_col)) {
    case_when(
      as.character(df[[sex_col]]) %in% c("0", "Female", "female", "F", "f") ~ "Female",
      as.character(df[[sex_col]]) %in% c("1", "Male", "male", "M", "m") ~ "Male",
      TRUE ~ as.character(df[[sex_col]])
    ) |> as.factor()
  } else {
    factor(NA)
  }

  df$age_num <- if (!is.na(age_col)) safe_num(df[[age_col]]) else NA_real_
  df$pmi_num <- if (!is.na(pmi_col)) safe_num(df[[pmi_col]]) else NA_real_
  df$rin_num <- if (!is.na(rin_col)) safe_num(df[[rin_col]]) else NA_real_

  if (!is.na(batch_col) && covar_is_usable_factor(df[[batch_col]], nrow(df))) {
    df$batch_factor <- as.factor(df[[batch_col]])
  } else {
    df$batch_factor <- factor(NA)
  }

  braak_col <- pick_existing(df, alias_sets$braak)
  if (!is.na(braak_col)) {
    df$braak_num <- safe_num(df[[braak_col]])
  } else {
    df$braak_num <- NA_real_
  }

  cerad_col <- pick_existing(df, alias_sets$cerad)
  if (!is.na(cerad_col)) {
    df$cerad_num <- safe_num(df[[cerad_col]])
  } else {
    df$cerad_num <- NA_real_
  }

  df$braak_num_std <- safe_z(df$braak_num)
  df$cerad_num_std <- safe_z(df$cerad_num)

  df
}

add_analysis_meta_covars <- function(meta_df, id_col, analysis_meta) {
  require_columns(meta_df, id_col, "add_analysis_meta_covars metadata")
  analysis_meta <- standardize_common_aliases(
    analysis_meta,
    context = "add_analysis_meta_covars analysis_meta",
    require_ids = TRUE
  )
  require_columns(analysis_meta, "individual_id", "add_analysis_meta_covars analysis_meta")

  covar_meta <- analysis_meta |> distinct(individual_id, .keep_all = TRUE)

  joined <- checked_left_join(
    meta_df,
    covar_meta,
    by = setNames("individual_id", id_col),
    label = "analysis metadata covariate",
    suffix = c("", ".analysis_meta")
  )

  joined$braak_num <- coalesce_alias_num(
    joined,
    c("braak_num", "braak_num.analysis_meta", "braak", "braak.analysis_meta", "braaksc", "braaksc.analysis_meta", "Braak")
  )
  joined$cerad_num <- coalesce_alias_num(
    joined,
    c("cerad_num", "cerad_num.analysis_meta", "cerad", "cerad.analysis_meta", "ceradsc", "ceradsc.analysis_meta", "CERAD")
  )
  joined$sex_label <- coalesce_alias_chr(
    joined,
    c("sex_label", "sex_label.analysis_meta", "sex", "sex.analysis_meta", "msex", "msex.analysis_meta")
  )
  joined
}

select_usable_covars <- function(meta_df, candidates, min_n = 20) {
  keep <- character(0)
  for (v in candidates) {
    if (!v %in% colnames(meta_df)) next
    if (is.numeric(meta_df[[v]]) || is.integer(meta_df[[v]])) {
      if (covar_is_usable_numeric(meta_df[[v]], min_n = min_n)) keep <- c(keep, v)
    } else {
      if (covar_is_usable_factor(meta_df[[v]], nrow(meta_df))) keep <- c(keep, v)
    }
  }
  unique(keep)
}

residualize_matrix <- function(mat, meta_df, sample_col, covars, min_n = 20) {
  mat <- as.matrix(mat)
  require_columns(meta_df, sample_col, "residualize_matrix metadata")
  meta_df <- meta_df[match(colnames(mat), meta_df[[sample_col]]), , drop = FALSE]
  if (!all(colnames(mat) == meta_df[[sample_col]])) {
    stop(
      "residualize_matrix could not align matrix columns to metadata sample column: ",
      sample_col,
      ". Available columns: ", available_cols_msg(meta_df),
      call. = FALSE
    )
  }

  covars <- intersect(covars, colnames(meta_df))
  if (length(covars) == 0) {
    warning("No usable covariates found; returning original matrix.", call. = FALSE)
    return(mat)
  }

  design_df <- meta_df |> select(all_of(covars)) |> mutate(across(where(is.character), as.factor))
  form <- as.formula(paste("y ~", paste(covars, collapse = " + ")))

  out <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)

  for (i in seq_len(nrow(mat))) {
    y <- as.numeric(mat[i, ])
    df <- bind_cols(tibble(y = y), design_df)
    keep <- complete.cases(df)

    if (sum(keep) < min_n || length(unique(y[keep])) < 3) next

    fit <- tryCatch(lm(form, data = df[keep, , drop = FALSE]), error = function(e) NULL)
    if (is.null(fit)) next

    adjusted <- rep(NA_real_, length(y))
    adjusted[keep] <- resid(fit) + mean(y[keep], na.rm = TRUE)
    out[i, ] <- adjusted
  }

  out
}

fit_stage_effects <- function(expr_mat, meta_df, sample_col, stage_col, stage_levels, modality_name, min_n = 10) {
  expr_mat <- as.matrix(expr_mat)
  require_columns(meta_df, c(sample_col, stage_col), "fit_stage_effects metadata")
  meta_df <- meta_df[match(colnames(expr_mat), meta_df[[sample_col]]), , drop = FALSE]
  if (!all(colnames(expr_mat) == meta_df[[sample_col]])) {
    stop(
      "fit_stage_effects could not align matrix columns to metadata sample column: ",
      sample_col,
      ". Available columns: ", available_cols_msg(meta_df),
      call. = FALSE
    )
  }

  fit_gene_contrast <- function(y, ref_stage, comp_stage) {
    df <- tibble(
      y = as.numeric(y),
      stage = factor(as.character(meta_df[[stage_col]]), levels = c(ref_stage, comp_stage))
    ) |>
      filter(!is.na(stage), is.finite(y))

    if (nrow(df) < min_n || length(unique(df$stage)) < 2) return(c(effect = NA_real_, p = NA_real_))

    fit <- tryCatch(lm(y ~ stage, data = df), error = function(e) NULL)
    if (is.null(fit)) return(c(effect = NA_real_, p = NA_real_))

    hit <- broom::tidy(fit) |> filter(term == paste0("stage", comp_stage))
    if (nrow(hit) != 1) return(c(effect = NA_real_, p = NA_real_))

    c(effect = hit$estimate, p = hit$p.value)
  }

  early <- t(apply(expr_mat, 1, fit_gene_contrast, ref_stage = stage_levels[1], comp_stage = stage_levels[2]))
  late  <- t(apply(expr_mat, 1, fit_gene_contrast, ref_stage = stage_levels[2], comp_stage = stage_levels[3]))
  total <- t(apply(expr_mat, 1, fit_gene_contrast, ref_stage = stage_levels[1], comp_stage = stage_levels[3]))

  tibble(
    gene = rownames(expr_mat),
    early_shift = early[, "effect"],
    early_p = early[, "p"],
    late_shift = late[, "effect"],
    late_p = late[, "p"],
    total_shift = total[, "effect"],
    total_p = total[, "p"]
  ) |>
    rename_with(~ paste0(.x, "__", modality_name),
                c(early_shift, early_p, late_shift, late_p, total_shift, total_p))
}

fit_adjusted_braak_beta <- function(mat, meta_df, sample_col, genes, covars,
                                    adjust_for_cerad = TRUE, min_n = 30) {
  mat <- as.matrix(mat)
  require_columns(meta_df, c(sample_col, "braak_num_std"), "fit_adjusted_braak_beta metadata")
  meta_df <- meta_df[match(colnames(mat), meta_df[[sample_col]]), , drop = FALSE]
  if (!all(colnames(mat) == meta_df[[sample_col]])) {
    stop(
      "fit_adjusted_braak_beta could not align matrix columns to metadata sample column: ",
      sample_col,
      ". Available columns: ", available_cols_msg(meta_df),
      call. = FALSE
    )
  }

  genes <- intersect(genes, rownames(mat))
  covars <- intersect(covars, colnames(meta_df))

  if (isTRUE(adjust_for_cerad) && "cerad_num_std" %in% colnames(meta_df)) {
    covars <- unique(c(covars, "cerad_num_std"))
  }

  covars <- covars[
    map_lgl(covars, function(v) {
      if (is.numeric(meta_df[[v]]) || is.integer(meta_df[[v]])) {
        covar_is_usable_numeric(meta_df[[v]], min_n = min_n)
      } else {
        covar_is_usable_factor(meta_df[[v]], nrow(meta_df))
      }
    })
  ]

  map_dfr(genes, function(g) {
    y <- safe_z(mat[g, ])

    df <- tibble(
      abundance = y,
      braak_num_std = meta_df$braak_num_std
    ) |>
      bind_cols(meta_df |> select(any_of(covars))) |>
      mutate(across(where(is.character), as.factor))

    keep <- complete.cases(df[, c("abundance", "braak_num_std", covars), drop = FALSE])

    if (sum(keep) < min_n || length(unique(df$braak_num_std[keep])) < 3) {
      return(tibble(gene = g, braak_beta = NA_real_, braak_p = NA_real_,
                    braak_n = sum(keep), adjusted_inverse_braak_beta = NA_real_))
    }

    form <- as.formula(paste("abundance ~", paste(c("braak_num_std", covars), collapse = " + ")))
    fit <- tryCatch(lm(form, data = df[keep, , drop = FALSE]), error = function(e) NULL)

    if (is.null(fit)) {
      return(tibble(gene = g, braak_beta = NA_real_, braak_p = NA_real_,
                    braak_n = sum(keep), adjusted_inverse_braak_beta = NA_real_))
    }

    hit <- broom::tidy(fit) |> filter(term == "braak_num_std")

    if (nrow(hit) != 1) {
      return(tibble(gene = g, braak_beta = NA_real_, braak_p = NA_real_,
                    braak_n = sum(keep), adjusted_inverse_braak_beta = NA_real_))
    }

    tibble(
      gene = g,
      braak_beta = hit$estimate,
      braak_p = hit$p.value,
      braak_n = sum(keep),
      adjusted_inverse_braak_beta = if_else(hit$estimate < 0, abs(hit$estimate), 0)
    )
  }) |>
    mutate(braak_padj = p.adjust(braak_p, method = "BH"))
}

paper_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "top",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

message("Loaded 01_utils.R")


fit_adjusted_cerad_beta <- function(mat, meta_df, sample_col, genes, covars, min_n = 30) {
  mat <- as.matrix(mat)
  require_columns(meta_df, c(sample_col, "cerad_num_std"), "fit_adjusted_cerad_beta metadata")
  meta_df <- meta_df[match(colnames(mat), meta_df[[sample_col]]), , drop = FALSE]
  if (!all(colnames(mat) == meta_df[[sample_col]])) {
    stop(
      "fit_adjusted_cerad_beta could not align matrix columns to metadata sample column: ",
      sample_col,
      ". Available columns: ", available_cols_msg(meta_df),
      call. = FALSE
    )
  }

  genes <- intersect(genes, rownames(mat))
  covars <- intersect(covars, colnames(meta_df))

  covars <- covars[
    purrr::map_lgl(covars, function(v) {
      if (is.numeric(meta_df[[v]]) || is.integer(meta_df[[v]])) {
        covar_is_usable_numeric(meta_df[[v]], min_n = min_n)
      } else {
        covar_is_usable_factor(meta_df[[v]], nrow(meta_df))
      }
    })
  ]

  purrr::map_dfr(genes, function(g) {
    y <- safe_z(mat[g, ])

    df <- tibble::tibble(
      abundance = y,
      cerad_num_std = meta_df$cerad_num_std
    ) |>
      dplyr::bind_cols(meta_df |> dplyr::select(dplyr::any_of(covars))) |>
      dplyr::mutate(dplyr::across(where(is.character), as.factor))

    keep <- stats::complete.cases(df[, c("abundance", "cerad_num_std", covars), drop = FALSE])

    if (sum(keep) < min_n || length(unique(df$cerad_num_std[keep])) < 3) {
      return(tibble::tibble(gene = g, cerad_beta = NA_real_, cerad_p = NA_real_,
                    cerad_n = sum(keep), adjusted_inverse_cerad_beta = NA_real_))
    }

    form <- stats::as.formula(paste("abundance ~", paste(c("cerad_num_std", covars), collapse = " + ")))
    fit <- tryCatch(stats::lm(form, data = df[keep, , drop = FALSE]), error = function(e) NULL)

    if (is.null(fit)) {
      return(tibble::tibble(gene = g, cerad_beta = NA_real_, cerad_p = NA_real_,
                    cerad_n = sum(keep), adjusted_inverse_cerad_beta = NA_real_))
    }

    hit <- broom::tidy(fit) |> dplyr::filter(term == "cerad_num_std")

    if (nrow(hit) != 1) {
      return(tibble::tibble(gene = g, cerad_beta = NA_real_, cerad_p = NA_real_,
                    cerad_n = sum(keep), adjusted_inverse_cerad_beta = NA_real_))
    }

    tibble::tibble(
      gene = g,
      cerad_beta = hit$estimate,
      cerad_p = hit$p.value,
      cerad_n = sum(keep),
      adjusted_inverse_cerad_beta = dplyr::if_else(hit$estimate < 0, abs(hit$estimate), 0)
    )
  }) |>
    dplyr::mutate(cerad_padj = p.adjust(cerad_p, method = "BH"))
}
