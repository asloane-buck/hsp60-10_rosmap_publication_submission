# ============================================================
# Run all supplemental figure scripts
# ============================================================

message("Starting supplemental figure pipeline.")

locate_supplemental_script_dir <- function() {
  this_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
    error = function(e) NA_character_
  )

  if (!is.na(this_file)) {
    return(dirname(this_file))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grepl(file_arg, args)])

  if (length(script_path) == 1 && !is.na(script_path) && nzchar(script_path)) {
    script_path <- normalizePath(script_path, mustWork = TRUE)
    return(dirname(script_path))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
    if (length(script_path) == 1 && !is.na(script_path) && nzchar(script_path)) {
      return(dirname(normalizePath(script_path, mustWork = TRUE)))
    }
  }

  candidate_dirs <- c(getwd(), file.path(getwd(), "R", "supplemental"))
  found <- candidate_dirs[file.exists(file.path(candidate_dirs, "00_supplemental_config.R"))]
  if (length(found) > 0) {
    return(normalizePath(found[[1]], mustWork = TRUE))
  }

  normalizePath(getwd(), mustWork = FALSE)
}

script_dir <- locate_supplemental_script_dir()

message("Script directory: ", script_dir)

# Source core pipeline files
source(file.path(script_dir, "00_supplemental_config.R"))
source(file.path(script_dir, "01_supplemental_load_inputs.R"))
source(file.path(script_dir, "02_supplemental_helper_functions.R"))

if ("hsp60_hsp10_interactor_inventory" %in% names(inputs) &&
    !"hsp60_hsp10_client_inventory" %in% names(inputs)) {
  inputs$hsp60_hsp10_client_inventory <- inputs$hsp60_hsp10_interactor_inventory
}
if ("hsp60_hsp10_client_inventory" %in% names(inputs) &&
    !"hsp60_hsp10_interactor_inventory" %in% names(inputs)) {
  inputs$hsp60_hsp10_interactor_inventory <- inputs$hsp60_hsp10_client_inventory
}

# Make key loaded inputs available globally for older supplemental scripts
if ("protein_meta" %in% names(inputs)) protein_meta <- inputs$protein_meta
if ("prot_meta" %in% names(inputs)) prot_meta <- inputs$prot_meta
if ("prot_meta_aligned" %in% names(inputs)) prot_meta_aligned <- inputs$prot_meta_aligned
if ("prot_mat" %in% names(inputs)) prot_mat <- inputs$prot_mat
if ("rna_meta" %in% names(inputs)) rna_meta <- inputs$rna_meta
if ("rna_mat" %in% names(inputs)) rna_mat <- inputs$rna_mat
if ("prot_mat_raw" %in% names(inputs)) prot_mat_raw <- inputs$prot_mat_raw

source(file.path(script_dir, "10_make_supplementary_figure_1_cohort_detection.R"))
source(file.path(script_dir, "11_make_supplementary_figure_2_matched_individual_sensitivity.R"))

# ============================================================
# Build base Hsp60/10 and mitochondrial background pathology tables
# needed by Supp Fig 4, Supp Fig 3, and Supp Fig 5
# ============================================================

if (!"regional_protein_screen" %in% names(inputs)) {
  stop("regional_protein_screen is missing from inputs; cannot build hsp_null_tbl/background_null_pool.", call. = FALSE)
}

if (!"hsp60_hsp10_interactor_inventory" %in% names(inputs)) {
  stop("hsp60_hsp10_interactor_inventory is missing from inputs.", call. = FALSE)
}

if (!"mito_background_tbl" %in% names(inputs)) {
  stop("mito_background_tbl is missing from inputs.", call. = FALSE)
}

# Cognition layer notes:
# The old pathology_vulnerability_score intentionally is not generated here.
# Supplementary Figure 3 now uses cognition_composite_score for the third
# abundance-matched mitochondrial null metric. Higher values should represent
# stronger coupling between lower Hsp60/10 client abundance and worse cognition.
#
# Important: regional_protein_screen usually contains only AD collapse and
# Braak/pathology columns. Cognition is therefore discovered from a separate
# gene-level cognition/candidate table when it is not already present in the
# regional table, then joined back by gene.

sf99_pick_first_col <- function(df, exact_candidates = character(), regex_candidates = character(), required = TRUE, label = "column") {
  exact_hit <- exact_candidates[exact_candidates %in% colnames(df)][1]
  if (!is.na(exact_hit)) return(exact_hit)

  for (pat in regex_candidates) {
    hits <- grep(pat, colnames(df), value = TRUE, ignore.case = TRUE)
    if (length(hits) > 0) return(hits[1])
  }

  if (isTRUE(required)) {
    stop(
      "Could not identify ", label, ".\n",
      "Tried exact candidates: ", paste(exact_candidates, collapse = ", "), "\n",
      "Tried regex candidates: ", paste(regex_candidates, collapse = ", "), "\n",
      "Available columns: ", paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }

  NA_character_
}

sf99_gene_col <- function(df, required = TRUE, label = "gene column") {
  sf99_pick_first_col(
    df,
    exact_candidates = c(
      "gene", "gene_symbol", "Gene", "Gene Symbol", "symbol", "SYMBOL",
      "hgnc_symbol", "HGNC", "display_gene", "curated_gene", "gene_name"
    ),
    regex_candidates = c("^gene$", "gene.*symbol", "symbol"),
    required = required,
    label = label
  )
}

sf99_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (sum(is.finite(x)) < 3 || stats::sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

sf99_orient_cognition_component <- function(x, col_name) {
  x <- suppressWarnings(as.numeric(x))
  cn <- tolower(col_name)

  # If the column is already explicitly oriented as inverse/worse/decline, keep it.
  if (grepl("inverse|worse|impair|declin|dement", cn) && !grepl("mmse", cn)) {
    return(x)
  }

  # MMSE is higher = better. Raw MMSE beta/rho columns are flipped so higher means worse cognition.
  # Columns already named inverse/lower/worse/decline are kept as-is.
  if (grepl("mmse", cn) && !grepl("inverse|lower|worse|impair|declin", cn)) {
    return(-x)
  }

  x
}

sf99_cognition_column_candidates <- list(
  composite_exact = c(
    "cognition_composite_score", "cognitive_composite_score",
    "cognition_association_score", "cognitive_association_score",
    "client_cognition_score", "cognition_score", "cognitive_score",
    "cognition_layer_score", "cognitive_layer_score",
    "cognition_priority_score", "cognitive_priority_score"
  ),
  composite_regex = c(
    "cognit.*composite", "cognit.*association.*score", "cognit.*priority.*score",
    "cognitive.*composite", "cognitive.*association.*score"
  ),
  final_dx_exact = c(
    "final_diagnosis_cognition_score", "final_diagnosis_association",
    "last_visit_diagnosis_association", "last_visit_diagnosis_beta",
    "final_diagnosis_beta", "cogdx_last_visit_beta", "last_visit_cogdx_beta",
    "final_dx_beta", "last_visit_dx_beta", "diagnosis_last_visit_beta",
    "final_diagnosis_magnitude", "last_visit_diagnosis_magnitude"
  ),
  final_dx_regex = c(
    "final.*diagnos.*(score|association|beta|rho|magnitude)",
    "last.*visit.*diagnos.*(score|association|beta|rho|magnitude)",
    "diagnos.*last.*visit.*(score|association|beta|rho|magnitude)",
    "cogdx.*(last|final).*(score|association|beta|rho|magnitude)",
    "(last|final).*cogdx.*(score|association|beta|rho|magnitude)"
  ),
  mmse_exact = c(
    "last_visit_mmse_association", "last_visit_mmse_beta",
    "last_mmse_beta", "mmse_last_visit_beta", "final_mmse_beta",
    "inverse_mmse_magnitude", "last_visit_inverse_mmse_magnitude",
    "mmse_association", "mmse_beta", "mmse_magnitude"
  ),
  mmse_regex = c(
    "last.*mmse.*(score|association|beta|rho|magnitude)",
    "mmse.*last.*(score|association|beta|rho|magnitude)",
    "inverse.*mmse.*(score|association|beta|rho|magnitude)",
    "mmse.*(score|association|beta|rho|magnitude)"
  ),
  death_dx_exact = c(
    "diagnosis_at_death_association", "diagnosis_at_death_beta",
    "death_diagnosis_association", "death_diagnosis_beta",
    "cogdx_death_beta", "final_cogdx_death_beta",
    "diagnosis_at_death_magnitude", "death_diagnosis_magnitude"
  ),
  death_dx_regex = c(
    "diagnos.*death.*(score|association|beta|rho|magnitude)",
    "death.*diagnos.*(score|association|beta|rho|magnitude)",
    "cogdx.*death.*(score|association|beta|rho|magnitude)",
    "death.*cogdx.*(score|association|beta|rho|magnitude)"
  )
)

sf99_make_cognition_lookup <- function(df, source_label = "unknown cognition source") {
  df <- tibble::as_tibble(df)
  gene_col <- sf99_gene_col(df, required = FALSE, label = paste(source_label, "gene column"))
  if (is.na(gene_col)) return(NULL)

  composite_col <- sf99_pick_first_col(
    df,
    exact_candidates = sf99_cognition_column_candidates$composite_exact,
    regex_candidates = sf99_cognition_column_candidates$composite_regex,
    required = FALSE,
    label = paste(source_label, "precomputed cognition composite")
  )

  final_dx_col <- sf99_pick_first_col(
    df,
    exact_candidates = sf99_cognition_column_candidates$final_dx_exact,
    regex_candidates = sf99_cognition_column_candidates$final_dx_regex,
    required = FALSE,
    label = paste(source_label, "final diagnosis / last-visit diagnosis component")
  )

  last_mmse_col <- sf99_pick_first_col(
    df,
    exact_candidates = sf99_cognition_column_candidates$mmse_exact,
    regex_candidates = sf99_cognition_column_candidates$mmse_regex,
    required = FALSE,
    label = paste(source_label, "last-visit MMSE component")
  )

  death_dx_col <- sf99_pick_first_col(
    df,
    exact_candidates = sf99_cognition_column_candidates$death_dx_exact,
    regex_candidates = sf99_cognition_column_candidates$death_dx_regex,
    required = FALSE,
    label = paste(source_label, "diagnosis-at-death component")
  )

  if (is.na(composite_col) && all(is.na(c(final_dx_col, last_mmse_col, death_dx_col)))) {
    return(NULL)
  }

  if (!is.na(composite_col)) {
    score_vec <- suppressWarnings(as.numeric(df[[composite_col]]))
    n_component_vec <- rep(NA_integer_, nrow(df))
    component_label <- composite_col
  } else {
    components <- list()
    finite_components <- list()

    if (!is.na(final_dx_col)) {
      v <- sf99_orient_cognition_component(df[[final_dx_col]], final_dx_col)
      components[[final_dx_col]] <- sf99_z(v)
      finite_components[[final_dx_col]] <- is.finite(v)
    }
    if (!is.na(last_mmse_col)) {
      v <- sf99_orient_cognition_component(df[[last_mmse_col]], last_mmse_col)
      components[[last_mmse_col]] <- sf99_z(v)
      finite_components[[last_mmse_col]] <- is.finite(v)
    }
    if (!is.na(death_dx_col)) {
      v <- sf99_orient_cognition_component(df[[death_dx_col]], death_dx_col)
      components[[death_dx_col]] <- sf99_z(v)
      finite_components[[death_dx_col]] <- is.finite(v)
    }

    score_vec <- rowSums(do.call(cbind, components), na.rm = TRUE)
    n_component_vec <- rowSums(do.call(cbind, finite_components), na.rm = TRUE)
    score_vec[n_component_vec == 0] <- NA_real_
    component_label <- paste(names(components), collapse = " + ")
  }

  out <- tibble::tibble(
    gene = clean_gene(df[[gene_col]]),
    cognition_composite_score = score_vec,
    cognition_metric_source = source_label,
    cognition_component_source = component_label,
    n_cognition_components = n_component_vec
  ) |>
    dplyr::filter(!is.na(.data$gene), nzchar(.data$gene), is.finite(.data$cognition_composite_score)) |>
    dplyr::group_by(.data$gene) |>
    dplyr::slice_max(order_by = abs(.data$cognition_composite_score), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  if (nrow(out) == 0) return(NULL)
  out
}

sf99_scan_project_for_cognition_tables <- function() {
  if (!exists("project_dir") || is.na(project_dir) || !dir.exists(project_dir)) {
    return(list())
  }

  candidate_files <- list.files(
    project_dir,
    pattern = "\\.(csv|tsv|txt|rds)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  # Keep this targeted so the pipeline does not inspect every large file.
  candidate_files <- candidate_files[
    grepl(
      "cogn|cognitive|cogdx|mmse|diagnosis|death|fig6|candidate|ranked|priorit|client",
      candidate_files,
      ignore.case = TRUE
    )
  ]

  # Prefer explicit cognition files, then candidate/prioritization tables.
  ord <- order(
    !grepl("cogn|cognitive|cogdx|mmse|diagnosis|death", candidate_files, ignore.case = TRUE),
    !grepl("fig6|candidate|ranked|priorit|client", candidate_files, ignore.case = TRUE),
    nchar(candidate_files)
  )
  candidate_files <- unique(candidate_files[ord])

  out <- list()
  for (path_i in head(candidate_files, 250)) {
    tbl_i <- tryCatch({
      ext <- tolower(tools::file_ext(path_i))
      if (ext == "rds") {
        obj <- readRDS(path_i)
        if (is.data.frame(obj) || tibble::is_tibble(obj)) tibble::as_tibble(obj) else NULL
      } else if (ext %in% c("tsv", "txt")) {
        readr::read_tsv(path_i, show_col_types = FALSE)
      } else {
        readr::read_csv(path_i, show_col_types = FALSE)
      }
    }, error = function(e) NULL)

    if (is.null(tbl_i) || !is.data.frame(tbl_i)) next

    lookup_i <- sf99_make_cognition_lookup(tbl_i, source_label = path_i)
    if (!is.null(lookup_i) && nrow(lookup_i) > 0) {
      message("Detected gene-level cognition table: ", path_i)
      message("  rows with usable cognition scores: ", nrow(lookup_i))
      out[[path_i]] <- lookup_i
    }
  }

  out
}

sf99_find_cognition_lookup <- function(regional_tbl, hsp_genes = character(), background_genes = character()) {
  lookup_candidates <- list()

  # 1. Try the regional table itself.
  regional_lookup <- sf99_make_cognition_lookup(regional_tbl, "regional_protein_screen")
  if (!is.null(regional_lookup)) {
    lookup_candidates[["regional_protein_screen"]] <- regional_lookup
  }

  # 2. Try already-loaded inputs.
  input_names <- names(inputs)[vapply(inputs, function(x) is.data.frame(x) || tibble::is_tibble(x), logical(1))]
  for (nm in input_names) {
    lookup_i <- sf99_make_cognition_lookup(inputs[[nm]], paste0("inputs$", nm))
    if (!is.null(lookup_i)) {
      lookup_candidates[[paste0("inputs$", nm)]] <- lookup_i
    }
  }

  # 3. Try data frames in the global environment.
  global_names <- ls(envir = .GlobalEnv)
  for (nm in global_names) {
    if (nm %in% c("inputs", "regional_dlpfc", "regional_base")) next
    obj <- tryCatch(get(nm, envir = .GlobalEnv), error = function(e) NULL)
    if (!is.data.frame(obj) && !tibble::is_tibble(obj)) next
    lookup_i <- sf99_make_cognition_lookup(obj, paste0("global$", nm))
    if (!is.null(lookup_i)) {
      lookup_candidates[[paste0("global$", nm)]] <- lookup_i
    }
  }

  # 4. Scan likely project-level CSV/RDS outputs, because main Figure 6/cognition
  # candidate files may not be part of supplemental candidate_csv_paths.
  project_lookups <- sf99_scan_project_for_cognition_tables()
  lookup_candidates <- c(lookup_candidates, project_lookups)

  if (length(lookup_candidates) == 0) {
    return(NULL)
  }

  target_genes <- unique(stats::na.omit(c(clean_gene(regional_tbl$gene_symbol), hsp_genes, background_genes)))

  candidate_summary <- purrr::imap_dfr(lookup_candidates, function(tbl, nm) {
    tibble::tibble(
      source = nm,
      n_rows = nrow(tbl),
      n_overlap_regional_or_background = sum(tbl$gene %in% target_genes),
      n_overlap_hsp = sum(tbl$gene %in% hsp_genes),
      n_overlap_background = sum(tbl$gene %in% background_genes)
    )
  }) |>
    dplyr::arrange(
      dplyr::desc(.data$n_overlap_hsp),
      dplyr::desc(.data$n_overlap_background),
      dplyr::desc(.data$n_overlap_regional_or_background),
      dplyr::desc(.data$n_rows)
    )

  readr::write_csv(candidate_summary, file.path(audits_dir, "COGNITION_TABLE_DISCOVERY_AUDIT.csv"))
  message("Cognition table discovery audit written to: ", file.path(audits_dir, "COGNITION_TABLE_DISCOVERY_AUDIT.csv"))
  print(candidate_summary)

  best_source <- candidate_summary$source[1]
  if (is.na(best_source) || candidate_summary$n_overlap_regional_or_background[1] == 0) {
    stop(
      "Cognition tables were detected, but none overlapped regional_protein_screen / Hsp60/10 / mitochondrial-background genes. ",
      "Check COGNITION_TABLE_DISCOVERY_AUDIT.csv.",
      call. = FALSE
    )
  }

  best_lookup <- lookup_candidates[[best_source]]
  message("Selected cognition table for supplemental null testing: ", best_source)
  message("  regional/background overlap: ", candidate_summary$n_overlap_regional_or_background[1])
  message("  Hsp60/10 overlap: ", candidate_summary$n_overlap_hsp[1])
  message("  background overlap: ", candidate_summary$n_overlap_background[1])

  best_lookup
}

regional_dlpfc <- inputs$regional_protein_screen |>
  tibble::as_tibble() |>
  dplyr::filter(.data$region_short == "DLPFC")

hsp_genes <- inputs$hsp60_hsp10_interactor_inventory |>
  tibble::as_tibble() |>
  dplyr::transmute(gene = clean_gene(.data$gene)) |>
  dplyr::filter(!is.na(.data$gene)) |>
  dplyr::pull(.data$gene) |>
  unique()

background_genes <- inputs$mito_background_tbl |>
  tibble::as_tibble() |>
  dplyr::transmute(gene = clean_gene(.data$gene)) |>
  dplyr::filter(!is.na(.data$gene)) |>
  dplyr::pull(.data$gene) |>
  unique()

cognition_lookup_tbl <- sf99_find_cognition_lookup(
  regional_tbl = regional_dlpfc,
  hsp_genes = hsp_genes,
  background_genes = background_genes
)

if (is.null(cognition_lookup_tbl) || nrow(cognition_lookup_tbl) == 0) {
  stop(
    "No usable gene-level cognition metric was found. The supplemental null analysis needs a table with a gene column plus either a cognition composite score, ",
    "or gene-level components for final diagnosis / last-visit diagnosis, last-visit MMSE, and/or diagnosis at death. ",
    "regional_protein_screen columns are: ", paste(colnames(regional_dlpfc), collapse = ", "),
    call. = FALSE
  )
}

readr::write_csv(cognition_lookup_tbl, file.path(audits_dir, "COGNITION_LOOKUP_SELECTED_FOR_SUPPLEMENTAL_NULL.csv"))

regional_dlpfc <- regional_dlpfc |>
  dplyr::mutate(gene = clean_gene(.data$gene_symbol)) |>
  dplyr::left_join(cognition_lookup_tbl, by = "gene")

n_missing_cognition <- sum(!is.finite(regional_dlpfc$cognition_composite_score))
if (n_missing_cognition > 0) {
  warning(
    "regional_dlpfc: ", n_missing_cognition,
    " proteins lack cognition_composite_score after joining the selected cognition table. ",
    "They will be excluded only from cognition-dependent null testing."
  )
}

regional_base <- regional_dlpfc |>
  dplyr::transmute(
    gene = clean_gene(.data$gene_symbol),
    protein_collapse_magnitude = suppressWarnings(as.numeric(.data$adjusted_collapse_magnitude)),
    inverse_braak_magnitude = suppressWarnings(as.numeric(.data$adjusted_inverse_braak_magnitude)),
    cognition_composite_score = suppressWarnings(as.numeric(.data$cognition_composite_score)),
    cognition_metric_source = .data$cognition_metric_source,
    cognition_component_source = .data$cognition_component_source,
    n_cognition_components = .data$n_cognition_components,
    matching_abundance = suppressWarnings(as.numeric(abs(.data$beta_ad))),
    abundance_bin = dplyr::ntile(.data$matching_abundance, 5),
    agora_nominated_target = 0
  ) |>
  dplyr::filter(
    !is.na(.data$gene),
    is.finite(.data$protein_collapse_magnitude),
    is.finite(.data$inverse_braak_magnitude),
    is.finite(.data$cognition_composite_score),
    is.finite(.data$matching_abundance),
    !is.na(.data$abundance_bin)
  ) |>
  dplyr::distinct(.data$gene, .keep_all = TRUE)

hsp_genes <- inputs$hsp60_hsp10_interactor_inventory |>
  tibble::as_tibble() |>
  dplyr::transmute(gene = clean_gene(.data$gene)) |>
  dplyr::filter(!is.na(.data$gene)) |>
  dplyr::pull(.data$gene) |>
  unique()

background_genes <- inputs$mito_background_tbl |>
  tibble::as_tibble() |>
  dplyr::transmute(gene = clean_gene(.data$gene)) |>
  dplyr::filter(!is.na(.data$gene)) |>
  dplyr::pull(.data$gene) |>
  unique()

hsp_null_tbl <- regional_base |>
  dplyr::filter(.data$gene %in% hsp_genes) |>
  dplyr::mutate(
    is_hsp60_10_client = TRUE,
    collapse_rank = dplyr::percent_rank(.data$protein_collapse_magnitude),
    inverse_braak_rank = dplyr::percent_rank(.data$inverse_braak_magnitude),
    cognition_rank = dplyr::percent_rank(.data$cognition_composite_score)
  )

background_null_pool <- regional_base |>
  dplyr::filter(.data$gene %in% background_genes, !.data$gene %in% hsp_genes) |>
  dplyr::mutate(
    is_hsp60_10_client = FALSE,
    collapse_rank = dplyr::percent_rank(.data$protein_collapse_magnitude),
    inverse_braak_rank = dplyr::percent_rank(.data$inverse_braak_magnitude),
    cognition_rank = dplyr::percent_rank(.data$cognition_composite_score)
  )

if (nrow(hsp_null_tbl) == 0) {
  stop("Built hsp_null_tbl has zero rows. Check gene overlap between regional_protein_screen and Hsp60/10 inventory.", call. = FALSE)
}

if (nrow(background_null_pool) == 0) {
  stop("Built background_null_pool has zero rows. Check gene overlap between regional_protein_screen and mito_background_tbl.", call. = FALSE)
}

inputs$hsp_null_tbl <- hsp_null_tbl
inputs$background_null_pool <- background_null_pool

readr::write_csv(hsp_null_tbl, file.path(audits_dir, "BASE_hsp_null_tbl_built_from_regional_DLPFC.csv"))
readr::write_csv(background_null_pool, file.path(audits_dir, "BASE_background_null_pool_built_from_regional_DLPFC.csv"))

message("Built base hsp_null_tbl/background_null_pool before Supp Fig 4:")
message("  hsp_null_tbl rows: ", nrow(hsp_null_tbl))
message("  background_null_pool rows: ", nrow(background_null_pool))

# Run this early because it computes the current adjusted Braak model
source(file.path(script_dir, "13_make_supplementary_figure_4_pathology_model_robustness.R"))



# Then standardize all downstream gene-level Braak values
# ============================================================
# Enforce current manuscript Braak metric for downstream supplements
# ============================================================

if (!exists("supfig4_outputs")) {
  stop("supfig4_outputs not found. Supp Fig 4 must run before downstream Braak-dependent supplemental figures.", call. = FALSE)
}

if (is.null(supfig4_outputs$plot_tbl) || !"gene" %in% colnames(supfig4_outputs$plot_tbl)) {
  stop("supfig4_outputs$plot_tbl is missing or lacks a gene column.", call. = FALSE)
}

current_braak_tbl <- supfig4_outputs$plot_tbl |>
  dplyr::select(
    gene,
    inverse_braak_beta_adjusted,
    braak_beta_adjusted,
    braak_p_adjusted,
    n_adjusted
  ) |>
  dplyr::mutate(gene = clean_gene(.data$gene)) |>
  dplyr::filter(
    !is.na(.data$gene),
    is.finite(.data$inverse_braak_beta_adjusted)
  ) |>
  dplyr::distinct(.data$gene, .keep_all = TRUE)

if (nrow(current_braak_tbl) == 0) {
  stop("Current adjusted Braak table is empty. Check Supp Fig 4 model fitting.", call. = FALSE)
}

apply_current_braak_metric <- function(tbl, label) {
  if (is.null(tbl)) {
    stop(label, " is NULL.", call. = FALSE)
  }

  tbl <- tibble::as_tibble(tbl)

  if (!"gene" %in% colnames(tbl)) {
    stop(
      label, " does not contain a gene column.\n",
      "Available columns: ", paste(colnames(tbl), collapse = ", "),
      call. = FALSE
    )
  }

  # Remove any previously patched adjusted-Braak columns before joining.
  # This prevents .x/.y suffixes if the script is rerun in the same session.
  tbl_clean <- tbl |>
    dplyr::select(
      -dplyr::any_of(c(
        "inverse_braak_beta_adjusted",
        "braak_beta_adjusted",
        "braak_p_adjusted",
        "n_adjusted",
        "braak_beta_cerad_adjusted",
        "braak_se_cerad_adjusted",
        "braak_p_cerad_adjusted",
        "n_cerad_adjusted",
        "inverse_braak_beta_cerad_adjusted",
        "inverse_braak_beta_adjusted.x",
        "inverse_braak_beta_adjusted.y",
        "braak_beta_adjusted.x",
        "braak_beta_adjusted.y",
        "braak_p_adjusted.x",
        "braak_p_adjusted.y",
        "n_adjusted.x",
        "n_adjusted.y",
        "braak_metric_source",
        "braak_metric_model"
      ))
    )

  tbl2 <- tbl_clean |>
    dplyr::mutate(gene = clean_gene(.data$gene)) |>
    dplyr::left_join(current_braak_tbl, by = "gene")

  if (!"inverse_braak_beta_adjusted" %in% colnames(tbl2)) {
    stop(
      label, " lacks inverse_braak_beta_adjusted after joining current_braak_tbl.\n",
      "Columns after join: ", paste(colnames(tbl2), collapse = ", "),
      call. = FALSE
    )
  }

  n_missing <- sum(!is.finite(tbl2$inverse_braak_beta_adjusted))

  if (n_missing > 0) {
    warning(
      label, ": ", n_missing,
      " genes lack current covariate-adjusted Braak values and will have NA for the current metric."
    )
  }

  tbl2 |>
    dplyr::mutate(
      inverse_braak_magnitude = .data$inverse_braak_beta_adjusted,
      braak_metric_source = "covariate_adjusted_no_CERAD",
      braak_metric_model = "protein_z ~ braak_z + age_z + sex + pmi_z"
    ) |>
    dplyr::filter(is.finite(.data$inverse_braak_magnitude))
}

hsp_tbl_for_patch <- if ("hsp_null_tbl" %in% names(inputs)) {
  inputs$hsp_null_tbl
} else if (exists("hsp_null_tbl", envir = .GlobalEnv)) {
  get("hsp_null_tbl", envir = .GlobalEnv)
} else {
  stop("Could not find hsp_null_tbl in inputs or global environment.", call. = FALSE)
}

background_tbl_for_patch <- if ("background_null_pool" %in% names(inputs)) {
  inputs$background_null_pool
} else if (exists("background_null_pool", envir = .GlobalEnv)) {
  get("background_null_pool", envir = .GlobalEnv)
} else {
  stop("Could not find background_null_pool in inputs or global environment.", call. = FALSE)
}

hsp_null_tbl <- apply_current_braak_metric(hsp_tbl_for_patch, "hsp_null_tbl")
background_null_pool <- apply_current_braak_metric(background_tbl_for_patch, "background_null_pool")

inputs$hsp_null_tbl <- hsp_null_tbl
inputs$background_null_pool <- background_null_pool

readr::write_csv(
  dplyr::bind_rows(
    hsp_null_tbl |> dplyr::mutate(source_table = "hsp_null_tbl"),
    background_null_pool |> dplyr::mutate(source_table = "background_null_pool")
  ) |>
    dplyr::select(
      source_table,
      gene,
      braak_metric_source,
      braak_metric_model,
      inverse_braak_magnitude,
      inverse_braak_beta_adjusted,
      braak_beta_adjusted,
      braak_p_adjusted,
      n_adjusted,
      dplyr::everything()
    ),
  file.path(audits_dir, "CURRENT_BRAAK_METRIC_ENFORCED_FOR_SUPPLEMENTAL_FIGURES.csv")
)

message("Current Braak metric enforced for downstream supplemental figures:")
message("  inverse_braak_magnitude := inverse_braak_beta_adjusted")
message("  model: protein_z ~ braak_z + age_z + sex + pmi_z")
message("  CERAD is not included.")

source(file.path(script_dir, "12_make_supplementary_figure_3_mitochondrial_specificity.R"))
source(file.path(script_dir, "15_make_supplementary_figure_6_regional_proteomics_validation.R"))
source(file.path(script_dir, "14_make_supplementary_figure_5_msbb_cross_cohort_validation.R"))

copy_manuscript_ready_supplemental_pdfs <- function() {
  ensure_dir(manuscript_ready_pdf_dir)

  stale_pdfs <- list.files(
    manuscript_ready_pdf_dir,
    pattern = "\\.pdf$",
    full.names = TRUE
  )
  if (length(stale_pdfs) > 0) {
    unlink(stale_pdfs)
  }

  final_pdf_names <- c(
    "Supplementary_Figure_1_cohort_detection.pdf",
    "Supplementary_Figure_2_matched_individual_sensitivity.pdf",
    "Supplementary_Figure_3_mitochondrial_specificity.pdf",
    "Supplementary_Figure_4_pathology_model_robustness.pdf",
    "Supplementary_Figure_5_msbb_cross_cohort_validation.pdf",
    "Supplementary_Figure_6_regional_proteomics_validation.pdf"
  )

  src_paths <- file.path(figures_dir, final_pdf_names)
  missing_paths <- src_paths[!file.exists(src_paths)]
  if (length(missing_paths) > 0) {
    stop(
      "Cannot populate manuscript-ready supplemental PDF folder; missing regenerated PDFs: ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
  }

  dest_paths <- file.path(manuscript_ready_pdf_dir, final_pdf_names)
  ok <- file.copy(src_paths, dest_paths, overwrite = TRUE)
  if (!all(ok)) {
    stop("Failed to copy one or more manuscript-ready supplemental PDFs.", call. = FALSE)
  }

  manifest <- tibble::tibble(
    supplemental_figure = paste0("Supplementary Figure ", seq_along(final_pdf_names)),
    pdf_path = dest_paths,
    size_bytes = file.info(dest_paths)$size
  )
  readr::write_csv(manifest, file.path(manuscript_ready_pdf_dir, "manifest.csv"))

  message("Manuscript-ready supplemental PDFs written to: ", manuscript_ready_pdf_dir)
  invisible(manifest)
}

manuscript_ready_pdf_manifest <- copy_manuscript_ready_supplemental_pdfs()
manuscript_ready_pdf_audit <- audit_manuscript_ready_supplemental_pdfs()
message(
  "Manuscript-ready supplemental PDF dimension audit written to: ",
  file.path(audits_dir, "manuscript_ready_supplemental_pdf_dimension_audit.csv")
)

message("Supplemental figure pipeline complete.")
message("Figures written to: ", figures_dir)
message("Panels written to: ", panels_dir)
message("Audits written to: ", audits_dir)
message("Tables written to: ", tables_dir)
message("Manuscript-ready supplemental PDFs written to: ", manuscript_ready_pdf_dir)




############################################################
## SUP FIG 5 MSBB DIAGNOSTIC AUDIT
############################################################

sf5_diag <- supfig5_outputs$msbb_tbl |>
  dplyr::mutate(
    group = as.character(group),
    is_hsp = group == "Hsp60/10 clients"
  )

cat("\nMSBB group sizes:\n")
print(sf5_diag |> dplyr::count(group))

cat("\nMSBB inverse beta summary:\n")
print(
  sf5_diag |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      n = sum(is.finite(msbb_inverse_braak_beta)),
      median = median(msbb_inverse_braak_beta, na.rm = TRUE),
      mean = mean(msbb_inverse_braak_beta, na.rm = TRUE),
      sd = sd(msbb_inverse_braak_beta, na.rm = TRUE),
      q25 = quantile(msbb_inverse_braak_beta, 0.25, na.rm = TRUE),
      q75 = quantile(msbb_inverse_braak_beta, 0.75, na.rm = TRUE),
      n_positive = sum(msbb_inverse_braak_beta > 0, na.rm = TRUE),
      pct_positive = mean(msbb_inverse_braak_beta > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    )
)

cat("\nWilcoxon tests:\n")

print(
  wilcox.test(
    msbb_inverse_braak_beta ~ group,
    data = sf5_diag |> dplyr::filter(is.finite(msbb_inverse_braak_beta)),
    exact = FALSE
  )
)

cat("\nOne-sided test: Hsp60/10 greater inverse beta than background\n")
print(
  wilcox.test(
    msbb_inverse_braak_beta ~ group,
    data = sf5_diag |> dplyr::filter(is.finite(msbb_inverse_braak_beta)),
    alternative = "less",
    exact = FALSE
  )
)

cat("\nRaw rho-based comparison, if available:\n")
if ("msbb_inverse_braak_rho" %in% colnames(sf5_diag)) {
  print(
    sf5_diag |>
      dplyr::group_by(group) |>
      dplyr::summarise(
        n = sum(is.finite(msbb_inverse_braak_rho)),
        median = median(msbb_inverse_braak_rho, na.rm = TRUE),
        mean = mean(msbb_inverse_braak_rho, na.rm = TRUE),
        pct_positive = mean(msbb_inverse_braak_rho > 0, na.rm = TRUE) * 100,
        .groups = "drop"
      )
  )
  
  print(
    wilcox.test(
      msbb_inverse_braak_rho ~ group,
      data = sf5_diag |> dplyr::filter(is.finite(msbb_inverse_braak_rho)),
      exact = FALSE
    )
  )
}

cat("\nCross-cohort Hsp60/10 overlap:\n")
print(
  supfig5_outputs$cross_tbl |>
    dplyr::filter(group == "Hsp60/10 clients") |>
    dplyr::summarise(
      n = dplyr::n(),
      n_rosmap_inverse = sum(rosmap_inverse_braak > 0, na.rm = TRUE),
      n_msbb_inverse_beta = sum(msbb_inverse_braak_beta > 0, na.rm = TRUE),
      n_inverse_both = sum(rosmap_inverse_braak > 0 & msbb_inverse_braak_beta > 0, na.rm = TRUE),
      pct_inverse_both = 100 * mean(rosmap_inverse_braak > 0 & msbb_inverse_braak_beta > 0, na.rm = TRUE),
      spearman_rho = cor(
        rosmap_inverse_braak,
        msbb_inverse_braak_beta,
        method = "spearman",
        use = "complete.obs"
      ),
      .groups = "drop"
    )
)

sf5_diag2 <- sf5_diag |>
  dplyr::mutate(
    group = factor(
      group,
      levels = c(
        "Hsp60/10 clients",
        "Non-client mitochondrial proteins"
      )
    )
  )

wilcox.test(
  msbb_inverse_braak_beta ~ group,
  data = sf5_diag2 |> dplyr::filter(is.finite(msbb_inverse_braak_beta)),
  alternative = "greater",
  exact = FALSE
)

stats::cor.test(
  supfig5_outputs$cross_tbl$rosmap_inverse_braak[
    supfig5_outputs$cross_tbl$group == "Hsp60/10 clients"
  ],
  supfig5_outputs$cross_tbl$msbb_inverse_braak_beta[
    supfig5_outputs$cross_tbl$group == "Hsp60/10 clients"
  ],
  method = "spearman",
  exact = FALSE
)
