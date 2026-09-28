options(stringsAsFactors = FALSE)

# =============================================================================
# Reviewer 2.2c: formal regional interaction analysis
#
# PRIMARY TEST:
#
# abundance ~ region * predictor
#             + region * cohort
#             + region * age
#             + region * sex
#             + region * PMI
#             + random intercept for participant
#
# predictor:
#   1) AD vs Control
#   2) continuous Braak stage I-VI
#
# Braak sensitivity:
#   categorical I-II / III-IV / V-VI
#
# Technical assay replicates are collapsed FIRST within participant x region.
# No separately estimated regional coefficients are subtracted to generate the
# primary interaction test.
# =============================================================================


# =============================================================================
# Packages
# =============================================================================

if (!requireNamespace("nlme", quietly = TRUE)) {
  stop(
    "Package 'nlme' is required but not installed.\n",
    "Run: install.packages('nlme')"
  )
}


# =============================================================================
# Paths
# =============================================================================

BASE <- Sys.getenv(
  "ROSMAP_REGIONAL_BASE",
  unset = ""
)

if (!nzchar(BASE)) {
  stop(
    "Environment variable ROSMAP_REGIONAL_BASE is required. ",
    "Set it to the local Brain Region Specificity data directory."
  )
}

BASE <- normalizePath(
  BASE,
  mustWork = TRUE
)

OUTDIR <- "outputs/reviewer_revisions/regional_formal_interaction"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

DLPFC_MATRIX <- file.path(
  BASE,
  "Raw Data",
  "n1086_residual_log2_batch.csv"
)

STG_MATRIX <- file.path(
  BASE,
  "Raw Data",
  "n278_residual_log2_batch.TCX.csv"
)

DLPFC_MAP <- file.path(
  BASE,
  "REGIONAL_FIRST_PASS_outputs",
  "AMP_DLPFC_sample_map_WITH_PHENO.csv"
)

STG_MAP <- file.path(
  BASE,
  "REGIONAL_FIRST_PASS_outputs",
  "AMP_STG_sample_map_WITH_PHENO.csv"
)

SOURCE_META <- file.path(
  BASE,
  "Metadata",
  "AMP-AD_DiverseCohorts_individual_metadata.csv"
)

CORRECTED_META <- file.path(
  "outputs/reviewer_revisions/regional_braak_fix",
  "regional_participant_metadata_corrected_braak.csv"
)

OLD_RESULTS <- file.path(
  BASE,
  "REGIONAL_COVARIATE_ADJUSTED_outputs",
  "AMP_covariate_adjusted_regional_screen_labeled_Hsp60_mito.csv"
)

required_files <- c(
  DLPFC_MATRIX,
  STG_MATRIX,
  DLPFC_MAP,
  STG_MAP,
  SOURCE_META,
  CORRECTED_META,
  OLD_RESULTS
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}


# =============================================================================
# Helpers
# =============================================================================

clean_id <- function(x) {

  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "nan", "<NA>")] <- NA

  sub("\\.0$", "", x)
}


clean_missing <- function(x) {

  x <- trimws(as.character(x))

  x[x %in% c(
    "",
    "NA",
    "NaN",
    "nan",
    "<NA>",
    "missing or unknown"
  )] <- NA

  x
}


get_accession <- function(x) {

  x <- trimws(as.character(x))

  ifelse(
    grepl("\\|", x),
    sub("^.*\\|", "", x),
    x
  )
}


normalize_logical <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(trimws(as.character(x)))

  out <- rep(NA, length(y))

  out[y %in% c("TRUE", "T", "1")] <- TRUE
  out[y %in% c("FALSE", "F", "0")] <- FALSE

  out
}


zscore <- function(x) {

  x <- suppressWarnings(as.numeric(x))

  s <- sd(x, na.rm = TRUE)

  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }

  as.numeric(scale(x))
}


# =============================================================================
# Convert matrix to accession-keyed numeric matrix
# =============================================================================

prepare_expression_matrix <- function(file, region_name) {

  x <- read.csv(
    file,
    check.names = FALSE
  )

  if (ncol(x) < 2) {
    stop(region_name, ": expression matrix has <2 columns.")
  }

  raw_id <- as.character(x[[1]])
  accession <- get_accession(raw_id)

  if (any(is.na(accession) | accession == "")) {
    stop(region_name, ": missing protein accession after parsing.")
  }

  dup <- unique(accession[duplicated(accession)])

  if (length(dup) > 0) {
    stop(
      region_name,
      ": duplicate accessions in expression matrix. First examples: ",
      paste(head(dup, 10), collapse = ", ")
    )
  }

  mat <- as.matrix(
    x[, -1, drop = FALSE]
  )

  storage.mode(mat) <- "numeric"

  rownames(mat) <- accession

  mapping <- data.frame(
    accession = accession,
    matrix_row_id = raw_id,
    region = region_name,
    stringsAsFactors = FALSE
  )

  list(
    matrix = mat,
    mapping = mapping
  )
}


# =============================================================================
# Collapse technical assays within participant x region
# =============================================================================

collapse_technical_replicates <- function(
  mat,
  smap,
  paired_ids,
  region_name
) {

  smap$individual_id <- clean_id(
    smap$individual_id
  )

  if (!"matrix_col" %in% names(smap)) {
    stop(region_name, ": sample map lacks matrix_col.")
  }

  if (!"trait_dataset" %in% names(smap)) {
    stop(region_name, ": sample map lacks trait_dataset.")
  }

  keep <- (
    !is.na(smap$individual_id) &
    smap$individual_id %in% paired_ids &
    smap$matrix_col %in% colnames(mat)
  )

  sm <- smap[
    keep,
    ,
    drop = FALSE
  ]

  if (nrow(sm) == 0) {
    stop(region_name, ": no paired assay columns were mapped.")
  }

  if (anyDuplicated(sm$matrix_col)) {
    dup <- unique(
      sm$matrix_col[duplicated(sm$matrix_col)]
    )

    stop(
      region_name,
      ": duplicated matrix_col values in map: ",
      paste(head(dup, 10), collapse = ", ")
    )
  }

  by_person <- split(
    sm$matrix_col,
    sm$individual_id
  )

  collapsed <- vapply(
    by_person,
    function(cols) {

      ans <- rowMeans(
        mat[, cols, drop = FALSE],
        na.rm = TRUE
      )

      ans[is.nan(ans)] <- NA_real_

      ans
    },
    FUN.VALUE = numeric(nrow(mat))
  )

  rownames(collapsed) <- rownames(mat)

  manifest <- data.frame(
    individual_id = names(by_person),
    region = region_name,
    n_technical_assays = lengths(by_person),
    trait_dataset = vapply(
      names(by_person),
      function(id) {

        vals <- sort(
          unique(
            as.character(
              sm$trait_dataset[
                sm$individual_id == id
              ]
            )
          )
        )

        paste(vals, collapse = ";")
      },
      FUN.VALUE = character(1)
    ),
    stringsAsFactors = FALSE
  )

  list(
    matrix = collapsed,
    manifest = manifest
  )
}


# =============================================================================
# Interaction-term extraction
# =============================================================================

find_interaction_term <- function(
  term_names,
  region_term,
  predictor_term
) {

  candidate1 <- paste0(
    region_term,
    ":",
    predictor_term
  )

  candidate2 <- paste0(
    predictor_term,
    ":",
    region_term
  )

  if (candidate1 %in% term_names) {
    return(candidate1)
  }

  if (candidate2 %in% term_names) {
    return(candidate2)
  }

  NA_character_
}


# =============================================================================
# Robust nlme fitting
#
# Primary optimizer: optim
# Automatic fallback: nlminb
#
# A returned fit accompanied by an optimizer/convergence warning is not
# accepted as clean. Such fits are automatically refit with nlminb.
# =============================================================================

fit_lme_with_fallback <- function(
  fixed,
  random,
  data,
  method
) {

  attempt <- function(
    optimizer,
    max_iter,
    max_eval
  ) {

    warnings_seen <- character()

    fit <- withCallingHandlers(

      tryCatch(
        nlme::lme(
          fixed = fixed,
          random = random,
          data = data,
          method = method,
          na.action = na.fail,
          control = nlme::lmeControl(
            opt = optimizer,
            msMaxIter = max_iter,
            msMaxEval = max_eval,
            returnObject = TRUE
          )
        ),
        error = function(e) e
      ),

      warning = function(w) {
        warnings_seen <<- c(
          warnings_seen,
          conditionMessage(w)
        )

        invokeRestart("muffleWarning")
      }
    )

    list(
      fit = fit,
      warnings = unique(warnings_seen)
    )
  }


  first <- attempt(
    optimizer = "optim",
    max_iter = 200,
    max_eval = 400
  )

  first_clean <- (
    !inherits(first$fit, "error") &&
    length(first$warnings) == 0
  )

  if (first_clean) {

    attr(
      first$fit,
      "optimizer_used"
    ) <- "optim"

    attr(
      first$fit,
      "fallback_used"
    ) <- FALSE

    attr(
      first$fit,
      "fallback_reason"
    ) <- NA_character_

    return(first$fit)
  }


  first_reason <- character()

  if (inherits(first$fit, "error")) {
    first_reason <- c(
      first_reason,
      paste0(
        "optim error: ",
        conditionMessage(first$fit)
      )
    )
  }

  if (length(first$warnings) > 0) {
    first_reason <- c(
      first_reason,
      paste0(
        "optim warning: ",
        paste(
          first$warnings,
          collapse = " | "
        )
      )
    )
  }


  second <- attempt(
    optimizer = "nlminb",
    max_iter = 1000,
    max_eval = 2000
  )

  second_clean <- (
    !inherits(second$fit, "error") &&
    length(second$warnings) == 0
  )

  if (second_clean) {

    attr(
      second$fit,
      "optimizer_used"
    ) <- "nlminb"

    attr(
      second$fit,
      "fallback_used"
    ) <- TRUE

    attr(
      second$fit,
      "fallback_reason"
    ) <- paste(
      first_reason,
      collapse = " | "
    )

    return(second$fit)
  }


  second_reason <- character()

  if (inherits(second$fit, "error")) {
    second_reason <- c(
      second_reason,
      paste0(
        "nlminb error: ",
        conditionMessage(second$fit)
      )
    )
  }

  if (length(second$warnings) > 0) {
    second_reason <- c(
      second_reason,
      paste0(
        "nlminb warning: ",
        paste(
          second$warnings,
          collapse = " | "
        )
      )
    )
  }

  simpleError(
    paste(
      c(
        first_reason,
        second_reason
      ),
      collapse = " || "
    )
  )
}


# =============================================================================
# Fit continuous/binary region x predictor mixed model
# =============================================================================

fit_interaction_model <- function(
  accession,
  dl_values,
  st_values,
  phenotype,
  predictor,
  min_people = 25
) {

  ids <- phenotype$individual_id

  wide <- phenotype

  wide$DLPFC <- dl_values[ids]
  wide$STG <- st_values[ids]

  long <- rbind(
    transform(
      wide,
      region = "DLPFC",
      abundance = DLPFC
    ),
    transform(
      wide,
      region = "STG",
      abundance = STG
    )
  )

  long$region <- factor(
    long$region,
    levels = c("DLPFC", "STG")
  )

  required <- c(
    "abundance",
    predictor,
    "region",
    "cohort",
    "age_z",
    "sex",
    "pmi_z",
    "individual_id"
  )

  long <- long[
    complete.cases(
      long[, required, drop = FALSE]
    ),
    ,
    drop = FALSE
  ]

  # Require both regions for every included participant.
  region_counts <- table(
    long$individual_id
  )

  complete_ids <- names(
    region_counts[region_counts == 2]
  )

  long <- long[
    long$individual_id %in% complete_ids,
    ,
    drop = FALSE
  ]

  n_people <- length(
    unique(long$individual_id)
  )

  if (n_people < min_people) {
    return(list(
      status = "insufficient_n",
      result = NULL,
      error = NA_character_
    ))
  }

  if (length(unique(long[[predictor]])) < 2) {
    return(list(
      status = "no_predictor_variation",
      result = NULL,
      error = NA_character_
    ))
  }

  long$individual_id <- factor(
    long$individual_id
  )

  long$cohort <- droplevels(
    factor(
      long$cohort,
      levels = c(
        "Mayo Clinic",
        "Emory"
      )
    )
  )

  long$sex <- droplevels(
    factor(long$sex)
  )

  # Explicit combined-region model.
  #
  # Region-specific nuisance effects are permitted by interacting region
  # with cohort, age, sex, and PMI.
  fixed_formula <- as.formula(
    paste0(
      "abundance ~ region * ",
      predictor,
      " + region * cohort",
      " + region * age_z",
      " + region * sex",
      " + region * pmi_z"
    )
  )

  fit <- fit_lme_with_fallback(
    fixed = fixed_formula,
    random = ~1 | individual_id,
    data = long,
    method = "REML"
  )

  if (inherits(fit, "error")) {
    return(list(
      status = "fit_error",
      result = NULL,
      error = conditionMessage(fit)
    ))
  }

  tt <- summary(fit)$tTable

  interaction_term <- find_interaction_term(
    rownames(tt),
    "regionSTG",
    predictor
  )

  if (is.na(interaction_term)) {
    return(list(
      status = "interaction_not_found",
      result = NULL,
      error = paste(
        "Terms:",
        paste(rownames(tt), collapse = "; ")
      )
    ))
  }

  beta <- unname(
    tt[interaction_term, "Value"]
  )

  se <- unname(
    tt[interaction_term, "Std.Error"]
  )

  df <- unname(
    tt[interaction_term, "DF"]
  )

  t_value <- unname(
    tt[interaction_term, "t-value"]
  )

  p_value <- unname(
    tt[interaction_term, "p-value"]
  )

  crit <- qt(
    0.975,
    df = df
  )

  # Reference-region predictor effect = DLPFC effect.
  if (!(predictor %in% rownames(tt))) {
    return(list(
      status = "main_predictor_not_found",
      result = NULL,
      error = paste(
        "Terms:",
        paste(rownames(tt), collapse = "; ")
      )
    ))
  }

  beta_dlpfc <- unname(
    tt[predictor, "Value"]
  )

  se_dlpfc <- unname(
    tt[predictor, "Std.Error"]
  )

  p_dlpfc <- unname(
    tt[predictor, "p-value"]
  )

  # STG slope = DLPFC main predictor effect + interaction.
  beta_stg <- beta_dlpfc + beta

  out <- data.frame(
    accession = accession,

    n_people = n_people,
    n_observations = nrow(long),

    optimizer_used = attr(
      fit,
      "optimizer_used"
    ),

    used_optimizer_fallback = isTRUE(
      attr(
        fit,
        "fallback_used"
      )
    ),

    fallback_reason = attr(
      fit,
      "fallback_reason"
    ),

    interaction_term = interaction_term,

    beta_region_x_predictor = beta,
    se_region_x_predictor = se,
    df_region_x_predictor = df,
    statistic_region_x_predictor = t_value,
    p_region_x_predictor = p_value,

    ci_low_region_x_predictor = beta - crit * se,
    ci_high_region_x_predictor = beta + crit * se,

    beta_predictor_DLPFC = beta_dlpfc,
    se_predictor_DLPFC = se_dlpfc,
    p_predictor_DLPFC = p_dlpfc,

    beta_predictor_STG = beta_stg,

    stringsAsFactors = FALSE
  )

  list(
    status = if (
      isTRUE(
        attr(
          fit,
          "fallback_used"
        )
      )
    ) {
      "ok_nlminb_fallback"
    } else {
      "ok"
    },
    result = out,
    error = NA_character_
  )
}


# =============================================================================
# Three-level categorical Braak sensitivity
#
# Formal global test:
#   full    = includes region : braak_bin3
#   reduced = removes region : braak_bin3
#
# Models use ML because fixed-effect structures are being compared.
# =============================================================================

fit_braak_bin3_interaction <- function(
  accession,
  dl_values,
  st_values,
  phenotype,
  min_people = 25
) {

  ids <- phenotype$individual_id

  wide <- phenotype
  wide$DLPFC <- dl_values[ids]
  wide$STG <- st_values[ids]

  long <- rbind(
    transform(
      wide,
      region = "DLPFC",
      abundance = DLPFC
    ),
    transform(
      wide,
      region = "STG",
      abundance = STG
    )
  )

  long$region <- factor(
    long$region,
    levels = c("DLPFC", "STG")
  )

  long$braak_bin3 <- factor(
    clean_missing(long$braak_bin3),
    levels = c(
      "I-II",
      "III-IV",
      "V-VI"
    )
  )

  required <- c(
    "abundance",
    "braak_bin3",
    "region",
    "cohort",
    "age_z",
    "sex",
    "pmi_z",
    "individual_id"
  )

  long <- long[
    complete.cases(
      long[, required, drop = FALSE]
    ),
    ,
    drop = FALSE
  ]

  region_counts <- table(
    long$individual_id
  )

  complete_ids <- names(
    region_counts[region_counts == 2]
  )

  long <- long[
    long$individual_id %in% complete_ids,
    ,
    drop = FALSE
  ]

  n_people <- length(
    unique(long$individual_id)
  )

  if (n_people < min_people) {
    return(list(
      status = "insufficient_n",
      result = NULL,
      error = NA_character_
    ))
  }

  long$individual_id <- factor(
    long$individual_id
  )

  long$cohort <- droplevels(
    factor(
      long$cohort,
      levels = c(
        "Mayo Clinic",
        "Emory"
      )
    )
  )

  long$sex <- droplevels(
    factor(long$sex)
  )

  long$braak_bin3 <- droplevels(
    long$braak_bin3
  )

  if (nlevels(long$braak_bin3) < 3) {
    return(list(
      status = "missing_braak_bin",
      result = NULL,
      error = NA_character_
    ))
  }

  full_formula <- abundance ~
    region * braak_bin3 +
    region * cohort +
    region * age_z +
    region * sex +
    region * pmi_z

  reduced_formula <- abundance ~
    region +
    braak_bin3 +
    region * cohort +
    region * age_z +
    region * sex +
    region * pmi_z

  full_fit <- fit_lme_with_fallback(
    fixed = full_formula,
    random = ~1 | individual_id,
    data = long,
    method = "ML"
  )

  if (inherits(full_fit, "error")) {
    return(list(
      status = "full_fit_error",
      result = NULL,
      error = conditionMessage(full_fit)
    ))
  }

  reduced_fit <- fit_lme_with_fallback(
    fixed = reduced_formula,
    random = ~1 | individual_id,
    data = long,
    method = "ML"
  )

  if (inherits(reduced_fit, "error")) {
    return(list(
      status = "reduced_fit_error",
      result = NULL,
      error = conditionMessage(reduced_fit)
    ))
  }

  cmp <- anova(
    reduced_fit,
    full_fit
  )

  global_p <- cmp$`p-value`[2]

  tt <- summary(full_fit)$tTable

  term_mid <- find_interaction_term(
    rownames(tt),
    "regionSTG",
    "braak_bin3III-IV"
  )

  term_high <- find_interaction_term(
    rownames(tt),
    "regionSTG",
    "braak_bin3V-VI"
  )

  extract_term <- function(term) {

    if (is.na(term)) {
      return(
        c(
          beta = NA_real_,
          se = NA_real_,
          df = NA_real_,
          statistic = NA_real_,
          p = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_
        )
      )
    }

    beta <- unname(tt[term, "Value"])
    se <- unname(tt[term, "Std.Error"])
    df <- unname(tt[term, "DF"])
    stat <- unname(tt[term, "t-value"])
    p <- unname(tt[term, "p-value"])

    crit <- qt(
      0.975,
      df = df
    )

    c(
      beta = beta,
      se = se,
      df = df,
      statistic = stat,
      p = p,
      ci_low = beta - crit * se,
      ci_high = beta + crit * se
    )
  }

  mid <- extract_term(term_mid)
  high <- extract_term(term_high)

  out <- data.frame(
    accession = accession,

    n_people = n_people,
    n_observations = nrow(long),

    full_optimizer_used = attr(
      full_fit,
      "optimizer_used"
    ),

    reduced_optimizer_used = attr(
      reduced_fit,
      "optimizer_used"
    ),

    used_optimizer_fallback = (
      isTRUE(
        attr(
          full_fit,
          "fallback_used"
        )
      ) ||
      isTRUE(
        attr(
          reduced_fit,
          "fallback_used"
        )
      )
    ),

    global_region_x_braak_bin3_p = global_p,

    beta_region_x_III_IV = unname(mid["beta"]),
    se_region_x_III_IV = unname(mid["se"]),
    p_region_x_III_IV = unname(mid["p"]),
    ci_low_region_x_III_IV = unname(mid["ci_low"]),
    ci_high_region_x_III_IV = unname(mid["ci_high"]),

    beta_region_x_V_VI = unname(high["beta"]),
    se_region_x_V_VI = unname(high["se"]),
    p_region_x_V_VI = unname(high["p"]),
    ci_low_region_x_V_VI = unname(high["ci_low"]),
    ci_high_region_x_V_VI = unname(high["ci_high"]),

    stringsAsFactors = FALSE
  )

  list(
    status = if (
      isTRUE(
        attr(
          full_fit,
          "fallback_used"
        )
      ) ||
      isTRUE(
        attr(
          reduced_fit,
          "fallback_used"
        )
      )
    ) {
      "ok_nlminb_fallback"
    } else {
      "ok"
    },
    result = out,
    error = NA_character_
  )
}


# =============================================================================
# Load participant metadata
# =============================================================================

cat("======================================================================\n")
cat("FORMAL PAIRED DLPFC/STG REGION x PREDICTOR ANALYSIS\n")
cat("======================================================================\n")

source_meta <- read.csv(
  SOURCE_META,
  check.names = FALSE
)

source_meta$individual_id <- clean_id(
  source_meta$individualID
)

corrected <- read.csv(
  CORRECTED_META,
  check.names = FALSE
)

corrected$individual_id <- clean_id(
  corrected$individual_id
)

if (anyDuplicated(source_meta$individual_id)) {
  stop("Source metadata contains duplicate participant IDs.")
}

if (anyDuplicated(corrected$individual_id)) {
  stop("Corrected metadata contains duplicate participant IDs.")
}

cohort_lookup <- source_meta[
  ,
  c(
    "individual_id",
    "cohort"
  ),
  drop = FALSE
]

phenotype <- merge(
  corrected,
  cohort_lookup,
  by = "individual_id",
  all.x = TRUE,
  validate = "one_to_one"
)

if (!"ad_outcome" %in% names(phenotype)) {
  stop("Corrected metadata lacks ad_outcome.")
}

for (col in c(
  "braak_stage_num",
  "braak_bin3",
  "age_death_num",
  "sex",
  "pmi_num"
)) {
  if (!col %in% names(phenotype)) {
    stop("Corrected metadata lacks required column: ", col)
  }
}

phenotype$braak_stage_num <- suppressWarnings(
  as.numeric(
    phenotype$braak_stage_num
  )
)

phenotype$age_death_num <- suppressWarnings(
  as.numeric(
    phenotype$age_death_num
  )
)

phenotype$pmi_num <- suppressWarnings(
  as.numeric(
    phenotype$pmi_num
  )
)

phenotype$sex <- clean_missing(
  phenotype$sex
)

phenotype$cohort <- clean_missing(
  phenotype$cohort
)

phenotype$ad_outcome <- clean_missing(
  phenotype$ad_outcome
)


# =============================================================================
# Load regional sample maps and identify paired participants
# =============================================================================

dl_map <- read.csv(
  DLPFC_MAP,
  check.names = FALSE
)

st_map <- read.csv(
  STG_MAP,
  check.names = FALSE
)

dl_map$individual_id <- clean_id(
  dl_map$individual_id
)

st_map$individual_id <- clean_id(
  st_map$individual_id
)

paired_ids <- intersect(
  unique(
    dl_map$individual_id[
      !is.na(dl_map$individual_id)
    ]
  ),
  unique(
    st_map$individual_id[
      !is.na(st_map$individual_id)
    ]
  )
)

cat(
  "Paired participants:",
  length(paired_ids),
  "\n"
)

if (length(paired_ids) != 215) {
  stop(
    "Expected 215 paired participants; found ",
    length(paired_ids)
  )
}


# =============================================================================
# Load expression matrices
# =============================================================================

cat("\nLoading DLPFC expression matrix...\n")

dl_obj <- prepare_expression_matrix(
  DLPFC_MATRIX,
  "DLPFC"
)

cat(
  "DLPFC:",
  nrow(dl_obj$matrix),
  "proteins x",
  ncol(dl_obj$matrix),
  "assays\n"
)

cat("\nLoading STG expression matrix...\n")

st_obj <- prepare_expression_matrix(
  STG_MATRIX,
  "STG"
)

cat(
  "STG:",
  nrow(st_obj$matrix),
  "proteins x",
  ncol(st_obj$matrix),
  "assays\n"
)


# =============================================================================
# Collapse technical replicate assays
# =============================================================================

cat("\nCollapsing technical replicates within participant x region...\n")

dl_collapsed <- collapse_technical_replicates(
  dl_obj$matrix,
  dl_map,
  paired_ids,
  "DLPFC"
)

st_collapsed <- collapse_technical_replicates(
  st_obj$matrix,
  st_map,
  paired_ids,
  "STG"
)

technical_manifest <- rbind(
  dl_collapsed$manifest,
  st_collapsed$manifest
)

write.csv(
  technical_manifest,
  file.path(
    OUTDIR,
    "technical_replicate_collapse_manifest.csv"
  ),
  row.names = FALSE
)

cat("\nTechnical assays per participant x region:\n")

print(
  table(
    technical_manifest$region,
    technical_manifest$n_technical_assays
  )
)


# =============================================================================
# Shared protein universe
# =============================================================================

common_accessions <- intersect(
  rownames(dl_collapsed$matrix),
  rownames(st_collapsed$matrix)
)

common_accessions <- sort(
  common_accessions
)

cat(
  "\nShared DLPFC/STG protein accessions:",
  length(common_accessions),
  "\n"
)

if (length(common_accessions) < 1000) {
  stop(
    "Unexpectedly small shared regional protein universe."
  )
}

DL <- dl_collapsed$matrix[
  common_accessions,
  ,
  drop = FALSE
]

ST <- st_collapsed$matrix[
  common_accessions,
  ,
  drop = FALSE
]


# =============================================================================
# Annotation
# =============================================================================

old <- read.csv(
  OLD_RESULTS,
  check.names = FALSE
)

required_annotation <- c(
  "gene_raw",
  "gene_symbol",
  "is_hsp60_client"
)

missing_annotation <- setdiff(
  required_annotation,
  names(old)
)

if (length(missing_annotation) > 0) {
  stop(
    "Old regional table lacks annotation columns: ",
    paste(
      missing_annotation,
      collapse = ", "
    )
  )
}

old$accession <- get_accession(
  old$gene_raw
)

annotation_cols <- intersect(
  c(
    "accession",
    "gene_symbol",
    "is_hsp60_client",
    "is_mitochondrial"
  ),
  names(old)
)

annotation_raw <- unique(
  old[, annotation_cols, drop = FALSE]
)

# Check that an accession does not have conflicting annotations.
annotation_split <- split(
  annotation_raw,
  annotation_raw$accession
)

annotation_conflicts <- vapply(
  annotation_split,
  function(x) {

    any(
      vapply(
        x[, setdiff(names(x), "accession"), drop = FALSE],
        function(v) {
          length(
            unique(
              v[!is.na(v)]
            )
          ) > 1
        },
        logical(1)
      )
    )
  },
  logical(1)
)

if (any(annotation_conflicts)) {
  stop(
    "Conflicting annotations for accessions: ",
    paste(
      head(
        names(annotation_conflicts)[
          annotation_conflicts
        ],
        10
      ),
      collapse = ", "
    )
  )
}

annotation <- do.call(
  rbind,
  lapply(
    annotation_split,
    function(x) x[1, , drop = FALSE]
  )
)

rownames(annotation) <- NULL

annotation$is_hsp60_client <- normalize_logical(
  annotation$is_hsp60_client
)

if ("is_mitochondrial" %in% names(annotation)) {
  annotation$is_mitochondrial <- normalize_logical(
    annotation$is_mitochondrial
  )
}

protein_manifest <- data.frame(
  accession = common_accessions,
  stringsAsFactors = FALSE
)

protein_manifest <- merge(
  protein_manifest,
  annotation,
  by = "accession",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  protein_manifest,
  file.path(
    OUTDIR,
    "shared_regional_protein_universe.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Paired participant phenotype manifest
# =============================================================================

paired_pheno <- phenotype[
  phenotype$individual_id %in% paired_ids,
  ,
  drop = FALSE
]

paired_pheno <- paired_pheno[
  order(
    paired_pheno$individual_id
  ),
  ,
  drop = FALSE
]

if (nrow(paired_pheno) != 215) {
  stop(
    "Expected 215 paired phenotype rows; found ",
    nrow(paired_pheno)
  )
}

# Standardized nuisance covariates for numerical stability.
paired_pheno$age_z <- zscore(
  paired_pheno$age_death_num
)

paired_pheno$pmi_z <- zscore(
  paired_pheno$pmi_num
)

paired_pheno$cohort <- factor(
  paired_pheno$cohort,
  levels = c(
    "Mayo Clinic",
    "Emory"
  )
)

paired_pheno$sex <- factor(
  paired_pheno$sex
)

# AD model eligibility.
paired_pheno$ad_eligible <- (
  paired_pheno$ad_outcome %in% c(
    "Control",
    "AD"
  ) &
  !is.na(paired_pheno$age_z) &
  !is.na(paired_pheno$sex) &
  !is.na(paired_pheno$pmi_z) &
  !is.na(paired_pheno$cohort)
)

paired_pheno$ad_status <- NA_real_

paired_pheno$ad_status[
  paired_pheno$ad_outcome == "Control"
] <- 0

paired_pheno$ad_status[
  paired_pheno$ad_outcome == "AD"
] <- 1

# Braak model eligibility.
paired_pheno$braak_eligible <- (
  paired_pheno$braak_stage_num %in% 1:6 &
  !is.na(paired_pheno$age_z) &
  !is.na(paired_pheno$sex) &
  !is.na(paired_pheno$pmi_z) &
  !is.na(paired_pheno$cohort)
)

ad_pheno <- paired_pheno[
  paired_pheno$ad_eligible,
  ,
  drop = FALSE
]

braak_pheno <- paired_pheno[
  paired_pheno$braak_eligible,
  ,
  drop = FALSE
]

cat(
  "\nAD paired eligible:",
  nrow(ad_pheno),
  "\n"
)

cat(
  "  Control:",
  sum(ad_pheno$ad_status == 0),
  "\n"
)

cat(
  "  AD:",
  sum(ad_pheno$ad_status == 1),
  "\n"
)

cat(
  "\nBraak paired eligible:",
  nrow(braak_pheno),
  "\n"
)

print(
  table(
    braak_pheno$braak_stage_num
  )
)

if (nrow(ad_pheno) != 138) {
  stop(
    "Expected 138 AD-eligible paired participants; found ",
    nrow(ad_pheno)
  )
}

if (sum(ad_pheno$ad_status == 0) != 34) {
  stop("Expected 34 paired Controls.")
}

if (sum(ad_pheno$ad_status == 1) != 104) {
  stop("Expected 104 paired AD participants.")
}

if (nrow(braak_pheno) != 134) {
  stop(
    "Expected 134 Braak-eligible paired participants; found ",
    nrow(braak_pheno)
  )
}

write.csv(
  paired_pheno,
  file.path(
    OUTDIR,
    "paired_participant_analysis_manifest.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Ensure matrix participant coverage
# =============================================================================

for (id in paired_pheno$individual_id) {

  if (!(id %in% colnames(DL))) {
    stop(
      "Paired participant missing from DLPFC collapsed matrix: ",
      id
    )
  }

  if (!(id %in% colnames(ST))) {
    stop(
      "Paired participant missing from STG collapsed matrix: ",
      id
    )
  }
}


# =============================================================================
# Regional scale audit
# =============================================================================

scale_audit <- data.frame(
  accession = common_accessions,

  mean_DLPFC = rowMeans(
    DL,
    na.rm = TRUE
  ),

  mean_STG = rowMeans(
    ST,
    na.rm = TRUE
  ),

  sd_DLPFC = apply(
    DL,
    1,
    sd,
    na.rm = TRUE
  ),

  sd_STG = apply(
    ST,
    1,
    sd,
    na.rm = TRUE
  ),

  stringsAsFactors = FALSE
)

scale_audit$sd_ratio_STG_over_DLPFC <- (
  scale_audit$sd_STG /
  scale_audit$sd_DLPFC
)

write.csv(
  scale_audit,
  file.path(
    OUTDIR,
    "regional_scale_audit.csv"
  ),
  row.names = FALSE
)

finite_ratio <- (
  is.finite(
    scale_audit$sd_ratio_STG_over_DLPFC
  )
)

cat("\nRegional residualized-abundance scale audit:\n")

cat(
  "  Median STG/DLPFC SD ratio:",
  median(
    scale_audit$sd_ratio_STG_over_DLPFC[
      finite_ratio
    ],
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "  IQR STG/DLPFC SD ratio:",
  paste(
    quantile(
      scale_audit$sd_ratio_STG_over_DLPFC[
        finite_ratio
      ],
      probs = c(0.25, 0.75),
      na.rm = TRUE
    ),
    collapse = " - "
  ),
  "\n"
)

cat(
  "  Spearman correlation of protein SDs:",
  cor(
    scale_audit$sd_DLPFC,
    scale_audit$sd_STG,
    method = "spearman",
    use = "complete.obs"
  ),
  "\n"
)


# =============================================================================
# Generic genome-wide runner
# =============================================================================

run_numeric_family <- function(
  predictor,
  phenotype_data,
  label
) {

  cat(
    "\n======================================================================\n"
  )

  cat(
    "RUNNING ",
    label,
    "\n",
    sep = ""
  )

  cat(
    "======================================================================\n"
  )

  results <- vector(
    "list",
    length(common_accessions)
  )

  status <- vector(
    "list",
    length(common_accessions)
  )

  for (i in seq_along(common_accessions)) {

    if (i %% 500 == 0) {
      cat(
        "  ",
        label,
        ": ",
        i,
        " / ",
        length(common_accessions),
        "\n",
        sep = ""
      )
    }

    accession <- common_accessions[i]

    ans <- fit_interaction_model(
      accession = accession,
      dl_values = DL[accession, ],
      st_values = ST[accession, ],
      phenotype = phenotype_data,
      predictor = predictor
    )

    status[[i]] <- data.frame(
      accession = accession,
      model_family = label,
      status = ans$status,
      error = ans$error,
      stringsAsFactors = FALSE
    )

    if (!is.null(ans$result)) {
      results[[i]] <- ans$result
    }
  }

  status <- do.call(
    rbind,
    status
  )

  ok <- !vapply(
    results,
    is.null,
    logical(1)
  )

  if (!any(ok)) {
    stop(
      "No successful models in family: ",
      label
    )
  }

  results <- do.call(
    rbind,
    results[ok]
  )

  results$p_region_x_predictor_fdr <- p.adjust(
    results$p_region_x_predictor,
    method = "BH"
  )

  results <- merge(
    results,
    annotation,
    by = "accession",
    all.x = TRUE,
    sort = FALSE
  )

  results <- results[
    order(
      results$p_region_x_predictor
    ),
    ,
    drop = FALSE
  ]

  list(
    results = results,
    status = status
  )
}


# =============================================================================
# AD interaction
# =============================================================================

ad_run <- run_numeric_family(
  predictor = "ad_status",
  phenotype_data = ad_pheno,
  label = "region_x_AD"
)

write.csv(
  ad_run$results,
  file.path(
    OUTDIR,
    "region_x_AD_interaction_all_proteins.csv"
  ),
  row.names = FALSE
)

write.csv(
  ad_run$status,
  file.path(
    OUTDIR,
    "region_x_AD_model_status.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Continuous Braak interaction
# =============================================================================

braak_run <- run_numeric_family(
  predictor = "braak_stage_num",
  phenotype_data = braak_pheno,
  label = "region_x_Braak_continuous"
)

write.csv(
  braak_run$results,
  file.path(
    OUTDIR,
    "region_x_Braak_continuous_interaction_all_proteins.csv"
  ),
  row.names = FALSE
)

write.csv(
  braak_run$status,
  file.path(
    OUTDIR,
    "region_x_Braak_continuous_model_status.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Three-bin Braak sensitivity
# =============================================================================

cat(
  "\n======================================================================\n"
)

cat(
  "RUNNING region_x_Braak_bin3 sensitivity\n"
)

cat(
  "======================================================================\n"
)

bin3_results <- vector(
  "list",
  length(common_accessions)
)

bin3_status <- vector(
  "list",
  length(common_accessions)
)

for (i in seq_along(common_accessions)) {

  if (i %% 500 == 0) {
    cat(
      "  Braak-bin3: ",
      i,
      " / ",
      length(common_accessions),
      "\n",
      sep = ""
    )
  }

  accession <- common_accessions[i]

  ans <- fit_braak_bin3_interaction(
    accession = accession,
    dl_values = DL[accession, ],
    st_values = ST[accession, ],
    phenotype = braak_pheno
  )

  bin3_status[[i]] <- data.frame(
    accession = accession,
    model_family = "region_x_Braak_bin3",
    status = ans$status,
    error = ans$error,
    stringsAsFactors = FALSE
  )

  if (!is.null(ans$result)) {
    bin3_results[[i]] <- ans$result
  }
}

bin3_status <- do.call(
  rbind,
  bin3_status
)

bin3_ok <- !vapply(
  bin3_results,
  is.null,
  logical(1)
)

if (!any(bin3_ok)) {
  stop(
    "No successful three-bin Braak interaction models."
  )
}

bin3_results <- do.call(
  rbind,
  bin3_results[bin3_ok]
)

bin3_results$global_region_x_braak_bin3_fdr <- p.adjust(
  bin3_results$global_region_x_braak_bin3_p,
  method = "BH"
)

bin3_results <- merge(
  bin3_results,
  annotation,
  by = "accession",
  all.x = TRUE,
  sort = FALSE
)

bin3_results <- bin3_results[
  order(
    bin3_results$global_region_x_braak_bin3_p
  ),
  ,
  drop = FALSE
]

write.csv(
  bin3_results,
  file.path(
    OUTDIR,
    "region_x_Braak_bin3_interaction_all_proteins.csv"
  ),
  row.names = FALSE
)

write.csv(
  bin3_status,
  file.path(
    OUTDIR,
    "region_x_Braak_bin3_model_status.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Hsp60/10-specific exports
#
# FDR columns remain genome-wide BH values. We do NOT recalculate a separate
# Hsp-only FDR here.
# =============================================================================

extract_hsp <- function(x) {

  if (!"is_hsp60_client" %in% names(x)) {
    stop(
      "Results table lacks is_hsp60_client annotation."
    )
  }

  x[
    x$is_hsp60_client %in% TRUE,
    ,
    drop = FALSE
  ]
}

ad_hsp <- extract_hsp(
  ad_run$results
)

braak_hsp <- extract_hsp(
  braak_run$results
)

bin3_hsp <- extract_hsp(
  bin3_results
)

write.csv(
  ad_hsp,
  file.path(
    OUTDIR,
    "region_x_AD_interaction_Hsp60_clients.csv"
  ),
  row.names = FALSE
)

write.csv(
  braak_hsp,
  file.path(
    OUTDIR,
    "region_x_Braak_continuous_interaction_Hsp60_clients.csv"
  ),
  row.names = FALSE
)

write.csv(
  bin3_hsp,
  file.path(
    OUTDIR,
    "region_x_Braak_bin3_interaction_Hsp60_clients.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Model-status audit
# =============================================================================

all_status <- rbind(
  ad_run$status,
  braak_run$status,
  bin3_status
)

write.csv(
  all_status,
  file.path(
    OUTDIR,
    "MASTER_regional_interaction_model_status.csv"
  ),
  row.names = FALSE
)

cat(
  "\n======================================================================\n"
)

cat(
  "MODEL STATUS\n"
)

cat(
  "======================================================================\n"
)

print(
  table(
    all_status$model_family,
    all_status$status
  )
)

actual_errors <- all_status[
  grepl(
    "error|not_found",
    all_status$status
  ),
  ,
  drop = FALSE
]


# =============================================================================
# Summary
# =============================================================================

cat(
  "\n======================================================================\n"
)

cat(
  "RESULT SUMMARY\n"
)

cat(
  "======================================================================\n"
)

cat(
  "Shared proteins:",
  length(common_accessions),
  "\n"
)

cat(
  "\nAD interaction:\n"
)

cat(
  "  Models:",
  nrow(ad_run$results),
  "\n"
)

cat(
  "  Genome-wide FDR < 0.05:",
  sum(
    ad_run$results$p_region_x_predictor_fdr < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "  Hsp clients tested:",
  nrow(ad_hsp),
  "\n"
)

cat(
  "  Hsp clients with genome-wide FDR < 0.05:",
  sum(
    ad_hsp$p_region_x_predictor_fdr < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nContinuous Braak interaction:\n"
)

cat(
  "  Models:",
  nrow(braak_run$results),
  "\n"
)

cat(
  "  Genome-wide FDR < 0.05:",
  sum(
    braak_run$results$p_region_x_predictor_fdr < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "  Hsp clients tested:",
  nrow(braak_hsp),
  "\n"
)

cat(
  "  Hsp clients with genome-wide FDR < 0.05:",
  sum(
    braak_hsp$p_region_x_predictor_fdr < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nThree-bin Braak interaction sensitivity:\n"
)

cat(
  "  Models:",
  nrow(bin3_results),
  "\n"
)

cat(
  "  Global interaction FDR < 0.05:",
  sum(
    bin3_results$global_region_x_braak_bin3_fdr < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "  Hsp clients tested:",
  nrow(bin3_hsp),
  "\n"
)

cat(
  "\nPRIMARY COEFFICIENT INTERPRETATION:\n"
)

cat(
  "  beta_region_x_predictor < 0:\n"
)

cat(
  "    predictor-abundance association is MORE NEGATIVE in STG than DLPFC.\n"
)

cat(
  "  beta_region_x_predictor > 0:\n"
)

cat(
  "    predictor-abundance association is MORE POSITIVE in STG than DLPFC.\n"
)

cat(
  "\nThe primary interaction is estimated directly in one combined-region\n"
)

cat(
  "mixed model; it is NOT obtained by subtracting independently estimated\n"
)

cat(
  "DLPFC and STG coefficients.\n"
)

cat(
  "\nOutputs written to:\n",
  OUTDIR,
  "\n",
  sep = ""
)


# =============================================================================
# Provenance
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

provenance <- data.frame(
  item = c(
    "analysis",
    "DLPFC_matrix",
    "STG_matrix",
    "DLPFC_map",
    "STG_map",
    "source_metadata",
    "corrected_braak_metadata",
    "paired_participants",
    "AD_eligible_participants",
    "Braak_eligible_participants",
    "shared_proteins"
  ),
  value = c(
    "Formal paired DLPFC/STG region x predictor mixed models",
    DLPFC_MATRIX,
    STG_MATRIX,
    DLPFC_MAP,
    STG_MAP,
    SOURCE_META,
    CORRECTED_META,
    length(paired_ids),
    nrow(ad_pheno),
    nrow(braak_pheno),
    length(common_accessions)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  provenance,
  file.path(
    OUTDIR,
    "analysis_provenance.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Hard warning for unexpected model failures
# =============================================================================

if (nrow(actual_errors) > 0) {

  write.csv(
    actual_errors,
    file.path(
      OUTDIR,
      "MODEL_ERRORS_REQUIRING_REVIEW.csv"
    ),
    row.names = FALSE
  )

  stop(
    "\nAnalysis outputs were written, but ",
    nrow(actual_errors),
    " model fits require review.\n",
    "See: ",
    file.path(
      OUTDIR,
      "MODEL_ERRORS_REQUIRING_REVIEW.csv"
    )
  )
}

cat(
  "\nFORMAL REGIONAL INTERACTION ANALYSIS COMPLETED WITH NO MODEL ERRORS.\n"
)
