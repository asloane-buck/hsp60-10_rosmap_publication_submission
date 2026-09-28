options(stringsAsFactors = FALSE)

if (!requireNamespace("nlme", quietly = TRUE)) {
  stop("Package 'nlme' is required.")
}

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

INDIR <- "outputs/reviewer_revisions/regional_formal_interaction"

OUTDIR <- file.path(
  INDIR,
  "validated_final"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

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

AD_FILE <- file.path(
  INDIR,
  "region_x_AD_interaction_all_proteins.csv"
)

BRAAK_FILE <- file.path(
  INDIR,
  "region_x_Braak_continuous_interaction_all_proteins.csv"
)

BIN3_FILE <- file.path(
  INDIR,
  "region_x_Braak_bin3_interaction_all_proteins.csv"
)

STATUS_FILE <- file.path(
  INDIR,
  "MASTER_regional_interaction_model_status.csv"
)

UNIVERSE_FILE <- file.path(
  INDIR,
  "shared_regional_protein_universe.csv"
)

BRAAK_RETRY_FILE <- file.path(
  INDIR,
  "continuous_Braak_flagged_model_retries.csv"
)


# =============================================================================
# Helpers
# =============================================================================

clean_id <- function(x) {

  x <- trimws(as.character(x))

  x[
    x %in% c(
      "",
      "NA",
      "NaN",
      "nan",
      "<NA>"
    )
  ] <- NA

  sub("\\.0$", "", x)
}


clean_missing <- function(x) {

  x <- trimws(as.character(x))

  x[
    x %in% c(
      "",
      "NA",
      "NaN",
      "nan",
      "<NA>",
      "missing or unknown"
    )
  ] <- NA

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


zscore <- function(x) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  s <- sd(
    x,
    na.rm = TRUE
  )

  if (!is.finite(s) || s == 0) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  as.numeric(
    scale(x)
  )
}


read_expression <- function(file) {

  x <- read.csv(
    file,
    check.names = FALSE
  )

  accession <- get_accession(
    x[[1]]
  )

  if (anyDuplicated(accession)) {
    stop(
      "Duplicate protein accessions in ",
      file
    )
  }

  mat <- as.matrix(
    x[, -1, drop = FALSE]
  )

  storage.mode(mat) <- "numeric"

  rownames(mat) <- accession

  mat
}


collapse_region <- function(
  mat,
  smap,
  paired_ids
) {

  smap$individual_id <- clean_id(
    smap$individual_id
  )

  sm <- smap[
    !is.na(smap$individual_id) &
      smap$individual_id %in% paired_ids &
      smap$matrix_col %in% colnames(mat),
    ,
    drop = FALSE
  ]

  by_person <- split(
    sm$matrix_col,
    sm$individual_id
  )

  ans <- vapply(
    by_person,
    function(cols) {

      y <- rowMeans(
        mat[, cols, drop = FALSE],
        na.rm = TRUE
      )

      y[is.nan(y)] <- NA_real_

      y
    },
    FUN.VALUE = numeric(
      nrow(mat)
    )
  )

  rownames(ans) <- rownames(mat)

  ans
}


make_long <- function(
  dl,
  st,
  phenotype
) {

  d <- phenotype

  d$DLPFC <- dl[
    d$individual_id
  ]

  d$STG <- st[
    d$individual_id
  ]

  long <- rbind(
    transform(
      d,
      region = "DLPFC",
      abundance = DLPFC
    ),
    transform(
      d,
      region = "STG",
      abundance = STG
    )
  )

  long$region <- factor(
    long$region,
    levels = c(
      "DLPFC",
      "STG"
    )
  )

  required <- c(
    "abundance",
    "ad_status",
    "region",
    "cohort",
    "age_z",
    "sex",
    "pmi_z",
    "individual_id"
  )

  long <- long[
    complete.cases(
      long[
        ,
        required,
        drop = FALSE
      ]
    ),
    ,
    drop = FALSE
  ]

  counts <- table(
    long$individual_id
  )

  complete_ids <- names(
    counts[counts == 2]
  )

  long <- long[
    long$individual_id %in%
      complete_ids,
    ,
    drop = FALSE
  ]

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

  long
}


# =============================================================================
# Load existing analysis
# =============================================================================

cat(
  "======================================================================\n"
)

cat(
  "FINALIZING VALIDATED REGIONAL INTERACTION RESULTS\n"
)

cat(
  "======================================================================\n"
)

ad <- read.csv(
  AD_FILE,
  check.names = FALSE
)

braak <- read.csv(
  BRAAK_FILE,
  check.names = FALSE
)

bin3 <- read.csv(
  BIN3_FILE,
  check.names = FALSE
)

status <- read.csv(
  STATUS_FILE,
  check.names = FALSE
)

universe <- read.csv(
  UNIVERSE_FILE,
  check.names = FALSE
)

braak_retry <- read.csv(
  BRAAK_RETRY_FILE,
  check.names = FALSE
)


# =============================================================================
# Validate expected pre-finalization state
# =============================================================================

q96b70_n <- sum(
  ad$accession == "Q96B70",
  na.rm = TRUE
)

valid_ad_handoff <-
  (
    nrow(ad) == 9004 &&
      q96b70_n == 0L
  ) ||
  (
    nrow(ad) == 9005 &&
      q96b70_n == 1L
  )

if (!valid_ad_handoff) {
  stop(
    "Unexpected AD pre-finalization state: ",
    nrow(ad),
    " rows and ",
    q96b70_n,
    " Q96B70 rows. Expected either ",
    "9004 rows with Q96B70 absent or ",
    "9005 rows with Q96B70 present exactly once."
  )
}

q96b70_already_present <-
  q96b70_n == 1L

if (nrow(braak) != 9002) {
  stop(
    "Expected 9002 original continuous Braak models; found ",
    nrow(braak)
  )
}

if (nrow(bin3) != 9001) {
  stop(
    "Expected 9001 Braak-bin3 models; found ",
    nrow(bin3)
  )
}

if (q96b70_already_present) {

  q96_existing <- ad[
    ad$accession == "Q96B70",
    ,
    drop = FALSE
  ]

  required_q96_cols <- c(
    "optimizer_used",
    "used_optimizer_fallback",
    "beta_region_x_predictor",
    "se_region_x_predictor",
    "p_region_x_predictor"
  )

  missing_q96_cols <- setdiff(
    required_q96_cols,
    names(q96_existing)
  )

  if (length(missing_q96_cols) > 0L) {
    stop(
      "Existing Q96B70 row is missing validation columns: ",
      paste(
        missing_q96_cols,
        collapse = ", "
      )
    )
  }

  if (
    as.character(
      q96_existing$optimizer_used[1]
    ) != "nlminb" ||
    !isTRUE(
      as.logical(
        q96_existing$used_optimizer_fallback[1]
      )
    )
  ) {
    stop(
      "Q96B70 is already present, but it is not annotated ",
      "as the validated nlminb fallback."
    )
  }

  cat(
    "\nQ96B70 already present as validated nlminb fallback; ",
    "will independently reconstruct it and verify agreement ",
    "without appending a duplicate.\n"
  )
}

q96fh0_retry <- braak_retry[
  braak_retry$accession == "Q96FH0",
  ,
  drop = FALSE
]

if (nrow(q96fh0_retry) != 1) {
  stop(
    "Expected exactly one Q96FH0 Braak retry row."
  )
}

if (
  q96fh0_retry$retry_status != "ok"
) {
  stop(
    "Q96FH0 nlminb retry was not clean."
  )
}


# =============================================================================
# Replace Q96FH0 continuous-Braak estimate with validated retry
# =============================================================================

hit <- which(
  braak$accession == "Q96FH0"
)

if (length(hit) != 1) {
  stop(
    "Expected exactly one Q96FH0 row in continuous Braak results."
  )
}

braak$beta_region_x_predictor[hit] <-
  q96fh0_retry$retry_beta

braak$se_region_x_predictor[hit] <-
  q96fh0_retry$retry_se

braak$p_region_x_predictor[hit] <-
  q96fh0_retry$retry_p

# Reconstruct t-statistic.
braak$statistic_region_x_predictor[hit] <-
  braak$beta_region_x_predictor[hit] /
  braak$se_region_x_predictor[hit]

# Recalculate 95% CI using stored model DF.
df_q96fh0 <-
  braak$df_region_x_predictor[hit]

crit_q96fh0 <- qt(
  0.975,
  df = df_q96fh0
)

braak$ci_low_region_x_predictor[hit] <-
  braak$beta_region_x_predictor[hit] -
  crit_q96fh0 *
  braak$se_region_x_predictor[hit]

braak$ci_high_region_x_predictor[hit] <-
  braak$beta_region_x_predictor[hit] +
  crit_q96fh0 *
  braak$se_region_x_predictor[hit]


# =============================================================================
# Refit Q96B70 AD using clean nlminb fallback
# =============================================================================

cat(
  "\nReconstructing Q96B70 AD model with nlminb fallback...\n"
)

source <- read.csv(
  SOURCE_META,
  check.names = FALSE
)

source$individual_id <- clean_id(
  source$individualID
)

corrected <- read.csv(
  CORRECTED_META,
  check.names = FALSE
)

corrected$individual_id <- clean_id(
  corrected$individual_id
)

cohort_lookup <- source[
  ,
  c(
    "individual_id",
    "cohort"
  ),
  drop = FALSE
]

pheno <- merge(
  corrected,
  cohort_lookup,
  by = "individual_id",
  all.x = TRUE,
  validate = "one_to_one"
)

pheno$age_death_num <- suppressWarnings(
  as.numeric(
    pheno$age_death_num
  )
)

pheno$pmi_num <- suppressWarnings(
  as.numeric(
    pheno$pmi_num
  )
)

pheno$cohort <- clean_missing(
  pheno$cohort
)

pheno$sex <- clean_missing(
  pheno$sex
)

pheno$ad_outcome <- clean_missing(
  pheno$ad_outcome
)

pheno$age_z <- zscore(
  pheno$age_death_num
)

pheno$pmi_z <- zscore(
  pheno$pmi_num
)

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

if (length(paired_ids) != 215) {
  stop(
    "Expected 215 paired participants."
  )
}

DL_raw <- read_expression(
  DLPFC_MATRIX
)

ST_raw <- read_expression(
  STG_MATRIX
)

DL <- collapse_region(
  DL_raw,
  dl_map,
  paired_ids
)

ST <- collapse_region(
  ST_raw,
  st_map,
  paired_ids
)

if (
  !("Q96B70" %in% rownames(DL)) ||
  !("Q96B70" %in% rownames(ST))
) {
  stop(
    "Q96B70 missing from a regional matrix."
  )
}

ad_pheno <- pheno[
  pheno$individual_id %in% paired_ids &
    pheno$ad_outcome %in% c(
      "Control",
      "AD"
    ) &
    !is.na(pheno$age_z) &
    !is.na(pheno$sex) &
    !is.na(pheno$pmi_z) &
    !is.na(pheno$cohort),
  ,
  drop = FALSE
]

ad_pheno$ad_status <- ifelse(
  ad_pheno$ad_outcome == "AD",
  1,
  0
)

q96_long <- make_long(
  DL["Q96B70", ],
  ST["Q96B70", ],
  ad_pheno
)

fixed_formula <- abundance ~
  region * ad_status +
  region * cohort +
  region * age_z +
  region * sex +
  region * pmi_z

mm <- model.matrix(
  fixed_formula,
  data = q96_long
)

if (
  qr(mm)$rank != ncol(mm)
) {
  stop(
    "Q96B70 fixed-effect design unexpectedly became rank deficient."
  )
}

warnings_seen <- character()

q96_fit <- withCallingHandlers(

  tryCatch(

    nlme::lme(
      fixed = fixed_formula,
      random = ~1 | individual_id,
      data = q96_long,
      method = "REML",
      na.action = na.fail,
      control = nlme::lmeControl(
        opt = "nlminb",
        msMaxIter = 1000,
        msMaxEval = 2000,
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

    invokeRestart(
      "muffleWarning"
    )
  }
)

if (inherits(q96_fit, "error")) {
  stop(
    "Q96B70 nlminb retry failed: ",
    conditionMessage(q96_fit)
  )
}

if (length(warnings_seen) > 0) {
  stop(
    "Q96B70 nlminb retry generated warning(s): ",
    paste(
      unique(warnings_seen),
      collapse = " | "
    )
  )
}

tt <- summary(
  q96_fit
)$tTable

interaction_candidates <- c(
  "regionSTG:ad_status",
  "ad_status:regionSTG"
)

interaction_term <-
  interaction_candidates[
    interaction_candidates %in%
      rownames(tt)
  ]

if (
  length(interaction_term) != 1
) {
  stop(
    "Could not uniquely identify Q96B70 AD interaction term."
  )
}

interaction_term <-
  interaction_term[1]

beta <-
  unname(
    tt[
      interaction_term,
      "Value"
    ]
  )

se <-
  unname(
    tt[
      interaction_term,
      "Std.Error"
    ]
  )

df <-
  unname(
    tt[
      interaction_term,
      "DF"
    ]
  )

stat <-
  unname(
    tt[
      interaction_term,
      "t-value"
    ]
  )

p <-
  unname(
    tt[
      interaction_term,
      "p-value"
    ]
  )

crit <- qt(
  0.975,
  df = df
)

if (
  !("ad_status" %in% rownames(tt))
) {
  stop(
    "Q96B70 AD main predictor term not found."
  )
}

beta_dlpfc <-
  unname(
    tt[
      "ad_status",
      "Value"
    ]
  )

se_dlpfc <-
  unname(
    tt[
      "ad_status",
      "Std.Error"
    ]
  )

p_dlpfc <-
  unname(
    tt[
      "ad_status",
      "p-value"
    ]
  )

beta_stg <-
  beta_dlpfc + beta

q96_annotation <- universe[
  universe$accession ==
    "Q96B70",
  ,
  drop = FALSE
]

if (nrow(q96_annotation) != 1) {
  stop(
    "Expected one Q96B70 annotation row."
  )
}

q96_row <- data.frame(
  accession = "Q96B70",

  n_people = length(
    unique(
      q96_long$individual_id
    )
  ),

  n_observations = nrow(
    q96_long
  ),

  interaction_term =
    interaction_term,

  beta_region_x_predictor =
    beta,

  se_region_x_predictor =
    se,

  df_region_x_predictor =
    df,

  statistic_region_x_predictor =
    stat,

  p_region_x_predictor =
    p,

  ci_low_region_x_predictor =
    beta - crit * se,

  ci_high_region_x_predictor =
    beta + crit * se,

  beta_predictor_DLPFC =
    beta_dlpfc,

  se_predictor_DLPFC =
    se_dlpfc,

  p_predictor_DLPFC =
    p_dlpfc,

  beta_predictor_STG =
    beta_stg,

  stringsAsFactors = FALSE
)

annotation_cols <- setdiff(
  names(ad),
  names(q96_row)
)

for (col in annotation_cols) {

  if (
    col ==
      "p_region_x_predictor_fdr"
  ) {
    next
  }

  if (col %in% names(q96_annotation)) {
    q96_row[[col]] <-
      q96_annotation[[col]][1]
  }
}

# Match result-table columns before adding.
missing_in_q96 <- setdiff(
  names(ad),
  names(q96_row)
)

for (col in missing_in_q96) {
  q96_row[[col]] <- NA
}

q96_row <- q96_row[
  ,
  names(ad),
  drop = FALSE
]

if (q96b70_already_present) {

  compare_cols <- c(
    "beta_region_x_predictor",
    "se_region_x_predictor",
    "p_region_x_predictor"
  )

  existing_values <- as.numeric(
    unlist(
      q96_existing[
        1,
        compare_cols,
        drop = FALSE
      ],
      use.names = FALSE
    )
  )

  reconstructed_values <- as.numeric(
    unlist(
      q96_row[
        1,
        compare_cols,
        drop = FALSE
      ],
      use.names = FALSE
    )
  )

  q96_differences <- abs(
    existing_values -
      reconstructed_values
  )

  names(q96_differences) <-
    compare_cols

  cat(
    "\nQ96B70 existing-vs-reconstructed differences:\n"
  )

  print(q96_differences)

  if (
    any(!is.finite(q96_differences)) ||
    max(q96_differences) > 1e-8
  ) {
    stop(
      "Existing Q96B70 validated row does not reproduce ",
      "the independent nlminb reconstruction."
    )
  }

  cat(
    "Existing Q96B70 row reproduces the independent ",
    "nlminb reconstruction; no duplicate appended.\n"
  )

} else {

  ad <- rbind(
    ad,
    q96_row
  )
}


# =============================================================================
# Recompute multiple-testing corrections
# =============================================================================

if (nrow(ad) != 9005) {
  stop(
    "Expected 9005 validated AD models; found ",
    nrow(ad)
  )
}

if (nrow(braak) != 9002) {
  stop(
    "Expected 9002 validated continuous-Braak models."
  )
}

ad$p_region_x_predictor_fdr <- p.adjust(
  ad$p_region_x_predictor,
  method = "BH"
)

braak$p_region_x_predictor_fdr <- p.adjust(
  braak$p_region_x_predictor,
  method = "BH"
)

# Recalculate bin3 FDR as a reproducibility check.
bin3$global_region_x_braak_bin3_fdr <- p.adjust(
  bin3$global_region_x_braak_bin3_p,
  method = "BH"
)

ad <- ad[
  order(
    ad$p_region_x_predictor
  ),
  ,
  drop = FALSE
]

braak <- braak[
  order(
    braak$p_region_x_predictor
  ),
  ,
  drop = FALSE
]

bin3 <- bin3[
  order(
    bin3$global_region_x_braak_bin3_p
  ),
  ,
  drop = FALSE
]


# =============================================================================
# Update status provenance
# =============================================================================

q96_status <- (
  status$accession == "Q96B70" &
  status$model_family == "region_x_AD"
)

if (sum(q96_status) != 1) {
  stop(
    "Expected one Q96B70 AD status row."
  )
}

status$status[
  q96_status
] <- "ok_nlminb_fallback"

status$error[
  q96_status
] <- NA

q96fh0_status <- (
  status$accession == "Q96FH0" &
  status$model_family ==
    "region_x_Braak_continuous"
)

if (sum(q96fh0_status) != 1) {
  stop(
    "Expected one Q96FH0 continuous-Braak status row."
  )
}

status$status[
  q96fh0_status
] <- "ok_nlminb_fallback"

status$error[
  q96fh0_status
] <- NA


# =============================================================================
# Hsp60/10 exports
# =============================================================================

to_logical <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  toupper(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "TRUE",
    "T",
    "1"
  )
}

ad_hsp <- ad[
  to_logical(
    ad$is_hsp60_client
  ),
  ,
  drop = FALSE
]

braak_hsp <- braak[
  to_logical(
    braak$is_hsp60_client
  ),
  ,
  drop = FALSE
]

bin3_hsp <- bin3[
  to_logical(
    bin3$is_hsp60_client
  ),
  ,
  drop = FALSE
]


# =============================================================================
# Write FINAL validated tables
# =============================================================================

write.csv(
  ad,
  file.path(
    OUTDIR,
    "FINAL_region_x_AD_interaction_all_proteins.csv"
  ),
  row.names = FALSE
)

write.csv(
  braak,
  file.path(
    OUTDIR,
    "FINAL_region_x_Braak_continuous_interaction_all_proteins.csv"
  ),
  row.names = FALSE
)

write.csv(
  bin3,
  file.path(
    OUTDIR,
    "FINAL_region_x_Braak_bin3_interaction_all_proteins.csv"
  ),
  row.names = FALSE
)

write.csv(
  ad_hsp,
  file.path(
    OUTDIR,
    "FINAL_region_x_AD_interaction_Hsp60_clients.csv"
  ),
  row.names = FALSE
)

write.csv(
  braak_hsp,
  file.path(
    OUTDIR,
    "FINAL_region_x_Braak_continuous_interaction_Hsp60_clients.csv"
  ),
  row.names = FALSE
)

write.csv(
  bin3_hsp,
  file.path(
    OUTDIR,
    "FINAL_region_x_Braak_bin3_interaction_Hsp60_clients.csv"
  ),
  row.names = FALSE
)

write.csv(
  status,
  file.path(
    OUTDIR,
    "FINAL_regional_interaction_model_status.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Validation summary
# =============================================================================

q96_final <- ad[
  ad$accession == "Q96B70",
  ,
  drop = FALSE
]

q96fh0_final <- braak[
  braak$accession == "Q96FH0",
  ,
  drop = FALSE
]

summary_table <- data.frame(
  metric = c(
    "AD_models",
    "AD_genomewide_FDR_lt_0.05",
    "AD_Hsp_tested",
    "AD_Hsp_genomewide_FDR_lt_0.05",

    "Braak_continuous_models",
    "Braak_continuous_genomewide_FDR_lt_0.05",
    "Braak_continuous_Hsp_tested",
    "Braak_continuous_Hsp_genomewide_FDR_lt_0.05",

    "Braak_bin3_models",
    "Braak_bin3_global_FDR_lt_0.05",
    "Braak_bin3_Hsp_tested"
  ),

  value = c(
    nrow(ad),

    sum(
      ad$p_region_x_predictor_fdr < 0.05,
      na.rm = TRUE
    ),

    nrow(ad_hsp),

    sum(
      ad_hsp$p_region_x_predictor_fdr < 0.05,
      na.rm = TRUE
    ),

    nrow(braak),

    sum(
      braak$p_region_x_predictor_fdr < 0.05,
      na.rm = TRUE
    ),

    nrow(braak_hsp),

    sum(
      braak_hsp$p_region_x_predictor_fdr < 0.05,
      na.rm = TRUE
    ),

    nrow(bin3),

    sum(
      bin3$global_region_x_braak_bin3_fdr < 0.05,
      na.rm = TRUE
    ),

    nrow(bin3_hsp)
  ),

  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(
    OUTDIR,
    "FINAL_regional_interaction_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# Hard validation
# =============================================================================

if (
  nrow(ad_hsp) != 300 ||
  nrow(braak_hsp) != 300 ||
  nrow(bin3_hsp) != 300
) {
  stop(
    "Expected 300 Hsp clients in each regional result family."
  )
}

if (
  q96_final$is_hsp60_client %in% TRUE
) {
  stop(
    "Q96B70 unexpectedly annotated as Hsp client."
  )
}

if (
  q96fh0_final$is_hsp60_client %in% TRUE
) {
  stop(
    "Q96FH0 unexpectedly annotated as Hsp client."
  )
}

remaining_fit_errors <- status[
  status$status %in% c(
    "fit_error",
    "error",
    "warning"
  ),
  ,
  drop = FALSE
]

if (nrow(remaining_fit_errors) > 0) {

  print(
    remaining_fit_errors
  )

  stop(
    "Unresolved fit errors/warnings remain."
  )
}


# =============================================================================
# Report
# =============================================================================

cat(
  "\n======================================================================\n"
)

cat(
  "FINAL VALIDATED REGIONAL RESULTS\n"
)

cat(
  "======================================================================\n"
)

print(
  summary_table,
  row.names = FALSE
)

cat(
  "\nQ96B70 / LENG9 validated AD fallback:\n"
)

print(
  q96_final[
    ,
    c(
      "accession",
      "gene_symbol",
      "n_people",
      "beta_region_x_predictor",
      "se_region_x_predictor",
      "p_region_x_predictor",
      "p_region_x_predictor_fdr"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat(
  "\nQ96FH0 validated Braak fallback:\n"
)

show_q96fh0 <- c(
  "accession",
  "gene_symbol",
  "is_hsp60_client",
  "n_people",
  "beta_region_x_predictor",
  "se_region_x_predictor",
  "p_region_x_predictor",
  "p_region_x_predictor_fdr"
)

show_q96fh0 <- show_q96fh0[
  show_q96fh0 %in%
    names(q96fh0_final)
]

print(
  q96fh0_final[
    ,
    show_q96fh0,
    drop = FALSE
  ],
  row.names = FALSE
)

cat(
  "\nFinal status table:\n"
)

print(
  table(
    status$model_family,
    status$status
  )
)

cat(
  "\nOutputs:\n",
  OUTDIR,
  "\n",
  sep = ""
)

cat(
  "\nFINAL REGIONAL INTERACTION TABLES VALIDATED.\n"
)
