from pathlib import Path
import shutil
import re

P02 = Path("R/02_load_data.R")
P03 = Path("R/03_build_adjusted_core_objects.R")

for p in (P02, P03):
    if not p.exists():
        raise SystemExit(f"Missing required production file: {p}")

s02 = P02.read_text()
s03 = P03.read_text()

# Refuse accidental double patching.
if "rna_excluded_ambiguous_batch_samples" in s02:
    raise SystemExit("R/02 already appears batch-patched.")

if "sequencing_batch_factor" in s03:
    raise SystemExit("R/03 already appears batch-patched.")


# =====================================================================
# R/02: centrally exclude the ambiguous RNA batch sample
# =====================================================================

anchor = "## RNA primary clinical-stage cohort remains the narrow cogdx 1/2/4"

if s02.count(anchor) != 1:
    raise SystemExit(
        f"Expected exactly one RNA cohort anchor in R/02; "
        f"found {s02.count(anchor)}"
    )

insert = '''## RNA primary clinical-stage cohort remains the narrow cogdx 1/2/4
##
## RNA technical-QC exclusion:
## 492_120515 has the ambiguous composite sequencing-batch annotation
## "0, 6, 7" and cannot be assigned to one technical batch.
rna_excluded_ambiguous_batch_samples <- c("492_120515")

if (!all(rna_excluded_ambiguous_batch_samples %in% analysis_meta$sample_id)) {
  stop(
    "Expected ambiguous-batch RNA sample was not present in analysis_meta.",
    call. = FALSE
  )
}
'''

s02 = s02.replace(anchor, insert, 1)

# Insert exclusion into the existing RNA metadata filter.
pattern = re.compile(
    r'(\.data\$sample_id\s+%in%\s+colnames\(rna_mat_raw\),\s*\n'
    r'\s*!is\.na\(\.data\$clinical_stage\))'
)

m = pattern.search(s02)

if m is None:
    raise SystemExit(
        "Could not uniquely locate RNA cohort filter in R/02."
    )

replacement = (
    m.group(1)
    + ",\n"
    + "    !.data$sample_id %in% rna_excluded_ambiguous_batch_samples"
)

s02 = s02[:m.start()] + replacement + s02[m.end():]

# Add hard cohort validation after RNA matrix/metadata alignment.
anchor2 = (
    "## Backward-compatible matrix names. No stage filtering is applied to protein."
)

if s02.count(anchor2) != 1:
    raise SystemExit(
        "Could not uniquely locate post-RNA-alignment anchor in R/02."
    )

validation = r'''## Hard validation of canonical RNA technical-QC cohort.
if ("492_120515" %in% rna_meta$sample_id ||
    "492_120515" %in% colnames(rna_mat_raw)) {
  stop(
    "Ambiguous sequencing-batch sample 492_120515 remains in RNA cohort.",
    call. = FALSE
  )
}

if (!"sequencing_batch" %in% colnames(rna_meta)) {
  stop(
    "RNA metadata is missing sequencing_batch.",
    call. = FALSE
  )
}

if (anyNA(rna_meta$sequencing_batch)) {
  stop(
    "Canonical RNA cohort contains missing sequencing_batch values.",
    call. = FALSE
  )
}

if (any(as.character(rna_meta$sequencing_batch) == "0, 6, 7")) {
  stop(
    'Canonical RNA cohort still contains ambiguous batch "0, 6, 7".',
    call. = FALSE
  )
}

if (ncol(rna_mat_raw) != 577L || nrow(rna_meta) != 577L) {
  stop(
    "Expected canonical RNA cohort n=577; found matrix n=",
    ncol(rna_mat_raw),
    " and metadata n=",
    nrow(rna_meta),
    ".",
    call. = FALSE
  )
}

expected_rna_stage_counts <- c(
  NCI = 200L,
  MCI = 158L,
  AD = 219L
)

observed_rna_stage_counts <- vapply(
  names(expected_rna_stage_counts),
  function(stage_name) {
    sum(as.character(rna_meta$clinical_stage) == stage_name)
  },
  integer(1)
)

if (!identical(
  unname(observed_rna_stage_counts),
  unname(expected_rna_stage_counts)
)) {
  stop(
    "Unexpected canonical RNA stage counts: ",
    paste(
      names(observed_rna_stage_counts),
      observed_rna_stage_counts,
      sep = "=",
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (dplyr::n_distinct(rna_meta$sequencing_batch) != 9L) {
  stop(
    "Expected 9 valid RNA sequencing batches after QC exclusion; found ",
    dplyr::n_distinct(rna_meta$sequencing_batch),
    ".",
    call. = FALSE
  )
}

'''

s02 = s02.replace(anchor2, validation + anchor2, 1)


# =====================================================================
# R/03: explicitly construct and require sequencing batch
# =====================================================================

# Require the source column.
old_required = '''  c("sample_id", "clinical_stage"),
  "RNA adjusted metadata"
'''

new_required = '''  c("sample_id", "clinical_stage", "sequencing_batch"),
  "RNA adjusted metadata"
'''

if s03.count(old_required) != 1:
    raise SystemExit(
        "Could not uniquely locate RNA require_columns block in R/03."
    )

s03 = s03.replace(old_required, new_required, 1)


# Create canonical factor before RNA metadata is reordered to the matrix.
align_anchor = '''rna_meta_adj <- rna_meta_adj[
  match(colnames(rna_mat_raw), rna_meta_adj$sample_id),
'''

if s03.count(align_anchor) != 1:
    raise SystemExit(
        "Could not uniquely locate RNA metadata alignment in R/03."
    )

batch_block = r'''if (anyNA(rna_meta_adj$sequencing_batch)) {
  stop(
    "RNA adjusted metadata contains missing sequencing_batch values.",
    call. = FALSE
  )
}

if (any(as.character(rna_meta_adj$sequencing_batch) == "0, 6, 7")) {
  stop(
    'RNA adjusted metadata contains excluded ambiguous batch "0, 6, 7".',
    call. = FALSE
  )
}

rna_meta_adj$sequencing_batch_factor <- factor(
  as.character(rna_meta_adj$sequencing_batch)
)

if (nlevels(rna_meta_adj$sequencing_batch_factor) != 9L) {
  stop(
    "Expected 9 RNA sequencing-batch levels; found ",
    nlevels(rna_meta_adj$sequencing_batch_factor),
    ".",
    call. = FALSE
  )
}

'''

s03 = s03.replace(
    align_anchor,
    batch_block + align_anchor,
    1
)


# Replace ONLY the RNA batch_factor candidate.
start = s03.find("rna_covars <- select_usable_covars(")
end = s03.find("protein_covars <- select_usable_covars(")

if start == -1 or end == -1 or end <= start:
    raise SystemExit(
        "Could not isolate RNA covariate-selection block in R/03."
    )

rna_block = s03[start:end]

if rna_block.count('"batch_factor"') != 1:
    raise SystemExit(
        "Expected exactly one batch_factor in RNA covariate block; "
        f"found {rna_block.count(chr(34) + 'batch_factor' + chr(34))}."
    )

rna_block_new = rna_block.replace(
    '"batch_factor"',
    '"sequencing_batch_factor"',
    1
)

s03 = s03[:start] + rna_block_new + s03[end:]


# Require every canonical RNA nuisance covariate.
protein_anchor = "protein_covars <- select_usable_covars("

required_rna = r'''required_rna_covars <- c(
  "age_num",
  "sex_factor",
  "pmi_num",
  "rin_num",
  "sequencing_batch_factor"
)

missing_required_rna_covars <- setdiff(
  required_rna_covars,
  rna_covars
)

if (length(missing_required_rna_covars) > 0L) {
  stop(
    "Required RNA nuisance covariates were not usable: ",
    paste(missing_required_rna_covars, collapse = ", "),
    call. = FALSE
  )
}

'''

if s03.count(protein_anchor) != 1:
    raise SystemExit(
        "Could not uniquely locate protein_covars block in R/03."
    )

s03 = s03.replace(
    protein_anchor,
    required_rna + protein_anchor,
    1
)


# =====================================================================
# Static validation before touching files
# =====================================================================

checks = {
    "R/02 exclusion": "rna_excluded_ambiguous_batch_samples" in s02,
    "R/02 n=577 guard": "Expected canonical RNA cohort n=577" in s02,
    "R/02 stage guard": "AD = 219L" in s02,
    "R/02 batch-level guard": "!= 9L" in s02,
    "R/03 sequencing factor":
        "rna_meta_adj$sequencing_batch_factor <- factor" in s03,
    "R/03 RNA covariate":
        '"sequencing_batch_factor"' in s03[start:],
    "R/03 required RNA covars":
        "required_rna_covars <- c(" in s03,
}

bad = [name for name, ok in checks.items() if not ok]

if bad:
    raise SystemExit(
        "Static validation failed before write: "
        + ", ".join(bad)
    )


# =====================================================================
# Back up current files and write
# =====================================================================

backup_dir = Path(
    "outputs/reviewer_revisions/"
    "RNA_batch_canonicalization/pre_patch_snapshot"
)

backup_dir.mkdir(
    parents=True,
    exist_ok=True
)

shutil.copy2(
    P02,
    backup_dir / "R_02_load_data_PRE_CANONICAL_BATCH_PATCH.R"
)

shutil.copy2(
    P03,
    backup_dir / "R_03_build_adjusted_core_objects_PRE_CANONICAL_BATCH_PATCH.R"
)

P02.write_text(s02)
P03.write_text(s03)

print("PATCH APPLIED")
print("")
print("R/02_load_data.R")
print("  - excludes 492_120515")
print("  - requires RNA n = 577")
print("  - requires NCI/MCI/AD = 200/158/219")
print("  - requires 9 valid sequencing batches")
print("")
print("R/03_build_adjusted_core_objects.R")
print("  - constructs sequencing_batch_factor explicitly")
print("  - replaces unusable RNA batch_factor candidate")
print("  - RNA adjustment is now:")
print("    age + sex + PMI + RIN + sequencing batch")
print("  - all five RNA nuisance covariates are mandatory")
