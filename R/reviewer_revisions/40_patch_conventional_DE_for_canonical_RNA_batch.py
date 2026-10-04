from pathlib import Path
import shutil

p24 = Path(
    "R/reviewer_revisions/"
    "24_run_conventional_stage_differential_FINAL.R"
)

p25 = Path(
    "R/reviewer_revisions/"
    "25_validate_conventional_stage_differential_FINAL.R"
)

for p in (p24, p25):
    if not p.exists():
        raise SystemExit(f"Missing file: {p}")

s24 = p24.read_text()
s25 = p25.read_text()

# ------------------------------------------------------------------
# Script 24: canonical RNA stage-count guard
# ------------------------------------------------------------------

old1 = (
    "expected_rna_counts <- "
    "c(NCI = 200L, MCI = 158L, AD = 220L)"
)

new1 = (
    "expected_rna_counts <- "
    "c(NCI = 200L, MCI = 158L, AD = 219L)"
)

if s24.count(old1) != 1:
    raise SystemExit(
        "Expected exactly one old RNA stage-count guard in script 24; "
        f"found {s24.count(old1)}."
    )

s24 = s24.replace(old1, new1, 1)


# ------------------------------------------------------------------
# Script 24: output validation/audit row
# ------------------------------------------------------------------

old2 = (
    '"RNA AD", as.character(observed_rna_counts["AD"]), '
    '"220", observed_rna_counts["AD"] == 220,'
)

new2 = (
    '"RNA AD", as.character(observed_rna_counts["AD"]), '
    '"219", observed_rna_counts["AD"] == 219,'
)

if s24.count(old2) != 1:
    raise SystemExit(
        "Expected exactly one old RNA AD audit row in script 24; "
        f"found {s24.count(old2)}."
    )

s24 = s24.replace(old2, new2, 1)


# ------------------------------------------------------------------
# Script 25: final validator expected RNA counts
# ------------------------------------------------------------------
#
# From the audited validator block:
#   "200",
#   "158",
#   "220",
#   "167",
#   "96",
#   "109",
#
# Change only that exact sequence so an unrelated 220 cannot be altered.
# ------------------------------------------------------------------

old3 = '''    "200",
    "158",
    "220",
    "167",
    "96",
    "109",'''

new3 = '''    "200",
    "158",
    "219",
    "167",
    "96",
    "109",'''

if s25.count(old3) != 1:
    raise SystemExit(
        "Expected exactly one old RNA/protein count block in script 25; "
        f"found {s25.count(old3)}."
    )

s25 = s25.replace(old3, new3, 1)


# ------------------------------------------------------------------
# Static safety checks
# ------------------------------------------------------------------

if "AD = 220L" in s24:
    raise SystemExit(
        "Stale AD = 220L remains in script 24."
    )

if (
    '"RNA AD", as.character(observed_rna_counts["AD"]), '
    '"220"' in s24
):
    raise SystemExit(
        "Stale RNA AD expected count remains in script 24."
    )

if old3 in s25:
    raise SystemExit(
        "Stale RNA count block remains in script 25."
    )


# ------------------------------------------------------------------
# Preserve immediate pre-patch files
# ------------------------------------------------------------------

backup = Path(
    "outputs/reviewer_revisions/"
    "RNA_batch_canonicalization/pre_patch_snapshot"
)

backup.mkdir(
    parents=True,
    exist_ok=True
)

shutil.copy2(
    p24,
    backup /
    "24_run_conventional_stage_differential_FINAL_PRE_BATCH_PATCH.R"
)

shutil.copy2(
    p25,
    backup /
    "25_validate_conventional_stage_differential_FINAL_PRE_BATCH_PATCH.R"
)

p24.write_text(s24)
p25.write_text(s25)

print("PATCH APPLIED")
print("Script 24:")
print("  RNA NCI/MCI/AD expected = 200/158/219")
print("  RNA AD audit expected = 219")
print("  DESeq2 still inherits production rna_covars")
print("Script 25:")
print("  final validator RNA AD expected = 219")
