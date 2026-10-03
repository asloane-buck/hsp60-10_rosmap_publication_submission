#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 60_run_final_reproducibility_gate.sh
#
# Final fresh-process reproducibility run for the ROSMAP revision
# analysis package. Run from repository root.
#
# This script:
#   1) records Git/code provenance;
#   2) runs the canonical main pipeline;
#   3) exports publication main figures if exporter is present;
#   4) runs the final 10,000-permutation R/92 sensitivity analysis;
#   5) runs the canonical supplemental pipeline;
#   6) reruns the major reviewer-revision analysis/validation chain;
#   7) records checksums, warnings, session information, and Git state.
#
# It does NOT git add/commit/reset/restore/clean anything.
# ============================================================

required_root_files=(
  "R/90_run_main_pipeline.R"
  "R/supplemental/90_run_supplemental_pipeline.R"
  "R/92_pre_submission_sensitivity_checks_PIPELINE_OBJECTS_FAST.R"
)

for f in "${required_root_files[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing $f"
    echo "Run this script from the repository root."
    exit 1
  fi
done

OUT="outputs/reviewer_revisions/final_reproducibility_run"
LOGDIR="$OUT/logs"
mkdir -p "$LOGDIR"

rm -f "$OUT/SUCCESS" "$OUT/FAILED"

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$OUT/start_utc.txt"
git rev-parse HEAD > "$OUT/git_HEAD_before.txt"
git status --short > "$OUT/git_status_before.txt"
git diff -- R > "$OUT/R_code_diff_before.patch"
shasum -a 256 "$OUT/R_code_diff_before.patch" \
  > "$OUT/R_code_diff_before.sha256"

R --version > "$OUT/R_version.txt" 2>&1 || true

printf "phase\tscript\tstatus\n" > "$OUT/phase_status.tsv"

CURRENT_PHASE="initialization"

on_error() {
  local code=$?
  printf "%s\n" "$CURRENT_PHASE" > "$OUT/FAILED"
  echo
  echo "============================================================"
  echo "FINAL REPRODUCIBILITY RUN FAILED"
  echo "Phase: $CURRENT_PHASE"
  echo "Exit status: $code"
  echo "============================================================"
  exit "$code"
}
trap on_error ERR

run_phase() {
  local label="$1"
  shift
  CURRENT_PHASE="$label"

  echo
  echo "============================================================"
  echo "PHASE: $label"
  echo "============================================================"

  local logfile="$LOGDIR/${label}.log"

  if "$@" 2>&1 | tee "$logfile"; then
    printf "%s\t%s\t0\n" "$label" "$*" \
      >> "$OUT/phase_status.tsv"
  else
    local code=${PIPESTATUS[0]}
    printf "%s\t%s\t%s\n" "$label" "$*" "$code" \
      >> "$OUT/phase_status.tsv"
    return "$code"
  fi
}

find_one_revision_script() {
  local prefix="$1"
  local matches=()

  while IFS= read -r path; do
    matches+=("$path")
  done < <(
    find R/reviewer_revisions \
      -maxdepth 1 \
      -type f \
      -name "${prefix}_*.R" \
      -print | sort
  )

  if [[ "${#matches[@]}" -eq 0 ]]; then
    echo "ERROR: no reviewer-revision script found for prefix ${prefix}" >&2
    return 1
  fi

  if [[ "${#matches[@]}" -gt 1 ]]; then
    echo "ERROR: multiple reviewer-revision scripts found for prefix ${prefix}:" >&2
    printf '  %s\n' "${matches[@]}" >&2
    echo "Refusing to guess which one is canonical." >&2
    return 1
  fi

  printf '%s\n' "${matches[0]}"
}

# ------------------------------------------------------------
# Static parse gate for canonical production scripts
# ------------------------------------------------------------

CURRENT_PHASE="parse_canonical_production"

Rscript - <<'RS' 2>&1 | tee "$LOGDIR/parse_canonical_production.log"
files <- c(
  "R/00_config.R",
  "R/01_utils.R",
  "R/02_load_data.R",
  "R/03_build_adjusted_core_objects.R",
  "R/04_build_pathway_sets_all_clients.R",
  "R/05_build_adjusted_all_client_tables.R",
  "R/06_build_cognition_objects.R",
  "R/07_build_matched_null_objects.R",
  "R/10_make_main_figure_1_pathway_remodeling.R",
  "R/11_make_main_figure_2_collapse_heterogeneity.R",
  "R/12_make_main_figure_3_pathology_coupling.R",
  "R/13_make_main_figure_4_cognition.R",
  "R/14_make_main_figure_5_matched_null_specificity.R",
  "R/15_make_main_figure_6_candidate_classification.R",
  "R/90_run_main_pipeline.R",
  "R/92_pre_submission_sensitivity_checks_PIPELINE_OBJECTS_FAST.R",
  "R/supplemental/90_run_supplemental_pipeline.R"
)

missing <- files[!file.exists(files)]
if (length(missing)) {
  stop(
    "Missing canonical production scripts: ",
    paste(missing, collapse = ", ")
  )
}

for (f in files) {
  parse(file = f)
  cat("PARSE OK:", f, "\n")
}
RS

printf "%s\t%s\t0\n" \
  "parse_canonical_production" \
  "canonical production scripts" \
  >> "$OUT/phase_status.tsv"

# ------------------------------------------------------------
# 1. Main production pipeline
# ------------------------------------------------------------

run_phase \
  "main_pipeline" \
  Rscript R/90_run_main_pipeline.R

if [[ -f "R/91_export_main_figures_publication_pdfs.R" ]]; then
  run_phase \
    "main_figure_export" \
    Rscript R/91_export_main_figures_publication_pdfs.R
else
  echo "R/91_export_main_figures_publication_pdfs.R not present; skipping exporter."
  printf "%s\t%s\t0\n" \
    "main_figure_export" \
    "SKIPPED_not_present" \
    >> "$OUT/phase_status.tsv"
fi

# ------------------------------------------------------------
# 2. Final archival R/92 sensitivity: 10,000 permutations
# ------------------------------------------------------------

CURRENT_PHASE="R92_sensitivity_10000"

env \
  SENSITIVITY_N_PERM=10000 \
  SENSITIVITY_N_PERM_CHECK1=10000 \
  SENSITIVITY_N_PERM_CHECK2=10000 \
  SENSITIVITY_WITHIN_STAGE_NULLS=FALSE \
  SENSITIVITY_SOURCE_UPSTREAM_IF_MISSING=TRUE \
  Rscript R/92_pre_submission_sensitivity_checks_PIPELINE_OBJECTS_FAST.R \
  2>&1 | tee "$LOGDIR/R92_sensitivity_10000.log"

printf "%s\t%s\t0\n" \
  "R92_sensitivity_10000" \
  "R/92_pre_submission_sensitivity_checks_PIPELINE_OBJECTS_FAST.R" \
  >> "$OUT/phase_status.tsv"

# ------------------------------------------------------------
# 3. Supplemental production pipeline
# ------------------------------------------------------------

run_phase \
  "supplemental_pipeline" \
  Rscript R/supplemental/90_run_supplemental_pipeline.R

# ------------------------------------------------------------
# 4. Reviewer-revision analysis/validation chain
# ------------------------------------------------------------

revision_scripts=(
  "R/reviewer_revisions/23_audit_conventional_stage_differential_warnings.R"
  "R/reviewer_revisions/24_run_conventional_stage_differential_FINAL.R"
  "R/reviewer_revisions/25_validate_conventional_stage_differential_FINAL.R"

  "R/reviewer_revisions/28_prepare_regional_joint_pathology_metadata.R"
  "R/reviewer_revisions/28_run_formal_regional_interaction.R"
  "R/reviewer_revisions/29_audit_regional_model_convergence.R"
  "R/reviewer_revisions/30_finalize_regional_interaction_results.R"
  "R/reviewer_revisions/31_audit_Hsp_regional_conclusion.R"

  "R/reviewer_revisions/32_audit_variancePartition_inputs.R"
  "R/reviewer_revisions/33_preflight_variancePartition_design.R"
  "R/reviewer_revisions/34_run_variancePartition_source_of_variation.R"
  "R/reviewer_revisions/35_audit_variancePartition_results.R"

  "R/reviewer_revisions/44_run_Hsp_pathway_PC1_sensitivity.R"
  "R/reviewer_revisions/45_validate_Hsp_pathway_PC1_sensitivity.R"

  "R/reviewer_revisions/46_build_alternative_mechanism_marker_manifest.R"
  "R/reviewer_revisions/47_run_alternative_mechanism_marker_stage_analysis.R"
  "R/reviewer_revisions/48_validate_alternative_mechanism_marker_stage_analysis.R"

  "R/reviewer_revisions/49_build_APOE_cohort_demographics.R"
  "R/reviewer_revisions/50_run_APOE_sensitivity.R"
  "R/reviewer_revisions/51_validate_APOE_sensitivity.R"

  "R/reviewer_revisions/52_build_demographic_tables.R"
  "R/reviewer_revisions/53_validate_demographic_tables.R"

  "R/reviewer_revisions/54_audit_reporting_requirements.R"
  "R/reviewer_revisions/56_resolve_reporting_requirements.R"
  "R/reviewer_revisions/57_validate_reporting_requirements_resolution.R"
  "R/reviewer_revisions/58_run_remaining_analysis_gates.R"
  "R/reviewer_revisions/59_validate_remaining_analysis_gates.R"
)

for script in "${revision_scripts[@]}"; do
  if [[ ! -f "$script" ]]; then
    echo "ERROR: missing reviewer-revision script: $script"
    exit 1
  fi

  label="reviewer_$(basename "$script" .R)"
  run_phase "$label" Rscript "$script"
done

# ------------------------------------------------------------
# 5. Final provenance and output inventory
# ------------------------------------------------------------

CURRENT_PHASE="final_provenance"

git rev-parse HEAD > "$OUT/git_HEAD_after.txt"
git status --short > "$OUT/git_status_after.txt"
git diff -- R > "$OUT/R_code_diff_after.patch"
shasum -a 256 "$OUT/R_code_diff_after.patch" \
  > "$OUT/R_code_diff_after.sha256"

if cmp -s \
  "$OUT/R_code_diff_before.patch" \
  "$OUT/R_code_diff_after.patch"; then
  echo "TRUE" > "$OUT/R_code_diff_unchanged.txt"
else
  echo "FALSE" > "$OUT/R_code_diff_unchanged.txt"
fi

Rscript - <<'RS'
out <- "outputs/reviewer_revisions/final_reproducibility_run"
writeLines(
  capture.output(sessionInfo()),
  file.path(out, "sessionInfo_final.txt")
)

pkgs <- installed.packages()
inventory <- data.frame(
  package = rownames(pkgs),
  version = pkgs[, "Version"],
  stringsAsFactors = FALSE
)
inventory <- inventory[order(inventory$package), ]
write.csv(
  inventory,
  file.path(out, "installed_packages.csv"),
  row.names = FALSE
)
RS

{
  find \
    outputs/main_figures/tables \
    outputs/supplemental_figures/audits \
    outputs/supplemental_figures/tables \
    outputs/reviewer_revisions/conventional_stage_differential_final \
    outputs/reviewer_revisions/regional_formal_interaction \
    outputs/reviewer_revisions/variancePartition_source_of_variation \
    outputs/reviewer_revisions/Hsp_pathway_PC1_sensitivity \
    outputs/reviewer_revisions/alternative_mechanism_marker_panel \
    outputs/reviewer_revisions/APOE_sensitivity \
    outputs/reviewer_revisions/demographic_tables \
    outputs/reviewer_revisions/reporting_requirements_resolution \
    outputs/reviewer_revisions/remaining_analysis_gates \
    -type f -name '*.csv' -print 2>/dev/null \
    | sort \
    | while IFS= read -r f; do
        shasum -a 256 "$f"
      done
} > "$OUT/scientific_csv_sha256.txt"

grep -RniE \
  'warning|deprecated' \
  "$LOGDIR" \
  > "$OUT/warning_lines.txt" \
  || true

wc -l < "$OUT/warning_lines.txt" \
  | tr -d ' ' \
  > "$OUT/warning_line_count.txt"

if [[ -f "renv.lock" ]]; then
  echo "TRUE" > "$OUT/renv_lock_present.txt"
  shasum -a 256 renv.lock > "$OUT/renv_lock.sha256"
else
  echo "FALSE" > "$OUT/renv_lock_present.txt"
fi

PERSONAL_USER="ashlynsloane"
MAC_PERSONAL_PATH="/Users/${PERSONAL_USER}"
LINUX_PERSONAL_PATH="/home/${PERSONAL_USER}"
BUCK_WORD="Buck"
INSTITUTE_WORD="Institute"
BUCK_PERSONAL_FRAGMENT="${BUCK_WORD} ${INSTITUTE_WORD}"

if git grep -n -E \
  -e "$MAC_PERSONAL_PATH" \
  -e "$LINUX_PERSONAL_PATH" \
  -e "$BUCK_PERSONAL_FRAGMENT" \
  -- \
  'R/**' \
  'README*' \
  'data/**' \
  ':!outputs/**' \
  > "$OUT/personal_path_scan.txt" 2>/dev/null; then
  echo "TRUE" > "$OUT/personal_path_hits_present.txt"
else
  : > "$OUT/personal_path_scan.txt"
  echo "FALSE" > "$OUT/personal_path_hits_present.txt"
fi

date -u +"%Y-%m-%dT%H:%M:%SZ" > "$OUT/end_utc.txt"
touch "$OUT/SUCCESS"

echo
echo "============================================================"
echo "FINAL REPRODUCIBILITY RUN COMPLETED"
echo "============================================================"
echo "Output: $OUT"
echo "Warning-line count: $(cat "$OUT/warning_line_count.txt")"
echo "R code diff unchanged during execution: $(cat "$OUT/R_code_diff_unchanged.txt")"
echo "renv.lock present: $(cat "$OUT/renv_lock_present.txt")"
echo "Tracked personal-path hits: $(cat "$OUT/personal_path_hits_present.txt")"
echo
echo "Next:"
echo "  Rscript R/reviewer_revisions/61_validate_final_reproducibility_gate.R"
