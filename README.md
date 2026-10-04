# Hsp60/10 mitochondrial chaperonin client analyses in ROSMAP

This repository contains R code used to generate the main and supplemental manuscript figures for the Hsp60/10 mitochondrial chaperonin client analysis in Alzheimer's disease using ROSMAP transcriptomic, proteomic, neuropathological, cognitive, and validation data.

## Data availability

Controlled-access ROSMAP data are not redistributed in this repository. Approved users should obtain the required transcriptomic, proteomic, clinical, neuropathological, and metadata files from the appropriate AD Knowledge Portal/Synapse study records and place them under `data/raw/` as described in `data/README.md`.

External annotation and validation files should be placed under `data/external/` or `data/supplemental_inputs/`. These files are ignored by Git by default unless licensing and lab/journal policy permit redistribution.

## Recommended entry points

From a fresh R session at the repository root, run the main analysis and standard figure outputs:

```r
source("R/90_run_main_pipeline.R")
```

To export final full-width main-figure PDF composites for manuscript submission:

```r
source("R/91_export_main_figures_publication_pdfs.R")
```

To run supplemental figures after the required local inputs/intermediate tables are available:

```r
source("R/supplemental/90_run_supplemental_pipeline.R")
```

## Script organization

Main pipeline scripts use numeric prefixes in run order:

- `00_config.R`: repository-relative paths, run flags, and output folders.
- `01_utils.R`-`07_build_matched_null_objects.R`: shared utilities and analysis objects.
- `10_make_main_figure_1_pathway_remodeling.R`-`15_make_main_figure_6_candidate_classification.R`: one script per main figure.
- `80_publication_export_helpers.R`: PDF export helpers.
- `90_run_main_pipeline.R`: complete main analysis and standard outputs.
- `91_export_main_figures_publication_pdfs.R`: final publication-size PDF export.

Supplemental scripts live under `R/supplemental/` and follow the same convention: configuration/load/helper scripts first, then one script per supplemental figure, then `90_run_supplemental_pipeline.R` as the wrapper.

## Configuration

Default paths are repository-relative. You can override key files with environment variables:

```bash
export HSP60_ROSMAP_PROJECT_DIR="/path/to/repo"
export HSP60_CLIENT_XLSX="/path/to/client_supplement.xlsx"
export MITOCARTA_XLS="/path/to/Human.MitoCarta3.0.xls"
export ROSMAP_CLINICAL_FILE="/path/to/ROSMAP_clinical.csv"
export AGORA_TARGET_FILE="/path/to/AGORA_nominated_targets.csv"
export MSBB_VALIDATION_FILE="/path/to/MSBB_validation_results.csv"
```

The clean expected local filenames are documented in `data/README.md` and `data/supplemental_inputs/README.md`.

## Software and environment

The analysis environment is recorded in `renv.lock` (R 4.5.2; Bioconductor 3.22). With the corresponding R version and the `renv` package installed, restore the recorded packages from the repository root:

```bash
Rscript -e 'renv::restore(lockfile = "renv.lock", prompt = FALSE)'
```

Figure export may additionally require system fonts and Ghostscript. Package restoration does not download controlled-access data or external annotation files.

## Reviewer revision analyses and validation

Reviewer analysis and audit scripts are in `R/reviewer_revisions/`. For a complete fresh-process run, after preparing the required local inputs:

```bash
bash R/reviewer_revisions/60_run_final_reproducibility_gate.sh
Rscript R/reviewer_revisions/61_validate_final_reproducibility_gate.R
```

The runner executes the main pipeline, publication figure export, 10,000-permutation sensitivity analysis, supplemental pipeline, and reviewer analysis/validation chain. It writes run records to `outputs/reviewer_revisions/final_reproducibility_run/`, including phase exit statuses, code provenance, scientific CSV checksums and session information. The separate validator checks the resulting run records. These records describe the run that produced them; subsequent changes to code or outputs require corresponding verification.

To validate existing remaining-analysis output tables without rerunning model fitting:

```bash
Rscript R/reviewer_revisions/59_validate_remaining_analysis_gates.R
```

Its required inputs are in `outputs/reviewer_revisions/remaining_analysis_gates/`. If those outputs need regeneration, run `R/reviewer_revisions/58_run_remaining_analysis_gates.R` before the validator. This analysis compares continuous and categorical Braak/CERAD models and audits Figure 1B, the matched cohort, BH correction and the use of side-by-side modality effects.

## Submission outputs

- Main publication figures: `outputs/main_figures/plots/manuscript_ready_pdf/`.
- Supplemental figure composites: `outputs/supplemental_figures/figures/`.
- Supplemental figure audit tables: `outputs/supplemental_figures/audits/`.
- Reconciled submission workbook: `supplementary_tables/Hsp60_10_Supplementary_Tables_REVISED.xlsx`.

The reconciled workbook contains Supplementary Tables 1–11. Historical exports and script backups are not substitutes for this workbook. Reviewer working outputs are ignored by Git by default; only individually reviewed aggregate results and provenance files should be selected for publication. Participant-level data and local input files remain excluded.

Figure 5 cognition resampling uses 306 scored clients and 603 scored non-client mitochondrial proteins from 609 detected background proteins. The verified saved-score result has observed mean 4.972222, null mean 3.194182, observed/null ratio 1.556650, and one-sided add-one empirical P = 1/10,001 (10,000 draws; seed 1405). Main Figure 5 and supplemental Figure S3 have separately saved null draws; use each figure's own source tables when checking its reported values.

## Code availability statement draft

Custom R analysis code used to generate the manuscript analyses and figures is available at this GitHub repository: https://github.com/asloane-buck/hsp60-10_rosmap_publication_submission and archived on Zenodo at https://doi.org/10.5281/zenodo.20851582. Controlled-access ROSMAP transcriptomic, proteomic, clinical, neuropathological, and metadata files are available through the applicable AD Knowledge Portal/Synapse access procedures and are not redistributed with this code repository.

## License

Code and repository materials authored by the project authors are released under the MIT License.

This license does not apply to third-party datasets or resources, including ROSMAP, AMP-AD, Synapse, Agora, STRING, MitoCarta, or other externally sourced resources, which remain subject to their own access terms and licenses.
