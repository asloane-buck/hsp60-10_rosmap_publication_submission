# Publication-ready script audit

## Current status

This repository has been sanitized and reorganized for a public code archive. Script names now follow a simpler publication convention:

- `R/90_run_main_pipeline.R` runs the main analysis pipeline and standard main-figure outputs.
- `R/91_export_main_figures_publication_pdfs.R` exports final publication-size main-figure PDFs.
- `R/supplemental/90_run_supplemental_pipeline.R` runs supplemental figures.
- Main figure scripts are named `10_make_main_figure_1_...R` through `15_make_main_figure_6_...R`.
- Supplemental figure scripts are named `10_make_supplementary_figure_1_...R` through `15_make_supplementary_figure_6_...R`.

## Sanitization changes

- Removed personal absolute path fallbacks from main and supplemental configuration.
- Replaced figure-specific hardcoded MitoCarta paths with configuration-based paths.
- Kept controlled-access ROSMAP/MSBB/clinical/proteomics inputs out of the repository.
- Added `data/README.md`, `data/raw/README.md`, `data/external/README.md`, and `data/supplemental_inputs/README.md` to document required local inputs.
- Added `.gitignore` rules to avoid committing raw data, generated tables, generated figures, workspaces, and credentials.
- Excluded generated artifacts such as `Rplots.pdf`.
- Kept the rewritten MSBB supplemental Figure 5 script and excluded the older duplicate version.

## Required local inputs

The code still expects controlled-access or external files to be present locally. These are intentionally not bundled. Key inputs include:

- ROSMAP RNA matrix, proteomics matrix, sample metadata, pathology/covariate metadata, and clinical/cognitive files.
- Hsp60/10 client supplement.
- Human MitoCarta file.
- AGORA nominated target table.
- Supplemental validation/intermediate inputs for MSBB, regional proteomics, matched-individual analyses, and final supplemental figure generation.

## Static checks performed

Prepared source files were scanned for obvious personal absolute local project paths, such as `/path/to/local/project/...`, and credential-like strings. No remaining personal local paths or obvious token/password strings were detected in tracked text files at the time of packaging.

## Not yet verified

The full R pipeline has not been executed in this environment because the controlled-access local inputs are not included. Before public release, run the checklist in `CODEX_REPRODUCIBILITY_CHECKLIST.md` from a clean local clone.

## Recommended public-release sequence

1. Put this repository in a private GitHub repo.
2. Add required local inputs on your machine only; do not commit them.
3. Run `source("R/90_run_main_pipeline.R")`.
4. Run `source("R/91_export_main_figures_publication_pdfs.R")`.
5. Run `source("R/supplemental/90_run_supplemental_pipeline.R")`.
6. Add `renv.lock` or `sessionInfo.txt`.
7. Ask PI/lab to approve data-use compliance.
8. Make public and archive a versioned release on Zenodo.
