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

## Software

This repository was prepared from R scripts using the following core packages: `tidyverse`, `janitor`, `readxl`, `broom`, `ggrepel`, `patchwork`, `scales`, `ggplot2`, `grid`, and optional export packages `ragg`, `svglite`, `ggtext`, and Ghostscript/font tools for PDF checks.

For final publication, add either:

- `renv.lock`, generated from the exact analysis environment; or
- a session-info file generated after a successful full run.

## Code availability statement draft

Custom R analysis code used to generate the manuscript analyses and figures is available at [GitHub repository URL] and archived at [Zenodo DOI]. Controlled-access ROSMAP transcriptomic, proteomic, clinical, neuropathological, and metadata files are available through the applicable AD Knowledge Portal/Synapse access procedures and are not redistributed with the code repository.

## License

Code and repository materials authored by the project authors are released under the MIT License.

This license does not apply to third-party datasets or resources, including ROSMAP, AMP-AD, Synapse, Agora, STRING, MitoCarta, or other externally sourced resources, which remain subject to their own access terms and licenses.
