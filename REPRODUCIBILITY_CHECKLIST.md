# Local reproducibility checklist

Use this checklist before making the repository public or archiving it on Zenodo.

## 1. Confirm no restricted data are tracked

Run from the repository root:

```bash
git status --short
git ls-files data/raw data/external data/supplemental_inputs
```

Expected: no controlled-access ROSMAP/MSBB/clinical/proteomics files are tracked.

## 2. Static scan for local paths and secrets

```bash
LOCAL_PROJECT_DIR="/path/to/local/project"
grep -RInE "${LOCAL_PROJECT_DIR}|SynapseToken|TOKEN|PASSWORD|SECRET|\.Renviron" . \
  --exclude-dir=.git \
  --exclude='*.pdf'
```

Expected: no personal paths, credentials, or tokens in tracked source files.

## 3. Run the main pipeline from a clean session

```r
rm(list = ls())
source("R/90_run_main_pipeline.R")
```

Expected: all six main figure scripts complete or fail only because a documented local input file is missing.

## 4. Run final main-figure PDF export

```r
rm(list = ls())
source("R/91_export_main_figures_publication_pdfs.R")
```

Expected: final PDFs are written under the configured plot output directory.

## 5. Run supplemental pipeline

```r
rm(list = ls())
source("R/supplemental/90_run_supplemental_pipeline.R")
```

Expected: Supplemental Figures 1-6 complete or fail only because a documented local input/intermediate file is missing.

## 6. Capture environment

Preferred:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

Alternative:

```r
sink("sessionInfo.txt")
sessionInfo()
sink()
```

## 7. Release workflow

1. Commit the sanitized source code only.
2. Push to a private GitHub repository.
3. Ask a lab member to run this checklist from a clean clone.
4. Make the repository public only after PI approval.
5. Create a GitHub release and archive that release with Zenodo.
