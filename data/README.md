# Data inputs

Controlled-access ROSMAP data and participant/sample-level files are **not** redistributed in this repository.

Place input files in the following locations after obtaining the required permissions/access:

```text
data/raw/metadata/Analysis_Meta_Merged.csv
data/raw/metadata/ROSMAP_clinical.csv
data/raw/metadata/ROSMAP_biospecimen_metadata.csv
data/raw/derived/ROSMAP_vst_gene_symbol_matrix_STAGE.csv
data/raw/proteomics/C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv
data/raw/proteomics/matched_metadata.csv
data/raw/proteomics/rosmap_50batch_specimen_metadata_for_batch_correction.csv
data/external/1-s2.0-S1355814523012117-MOESM4_ESM.xlsx
data/external/Human.MitoCarta3.0.xls
data/external/AGORA_nominated_targets.csv
```

Accepted alternate filenames are handled in several loader functions for backwards compatibility, but the names above are the clean expected names for a fresh reproducibility run.

Additional supplemental inputs are documented in `data/supplemental_inputs/README.md`.

## Access notes

- ROSMAP transcriptomic, proteomic, clinical, neuropathological, and metadata files should be downloaded by approved users from the AD Knowledge Portal/Synapse and kept local.
- MSBB or other validation cohort files should be kept local if their data-use terms do not allow redistribution.
- Public/external annotation files may be included only if their license and your lab/journal policy permit redistribution. Otherwise, keep them local and cite the source in the manuscript/repository documentation. The Hsp60/10 client inventory citation and expected filename are listed in `data/external/README.md`.
- Do not commit subject-level files, sample-level files, Synapse downloads, credentials, `.Renviron`, or generated outputs.
