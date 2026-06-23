# Supplemental inputs

Supplemental scripts can use outputs from the main pipeline and/or local intermediate files. These files are ignored by Git by default.

Expected clean filenames, where applicable:

```text
data/supplemental_inputs/main_analysis_workspace.RData
data/supplemental_inputs/robustness_workspace.RData
data/supplemental_inputs/ROSMAP_vsd.rds
data/supplemental_inputs/Hsp60_10_all_clients.csv
data/supplemental_inputs/Broad_MitoCarta_non_Hsp60_10.csv
data/supplemental_inputs/TCA_pyruvate_metabolism_non_Hsp60_10.csv
data/supplemental_inputs/AMP_covariate_adjusted_regional_screen_labeled_Hsp60_mito.csv
data/supplemental_inputs/AMP_covariate_adjusted_Braak_models_by_region.csv
data/supplemental_inputs/AMP_covariate_adjusted_Hsp60_clients_DLPFC_vs_STG.csv
data/supplemental_inputs/MSBB_covariate_adjusted_braak_and_collapse.csv
data/supplemental_inputs/MSBB_standardized_braak_collapse_percentiles.csv
data/supplemental_inputs/MSBB_covariate_adjusted_braak_collapse_percentiles.csv
data/supplemental_inputs/COMBINED_ROSMAP_MSBB_validation_results.csv
data/supplemental_inputs/MSBB_braak_effects.csv
data/supplemental_inputs/MSBB_collapse.csv
```

Some supplemental panels can also find intermediate objects loaded from the optional `.RData` workspaces. Do not commit those workspaces unless your PI/lab has explicitly approved redistribution and all data-use terms permit it.
