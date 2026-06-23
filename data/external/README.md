# External annotation files

This directory is intentionally kept empty in the public repository. Place public/non-controlled annotation inputs here after downloading them from the original sources. These files are ignored by default to avoid accidental redistribution of third-party datasets.

## Required files

### Hsp60/10 client inventory

Expected path:

```text
data/external/1-s2.0-S1355814523012117-MOESM4_ESM.xlsx
```

Source/citation:

Bie AS, Cömert C, Körner R, Corydon TJ, Palmfeldt J, Hipp MS, Hartl FU, Bross P. An inventory of interactors of the human HSP60/HSP10 chaperonin in the mitochondrial matrix space. *Cell Stress Chaperones*. 2020;25(3):407-416. doi:10.1007/s12192-020-01080-6.

Repository note:

The analysis uses Supplementary Table S1 from the source publication as the Hsp60/10 client inventory. Download the supplementary Excel file from the publisher/source page and save it using the expected filename above. Do not commit the spreadsheet unless redistribution is explicitly permitted by the source license and approved by the authors/lab.

### MitoCarta 3.0

Expected path:

```text
data/external/Human.MitoCarta3.0.xls
```

Source:

Broad Institute MitoCarta 3.0 download page.

### AGORA nominated targets

Expected path:

```text
data/external/AGORA_nominated_targets.csv
```

Source:

Agora AD Knowledge Portal nominated targets page. Download the table from Agora and save it using the expected filename above.

## Do not commit

Do not commit controlled-access human omics, clinical, neuropathology, sample-level metadata, Synapse downloads, credentials, `.Renviron`, generated outputs, or files whose redistribution terms are unclear.
