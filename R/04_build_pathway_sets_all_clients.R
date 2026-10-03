############################################################
## 04_build_pathway_sets_all_clients.R
## Build all-client and comparator pathway gene sets.
############################################################

require_objects(
  c("cfg", "prot_mat", "rna_mat", "prot_meta_adj", "rna_meta_adj"),
  context = "04_build_pathway_sets_all_clients.R"
)

read_hsp60_client_inventory <- function() {
  if (!file.exists(cfg$hsp60_client_file)) {
    stop(
      "Cannot find Hsp60/10 client Excel file:\n",
      cfg$hsp60_client_file,
      "\nUpdate cfg$hsp60_client_file in 00_config.R.",
      call. = FALSE
    )
  }

  raw <- tryCatch(
    readxl::read_excel(
      cfg$hsp60_client_file,
      sheet = cfg$hsp60_client_sheet,
      skip = cfg$hsp60_client_skip
    ),
    error = function(e) {
      stop(
        "Could not read Hsp60/10 client Excel file: ", cfg$hsp60_client_file,
        "\nOriginal error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  gene_col <- c("Gene name", "gene", "Gene", "SYMBOL", "symbol", "gene_symbol")[
    c("Gene name", "gene", "Gene", "SYMBOL", "symbol", "gene_symbol") %in% colnames(raw)
  ][1]

  if (is.na(gene_col)) {
    stop(
      "Could not find a gene column in Hsp60/10 client supplement. Available columns: ",
      paste(colnames(raw), collapse = ", "),
      call. = FALSE
    )
  }

  raw |>
    mutate(gene = canonical_gene_symbol(.data[[gene_col]])) |>
    filter(
      !is.na(gene),
      gene != "",
      !gene %in% c("N/A", "NA", "HSPD1", "HSPE1", "HSPE1-MOB4")
    ) |>
    distinct(gene, .keep_all = TRUE)
}

hsp60_10_client_inventory <- read_hsp60_client_inventory()
all_hsp60_10_clients <- hsp60_10_client_inventory$gene

top19_hsp60_10_occupancy_genes <- canonical_gene_symbol(c(
  "ATP5B", "MDH2", "ATP5A1", "TRAP1", "HSPA9",
  "GOT2", "PRDX3", "SHMT2", "HSD17B10", "TUFM",
  "C1QBP", "MT-CO2", "SSBP1", "ECH1", "IDH3A",
  "OAT", "ECHS1", "MRPS23", "ETFB"
))

read_mitocarta_sets <- function() {
  if (!file.exists(cfg$mitocarta_xls)) {
    warning("MitoCarta file not found. Comparator pathway sets will use fallbacks only:\n", cfg$mitocarta_xls, call. = FALSE)
    return(NULL)
  }

  mito_raw <- tryCatch(
    readxl::read_excel(
      cfg$mitocarta_xls,
      sheet = "A Human MitoCarta3.0",
      col_types = "text",
      na = c("", "NA")
    ),
    error = function(e) {
      warning(
        "Could not read MitoCarta file. Comparator pathway sets will use fallbacks only:\n",
        cfg$mitocarta_xls,
        "\nOriginal error: ", conditionMessage(e),
        call. = FALSE
      )
      return(NULL)
    }
  )
  if (is.null(mito_raw)) return(NULL)

  if (!"Symbol" %in% colnames(mito_raw)) {
    warning("MitoCarta file missing Symbol column. Comparator pathway sets will use fallbacks only.", call. = FALSE)
    return(NULL)
  }

  mito_raw |>
    mutate(
      gene = canonical_gene_symbol(Symbol),
      pathway_text = stringr::str_to_lower(
        paste(
          `MitoCarta3.0_SubMitoLocalization`,
          `MitoCarta3.0_MitoPathways`,
          Description,
          sep = " ; "
        )
      )
    ) |>
    filter(!is.na(gene), gene != "") |>
    distinct(gene, .keep_all = TRUE)
}

mitocarta_tbl <- read_mitocarta_sets()

make_mito_set <- function(pattern) {
  if (is.null(mitocarta_tbl)) return(character(0))
  mitocarta_tbl |>
    filter(str_detect(pathway_text, regex(pattern, ignore_case = TRUE))) |>
    pull(gene) |>
    clean_gene_symbols()
}

mitocarta_genes <- if (is.null(mitocarta_tbl)) character(0) else clean_gene_symbols(mitocarta_tbl$gene)
broad_mito_non_hsp60_10 <- setdiff(mitocarta_genes, all_hsp60_10_clients)

pathway_gene_sets <- list(
  Hsp60_10_all_clients = all_hsp60_10_clients,
  Hsp60_10_top19_occupancy = top19_hsp60_10_occupancy_genes,
  Broad_MitoCarta_non_Hsp60_10 = broad_mito_non_hsp60_10,
  Mito_translation_non_Hsp60_10 = setdiff(make_mito_set("translation|ribosom|mitoribosom|trna|rrna|rna processing|mitochondrial gene expression"), all_hsp60_10_clients),
  OXPHOS_ETC_non_Hsp60_10 = setdiff(make_mito_set("oxidative phosphorylation|oxphos|respiratory chain|electron transport|complex i|complex ii|complex iii|complex iv|complex v|atp synth"), all_hsp60_10_clients),
  TCA_pyruvate_metabolism_non_Hsp60_10 = setdiff(make_mito_set("tca|citric acid|krebs|pyruvate|pdh|tricarboxylic|dehydrogenase|acetyl"), all_hsp60_10_clients),
  FAO_metabolism_non_Hsp60_10 = setdiff(make_mito_set("fatty acid|beta oxidation|β-oxidation|acyl|carnitine|lipid"), all_hsp60_10_clients),
  Mito_protein_quality_control_non_Hsp60_10 = setdiff(make_mito_set("protein import|protein sorting|protein folding|protease|quality control|chaperone|unfolded|aaa"), all_hsp60_10_clients),
  Proteasome_core_non_Hsp60_10 = setdiff(canonical_gene_symbol(c(
    "PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7","PSMA8",
    "PSMB1","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMB10","PSMB11",
    "PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6",
    "PSMD1","PSMD2","PSMD3","PSMD4","PSMD5","PSMD6","PSMD7","PSMD8","PSMD9","PSMD10","PSMD11","PSMD12","PSMD13","PSMD14",
    "PSME1","PSME2","PSME3","PSME4","PSMF1","POMP","ADRM1","UCHL5","USP14"
  )), all_hsp60_10_clients),
  Lysosome_core_non_Hsp60_10 = setdiff(canonical_gene_symbol(c(
    "LAMP1","LAMP2","LAMP3","CTSA","CTSB","CTSC","CTSD","CTSF","CTSH","CTSK","CTSL","CTSO","CTSS","CTSV","CTSZ",
    "GAA","GBA1","GLB1","HEXA","HEXB","IDUA","IDS","ARSA","ARSB","GALC","MAN2B1","NAGLU","SGSH","GUSB","FUCA1","NEU1",
    "NPC1","NPC2","ATP6V0A1","ATP6V0A2","ATP6V0A4","ATP6V0B","ATP6V0C","ATP6V0D1","ATP6V0D2","ATP6V0E1","ATP6V0E2",
    "ATP6V1A","ATP6V1B1","ATP6V1B2","ATP6V1C1","ATP6V1C2","ATP6V1D","ATP6V1E1","ATP6V1E2","ATP6V1F","ATP6V1G1","ATP6V1G2","ATP6V1G3","ATP6V1H",
    "TFEB","TFE3","MITF","MCOLN1","RAB7A","RAB5A","RAB5B","RAB5C","RAB9A","RAB11A"
  )), all_hsp60_10_clients),
  UPRmt_core_non_Hsp60_10 = setdiff(canonical_gene_symbol(c(
    "ATF5","ATF4","DDIT3","DDIT4","CEBPB","DELE1","EIF2AK1","OMA1",
    "HSPD1","HSPE1","HSPA9","DNAJA3","LONP1","CLPP","CLPX","HTRA2","YME1L1","AFG3L2","SPG7"
  )), all_hsp60_10_clients)
)

pathway_gene_sets <- purrr::map(pathway_gene_sets, clean_gene_symbols)

prot_scores <- compute_pathway_scores(
  prot_mat,
  pathway_gene_sets,
  sample_col = "SampleID"
) |>
  checked_left_join(
    prot_meta_adj |>
      dplyr::select(
        SampleID,
        IndividualID,
        cogdx_num,
        clinical_stage,
        clinical_stage_broad,
        protein_legacy_dx
      ),
    by = "SampleID",
    label = "protein pathway scores to canonical clinical metadata"
  )

rna_scores <- compute_pathway_scores(
  rna_mat,
  pathway_gene_sets,
  sample_col = "sample_id"
) |>
  checked_left_join(
    rna_meta_adj |>
      dplyr::select(
        sample_id,
        individual_id,
        cogdx_num,
        clinical_stage,
        clinical_stage_broad
      ),
    by = "sample_id",
    label = "RNA pathway scores to canonical clinical metadata"
  )

pathway_overlap_audit <- imap_dfr(pathway_gene_sets, function(genes, pathway) {
  is_hsp_set <- str_detect(pathway, regex("^Hsp60_10", ignore_case = TRUE))
  tibble(
    pathway = pathway,
    is_hsp60_10_set = is_hsp_set,
    n_genes = length(genes),
    n_overlap_hsp60_10_clients = length(intersect(genes, all_hsp60_10_clients)),
    overlap_hsp60_10_clients = paste(intersect(genes, all_hsp60_10_clients), collapse = ";")
  )
})

bad <- pathway_overlap_audit |> filter(!is_hsp60_10_set, n_overlap_hsp60_10_clients > 0)
if (nrow(bad) > 0) {
  print(bad)
  stop("Comparator sets still overlap with Hsp60/10 clients.", call. = FALSE)
}

pathway_detection_audit <- imap_dfr(pathway_gene_sets, function(genes, pathway) {
  tibble(
    pathway = pathway,
    n_total = length(genes),
    n_detected_protein = length(intersect(genes, rownames(prot_mat))),
    n_detected_rna = length(intersect(genes, rownames(rna_mat))),
    n_detected_both = length(intersect(intersect(genes, rownames(prot_mat)), rownames(rna_mat)))
  )
})

write_tbl(hsp60_10_client_inventory, "hsp60_10_client_inventory_from_supplement")
write_tbl(pathway_overlap_audit, "pathway_gene_sets_overlap_audit")
write_tbl(pathway_detection_audit, "pathway_gene_sets_detection_audit")
save_obj(hsp60_10_client_inventory, "hsp60_10_client_inventory")
save_obj(all_hsp60_10_clients, "all_hsp60_10_clients")
save_obj(top19_hsp60_10_occupancy_genes, "top19_hsp60_10_occupancy_genes")
save_obj(pathway_gene_sets, "pathway_gene_sets")
save_obj(prot_scores, "prot_scores")
save_obj(rna_scores, "rna_scores")

saveRDS(pathway_gene_sets, file.path(cfg$pathway_dir, "pathway_gene_sets_all_clients_NO_Hsp60_10_overlap.rds"))

message("Loaded 04_build_pathway_sets_all_clients.R")
