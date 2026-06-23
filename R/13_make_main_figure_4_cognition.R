############################################################
## 13_make_main_figure_4_cognition.R
##
## Integrated Figure 4:
##   A. Network-level Hsp60/10 client score by cognitive diagnosis
##   B. Client-level cognition-prioritized proteins across outcomes
##   C. MitoCarta class breakdown by cognition support tier
##
## This replaces:
##   - old 13_main_fig4_cognition.R
##   - exploratory 14_client_level_cognition_prioritization.R
##   - exploratory 16_plot_client_level_cognition_MitoCarta_clean.R
##
## Run after:
##   source("00_config.R")
##   source("01_utils.R")
##   source("06_build_cognition_objects.R")
##
## Requires:
##   - cfg
##   - cognition_model_df
##   - prot_mat
##   - all_hsp60_10_client_tbl
##
## Outputs:
##   - Integrated Figure 4 PDF/PNG/SVG
##   - aggregate cognition panel source data
##   - network-level cognition forest table
##   - client-level cognition model table
##   - client-level cognition gene summary
##   - MitoCarta class breakdown table
##   - figure source data
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(scales)
  library(grid)
})

############################################################
## 0. Setup
############################################################

if (!exists("cfg")) {
  stop("cfg does not exist. Run source('00_config.R') first.", call. = FALSE)
}

require_objects(
  c("cfg", "cognition_model_df", "prot_mat", "all_hsp60_10_client_tbl"),
  context = "13_make_main_figure_4_cognition.R"
)

require_columns(
  cognition_model_df,
  c(
    "SampleID",
    "individualID",
    "Hsp60_client_score_z",
    "cogdx",
    "dcfdx_lv",
    "mmse_last_valid",
    "age_death",
    "sex",
    "educ",
    "pmi",
    "Braak",
    "CERAD"
  ),
  "Figure 4 integrated cognition_model_df"
)

out_dir <- file.path(cfg$plot_dir, "main_fig4_cognition")
tab_dir <- file.path(cfg$table_dir, "main_fig4_cognition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

## Same MitoCarta path used in Figure 6 / cognition-MitoCarta script.
mitocarta_path <- cfg$mitocarta_xls

if (!file.exists(mitocarta_path)) {
  stop(
    "MitoCarta file not found:\n",
    mitocarta_path,
    "\nCheck the path and filename.",
    call. = FALSE
  )
}

top_n_genes <- 25

############################################################
## 1. Helper functions
############################################################

clean_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

clean_col <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    tolower()
}

save_fig4_plot <- function(plot, name, width = 15.8, height = 10.2, dpi = 600) {
  pdf_out <- file.path(out_dir, paste0(name, ".pdf"))
  png_out <- file.path(out_dir, paste0(name, ".png"))
  svg_out <- file.path(out_dir, paste0(name, ".svg"))
  
  ggsave(
    pdf_out,
    plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf,
    bg = "white",
    limitsize = FALSE,
    useDingbats = FALSE
  )
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_out,
      plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = ragg::agg_png,
      bg = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      png_out,
      plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      svg_out,
      plot,
      width = width,
      height = height,
      units = "in",
      device = svglite::svglite,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  message("Saved PDF: ", pdf_out)
  message("Saved PNG: ", png_out)
  if (file.exists(svg_out)) message("Saved SVG: ", svg_out)
  
  invisible(plot)
}

write_fig4_table <- function(x, name) {
  out <- file.path(tab_dir, paste0(name, ".csv"))
  readr::write_csv(x, out)
  message("Wrote table: ", out)
  invisible(x)
}

theme_fig4 <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2.0, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.7, hjust = 0, color = "grey35"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 0.7),
      axis.line = element_line(linewidth = 0.42),
      axis.ticks = element_line(linewidth = 0.32),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(7, 9, 7, 7)
    )
}

############################################################
## 2. Prepare cognition dataframe
############################################################

plot_df <- cognition_model_df %>%
  mutate(
    SampleID = as.character(SampleID),
    individualID = as.character(individualID),
    
    cogdx_severity = case_when(
      cogdx == 1 ~ 0,
      cogdx %in% c(2, 3) ~ 1,
      cogdx %in% c(4, 5) ~ 2,
      TRUE ~ NA_real_
    ),
    
    dcfdx_lv_severity = case_when(
      dcfdx_lv == 1 ~ 0,
      dcfdx_lv %in% c(2, 3) ~ 1,
      dcfdx_lv %in% c(4, 5) ~ 2,
      TRUE ~ NA_real_
    ),
    
    ## Reoriented so higher = better cognition.
    cogdx_better = case_when(
      cogdx_severity == 0 ~ 2,
      cogdx_severity == 1 ~ 1,
      cogdx_severity == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    
    dcfdx_lv_better = case_when(
      dcfdx_lv_severity == 0 ~ 2,
      dcfdx_lv_severity == 1 ~ 1,
      dcfdx_lv_severity == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    
    cogdx_better_z = as.numeric(scale(cogdx_better)),
    dcfdx_lv_better_z = as.numeric(scale(dcfdx_lv_better)),
    mmse_last_valid_z = as.numeric(scale(mmse_last_valid)),
    
    cog_group = case_when(
      cogdx == 1 ~ "NCI",
      cogdx %in% c(2, 3) ~ "MCI",
      cogdx %in% c(4, 5) ~ "AD",
      TRUE ~ NA_character_
    ),
    cog_group = factor(cog_group, levels = c("NCI", "MCI", "AD")),
    sex = as.factor(sex)
  )

############################################################
## 3. Panel A: network-level client score by cognition group
############################################################

diag_df <- plot_df %>%
  filter(
    !is.na(cog_group),
    !is.na(Hsp60_client_score_z),
    !is.na(age_death),
    !is.na(sex),
    !is.na(educ),
    !is.na(pmi),
    !is.na(Braak),
    !is.na(CERAD)
  ) %>%
  mutate(
    cog_group_num = case_when(
      cog_group == "NCI" ~ 0,
      cog_group == "MCI" ~ 1,
      cog_group == "AD" ~ 2,
      TRUE ~ NA_real_
    )
  )

fit_hsp_adjust <- lm(
  Hsp60_client_score_z ~ age_death + sex + educ + pmi + Braak + CERAD,
  data = diag_df
)

diag_df <- diag_df %>%
  mutate(
    adjusted_hsp_score = resid(fit_hsp_adjust) + mean(Hsp60_client_score_z, na.rm = TRUE)
  )

fit_group_trend <- lm(
  Hsp60_client_score_z ~ cog_group_num + age_death + sex + educ + pmi + Braak + CERAD,
  data = diag_df
)

trend_stats <- broom::tidy(fit_group_trend, conf.int = TRUE) %>%
  filter(term == "cog_group_num")

trend_label <- paste0(
  "Adjusted trend = ",
  round(trend_stats$estimate, 2),
  " per step\n",
  ifelse(
    trend_stats$p.value < 0.001,
    "p < 0.001",
    paste0("p = ", scales::pvalue(trend_stats$p.value, accuracy = 0.001))
  )
)

group_stats <- diag_df %>%
  group_by(cog_group) %>%
  summarise(
    n = n(),
    mean_adj = mean(adjusted_hsp_score, na.rm = TRUE),
    sem_adj = sd(adjusted_hsp_score, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

y_min_A <- min(diag_df$adjusted_hsp_score, na.rm = TRUE)
y_max_A <- max(diag_df$adjusted_hsp_score, na.rm = TRUE)

pA <- ggplot(diag_df, aes(x = cog_group, y = adjusted_hsp_score)) +
  geom_violin(
    width = 0.78,
    fill = "grey93",
    color = NA,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.32,
    outlier.shape = NA,
    linewidth = 0.58,
    fill = "white"
  ) +
  geom_jitter(
    width = 0.095,
    alpha = 0.45,
    size = 1.15
  ) +
  geom_point(
    data = group_stats,
    aes(x = cog_group, y = mean_adj),
    inherit.aes = FALSE,
    size = 2.8,
    shape = 23,
    fill = "white",
    stroke = 0.9
  ) +
  geom_errorbar(
    data = group_stats,
    aes(
      x = cog_group,
      ymin = mean_adj - sem_adj,
      ymax = mean_adj + sem_adj
    ),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.62
  ) +
  geom_text(
    data = group_stats,
    aes(
      x = cog_group,
      y = y_min_A - 0.22,
      label = paste0("n = ", n)
    ),
    inherit.aes = FALSE,
    size = 2.8
  ) +
  annotate(
    "text",
    x = 1.88,
    y = y_max_A + 0.24,
    label = trend_label,
    size = 2.18,
    fontface = "bold",
    lineheight = 0.90
  ) +
  coord_cartesian(
    ylim = c(y_min_A - 0.42, y_max_A + 0.48),
    clip = "off"
  ) +
  labs(
    title = "A. Client network abundance declines with cognitive impairment",
    subtitle = "Adjusted Hsp60/10 score across final cognitive diagnosis",
    x = "Final cognitive diagnosis\nat death",
    y = "Adjusted Hsp60/10 client score"
  ) +
  theme_fig4(base_size = 10.5) +
  theme(
    axis.title.x = element_text(size = 9.5, face = "bold", lineheight = 0.92),
    axis.title.y = element_text(margin = margin(r = 0), face = "bold"),
    plot.margin = margin(10, 12, 10, 8)
  )

############################################################
## 4. Network-level forest table for source data only
############################################################

run_network_model <- function(outcome, label, data) {
  formula_use <- as.formula(
    paste0(
      outcome,
      " ~ Hsp60_client_score_z + age_death + sex + educ + pmi + Braak + CERAD"
    )
  )
  
  fit <- lm(formula_use, data = data)
  
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term == "Hsp60_client_score_z") %>%
    mutate(
      outcome = label,
      n = nobs(fit),
      r_squared = summary(fit)$r.squared,
      adj_r_squared = summary(fit)$adj.r.squared
    ) %>%
    select(outcome, n, estimate, std.error, statistic, p.value, conf.low, conf.high, r_squared, adj_r_squared)
}

network_cognition_tbl <- bind_rows(
  run_network_model("cogdx_better_z", "Final cognitive diagnosis", plot_df),
  run_network_model("dcfdx_lv_better_z", "Last-valid cognitive diagnosis", plot_df),
  run_network_model("mmse_last_valid_z", "Last-valid MMSE", plot_df)
)

############################################################
## 5. Client-level cognition models
############################################################

all_hsp_clients <- all_hsp60_10_client_tbl %>%
  transmute(gene = clean_gene(gene)) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene)

protein_mat <- prot_mat
rownames(protein_mat) <- clean_gene(rownames(protein_mat))
colnames(protein_mat) <- as.character(colnames(protein_mat))

detected_clients <- intersect(
  rownames(protein_mat),
  all_hsp_clients$gene
)

missing_clients <- setdiff(
  all_hsp_clients$gene,
  detected_clients
)

message("\nFigure 4 client-level cognition models")
message("Curated Hsp60/10 clients: ", nrow(all_hsp_clients))
message("Detected clients in prot_mat: ", length(detected_clients))
message("Missing clients: ", length(missing_clients))

if (length(detected_clients) < 5) {
  stop("Too few detected Hsp60/10 clients in prot_mat.", call. = FALSE)
}

client_mat <- protein_mat[detected_clients, , drop = FALSE]
client_mat <- apply(client_mat, 2, as.numeric)
rownames(client_mat) <- detected_clients

client_mat_z <- t(scale(t(client_mat)))

client_mat_z <- client_mat_z[
  rowSums(!is.na(client_mat_z)) > 0,
  ,
  drop = FALSE
]

client_long <- as.data.frame(client_mat_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "SampleID",
    values_to = "client_abundance_z"
  ) %>%
  mutate(
    gene = clean_gene(gene),
    SampleID = as.character(SampleID),
    client_abundance_z = as.numeric(client_abundance_z)
  )

cog_df <- plot_df %>%
  select(
    SampleID,
    individualID,
    cogdx,
    dcfdx_lv,
    mmse_last_valid,
    cogdx_better_z,
    dcfdx_lv_better_z,
    mmse_last_valid_z,
    age_death,
    sex,
    educ,
    pmi,
    Braak,
    CERAD
  )

analysis_long <- client_long %>%
  checked_left_join(
    cog_df,
    by = "SampleID",
    label = "client-level abundance to cognition metadata"
  )

run_client_model <- function(data, outcome, outcome_label) {
  
  model_df <- data %>%
    select(
      gene,
      client_abundance_z,
      all_of(outcome),
      age_death,
      sex,
      educ,
      pmi,
      Braak,
      CERAD
    ) %>%
    drop_na()
  
  message(
    "Running ", outcome_label,
    ": n rows = ", nrow(model_df),
    "; n genes = ", n_distinct(model_df$gene)
  )
  
  model_df %>%
    group_by(gene) %>%
    group_modify(~ {
      
      dat <- .x
      
      if (
        nrow(dat) < 40 ||
        sd(dat$client_abundance_z, na.rm = TRUE) == 0 ||
        sd(dat[[outcome]], na.rm = TRUE) == 0
      ) {
        return(
          tibble(
            outcome = outcome_label,
            n = nrow(dat),
            estimate = NA_real_,
            std.error = NA_real_,
            statistic = NA_real_,
            p.value = NA_real_,
            conf.low = NA_real_,
            conf.high = NA_real_,
            r_squared = NA_real_,
            adj_r_squared = NA_real_,
            model_status = "skipped_low_n_or_no_variance"
          )
        )
      }
      
      form <- as.formula(
        paste0(
          outcome,
          " ~ client_abundance_z + age_death + sex + educ + pmi + Braak + CERAD"
        )
      )
      
      fit <- tryCatch(
        lm(form, data = dat),
        error = function(e) e
      )
      
      if (inherits(fit, "error")) {
        return(
          tibble(
            outcome = outcome_label,
            n = nrow(dat),
            estimate = NA_real_,
            std.error = NA_real_,
            statistic = NA_real_,
            p.value = NA_real_,
            conf.low = NA_real_,
            conf.high = NA_real_,
            r_squared = NA_real_,
            adj_r_squared = NA_real_,
            model_status = paste0("model_error: ", fit$message)
          )
        )
      }
      
      tidy_fit <- broom::tidy(fit, conf.int = TRUE) %>%
        filter(term == "client_abundance_z")
      
      if (nrow(tidy_fit) != 1) {
        return(
          tibble(
            outcome = outcome_label,
            n = nrow(dat),
            estimate = NA_real_,
            std.error = NA_real_,
            statistic = NA_real_,
            p.value = NA_real_,
            conf.low = NA_real_,
            conf.high = NA_real_,
            r_squared = summary(fit)$r.squared,
            adj_r_squared = summary(fit)$adj.r.squared,
            model_status = "missing_client_term"
          )
        )
      }
      
      tidy_fit %>%
        transmute(
          outcome = outcome_label,
          n = nobs(fit),
          estimate,
          std.error,
          statistic,
          p.value,
          conf.low,
          conf.high,
          r_squared = summary(fit)$r.squared,
          adj_r_squared = summary(fit)$adj.r.squared,
          model_status = "ok"
        )
    }) %>%
    ungroup()
}

client_cognition_all <- bind_rows(
  run_client_model(
    analysis_long,
    outcome = "mmse_last_valid_z",
    outcome_label = "Last-valid MMSE"
  ),
  run_client_model(
    analysis_long,
    outcome = "dcfdx_lv_better_z",
    outcome_label = "Last-valid cognitive diagnosis"
  ),
  run_client_model(
    analysis_long,
    outcome = "cogdx_better_z",
    outcome_label = "Final cognitive diagnosis"
  )
) %>%
  group_by(outcome) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      estimate > 0 ~ "higher abundance = better cognition",
      estimate < 0 ~ "higher abundance = worse cognition",
      TRUE ~ NA_character_
    ),
    cognition_associated_nominal = model_status == "ok" & p.value < 0.05 & estimate > 0,
    cognition_associated_fdr = model_status == "ok" & q.value < 0.05 & estimate > 0
  ) %>%
  arrange(outcome, q.value, desc(estimate))

client_cognition_summary <- client_cognition_all %>%
  filter(model_status == "ok") %>%
  group_by(gene) %>%
  summarise(
    n_outcomes_tested = n(),
    
    n_nominal_better_cognition = sum(cognition_associated_nominal, na.rm = TRUE),
    n_fdr_better_cognition = sum(cognition_associated_fdr, na.rm = TRUE),
    
    best_p = min(p.value, na.rm = TRUE),
    best_q = min(q.value, na.rm = TRUE),
    
    mean_beta = mean(estimate, na.rm = TRUE),
    median_beta = median(estimate, na.rm = TRUE),
    
    mmse_beta = estimate[outcome == "Last-valid MMSE"][1],
    mmse_p = p.value[outcome == "Last-valid MMSE"][1],
    mmse_q = q.value[outcome == "Last-valid MMSE"][1],
    
    last_dx_beta = estimate[outcome == "Last-valid cognitive diagnosis"][1],
    last_dx_p = p.value[outcome == "Last-valid cognitive diagnosis"][1],
    last_dx_q = q.value[outcome == "Last-valid cognitive diagnosis"][1],
    
    final_dx_beta = estimate[outcome == "Final cognitive diagnosis"][1],
    final_dx_p = p.value[outcome == "Final cognitive diagnosis"][1],
    final_dx_q = q.value[outcome == "Final cognitive diagnosis"][1],
    
    .groups = "drop"
  ) %>%
  mutate(
    cognition_priority_score =
      n_fdr_better_cognition * 3 +
      n_nominal_better_cognition +
      if_else(!is.na(mmse_beta) & mmse_beta > 0, 0.5, 0),
    cognition_priority_tier = case_when(
      n_fdr_better_cognition >= 2 ~ "FDR-supported across multiple cognition outcomes",
      n_fdr_better_cognition == 1 ~ "FDR-supported in one cognition outcome",
      n_nominal_better_cognition >= 2 ~ "Nominal across multiple cognition outcomes",
      n_nominal_better_cognition == 1 ~ "Nominal in one cognition outcome",
      TRUE ~ "Not cognition-prioritized"
    )
  ) %>%
  arrange(
    desc(cognition_priority_score),
    best_q,
    desc(mean_beta)
  )

############################################################
## 6. Read MitoCarta and assign functional classes
############################################################

sheets <- readxl::excel_sheets(mitocarta_path)

read_mitocarta_sheet <- function(sh) {
  raw <- readxl::read_excel(mitocarta_path, sheet = sh)
  names(raw) <- clean_col(names(raw))
  raw
}

candidate_sheets <- purrr::map(sheets, function(sh) {
  tmp <- tryCatch(read_mitocarta_sheet(sh), error = function(e) NULL)
  if (is.null(tmp)) return(NULL)
  
  cols <- colnames(tmp)
  
  gene_col <- intersect(
    c("symbol", "gene_symbol", "genesymbol", "hgnc_symbol", "gene", "gene_name"),
    cols
  )[1]
  
  pathway_col <- intersect(
    c(
      "mitopathways",
      "mito_pathways",
      "mitocarta3_0_mitopathways",
      "mitocarta3_mitopathways",
      "mitocarta_pathways",
      "mitocarta_pathway",
      "pathways",
      "pathway",
      "category",
      "categories"
    ),
    cols
  )[1]
  
  list(
    sheet = sh,
    table = tmp,
    gene_col = gene_col,
    pathway_col = pathway_col,
    n_rows = nrow(tmp),
    n_cols = ncol(tmp)
  )
})

candidate_sheets <- candidate_sheets[!purrr::map_lgl(candidate_sheets, is.null)]

sheet_choice_tbl <- purrr::map_dfr(candidate_sheets, function(x) {
  tibble(
    sheet = x$sheet,
    n_rows = x$n_rows,
    n_cols = x$n_cols,
    gene_col = ifelse(is.na(x$gene_col), NA_character_, x$gene_col),
    pathway_col = ifelse(is.na(x$pathway_col), NA_character_, x$pathway_col),
    score = as.integer(!is.na(x$gene_col)) * 10 +
      as.integer(!is.na(x$pathway_col)) * 5 +
      log10(max(x$n_rows, 1))
  )
}) %>%
  arrange(desc(score), desc(n_rows))

chosen_sheet <- sheet_choice_tbl$sheet[1]
chosen <- candidate_sheets[[which(purrr::map_chr(candidate_sheets, "sheet") == chosen_sheet)[1]]]

if (is.na(chosen$gene_col)) {
  stop("Could not identify a gene symbol column in MitoCarta file.", call. = FALSE)
}

mito_raw <- chosen$table
gene_col <- chosen$gene_col
pathway_col <- chosen$pathway_col

desc_col <- intersect(
  c("description", "protein_description", "name", "gene_description"),
  colnames(mito_raw)
)[1]

sub_loc_col <- intersect(
  c(
    "sub_mito_localization",
    "submito_localization",
    "mitocarta3_0_submito_localization",
    "submitochondrial_location",
    "localization"
  ),
  colnames(mito_raw)
)[1]

mitocarta_tbl <- mito_raw %>%
  transmute(
    gene = clean_gene(.data[[gene_col]]),
    mitocarta_sheet = chosen_sheet,
    mitocarta_pathways_raw = if (!is.na(pathway_col)) as.character(.data[[pathway_col]]) else NA_character_,
    mitocarta_description = if (!is.na(desc_col)) as.character(.data[[desc_col]]) else NA_character_,
    mitocarta_submito = if (!is.na(sub_loc_col)) as.character(.data[[sub_loc_col]]) else NA_character_
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene, .keep_all = TRUE)

classify_mitocarta <- function(gene, pathways, description = NA_character_) {
  g <- clean_gene(gene)
  x <- paste(pathways, description, sep = " | ")
  x <- stringr::str_to_lower(x)
  x[is.na(x)] <- ""
  
  dplyr::case_when(
    stringr::str_detect(x, "ribosom|translation|mitoribosom|mitochondrial central dogma|mrna|trna|rrna|rna metabolism|transcription") |
      stringr::str_detect(g, "^MRPL|^MRPS") |
      g %in% c("TUFM", "TSFM", "GFM1", "GFM2", "DAP3", "LRPPRC", "SLIRP", "PTCD3") ~
      "Mito translation / ribosome",
    
    stringr::str_detect(x, "tca|tricarboxylic|pyruvate|carbohydrate|ketoglutarate|dehydrogenase|nad|redox|fe-s|heme|iron-sulfur|organic acid") |
      g %in% c(
        "PDHA1", "PDHB", "PDHX", "DLAT", "DLD", "DLST", "OGDH",
        "IDH3A", "IDH3B", "IDH3G", "IDH2", "MDH2", "CS", "FH", "ACO2",
        "SUCLA2", "SUCLG1", "SUCLG2", "GOT2", "SHMT2", "GLUD1", "GLUD2",
        "OAT", "ALDH2", "ME2"
      ) ~
      "TCA / pyruvate / redox metabolism",
    
    stringr::str_detect(x, "oxidative phosphorylation|oxphos|respiratory chain|electron transport|complex i|complex ii|complex iii|complex iv|complex v|atp synthase") |
      stringr::str_detect(g, "^NDUF|^UQCR|^COX|^SDH|^ATP5|^MT-") ~
      "OXPHOS / ETC",
    
    stringr::str_detect(x, "fatty acid|lipid|beta.oxidation|β.oxidation|acyl|carnitine") |
      g %in% c("ECHS1", "ECH1", "HADHA", "HADHB", "ACAA2", "ACADVL", "ETFA", "ETFB", "ETFDH") ~
      "FAO / lipid metabolism",
    
    stringr::str_detect(x, "protein import|protein sorting|protein homeostasis|proteostasis|protein folding|chaperone|protease|quality control|stress|detoxification|antioxidant|ros") |
      g %in% c("HSPA9", "TRAP1", "CLPX", "CLPP", "LONP1", "AFG3L2", "SPG7", "PHB", "PHB2", "PRDX3", "SOD2", "TXN2") ~
      "Proteostasis / stress",
    
    stringr::str_detect(x, "transport|carrier|membrane|import|export|transloc|dynamics|fission|fusion|signaling|apoptosis|calcium") ~
      "Transport / membrane / signaling",
    
    stringr::str_detect(x, "metabolism|amino acid|nucleotide|one-carbon|folate|cofactor|vitamin|porphyrin|ubiquinone|coq") ~
      "Other mitochondrial metabolism",
    
    x != "" ~ "Other MitoCarta-annotated clients",
    TRUE ~ "No MitoCarta annotation"
  )
}

functional_levels <- c(
  "Mito translation / ribosome",
  "TCA / pyruvate / redox metabolism",
  "OXPHOS / ETC",
  "FAO / lipid metabolism",
  "Proteostasis / stress",
  "Transport / membrane / signaling",
  "Other mitochondrial metabolism",
  "Other MitoCarta-annotated clients",
  "No MitoCarta annotation"
)

functional_pal <- c(
  "Mito translation / ribosome" = "#7A2E55",
  "TCA / pyruvate / redox metabolism" = "#C8962E",
  "OXPHOS / ETC" = "#3B6EA8",
  "FAO / lipid metabolism" = "#7A8F2A",
  "Proteostasis / stress" = "#4F8A5B",
  "Transport / membrane / signaling" = "#7B6EA8",
  "Other mitochondrial metabolism" = "#8E6C8A",
  "Other MitoCarta-annotated clients" = "grey58",
  "No MitoCarta annotation" = "grey82"
)

functional_legend_labels <- c(
  "Mito translation / ribosome" = "Translation / ribosome",
  "TCA / pyruvate / redox metabolism" = "TCA / pyruvate / redox",
  "OXPHOS / ETC" = "OXPHOS / ETC",
  "FAO / lipid metabolism" = "FAO / lipid",
  "Proteostasis / stress" = "Proteostasis / stress",
  "Transport / membrane / signaling" = "Transport / membrane / signaling",
  "Other mitochondrial metabolism" = "Other mito metabolism",
  "Other MitoCarta-annotated clients" = "Other MitoCarta",
  "No MitoCarta annotation" = "No annotation"
)

functional_display_map <- c(
  "Mito translation / ribosome" = "Translation /\nribosome",
  "TCA / pyruvate / redox metabolism" = "TCA /\nredox",
  "OXPHOS / ETC" = "OXPHOS /\nETC",
  "FAO / lipid metabolism" = "FAO /\nlipid",
  "Proteostasis / stress" = "Proteostasis",
  "Transport / membrane / signaling" = "Transport /\nsignaling",
  "Other mitochondrial metabolism" = "Other mito",
  "Other MitoCarta-annotated clients" = "Other\nMitoCarta",
  "No MitoCarta annotation" = "No\nannotation"
)

functional_display_levels <- unname(functional_display_map[functional_levels])

class_tbl <- mitocarta_tbl %>%
  mutate(
    functional_class = classify_mitocarta(
      gene,
      mitocarta_pathways_raw,
      mitocarta_description
    ),
    functional_class = factor(functional_class, levels = functional_levels),
    functional_class_display = recode(
      as.character(functional_class),
      !!!functional_display_map
    ),
    functional_class_display = factor(
      functional_class_display,
      levels = functional_display_levels
    )
  )

############################################################
## 7. Join cognition summary to MitoCarta class annotations
############################################################

rank_tbl_all_unfiltered <- client_cognition_summary %>%
  left_join(
    class_tbl %>%
      select(
        gene,
        functional_class,
        functional_class_display,
        mitocarta_pathways_raw,
        mitocarta_description,
        mitocarta_submito
      ),
    by = "gene"
  ) %>%
  mutate(
    functional_class = if_else(
      is.na(as.character(functional_class)),
      "No MitoCarta annotation",
      as.character(functional_class)
    ),
    functional_class = factor(functional_class, levels = functional_levels),
    
    functional_class_display = if_else(
      is.na(as.character(functional_class_display)),
      functional_display_map["No MitoCarta annotation"],
      as.character(functional_class_display)
    ),
    functional_class_display = factor(
      functional_class_display,
      levels = functional_display_levels
    ),
    
    n_fdr_better_cognition = replace_na(n_fdr_better_cognition, 0),
    n_nominal_better_cognition = replace_na(n_nominal_better_cognition, 0),
    cognition_priority_score = replace_na(cognition_priority_score, 0),
    mean_beta = replace_na(mean_beta, -Inf),
    best_q = replace_na(best_q, Inf)
  )

rank_tbl_top_pool <- rank_tbl_all_unfiltered %>%
  filter(n_fdr_better_cognition >= 2)

if (nrow(rank_tbl_top_pool) == 0) {
  stop("No genes have n_fdr_better_cognition >= 2.", call. = FALSE)
}

top_genes <- rank_tbl_top_pool %>%
  arrange(
    desc(n_fdr_better_cognition),
    desc(cognition_priority_score),
    desc(mean_beta),
    best_q
  ) %>%
  slice_head(n = top_n_genes) %>%
  pull(gene)

display_tbl <- rank_tbl_top_pool %>%
  filter(gene %in% top_genes) %>%
  arrange(
    functional_class,
    desc(n_fdr_better_cognition),
    desc(cognition_priority_score),
    desc(mean_beta),
    best_q
  )

gene_order <- display_tbl$gene

############################################################
## 8. Panel B: client-level cognition dot plot
############################################################

outcome_map <- c(
  "Last-valid MMSE" = "MMSE",
  "Last-valid cognitive diagnosis" = "Last dx",
  "Final cognitive diagnosis" = "Final dx"
)

outcome_display_levels <- unname(outcome_map)
consistency_col_label <- "FDR outcomes"
plot_col_levels <- c(outcome_display_levels, consistency_col_label)

client_dot_df <- client_cognition_all %>%
  filter(
    gene %in% top_genes,
    model_status == "ok",
    outcome %in% names(outcome_map)
  ) %>%
  left_join(
    display_tbl %>%
      select(
        gene,
        functional_class,
        functional_class_display,
        n_fdr_better_cognition,
        n_nominal_better_cognition,
        cognition_priority_score,
        mean_beta,
        best_q
      ),
    by = "gene"
  ) %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    functional_class = factor(functional_class, levels = functional_levels),
    functional_class_display = factor(functional_class_display, levels = functional_display_levels),
    outcome_display = recode(outcome, !!!outcome_map),
    plot_col = factor(outcome_display, levels = plot_col_levels),
    pt_alpha = case_when(
      q.value < 0.05 & estimate > 0 ~ 0.95,
      p.value < 0.05 & estimate > 0 ~ 0.65,
      TRUE ~ 0.25
    ),
    neglog10_q = -log10(q.value),
    neglog10_q = if_else(is.infinite(neglog10_q), NA_real_, neglog10_q),
    neglog10_q = pmin(neglog10_q, 8)
  )

consistency_df <- display_tbl %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    functional_class = factor(functional_class, levels = functional_levels),
    functional_class_display = factor(functional_class_display, levels = functional_display_levels),
    consistency_col = factor(consistency_col_label, levels = plot_col_levels),
    consistency_label = as.character(n_fdr_better_cognition)
  )

legend_key_df <- tibble(
  plot_col = factor(plot_col_levels[1], levels = plot_col_levels),
  gene = factor(gene_order[1], levels = rev(gene_order)),
  functional_class = factor(functional_levels, levels = functional_levels),
  functional_class_display = factor(functional_display_levels[1], levels = functional_display_levels)
)

pB <- ggplot() +
  geom_tile(
    data = consistency_df,
    aes(x = consistency_col, y = gene),
    width = 0.82,
    height = 0.78,
    fill = "grey96",
    color = "grey85",
    linewidth = 0.35
  ) +
  geom_text(
    data = consistency_df,
    aes(x = consistency_col, y = gene, label = consistency_label),
    size = 2.75,
    fontface = "bold"
  ) +
  geom_point(
    data = legend_key_df,
    aes(x = plot_col, y = gene, fill = functional_class),
    inherit.aes = FALSE,
    shape = 22,
    size = 2.8,
    alpha = 0,
    color = NA,
    na.rm = TRUE,
    show.legend = TRUE
  ) +
  geom_point(
    data = client_dot_df,
    aes(
      x = plot_col,
      y = gene,
      size = neglog10_q,
      fill = functional_class,
      alpha = pt_alpha
    ),
    shape = 21,
    color = "black",
    stroke = 0.32
  ) +
  facet_grid(
    functional_class_display ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_alpha_identity() +
  scale_fill_manual(
    values = functional_pal,
    breaks = rev(functional_levels),
    drop = FALSE,
    name = "MitoCarta class",
    labels = function(x) stringr::str_wrap(unname(functional_legend_labels[x]), width = 22)
  ) +
  scale_size_continuous(
    range = c(1.8, 6.3),
    breaks = c(1.3, 2, 3),
    labels = c("0.05", "0.01", "0.001"),
    name = expression(-log[10](FDR))
  ) +
  scale_x_discrete(
    limits = plot_col_levels,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "B. Cognition-prioritized clients show consistent effects across outcomes",
    subtitle = paste0(
      "Top ", top_n_genes,
      " clients; dot size = -log10(FDR); right column = number of FDR-supported outcomes"
    ),
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_fig4(base_size = 9.7) +
  theme(
    axis.text.x = element_text(size = 7.4, angle = 18, hjust = 1, vjust = 1, face = "bold", lineheight = 0.92),
    axis.text.y = element_text(size = 8.0),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey96", color = NA),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 7.0,
      lineheight = 0.84,
      margin = margin(r = 0.6, l = 0.6)
    ),
    panel.spacing.y = unit(0.24, "lines"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 1, r = 1, b = 1, l = 1),
    legend.spacing.y = unit(1, "pt"),
    legend.spacing.x = unit(2, "pt"),
    legend.title = element_text(face = "bold", size = 6.8),
    legend.text = element_text(size = 5.3, lineheight = 0.82),
    legend.key.size = unit(0.17, "cm"),
    plot.margin = margin(7, 8, 7, 5)
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      ncol = 3,
      byrow = TRUE,
      override.aes = list(shape = 22, size = 3.0, alpha = 1, color = "grey35", linewidth = 0.2)
    ),
    size = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(fill = "grey70", alpha = 1)
    )
  )

############################################################
## 9. Panel C: functional-class breakdown by FDR tier
## Uses unfiltered all detected client genes.
############################################################

functional_class_short_map <- c(
  "Mito translation / ribosome" = "Translation /\nribosome",
  "TCA / pyruvate / redox metabolism" = "TCA / pyruvate /\nredox",
  "OXPHOS / ETC" = "OXPHOS /\nETC",
  "FAO / lipid metabolism" = "FAO /\nlipid",
  "Proteostasis / stress" = "Proteostasis /\nstress",
  "Transport / membrane / signaling" = "Transport /\nmembrane /\nsignaling",
  "Other mitochondrial metabolism" = "Other mito\nmetabolism",
  "Other MitoCarta-annotated clients" = "Other\nMitoCarta",
  "No MitoCarta annotation" = "No MitoCarta\nannotation"
)

functional_class_short_levels <- unname(functional_class_short_map[functional_levels])

class_breakdown_tbl <- rank_tbl_all_unfiltered %>%
  mutate(
    cognition_support_tier = case_when(
      n_fdr_better_cognition == 3 ~ "FDR in 3/3",
      n_fdr_better_cognition == 2 ~ "FDR in 2/3",
      n_fdr_better_cognition == 1 ~ "FDR in 1/3",
      TRUE ~ "Not FDR-supported"
    ),
    cognition_support_tier = factor(
      cognition_support_tier,
      levels = c("FDR in 3/3", "FDR in 2/3", "FDR in 1/3", "Not FDR-supported")
    ),
    functional_class = factor(functional_class, levels = functional_levels),
    functional_class_short = recode(as.character(functional_class), !!!functional_class_short_map),
    functional_class_short = factor(
      functional_class_short,
      levels = rev(functional_class_short_levels)
    )
  ) %>%
  count(cognition_support_tier, functional_class, functional_class_short, name = "n") %>%
  group_by(cognition_support_tier) %>%
  mutate(
    tier_total = sum(n),
    fraction_within_tier = n / tier_total
  ) %>%
  ungroup()

pC <- ggplot(
  class_breakdown_tbl,
  aes(x = cognition_support_tier, y = functional_class_short)
) +
  geom_point(
    aes(size = n, color = functional_class),
    alpha = 0.88
  ) +
  geom_text(
    aes(label = n),
    size = 2.35,
    fontface = "bold",
    color = "white"
  ) +
  scale_color_manual(
    values = functional_pal,
    drop = FALSE,
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(2.2, 7.2),
    guide = "none"
  ) +
  labs(
    title = "C. Cognition-supported client classes",
    subtitle = "Counts show all detected clients by MitoCarta class \nand number of FDR-supported cognition outcomes",
    x = NULL,
    y = NULL
  ) +
  theme_fig4(base_size = 9.6) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 6.4,
      lineheight = 0.88
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 6.25,
      lineheight = 0.84
    ),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25),
    plot.margin = margin(7, 9, 10, 7)
  )

############################################################
## 10. Assemble integrated Figure 4
############################################################

left_col <- pA / pC +
  plot_layout(heights = c(1.08, 0.92))

fig4_integrated <- left_col | pB

fig4_integrated <- fig4_integrated +
  plot_layout(widths = c(1.34, 1.18), guides = "collect") +
  plot_annotation(
    title = "Hsp60/10 Client Network Abundance Tracks Cognitive Preservation",
    subtitle = "Network-level decline is accompanied by cognition-prioritized mitochondrial client classes in ROSMAP",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17.0, hjust = 0.5),
      plot.subtitle = element_text(size = 10.6, hjust = 0.5, color = "grey35"),
      plot.margin = margin(8, 10, 8, 10),
      legend.position = "bottom"
    )
  ) &
  theme(
    legend.position = "bottom"
  )

fig4_integrated

############################################################
## 11. Save outputs
############################################################

save_fig4_plot(
  fig4_integrated,
  "Figure4_Hsp60_10_cognition",
  width = 16.2,
  height = 10.4,
  dpi = 600
)

write_fig4_table(diag_df, "Fig4A_network_diagnosis_panel_dataframe")
write_fig4_table(group_stats, "Fig4A_network_diagnosis_group_summary")
write_fig4_table(network_cognition_tbl, "Fig4_network_level_cognition_models_source_data")
write_fig4_table(client_cognition_all, "Fig4_client_level_cognition_all_models")
write_fig4_table(client_cognition_summary, "Fig4_client_level_cognition_summary_by_gene")
write_fig4_table(rank_tbl_all_unfiltered, "Fig4_client_level_cognition_summary_with_MitoCarta_classes")
write_fig4_table(client_dot_df, "Fig4B_client_level_dotplot_source_data")
write_fig4_table(class_breakdown_tbl, "Fig4C_MitoCarta_class_breakdown_by_cognition_FDR_tier")
write_fig4_table(sheet_choice_tbl, "Fig4_MitoCarta_sheet_choice_audit")
write_fig4_table(mitocarta_tbl, "Fig4_MitoCarta_join_table_cleaned")
write_fig4_table(
  tibble(
    n_curated_hsp60_10_clients = nrow(all_hsp_clients),
    n_detected_clients_in_prot_mat = length(detected_clients),
    n_missing_clients = length(missing_clients),
    n_top_genes_displayed_panel_B = length(top_genes),
    mitocarta_sheet_used = chosen_sheet,
    mitocarta_gene_col_used = gene_col,
    mitocarta_pathway_col_used = ifelse(is.na(pathway_col), "NONE", pathway_col)
  ),
  "Fig4_integrated_audit"
)

message("\n============================================================")
message("Integrated Figure 4 complete.")
message("Curated Hsp60/10 clients: ", nrow(all_hsp_clients))
message("Detected clients in prot_mat: ", length(detected_clients))
message("Displayed genes in Panel B: ", length(top_genes))
message("MitoCarta sheet used: ", chosen_sheet)
message("Outputs written to:")
message(out_dir)
message(tab_dir)
message("============================================================")
