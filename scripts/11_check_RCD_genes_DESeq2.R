# ============================================================
# RCD GENE ANALYSIS ACROSS RETAINED DESeq2 RESULTS
# Marchantia polymorpha
#
# GENE SET:
#   metadata/RCD_Marchantia_genes.tsv
#
# Main DEG definition:
#   padj < 0.05 & |log2FoldChange| >= 1
#
# IMPORTANT:
#   Input files are the FULL DESeq2/apeglm result files.
#
# ============================================================


# ============================================================
# LOAD PACKAGES
# ============================================================

library(dplyr)
library(tidyr)


# ============================================================
# PROJECT PATH
# ============================================================

PROJECT_DIR <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD"

setwd(PROJECT_DIR)


# ============================================================
# RCD GENE LIST
#
# Marchantia-only RCD gene set
# ============================================================

GENE_LIST <- file.path(
  PROJECT_DIR,
  "metadata/RCD_Marchantia_genes.tsv"
)


# ============================================================
# OUTPUT DIRECTORY
#
# Separate from the Arabidopsis-orthology gene set results
# ============================================================

OUT_DIR <- file.path(
  PROJECT_DIR,
  "results/Marchantia_polymorpha/DESeq2/RCD_gene_check_Marchantia"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# DESeq2 RESULT FILES
#
# These are the 9 retained comparisons.
#
# IMPORTANT:
# Use the FULL DESeq2/apeglm result files.
# ============================================================

files <- c(
  
  # ----------------------------------------------------------
  # Tan et al. 2023 — D2
  # ----------------------------------------------------------
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Tan2023_D2/cold_vs_control_apeglm.csv"
  ),
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Tan2023_D2/darkness_vs_control_apeglm.csv"
  ),
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Tan2023_D2/heat_vs_control_apeglm.csv"
  ),
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Tan2023_D2/high_light_vs_control_apeglm.csv"
  ),
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Tan2023_D2/nitrogen_deficiency_vs_control_apeglm.csv"
  ),
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Tan2023_D2/osmotic_vs_control_apeglm.csv"
  ),
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Tan2023_D2/salt_vs_control_apeglm.csv"
  ),
  
  
  # ----------------------------------------------------------
  # Schroder et al. 2023
  # ----------------------------------------------------------
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Schroder2023/high_light_vs_control_apeglm.csv"
  ),
  
  
  # ----------------------------------------------------------
  # Grenz et al. 2025
  # ----------------------------------------------------------
  
  file.path(
    PROJECT_DIR,
    "results/Marchantia_polymorpha/DESeq2/Grenz2025/infected_vs_control_apeglm.csv"
  )
  
)


# ============================================================
# CHECK FILES
# ============================================================

missing_files <- files[
  !file.exists(files)
]

if (length(missing_files) > 0) {
  
  cat("\nERROR: The following DESeq2 files are missing:\n\n")
  
  print(missing_files)
  
  stop(
    "\nMissing DESeq2 result files. Analysis stopped."
  )
  
}


cat("\n")
cat("============================================================\n")
cat("ALL DESeq2 RESULT FILES FOUND\n")
cat("============================================================\n")

cat(
  "\nNumber of DESeq2 files:",
  length(files),
  "\n\n"
)


# ============================================================
# CHECK GENE LIST FILE
# ============================================================

if (!file.exists(GENE_LIST)) {
  
  stop(
    paste(
      "\nERROR: PCD gene list not found:",
      GENE_LIST
    )
  )
  
}


# ============================================================
# LOAD PCD GENE LIST
# ============================================================

rcd <- read.delim(
  GENE_LIST,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# CHECK RCD GENE LIST
# ============================================================

if (!"gene_id" %in% colnames(rcd)) {
  
  stop(
    "ERROR: RCD gene list must contain a column called 'gene_id'."
  )
  
}


# ============================================================
# CLEAN RCD GENE IDs
# ============================================================

rcd_genes <- unique(
  trimws(rcd$gene_id)
)

rcd_genes <- rcd_genes[
  !is.na(rcd_genes) &
    rcd_genes != ""
]


cat(
  "Number of unique RCD genes in input list: ",
  length(rcd_genes),
  "\n",
  sep = ""
)


# ============================================================
# CHECK FOR DUPLICATE IDs IN ORIGINAL FILE
# ============================================================

duplicate_rcd_genes <- rcd$gene_id[
  duplicated(trimws(rcd$gene_id))
]

duplicate_rcd_genes <- unique(
  duplicate_rcd_genes[
    !is.na(duplicate_rcd_genes) &
      trimws(duplicate_rcd_genes) != ""
  ]
)

if (length(duplicate_rcd_genes) > 0) {
  
  cat("\nWARNING:\n")
  cat(
    "Duplicate gene IDs were present in the input gene list.\n"
  )
  
  cat(
    "Duplicates:\n",
    paste(
      duplicate_rcd_genes,
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "Duplicates were collapsed to unique gene IDs for analysis.\n"
  )
  
}


# ============================================================
# PROCESS EACH DESeq2 FILE
# ============================================================

results_list <- list()


for (f in files) {
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("Processing:\n")
  cat(basename(f), "\n")
  cat("------------------------------------------------------------\n")
  
  
  # ==========================================================
  # READ FULL DESeq2 RESULT
  # ==========================================================
  
  df <- read.csv(
    f,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  # ==========================================================
  # FIRST COLUMN = GENE ID
  # ==========================================================
  
  colnames(df)[1] <- "gene_id"
  
  df$gene_id <- trimws(
    df$gene_id
  )
  
  
  # ==========================================================
  # CHECK REQUIRED COLUMNS
  # ==========================================================
  
  required_columns <- c(
    "log2FoldChange",
    "padj"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste(
        "\nMissing required columns in:",
        f,
        "\nMissing:",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
    
  }
  
  
  # ==========================================================
  # STUDY NAME
  # ==========================================================
  
  study <- basename(
    dirname(f)
  )
  
  
  # ==========================================================
  # COMPARISON NAME
  # ==========================================================
  
  comparison <- tools::file_path_sans_ext(
    basename(f)
  )
  
  
  # ==========================================================
  # COMPARISON ID
  # ==========================================================
  
  comparison_id <- paste(
    study,
    comparison,
    sep = "__"
  )
  
  
  # ==========================================================
  # EXTRACT PCD GENES
  # ==========================================================
  
  matched <- df %>%
    
    filter(
      gene_id %in% rcd_genes
    ) %>%
    
    mutate(
      
      study = study,
      
      comparison = comparison,
      
      comparison_id = comparison_id
      
    )
  
  
  # ==========================================================
  # CHECK NUMBER OF PCD GENES FOUND
  # ==========================================================
  
  cat(
    "RCD genes found:",
    nrow(matched),
    "\n"
  )
  
  
  # ==========================================================
  # CHECK FOR DUPLICATE RCD GENE IDs
  # ==========================================================
  
  duplicated_genes <- matched$gene_id[
    duplicated(matched$gene_id)
  ]
  
  if (length(duplicated_genes) > 0) {
    
    warning(
      paste(
        "Duplicate RCD gene IDs detected in:",
        basename(f),
        "\nGenes:",
        paste(
          unique(duplicated_genes),
          collapse = ", "
        )
      )
    )
    
  }
  
  
  # ==========================================================
  # ADD DEG CLASSIFICATION
  # ==========================================================
  
  matched <- matched %>%
    
    mutate(
      
      # ------------------------------------------------------
      # Valid statistics
      # ------------------------------------------------------
      
      valid_stats =
        !is.na(padj) &
        !is.na(log2FoldChange),
      
      
      # ------------------------------------------------------
      # MAIN DEG DEFINITION
      #
      # padj < 0.05
      # |log2FC| >= 1
      # ------------------------------------------------------
      
      DEG =
        valid_stats &
        padj < 0.05 &
        abs(log2FoldChange) >= 1,
      
      
      # ------------------------------------------------------
      # Direction
      # ------------------------------------------------------
      
      direction = case_when(
        
        !valid_stats ~ "NA",
        
        !DEG ~ "Not_DEG",
        
        log2FoldChange > 0 ~ "Up",
        
        log2FoldChange < 0 ~ "Down",
        
        TRUE ~ "Not_DEG"
        
      )
      
    )
  
  
  # ==========================================================
  # PRINT RCD DEG COUNTS
  # ==========================================================
  
  cat(
    "RCD genes with valid statistics:",
    sum(matched$valid_stats),
    "\n"
  )
  
  cat(
    "RCD DEGs:",
    sum(matched$DEG),
    "\n"
  )
  
  cat(
    "  Up:",
    sum(
      matched$DEG &
        matched$direction == "Up"
    ),
    "\n"
  )
  
  cat(
    "  Down:",
    sum(
      matched$DEG &
        matched$direction == "Down"
    ),
    "\n"
  )
  
  
  # ==========================================================
  # STORE RESULT
  # ==========================================================
  
  results_list[[paste(study, comparison, sep = "__")]] <- matched
  
}


# ============================================================
# COMBINE RESULTS
# ============================================================

pcd_deseq2 <- bind_rows(
  results_list
)


# ============================================================
# CHECK COMBINED RESULT
# ============================================================

if (nrow(pcd_deseq2) == 0) {
  
  stop(
    "\nERROR: No PCD genes were found in any DESeq2 result."
  )
  
}


cat("\n")
cat("============================================================\n")
cat("COMBINED PCD RESULTS\n")
cat("============================================================\n")

cat(
  "Total PCD gene-DESeq2 matches:",
  nrow(pcd_deseq2),
  "\n"
)


# ============================================================
# REORDER IMPORTANT COLUMNS
# ============================================================

pcd_deseq2 <- pcd_deseq2 %>%
  
  select(
    gene_id,
    study,
    comparison,
    comparison_id,
    log2FoldChange,
    padj,
    valid_stats,
    DEG,
    direction,
    everything()
  )


# ============================================================
# SAVE COMPLETE PCD DESeq2 RESULTS
#
# Contains actual DESeq2 statistics for every PCD gene
# found in every comparison.
# ============================================================

write.table(
  pcd_deseq2,
  file.path(
    OUT_DIR,
    "PCD_Marchantia_genes_all_DESeq2_results.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# CREATE COMPLETE PCD × COMPARISON GRID
# ============================================================

comparison_ids <- unique(
  pcd_deseq2$comparison_id
)


pcd_grid <- expand.grid(
  gene_id = pcd_genes,
  comparison_id = comparison_ids,
  stringsAsFactors = FALSE
)


# ============================================================
# ADD DEG INFORMATION
# ============================================================

deg_table <- pcd_grid %>%
  
  left_join(
    
    pcd_deseq2 %>%
      
      select(
        gene_id,
        comparison_id,
        log2FoldChange,
        padj,
        valid_stats,
        DEG,
        direction
      ),
    
    by = c(
      "gene_id",
      "comparison_id"
    )
    
  ) %>%
  
  mutate(
    
    # --------------------------------------------------------
    # Gene found in a particular DESeq2 result
    # --------------------------------------------------------
    
    found_in_DESeq2 =
      !is.na(valid_stats),
    
    
    # --------------------------------------------------------
    # IMPORTANT:
    # Keep DEG = FALSE only when the gene was actually found
    # and had valid statistics.
    #
    # Genes not found remain NA here.
    # --------------------------------------------------------
    
    DEG = case_when(
      
      !found_in_DESeq2 ~ NA,
      
      TRUE ~ DEG
      
    ),
    
    
    # --------------------------------------------------------
    # Direction
    # --------------------------------------------------------
    
    direction = case_when(
      
      !found_in_DESeq2 ~ "Not_found",
      
      TRUE ~ direction
      
    )
    
  )


# ============================================================
# SAVE LONG PCD MATRIX
# ============================================================

write.table(
  deg_table,
  file.path(
    OUT_DIR,
    "PCD_Marchantia_gene_DEG_long_LFC1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# PCD × COMPARISON DEG MATRIX
#
# 1    = DEG
# 0    = found but not DEG
# NA   = gene not found in that DESeq2 result
# ============================================================

deg_wide <- deg_table %>%
  
  select(
    gene_id,
    comparison_id,
    DEG
  ) %>%
  
  mutate(
    DEG = as.integer(DEG)
  ) %>%
  
  pivot_wider(
    names_from = comparison_id,
    values_from = DEG,
    values_fill = list(DEG = 0)
  )


write.table(
  deg_wide,
  file.path(
    OUT_DIR,
    "PCD_Marchantia_gene_DEG_matrix_LFC1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# PCD × COMPARISON DIRECTION MATRIX
#
# Up
# Down
# Not_DEG
# Not_found
# ============================================================

direction_wide <- deg_table %>%
  
  select(
    gene_id,
    comparison_id,
    direction
  ) %>%
  
  pivot_wider(
    names_from = comparison_id,
    values_from = direction,
    values_fill = "Not_found"
  )


write.table(
  direction_wide,
  file.path(
    OUT_DIR,
    "PCD_Marchantia_gene_direction_matrix_LFC1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# PCD GENE SUMMARY
# ============================================================

gene_summary <- deg_table %>%
  
  group_by(
    gene_id
  ) %>%
  
  summarise(
    
    n_comparisons =
      n(),
    
    n_found =
      sum(found_in_DESeq2),
    
    n_DEGs =
      sum(
        DEG,
        na.rm = TRUE
      ),
    
    n_upregulated =
      sum(
        DEG &
          direction == "Up",
        na.rm = TRUE
      ),
    
    n_downregulated =
      sum(
        DEG &
          direction == "Down",
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(
    desc(n_DEGs),
    desc(n_found),
    gene_id
  )


write.table(
  gene_summary,
  file.path(
    OUT_DIR,
    "PCD_Marchantia_gene_DEG_summary_LFC1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# SUMMARY BY COMPARISON
# ============================================================

comparison_summary <- deg_table %>%
  
  mutate(
    
    study = sub(
      "__.*$",
      "",
      comparison_id
    ),
    
    comparison = sub(
      "^[^_]+__",
      "",
      comparison_id
    )
    
  ) %>%
  
  group_by(
    study,
    comparison
  ) %>%
  
  summarise(
    
    PCD_genes_in_input =
      length(pcd_genes),
    
    PCD_genes_found =
      sum(found_in_DESeq2),
    
    PCD_genes_not_found =
      sum(!found_in_DESeq2),
    
    PCD_DEGs =
      sum(
        DEG,
        na.rm = TRUE
      ),
    
    PCD_upregulated =
      sum(
        DEG &
          direction == "Up",
        na.rm = TRUE
      ),
    
    PCD_downregulated =
      sum(
        DEG &
          direction == "Down",
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )


write.table(
  comparison_summary,
  file.path(
    OUT_DIR,
    "PCD_Marchantia_DEG_summary_by_comparison_LFC1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# PCD GENES MISSING FROM ALL DESeq2 RESULTS
# ============================================================

genes_found <- unique(
  pcd_deseq2$gene_id
)


genes_missing <- setdiff(
  pcd_genes,
  genes_found
)


write.table(
  data.frame(
    gene_id = genes_missing
  ),
  file.path(
    OUT_DIR,
    "PCD_Marchantia_genes_not_found_in_DESeq2.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("MARCHANTIA PCD GENE / DESeq2 ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat(
  "PCD genes in input list: ",
  length(pcd_genes),
  "\n",
  sep = ""
)

cat(
  "DESeq2 files analyzed: ",
  length(files),
  "\n",
  sep = ""
)

cat(
  "PCD genes found in >=1 DESeq2 result: ",
  length(genes_found),
  "\n",
  sep = ""
)

cat(
  "PCD genes missing from ALL results: ",
  length(genes_missing),
  "\n",
  sep = ""
)

cat(
  "Total PCD gene-DESeq2 matches: ",
  nrow(pcd_deseq2),
  "\n",
  sep = ""
)

cat(
  "Total PCD DEGs across comparisons: ",
  sum(
    pcd_deseq2$DEG,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Total PCD upregulated calls: ",
  sum(
    pcd_deseq2$DEG &
      pcd_deseq2$direction == "Up",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Total PCD downregulated calls: ",
  sum(
    pcd_deseq2$DEG &
      pcd_deseq2$direction == "Down",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat("\n")
cat("MAIN DEG DEFINITION:\n")
cat("  padj < 0.05 & |log2FoldChange| >= 1\n")


cat("\n")
cat("PCD GENE LIST:\n")
cat(
  "  metadata/PCD_Marchantia_genes.tsv\n"
)


cat("\n")
cat("OUTPUT DIRECTORY:\n")
cat(
  OUT_DIR,
  "\n"
)


cat("\n")
cat("CREATED FILES:\n")

cat(
  "  PCD_Marchantia_genes_all_DESeq2_results.tsv\n"
)

cat(
  "  PCD_Marchantia_gene_DEG_long_LFC1.tsv\n"
)

cat(
  "  PCD_Marchantia_gene_DEG_matrix_LFC1.tsv\n"
)

cat(
  "  PCD_Marchantia_gene_direction_matrix_LFC1.tsv\n"
)

cat(
  "  PCD_Marchantia_gene_DEG_summary_LFC1.tsv\n"
)

cat(
  "  PCD_Marchantia_DEG_summary_by_comparison_LFC1.tsv\n"
)

cat(
  "  PCD_Marchantia_genes_not_found_in_DESeq2.tsv\n"
)


cat("\n")
cat("============================================================\n")

