# ============================================================
# Conserved Gene Matrix from 9 DESeq2 Results
# Keep genes present in ALL 9 result files
# ============================================================

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base_dir <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/results/Marchantia_polymorpha/DESeq2"

output_file <- file.path(base_dir, "conserved_DEG_matrix_all_genes.csv")
log_file    <- file.path(base_dir, "conserved_DEG_matrix_all_genes.log")


# ------------------------------------------------------------
# Start log
# ------------------------------------------------------------

log_con <- file(log_file, open = "wt")

log_message <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_con)
}


# ------------------------------------------------------------
# DESeq2 result files
# ------------------------------------------------------------

files <- c(
  Grenz2025_infected =
    file.path(base_dir, "Grenz2025/infected_vs_control_apeglm.csv"),
  
  Schroder2023_high_light =
    file.path(base_dir, "Schroder2023/high_light_vs_control_apeglm.csv"),
  
  Tan2023_D2_salt =
    file.path(base_dir, "Tan2023/salt_vs_control_apeglm.csv"),
  
  Tan2023_D2_nitrogen_deficiency =
    file.path(base_dir, "Tan2023/nitrogen_deficiency_vs_control_apeglm.csv"),
  
  Tan2023_D2_osmotic =
    file.path(base_dir, "Tan2023/osmotic_vs_control_apeglm.csv"),
  
  Tan2023_D2_high_light =
    file.path(base_dir, "Tan2023/high_light_vs_control_apeglm.csv"),
  
  Tan2023_D2_heat =
    file.path(base_dir, "Tan2023/heat_vs_control_apeglm.csv"),
  
  Tan2023_D2_darkness =
    file.path(base_dir, "Tan2023/darkness_vs_control_apeglm.csv"),
  
  Tan2023_D2_cold =
    file.path(base_dir, "Tan2023/cold_vs_control_apeglm.csv")
)


# ------------------------------------------------------------
# Read all DESeq2 results
# ------------------------------------------------------------

log_message("============================================================")
log_message("Conserved Gene Matrix")
log_message("============================================================")
log_message("")
log_message("Number of DESeq2 result files: ", length(files))
log_message("")

results <- list()

for (i in seq_along(files)) {
  
  name <- names(files)[i]
  
  log_message("Reading: ", name)
  
  df <- read.csv(
    files[i],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  results[[name]] <- df
  
  log_message("  Genes: ", nrow(df))
  log_message("  Columns: ", ncol(df))
  log_message("")
}



# ------------------------------------------------------------
# Filter significant DEGs
# ------------------------------------------------------------

# ------------------------------------------------------------
# Filter significant DEGs
# padj < 0.05 AND |log2FoldChange| >= 1
# ------------------------------------------------------------

for (name in names(results)) {
  
  results[[name]] <- results[[name]][
    !is.na(results[[name]]$padj) &
      results[[name]]$padj < 0.05 &
      !is.na(results[[name]]$log2FoldChange) &
      abs(results[[name]]$log2FoldChange) >= 1,
  ]
  
  log_message(
    name,
    " significant DEGs (padj < 0.05, |LFC| >= 1): ",
    nrow(results[[name]])
  )
}

log_message("")

# ------------------------------------------------------------
# Identify gene ID column
# ------------------------------------------------------------

gene_col <- "gene_id"


# ------------------------------------------------------------
# Find genes present in ALL 9 files
# ------------------------------------------------------------

gene_lists <- lapply(results, function(x) x[[gene_col]])

conserved_genes <- Reduce(intersect, gene_lists)


log_message("============================================================")
log_message("Gene intersection")
log_message("============================================================")
log_message("")
log_message("Genes present in all result files: ", length(conserved_genes))
log_message("")


# ------------------------------------------------------------
# Create conserved matrix
# ------------------------------------------------------------

conserved_tables <- list()

for (name in names(results)) {
  
  df <- results[[name]]
  
  df <- df[df[[gene_col]] %in% conserved_genes, ]
  
  # Keep the same original DESeq2 columns,
  # but prefix them with the comparison name
  colnames(df)[colnames(df) != gene_col] <-
    paste0(name, "_", colnames(df)[colnames(df) != gene_col])
  
  conserved_tables[[name]] <- df
}


# ------------------------------------------------------------
# Merge all comparisons
# ------------------------------------------------------------

conserved_matrix <- conserved_tables[[1]]

for (i in 2:length(conserved_tables)) {
  
  conserved_matrix <- merge(
    conserved_matrix,
    conserved_tables[[i]],
    by = gene_col,
    all = TRUE,
    sort = FALSE
  )
}


# ------------------------------------------------------------
# Put gene_id first and sort by gene ID
# ------------------------------------------------------------

conserved_matrix <- conserved_matrix[
  order(conserved_matrix[[gene_col]]),
]

conserved_matrix <- conserved_matrix[
  ,
  c(gene_col, setdiff(colnames(conserved_matrix), gene_col))
]


# ------------------------------------------------------------
# Save matrix
# ------------------------------------------------------------

write.csv(
  conserved_matrix,
  output_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

log_message("============================================================")
log_message("FINAL SUMMARY")
log_message("============================================================")
log_message("")
log_message("Input files: ", length(files))
log_message("Genes present in ALL files: ", length(conserved_genes))
log_message("Final matrix rows: ", nrow(conserved_matrix))
log_message("Final matrix columns: ", ncol(conserved_matrix))
log_message("")
log_message("Output matrix:")
log_message(output_file)
log_message("")
log_message("Log file:")
log_message(log_file)
log_message("")
log_message("Completed successfully.")
log_message("============================================================")


# ------------------------------------------------------------
# Close log
# ------------------------------------------------------------

close(log_con)
