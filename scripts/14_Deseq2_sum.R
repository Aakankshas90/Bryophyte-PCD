###############################################################
# Summarize all Marchantia DESeq2 comparisons
###############################################################

library(dplyr)

project_dir <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD"
setwd(project_dir)

de_dir <- "results/Marchantia_polymorpha/DESeq2/New Folder"


###############################################################
# Define all DESeq2 result files
###############################################################

files <- c(
  "Tan2023_D2/cold_vs_control.csv",
  "Tan2023_D2/darkness_vs_control.csv",
  "Tan2023_D2/heat_vs_control.csv",
  "Tan2023_D2/high_light_vs_control.csv",
  "Tan2023_D2/nitrogen_deficiency_vs_control.csv",
  "Tan2023_D2/osmotic_vs_control.csv",
  "Tan2023_D2/salt_vs_control.csv",
  "Schroder2023/high_light_vs_control.csv",
  "Grenz2025/infected_vs_control.csv",
  "calcium_deficient_vs_control.csv"
)


###############################################################
# Summarize
###############################################################

summary_list <- list()

for (f in files) {
  
  x <- read.csv(
    file.path(de_dir, f),
    row.names = 1,
    check.names = FALSE
  )
  
  sig <- !is.na(x$padj) & x$padj < 0.05
  
  summary_list[[f]] <- data.frame(
    comparison = gsub(
      "_vs_control.csv",
      "",
      basename(f)
    ),
    study = ifelse(
      grepl("/", f),
      dirname(f),
      "Leong2023"
    ),
    total_genes = nrow(x),
    significant = sum(sig),
    upregulated = sum(
      sig & x$log2FoldChange > 0
    ),
    downregulated = sum(
      sig & x$log2FoldChange < 0
    )
  )
}


###############################################################
# Combine results
###############################################################

summary_table <- bind_rows(summary_list)

print(summary_table)


###############################################################
# Save
###############################################################

write.table(
  summary_table,
  file.path(
    de_dir,
    "DESeq2_summary_all_comparisons.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nDESeq2 summary completed.\n")

