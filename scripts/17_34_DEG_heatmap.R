# ============================================================
# Conserved DEG analysis: 34 genes x 9 comparisons
# ============================================================

library(dplyr)
library(pheatmap)

# ------------------------------------------------------------
# 1. Input file
# ------------------------------------------------------------

input_file <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/results/Marchantia_polymorpha/DESeq2/34_gene/34_DEG_matrix.csv"

df <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Extract the nine log2FoldChange columns
# ------------------------------------------------------------

lfc_cols <- c(
  "Grenz2025_infected_log2FoldChange",
  "Schroder2023_high_light_log2FoldChange",
  "Tan2023_D2_salt_log2FoldChange",
  "Tan2023_D2_nitrogen_deficiency_log2FoldChange",
  "Tan2023_D2_osmotic_log2FoldChange",
  "Tan2023_D2_high_light_log2FoldChange",
  "Tan2023_D2_heat_log2FoldChange",
  "Tan2023_D2_darkness_log2FoldChange",
  "Tan2023_D2_cold_log2FoldChange"
)

lfc_df <- df %>%
  select(gene_id, all_of(lfc_cols))

# ------------------------------------------------------------
# 3. Create log2FC matrix
# ------------------------------------------------------------

lfc_mat <- as.data.frame(lfc_df)

rownames(lfc_mat) <- lfc_mat$gene_id
lfc_mat$gene_id <- NULL

lfc_mat <- as.matrix(lfc_mat)

mode(lfc_mat) <- "numeric"

# ------------------------------------------------------------
# 4. Rename columns for the heatmap
# ------------------------------------------------------------

colnames(lfc_mat) <- c(
  "Infection",
  "High light\n(Schroder)",
  "Salt",
  "N deficiency",
  "Osmotic",
  "High light\n(Tan)",
  "Heat",
  "Darkness",
  "Cold"
)

# ------------------------------------------------------------
# 5. Generate clustered heatmap
# ------------------------------------------------------------

heatmap_colors <- colorRampPalette(
  c("blue", "white", "red")
)(100)

max_abs_lfc <- max(abs(lfc_mat), na.rm = TRUE)

pdf(
  "/Users/aakanksha/Desktop/github/Bryophyte-PCD/results/Marchantia_polymorpha/DESeq2/34_gene/34_gene_log2FC_clustered_heatmap.pdf",
  width = 10,
  height = 12
)

pheatmap(
  lfc_mat,
  color = heatmap_colors,
  breaks = seq(-max_abs_lfc, max_abs_lfc, length.out = 101),
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 9,
  angle_col = 45,
  main = "Conserved DEGs across nine stress comparisons"
)

dev.off()

# ------------------------------------------------------------
# 6. Direction of each response
# ------------------------------------------------------------

direction_mat <- ifelse(
  lfc_mat > 0,
  "Up",
  ifelse(lfc_mat < 0, "Down", "Zero")
)

direction_pattern <- apply(
  direction_mat,
  1,
  paste,
  collapse = ";"
)

# ------------------------------------------------------------
# 7. Count up/down responses across all nine comparisons
# ------------------------------------------------------------

total_up <- rowSums(lfc_mat > 0)
total_down <- rowSums(lfc_mat < 0)

# ------------------------------------------------------------
# 8. Identify consistently up/down/mixed genes
# ------------------------------------------------------------

overall_pattern <- ifelse(
  total_up == 9,
  "Consistently up",
  ifelse(
    total_down == 9,
    "Consistently down",
    "Mixed direction"
  )
)

# ------------------------------------------------------------
# 9. Abiotic stress responses
# ------------------------------------------------------------

abiotic_mat <- lfc_mat[, -1, drop = FALSE]

abiotic_up <- rowSums(abiotic_mat > 0)
abiotic_down <- rowSums(abiotic_mat < 0)

abiotic_predominant <- ifelse(
  abiotic_up > abiotic_down,
  "Up",
  ifelse(
    abiotic_down > abiotic_up,
    "Down",
    "Mixed"
  )
)

# ------------------------------------------------------------
# 10. Infection direction
# ------------------------------------------------------------

infection_direction <- ifelse(
  lfc_mat[, "Infection"] > 0,
  "Up",
  ifelse(
    lfc_mat[, "Infection"] < 0,
    "Down",
    "Zero"
  )
)

# ------------------------------------------------------------
# 11. Compare infection with predominant abiotic direction
# ------------------------------------------------------------

infection_vs_abiotic <- ifelse(
  abiotic_predominant == "Mixed",
  "Mixed abiotic response",
  ifelse(
    infection_direction == abiotic_predominant,
    "Same direction",
    "Opposite direction"
  )
)

# ------------------------------------------------------------
# 12. Create classification table
# ------------------------------------------------------------

classification <- data.frame(
  Gene = rownames(lfc_mat),
  Overall_pattern = overall_pattern,
  Direction_pattern = direction_pattern,
  Total_up = total_up,
  Total_down = total_down,
  Abiotic_up = abiotic_up,
  Abiotic_down = abiotic_down,
  Abiotic_predominant_direction = abiotic_predominant,
  Infection_direction = infection_direction,
  Infection_vs_abiotic = infection_vs_abiotic,
  stringsAsFactors = FALSE
)

# Add the nine log2FC values
classification <- cbind(
  classification,
  as.data.frame(lfc_mat)
)

# ------------------------------------------------------------
# 13. Save classification table
# ------------------------------------------------------------

write.csv(
  classification,
  "/Users/aakanksha/Desktop/github/Bryophyte-PCD/results/Marchantia_polymorpha/DESeq2/34_gene/34_gene_direction_classification.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Print summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("34-GENE CONSERVED DEG ANALYSIS\n")
cat("============================================\n")

cat("\nTotal genes:", nrow(lfc_mat), "\n")

cat(
  "\nConsistently UP across all 9:",
  sum(overall_pattern == "Consistently up"),
  "\n"
)

cat(
  "Consistently DOWN across all 9:",
  sum(overall_pattern == "Consistently down"),
  "\n"
)

cat(
  "Mixed direction:",
  sum(overall_pattern == "Mixed direction"),
  "\n"
)

cat(
  "\nInfection opposite to predominant abiotic direction:",
  sum(infection_vs_abiotic == "Opposite direction"),
  "\n"
)

cat(
  "Same direction:",
  sum(infection_vs_abiotic == "Same direction"),
  "\n"
)

cat(
  "Mixed abiotic response:",
  sum(infection_vs_abiotic == "Mixed abiotic response"),
  "\n"
)

cat("\n============================================\n")
cat("Genes consistently UP:\n")
print(
  classification$Gene[
    classification$Overall_pattern == "Consistently up"
  ]
)

cat("\nGenes consistently DOWN:\n")
print(
  classification$Gene[
    classification$Overall_pattern == "Consistently down"
  ]
)

cat("\nGenes with infection opposite to predominant abiotic direction:\n")
print(
  classification$Gene[
    classification$Infection_vs_abiotic == "Opposite direction"
  ]
)

cat("\n============================================\n")
cat("Output files:\n")
cat("34_gene_log2FC_clustered_heatmap.pdf\n")
cat("34_gene_direction_classification.csv\n")
cat("============================================\n")

