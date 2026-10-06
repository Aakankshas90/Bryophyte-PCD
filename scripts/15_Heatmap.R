# ============================================================
# Marchantia polymorpha
#
# Current DESeq2 results
#
# INCLUDED:
#   Grenz2025 infected
#   Schroder2023 high light
#   Tan2023_D2: 7 stresses
#
#
# DEG definition:
#   padj < 0.05
#   AND |log2FoldChange| >= 1
#
# Gene sets:
#   Marchantia PCD genes
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. PROJECT
# ------------------------------------------------------------

PROJECT_DIR <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD"

DEG_DIR <- file.path(
  PROJECT_DIR,
  "results/Marchantia_polymorpha/DESeq2"
)

# ------------------------------------------------------------
# 2. CURRENT DESeq2 FILES
# ------------------------------------------------------------

files <- c(
  
  file.path(DEG_DIR,
            "Tan2023/cold_vs_control_apeglm.csv"),
  
  file.path(DEG_DIR,
            "Tan2023/darkness_vs_control_apeglm.csv"),
  
  file.path(DEG_DIR,
            "Tan2023/heat_vs_control_apeglm.csv"),
  
  file.path(DEG_DIR,
            "Tan2023/high_light_vs_control_apeglm.csv"),
  
  file.path(DEG_DIR,
            "Tan2023/nitrogen_deficiency_vs_control_apeglm.csv"),
  
  file.path(DEG_DIR,
            "Tan2023/osmotic_vs_control_apeglm.csv"),
  
  file.path(DEG_DIR,
            "Tan2023/salt_vs_control_apeglm.csv"),
  
  file.path(DEG_DIR,
            "Schroder2023/high_light_vs_control_apeglm.csv"),
  
  file.path(DEG_DIR,
            "Grenz2025/infected_vs_control_apeglm.csv")
)

# ------------------------------------------------------------
# Check files
# ------------------------------------------------------------

if (any(!file.exists(files))) {
  
  cat("\nMissing files:\n")
  
  print(files[!file.exists(files)])
  
  stop("Analysis stopped.")
  
}

cat("\nAll 9 DESeq2 files found.\n")


# ============================================================
# 3. FUNCTION
# ============================================================

run_pcd <- function(
    gene_file,
    output_dir,
    title_text
) {
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ----------------------------------------------------------
  # Read PCD genes
  # ----------------------------------------------------------
  
  genes <- read.delim(
    gene_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    pull(gene_id) %>%
    trimws() %>%
    unique()
  
  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]
  
  cat("\n============================================\n")
  cat(title_text, "\n")
  cat("PCD genes:", length(genes), "\n")
  cat("============================================\n")
  
  # ----------------------------------------------------------
  # Read the 9 DESeq2 files
  # ----------------------------------------------------------
  
  all_results <- lapply(files, function(f) {
    
    x <- read.csv(
      f,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    # First column is gene ID
    colnames(x)[1] <- "gene_id"
    
    x$gene_id <- trimws(x$gene_id)
    
    # Study
    x$study <- basename(dirname(f))
    
    # Comparison
    x$comparison <- tools::file_path_sans_ext(
      basename(f)
    )
    
    # Keep only PCD genes
    x %>%
      filter(gene_id %in% genes) %>%
      select(
        gene_id,
        study,
        comparison,
        baseMean,
        log2FoldChange,
        lfcSE,
        pvalue,
        padj
      )
    
  })
  
  df <- bind_rows(all_results)
  
  # ----------------------------------------------------------
  # Classify DEGs
  # ----------------------------------------------------------
  
  df <- df %>%
    mutate(
      
      DEG =
        !is.na(padj) &
        !is.na(log2FoldChange) &
        padj < 0.05 &
        abs(log2FoldChange) >= 1,
      
      response = case_when(
        
        is.na(padj) |
          is.na(log2FoldChange) ~ "NA",
        
        DEG &
          log2FoldChange > 0 ~ "UP",
        
        DEG &
          log2FoldChange < 0 ~ "DOWN",
        
        TRUE ~ "NS"
        
      ),
      
      stress = str_remove(
        comparison,
        "_apeglm$"
      ),
      
      stress = str_remove(
        stress,
        "_vs_control$"
      )
      
    )
  
  # ----------------------------------------------------------
  # Save complete table
  # ----------------------------------------------------------
  
  write.table(
    df,
    file.path(
      output_dir,
      "PCD_DE_complete_9_comparisons_LFC1.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # Significant PCD DEGs
  # ----------------------------------------------------------
  
  significant <- df %>%
    filter(DEG) %>%
    arrange(padj)
  
  write.table(
    significant,
    file.path(
      output_dir,
      "PCD_DE_significant_padj0.05_LFC1.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # Summary by comparison
  # ----------------------------------------------------------
  
  comparison_summary <- df %>%
    group_by(
      study,
      stress
    ) %>%
    summarise(
      
      total_PCD_genes = length(genes),
      
      PCD_genes_found = n(),
      
      significant = sum(
        DEG,
        na.rm = TRUE
      ),
      
      upregulated = sum(
        DEG &
          log2FoldChange > 0,
        na.rm = TRUE
      ),
      
      downregulated = sum(
        DEG &
          log2FoldChange < 0,
        na.rm = TRUE
      ),
      
      .groups = "drop"
      
    )
  
  write.table(
    comparison_summary,
    file.path(
      output_dir,
      "PCD_DE_summary_by_comparison_LFC1.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # Summary by gene
  # ----------------------------------------------------------
  
  gene_summary <- df %>%
    group_by(gene_id) %>%
    summarise(
      
      comparisons_found = n(),
      
      significant_comparisons =
        sum(
          padj < 0.05,
          na.rm = TRUE
        ),
      
      DEG_comparisons =
        sum(
          DEG,
          na.rm = TRUE
        ),
      
      upregulated_comparisons =
        sum(
          DEG &
            log2FoldChange > 0,
          na.rm = TRUE
        ),
      
      downregulated_comparisons =
        sum(
          DEG &
            log2FoldChange < 0,
          na.rm = TRUE
        ),
      
      .groups = "drop"
      
    ) %>%
    arrange(
      desc(DEG_comparisons),
      desc(significant_comparisons)
    )
  
  write.table(
    gene_summary,
    file.path(
      output_dir,
      "PCD_DE_summary_by_gene_LFC1.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # ==========================================================
  # HEATMAP
  # ==========================================================
  
  comparison_order <- c(
    "Grenz2025\nInfected",
    "Schroder2023\nHigh light",
    "Tan2023\nCold",
    "Tan2023\nDarkness",
    "Tan2023\nHeat",
    "Tan2023\nHigh light",
    "Tan2023\nNitrogen deficiency",
    "Tan2023\nOsmotic",
    "Tan2023\nSalt"
  )
  
  # ----------------------------------------------------------
  # Create heatmap labels directly from study + comparison
  # ----------------------------------------------------------
  
  df <- df %>%
    mutate(
      comparison_label = case_when(
        
        study == "Grenz2025" &
          comparison == "infected_vs_control_apeglm" ~
          "Grenz2025\nInfected",
        
        study == "Schroder2023" &
          comparison == "high_light_vs_control_apeglm" ~
          "Schroder2023\nHigh light",
        
        study == "Tan2023" &
          comparison == "cold_vs_control_apeglm" ~
          "Tan2023\nCold",
        
        study == "Tan2023" &
          comparison == "darkness_vs_control_apeglm" ~
          "Tan2023\nDarkness",
        
        study == "Tan2023" &
          comparison == "heat_vs_control_apeglm" ~
          "Tan2023\nHeat",
        
        study == "Tan2023" &
          comparison == "high_light_vs_control_apeglm" ~
          "Tan2023\nHigh light",
        
        study == "Tan2023" &
          comparison == "nitrogen_deficiency_vs_control_apeglm" ~
          "Tan2023\nNitrogen deficiency",
        
        study == "Tan2023" &
          comparison == "osmotic_vs_control_apeglm" ~
          "Tan2023\nOsmotic",
        
        study == "Tan2023" &
          comparison == "salt_vs_control_apeglm" ~
          "Tan2023\nSalt",
        
        TRUE ~ NA_character_
      )
    )
  
  # ----------------------------------------------------------
  # IMPORTANT CHECK
  # ----------------------------------------------------------
  
  cat("\nHeatmap comparisons detected:\n")
  print(
    df %>%
      count(
        study,
        comparison,
        comparison_label
      )
  )
  
  # Stop if anything failed to receive a label
  if (any(is.na(df$comparison_label))) {
    
    cat("\nERROR: Some comparisons did not receive a heatmap label:\n")
    
    print(
      df %>%
        filter(is.na(comparison_label)) %>%
        distinct(study, comparison)
    )
    
    stop(
      "Heatmap labeling failed. Check the comparison names above."
    )
  }
  
  # ----------------------------------------------------------
  # Set factor order
  # ----------------------------------------------------------
  
  df$comparison_label <- factor(
    df$comparison_label,
    levels = comparison_order
  )
  
  # ----------------------------------------------------------
  # Gene order
  # ----------------------------------------------------------
  
  gene_order <- df %>%
    group_by(gene_id) %>%
    summarise(
      n_sig = sum(
        DEG,
        na.rm = TRUE
      ),
      mean_abs_LFC = mean(
        abs(log2FoldChange),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    arrange(
      desc(n_sig),
      desc(mean_abs_LFC)
    ) %>%
    pull(gene_id)
  
  # ----------------------------------------------------------
  # Heatmap data
  # ----------------------------------------------------------
  
  heatmap_data <- df %>%
    select(
      gene_id,
      comparison_label,
      log2FoldChange,
      padj,
      DEG            # <- carry through the correctly-thresholded flag
    ) %>%
    mutate(
      gene_id = factor(gene_id, levels = rev(gene_order)),
      
      significance = case_when(
        !DEG ~ "",
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE ~ ""
      ),
      
      log2FC_plot = pmax(pmin(log2FoldChange, 4), -4)
    )
  
  cat("\n============================================\n")
  cat("HEATMAP SIGNIFICANCE CHECK\n")
  cat("============================================\n")
  
  cat(
    "\nTotal significant gene-condition pairs:",
    sum(heatmap_data$DEG, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "One star (padj < 0.05):",
    sum(
      heatmap_data$DEG &
        heatmap_data$padj < 0.05 &
        heatmap_data$padj >= 0.01,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Two stars (padj < 0.01):",
    sum(
      heatmap_data$DEG &
        heatmap_data$padj < 0.01 &
        heatmap_data$padj >= 0.001,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Three stars (padj < 0.001):",
    sum(
      heatmap_data$DEG &
        heatmap_data$padj < 0.001,
      na.rm = TRUE
    ),
    "\n"
  )
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  # Gene ID to gene symbol mapping
  gene_symbols <- c(
    "Mp1g17170" = "MpTHIO",
    "Mp1g28640" = "MpAIG2",
    "Mp3g20340" = "MpRBOHA",
    "Mp7g17040" = "MpMCA-I",
    "Mp3g05030" = "MpATG1",
    "Mp6g18570" = "MpATG2",
    "Mp1g12840" = "MpATG5",
    "Mp2g07850" = "MpATG7",
    "Mp1g25570" = "MpATG10",
    "Mp2g25440" = "MpATG12",
    "Mp7g03210" = "MpATG13",
    "Mp4g19650" = "MpZOU1",
    "Mp8g04130" = "MpZOU2",
    "Mp4g04910" = "MpICE1"
  )
  
  heatmap_data <- heatmap_data %>%
    mutate(
      gene_id_chr = as.character(gene_id),
      gene_label = if_else(
        gene_id_chr %in% names(gene_symbols),
        paste0(
          unname(gene_symbols[gene_id_chr]),
          " (",
          gene_id_chr,
          ")"
        ),
        gene_id_chr
      )
    )
  
  
  p <- ggplot(
    heatmap_data,
    aes(
      x = comparison_label,
      y = gene_label,
      fill = log2FC_plot
    )
  ) +
    
    geom_tile(
      color = "white",
      linewidth = 0.3
    ) +
    
    geom_text(
      aes(
        label = significance
      ),
      size = 3
    ) +
    
    scale_fill_gradient2(
      name = "log2 fold change",
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-4, 4),
      oob = scales::squish
    ) +
    
    labs(
      title = title_text,
      subtitle =
        "* padj < 0.05, ** padj < 0.01, *** padj < 0.001",
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal(
      base_size = 10
    ) +
    
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 0.2,
        vjust = 0.1
      ),
      axis.text.y = element_text(
        size = 8
      ),
      panel.grid = element_blank(),
      plot.title = element_text(
        face = "bold"
      )
    )
  
  # ----------------------------------------------------------
  # Save heatmap
  # ----------------------------------------------------------
  
 # ggsave(
  #  file.path(
   #   output_dir,
    #  "PCD_DE_heatmap_9_comparisons_LFC1.pdf"
  #  ),
  #  p,
  #  width = 11,
  #  height = 8
  #)
  
  ggsave(
    file.path(
      output_dir,
      "PCD_DE_heatmap_9_comparisons_LFC1.png"
    ),
    p,
    width = 11,
    height = 8,
    dpi = 300
  )
  
  # ----------------------------------------------------------
  # Final console summary
  # ----------------------------------------------------------
  
  cat("\nCompleted:", title_text, "\n")
  cat("PCD genes:", length(genes), "\n")
  cat("PCD gene-condition observations:", nrow(df), "\n")
  cat("Significant PCD DEGs:", nrow(significant), "\n")
  cat(
    "UP:",
    sum(significant$log2FoldChange > 0),
    "\n"
  )
  cat(
    "DOWN:",
    sum(significant$log2FoldChange < 0),
    "\n"
  )
  
  invisible(df)
}


# ============================================================
# 4. MARCHANTIA-ONLY PCD GENES
# ============================================================

marchantia <- run_pcd(
  
  gene_file = file.path(
    PROJECT_DIR,
    "metadata/PCD_Marchantia_genes.tsv"
  ),
  
  output_dir = file.path(
    DEG_DIR,
    "PCD_gene_check"
  ),
  
  title_text =
    "RCD-associated genes"
)


# ============================================================
# 6. DONE
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("ITERATION 3 — COMPLETE\n")
cat("============================================================\n")

cat("\nCurrent analysis:\n")
cat("  9 comparisons\n")
cat("  padj < 0.05\n")
cat("  |log2FoldChange| >= 1\n")

cat("\nOutputs:\n")
cat(
  "  PCD_gene_check/\n"
)

cat("\n============================================================\n")

