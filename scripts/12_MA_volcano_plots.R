############################################################
# 13_MA_volcano_plots.R
#
# MA plots and volcano plots for DESeq2 results
#
# Designed for running directly in RStudio
############################################################

library(ggplot2)
library(dplyr)


############################################################
# PATHS
############################################################

PROJECT_DIR <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD"

INPUT_DIR <- file.path(
  PROJECT_DIR,
  "results/Marchantia_polymorpha/DESeq2"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results/Marchantia_polymorpha/DESeq2/plots"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# PLOT FUNCTION
############################################################

make_plots <- function(
    result_file,
    comparison_name,
    output_prefix
) {
  
  cat("\n========================================\n")
  cat(comparison_name, "\n")
  cat("========================================\n")
  
  
  ##########################################################
  # Read DESeq2 results
  ##########################################################
  
  res <- read.csv(
    result_file,
    row.names = 1,
    check.names = FALSE
  )
  
  
  ##########################################################
  # Check required columns
  ##########################################################
  
  required_cols <- c(
    "baseMean",
    "log2FoldChange",
    "padj"
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(res)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      paste(
        "Missing columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  
  ##########################################################
  # Clean data
  ##########################################################
  
  res <- res %>%
    mutate(
      gene = rownames(res)
    )
  
  
  ##########################################################
  # Significance categories
  ##########################################################
  
  res <- res %>%
    mutate(
      significance = case_when(
        
        !is.na(padj) &
          padj < 0.05 &
          log2FoldChange >= 1
        ~ "Upregulated",
        
        !is.na(padj) &
          padj < 0.05 &
          log2FoldChange <= -1
        ~ "Downregulated",
        
        TRUE
        ~ "Not significant"
      )
    )
  
  
  ##########################################################
  # Print summary
  ##########################################################
  
  cat(
    "Total genes:",
    nrow(res),
    "\n"
  )
  
  cat(
    "Significant:",
    sum(
      res$padj < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Upregulated:",
    sum(
      res$padj < 0.05 &
        res$log2FoldChange >= 1,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Downregulated:",
    sum(
      res$padj < 0.05 &
        res$log2FoldChange <= -1,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  ##########################################################
  # MA PLOT
  ##########################################################
  
  ma_data <- res %>%
    filter(
      !is.na(baseMean),
      !is.na(log2FoldChange),
      baseMean > 0
    )
  
  
  ma_plot <- ggplot(
    ma_data,
    aes(
      x = log10(baseMean + 1),
      y = log2FoldChange
    )
  ) +
    
    geom_point(
      aes(color = significance),
      size = 1,
      alpha = 0.5
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = c(-1, 1),
      linetype = "dotted"
    ) +
    
    scale_color_manual(
      values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not significant" = "grey70"
      )
    ) +
    
    labs(
      title = paste0(
        "MA plot: ",
        comparison_name
      ),
      x = "log10(baseMean + 1)",
      y = "log2 fold change",
      color = "Category"
    ) +
    
    theme_classic(base_size = 14)
  
  
  ##########################################################
  # Save MA plot
  ##########################################################
  
  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      paste0(output_prefix, "_MA.png")
    ),
    plot = ma_plot,
    width = 8,
    height = 7,
    dpi = 300
  )
  
  
  ##########################################################
  # VOLCANO PLOT
  ##########################################################
  
  volcano_data <- res %>%
    filter(
      !is.na(log2FoldChange),
      !is.na(padj)
    ) %>%
    
    mutate(
      neg_log10_padj =
        -log10(pmax(padj, 1e-300))
    )
  
  
  volcano_plot <- ggplot(
    volcano_data,
    aes(
      x = log2FoldChange,
      y = neg_log10_padj
    )
  ) +
    
    geom_point(
      aes(color = significance),
      size = 1,
      alpha = 0.5
    ) +
    
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    
    scale_color_manual(
      values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not significant" = "grey70"
      )
    ) +
    
    labs(
      title = paste0(
        "Volcano plot: ",
        comparison_name
      ),
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = "Category"
    ) +
    
    theme_classic(base_size = 14)
  
  
  ##########################################################
  # Save volcano plot
  ##########################################################
  
  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      paste0(output_prefix, "_volcano.png")
    ),
    plot = volcano_plot,
    width = 8,
    height = 7,
    dpi = 300
  )
  
  
  ##########################################################
  # Also save PDF versions
  ##########################################################
  
  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      paste0(output_prefix, "_MA.pdf")
    ),
    plot = ma_plot,
    width = 8,
    height = 7
  )
  
  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      paste0(output_prefix, "_volcano.pdf")
    ),
    plot = volcano_plot,
    width = 8,
    height = 7
  )
  
  
  cat("Plots saved.\n")
}


############################################################
# TAN2023 D2
############################################################

tan_dir <- file.path(
  INPUT_DIR,
  "Tan2023_D2"
)

tan_stresses <- c(
  "cold",
  "darkness",
  "heat",
  "high_light",
  "osmotic",
  "nitrogen_deficiency",
  "salt"
)


for (stress in tan_stresses) {
  
  result_file <- file.path(
    tan_dir,
    paste0(
      stress,
      "_vs_control_apeglm.csv"
    )
  )
  
  if (file.exists(result_file)) {
    
    make_plots(
      result_file = result_file,
      
      comparison_name = paste(
        "Tan2023_D2:",
        stress,
        "vs control"
      ),
      
      output_prefix = paste0(
        "Tan2023_D2_",
        stress,
        "_vs_control"
      )
    )
    
  } else {
    
    cat(
      "\nWARNING: File not found:",
      result_file,
      "\n"
    )
  }
}


############################################################
# SCHRODER 2023
############################################################

schroder_file <- file.path(
  INPUT_DIR,
  "Schroder2023",
  "high_light_vs_control_apeglm.csv"
)

if (file.exists(schroder_file)) {
  
  make_plots(
    result_file = schroder_file,
    
    comparison_name =
      "Schroder2023: high light vs control",
    
    output_prefix =
      "Schroder2023_high_light_vs_control"
  )
}


############################################################
# GRENZ 2025
############################################################

grenz_file <- file.path(
  INPUT_DIR,
  "Grenz2025",
  "infected_vs_control_apeglm.csv"
)

if (file.exists(grenz_file)) {
  
  make_plots(
    result_file = grenz_file,
    
    comparison_name =
      "Grenz2025: infected vs control",
    
    output_prefix =
      "Grenz2025_infected_vs_control"
  )
}


############################################################
# LEONG 2022/2023
############################################################

leong_file <- file.path(
  INPUT_DIR,
  "Leong2022",
  "calcium_deficiency_vs_control_apeglm.csv"
)

if (file.exists(leong_file)) {
  
  make_plots(
    result_file = leong_file,
    
    comparison_name =
      "Leong2022: calcium deficiency vs control",
    
    output_prefix =
      "Leong2022_calcium_deficiency_vs_control"
  )
}


############################################################
# DONE
############################################################

cat("\n")
cat("========================================\n")
cat("ALL PLOTS COMPLETE\n")
cat("========================================\n")

cat(
  "Output directory:\n",
  OUTPUT_DIR,
  "\n"
)

