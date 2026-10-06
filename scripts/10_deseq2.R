############################################################
# DESeq2 analysis — individual Marchantia studies
#
# Input:
#   Combined tximport counts.tsv
#   Study-specific metadata files
#
# Analyses:
#   Tan2023_D2      : 7 stresses vs control
#   Schroder2023    : high_light vs control
#   Grenz2025       : infected vs control
#   Leong2022       : calcium_deficiency vs control
#
# Filtering:
#   Gene must have >=10 counts in at least 50% of samples
#
# DESeq2:
#   design = ~ condition
#   control = reference
#   LFC shrinkage = apeglm
############################################################


library(DESeq2)
library(apeglm)


############################################################
# 1. PATHS
############################################################

counts_file <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/results/Marchantia_polymorpha/tximport/combined/counts.tsv"

tan_metadata_file <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/metadata/Tan2023_D2_samplesheet.tsv"

schroder_metadata_file <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/metadata/Schroder2023_samplesheet.tsv"

leong_metadata_file <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/metadata/Leong2022_samplesheet.tsv"

grenz_metadata_file <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/metadata/Grenz2025_samplesheet.tsv"


############################################################
# 2. OUTPUT DIRECTORIES
############################################################

output_base <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/results/Marchantia_polymorpha/DESeq2"

dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

dir.create(
  file.path(output_base, "Tan2023_D2"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(output_base, "Schroder2023"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(output_base, "Leong2022"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(output_base, "Grenz2025"),
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. LOAD COMBINED COUNTS
############################################################

cat("\nLoading combined count matrix...\n")

counts <- read.delim(
  counts_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

counts <- as.matrix(counts)

cat(
  "Combined counts:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "samples\n\n"
)


############################################################
# 4. LOAD METADATA
############################################################

tan_metadata <- read.delim(
  tan_metadata_file,
  row.names = 1,
  check.names = FALSE
)

schroder_metadata <- read.delim(
  schroder_metadata_file,
  row.names = 1,
  check.names = FALSE
)

leong_metadata <- read.delim(
  leong_metadata_file,
  row.names = 1,
  check.names = FALSE
)

grenz_metadata <- read.delim(
  grenz_metadata_file,
  row.names = 1,
  check.names = FALSE
)


############################################################
# 5. FUNCTION TO RUN ONE DESEQ2 COMPARISON
############################################################

run_deseq2 <- function(
    counts,
    metadata,
    study_name,
    stress,
    output_dir
) {
  
  cat("\n")
  cat("========================================\n")
  cat(study_name, ":", stress, "vs control\n")
  cat("========================================\n")
  
  
  ########################################################
  # Check metadata
  ########################################################
  
  if (!"condition" %in% colnames(metadata)) {
    stop(
      paste(
        "Metadata for",
        study_name,
        "does not contain a 'condition' column."
      )
    )
  }
  
  
  ########################################################
  # Make sure metadata samples exist in count matrix
  ########################################################
  
  missing_samples <- setdiff(
    rownames(metadata),
    colnames(counts)
  )
  
  if (length(missing_samples) > 0) {
    
    stop(
      paste(
        "These metadata samples are missing from counts:",
        paste(missing_samples, collapse = ", ")
      )
    )
  }
  
  
  ########################################################
  # Select only this study's samples
  ########################################################
  
  counts_study <- counts[, rownames(metadata), drop = FALSE]
  
  
  ########################################################
  # Select only control + requested stress
  ########################################################
  
  keep_samples <- metadata$condition %in%
    c("control", stress)
  
  metadata_sub <- metadata[keep_samples, , drop = FALSE]
  
  counts_sub <- counts_study[
    ,
    rownames(metadata_sub),
    drop = FALSE
  ]
  
  
  ########################################################
  # IMPORTANT CHECK
  ########################################################
  
  cat("\nSamples included:\n")
  
  print(
    data.frame(
      sample = rownames(metadata_sub),
      condition = metadata_sub$condition
    )
  )
  
  cat("\nCondition counts:\n")
  
  print(table(metadata_sub$condition))
  
  
  ########################################################
  # Check that both groups are present
  ########################################################
  
  if (
    !all(
      c("control", stress) %in%
      unique(as.character(metadata_sub$condition))
    )
  ) {
    
    stop(
      paste(
        "\nERROR:",
        study_name,
        stress,
        "does not contain both control and",
        stress,
        "samples."
      )
    )
  }
  
  
  ########################################################
  # Set condition factor
  ########################################################
  
  metadata_sub$condition <- factor(
    metadata_sub$condition,
    levels = c("control", stress)
  )
  
  
  ########################################################
  # Check sample order
  ########################################################
  
  stopifnot(
    identical(
      colnames(counts_sub),
      rownames(metadata_sub)
    )
  )
  
  
  ########################################################
  # Create DESeq2 object
  ########################################################
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData = metadata_sub,
    design = ~ condition
  )
  
  
  ########################################################
  # Low-count filtering
  #
  # Keep genes with >=10 counts in at least 50% of samples
  ########################################################
  
  keep_genes <- rowSums(
    counts(dds) >= 10
  ) >= ceiling(
    0.5 * ncol(dds)
  )
  
  dds <- dds[keep_genes, ]
  
  
  cat(
    "\nGenes retained after filtering:",
    nrow(dds),
    "\n"
  )
  
  
  ########################################################
  # Run DESeq2
  ########################################################
  
  dds <- DESeq(dds)
  
  
  ########################################################
  # Get coefficient
  ########################################################
  
  results_names <- resultsNames(dds)
  
  print(results_names)
  
  
  coef_name <- results_names[
    grep(
      paste0(
        "^condition_",
        stress,
        "_vs_control$"
      ),
      results_names
    )
  ]
  
  
  if (length(coef_name) != 1) {
    
    stop(
      paste(
        "Could not uniquely identify coefficient for",
        stress
      )
    )
  }
  
  
  ########################################################
  # Standard DESeq2 results
  ########################################################
  
  res <- results(
    dds,
    contrast = c(
      "condition",
      stress,
      "control"
    )
  )
  
  
  ########################################################
  # LFC shrinkage with apeglm
  ########################################################
  
  res_shrunk <- lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )
  
  
  ########################################################
  # Print DEG summary
  ########################################################
  
  cat(
    "padj < 0.05:",
    sum(
      res_shrunk$padj < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "padj < 0.05 & |LFC| >= 0.58:",
    sum(
      res_shrunk$padj < 0.05 &
        abs(res_shrunk$log2FoldChange) >= 0.58,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "padj < 0.05 & |LFC| >= 1:",
    sum(
      res_shrunk$padj < 0.05 &
        abs(res_shrunk$log2FoldChange) >= 1,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "padj < 0.01 & |LFC| >= 1:",
    sum(
      res_shrunk$padj < 0.01 &
        abs(res_shrunk$log2FoldChange) >= 1,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  ########################################################
  # Save standard results
  ########################################################
  
  write.csv(
    as.data.frame(res),
    file = file.path(
      output_dir,
      paste0(
        stress,
        "_vs_control_DESeq2.csv"
      )
    )
  )
  
  
  ########################################################
  # Save apeglm-shrunk results
  ########################################################
  
  write.csv(
    as.data.frame(res_shrunk),
    file = file.path(
      output_dir,
      paste0(
        stress,
        "_vs_control_apeglm.csv"
      )
    )
  )
  
  
  ########################################################
  # Return objects
  ########################################################
  
  return(
    list(
      dds = dds,
      results = res,
      shrunk = res_shrunk,
      metadata = metadata_sub
    )
  )
}


############################################################
# 6. TAN2023_D2
############################################################

tan_stresses <- c(
  "cold",
  "darkness",
  "heat",
  "high_light",
  "osmotic",
  "nitrogen_deficiency",
  "salt"
)

tan_results <- list()

for (stress in tan_stresses) {
  
  tan_results[[stress]] <- run_deseq2(
    counts = counts,
    metadata = tan_metadata,
    study_name = "Tan2023_D2",
    stress = stress,
    output_dir = file.path(
      output_base,
      "Tan2023_D2"
    )
  )
}


############################################################
# 7. SCHRODER2023
############################################################

schroder_results <- run_deseq2(
  counts = counts,
  metadata = schroder_metadata,
  study_name = "Schroder2023",
  stress = "high_light",
  output_dir = file.path(
    output_base,
    "Schroder2023"
  )
)


############################################################
# 8. GRENZ2025
############################################################

grenz_results <- run_deseq2(
  counts = counts,
  metadata = grenz_metadata,
  study_name = "Grenz2025",
  stress = "infected",
  output_dir = file.path(
    output_base,
    "Grenz2025"
  )
)


############################################################
# 9. LEONG2022
############################################################

leong_results <- run_deseq2(
  counts = counts,
  metadata = leong_metadata,
  study_name = "Leong2022",
  stress = "calcium_deficiency",
  output_dir = file.path(
    output_base,
    "Leong2022"
  )
)


############################################################
# DONE
############################################################

cat("\n")
cat("========================================\n")
cat("ALL DESEQ2 ANALYSES COMPLETED\n")
cat("========================================\n")