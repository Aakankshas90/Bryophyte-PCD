############################################################
# WGCNA ANALYSIS OF BATCH-CORRECTED VST RNA-SEQ DATA
#
# Input:
# vst_batchcorrected_count.csv
############################################################


############################################################
# 0. INSTALL AND LOAD REQUIRED PACKAGES
############################################################

# Uncomment if packages are not installed

# install.packages(c(
#   "tidyverse",
#   "magrittr",
#   "WGCNA"
# ))

library(tidyverse)
library(magrittr)
library(WGCNA)

options(stringsAsFactors = FALSE)

# Allow multi-threading
allowWGCNAThreads()


############################################################
# 1. SET WORKING DIRECTORY
############################################################

setwd("C:/bg/marchantia")

# Create output directory

output_dir <- "WGCNA_results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


############################################################
# 2. LOAD EXPRESSION DATA
############################################################

data <- read.csv(
  "vst_batchcorrected_count.csv",
  check.names = FALSE
)

# Check structure

str(data)

# Rename first column as Geneid

names(data)[1] <- "Geneid"

# Remove duplicate genes if present

data <- data %>%
  distinct(Geneid, .keep_all = TRUE)

# Get sample columns

col_sel <- names(data)[-1]

cat("Number of genes:", nrow(data), "\n")
cat("Number of samples:", length(col_sel), "\n")


############################################################
# 3. SAMPLE DISTRIBUTION / OUTLIER CHECK
############################################################

mdata <- data %>%
  tidyr::pivot_longer(
    cols = all_of(col_sel),
    names_to = "name",
    values_to = "value"
  ) %>%
  mutate(
    group = sub("_.*", "", name)
  )

# Violin plot

p <- mdata %>%
  ggplot(
    aes(
      x = name,
      y = value
    )
  ) +
  geom_violin(
    fill = "lightblue"
  ) +
  geom_point(
    alpha = 0.15,
    size = 0.3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  ) +
  labs(
    title = "Sample Expression Distributions",
    x = "Samples",
    y = "VST Expression"
  ) +
  facet_grid(
    cols = vars(group),
    scales = "free_x",
    space = "free_x"
  )

print(p)

ggsave(
  filename = file.path(
    output_dir,
    "Sample_expression_distribution.pdf"
  ),
  plot = p,
  width = 18,
  height = 8
)


############################################################
# 4. PREPARE EXPRESSION MATRIX FOR WGCNA
#
# WGCNA requires:
# Rows = samples
# Columns = genes
############################################################

de_input <- as.matrix(
  data[, -1]
)

rownames(de_input) <- data$Geneid

# Ensure numeric matrix

storage.mode(de_input) <- "numeric"

# Transpose:
# Rows = samples
# Columns = genes

input_mat <- t(de_input)

# Check dimensions

cat(
  "WGCNA matrix dimensions:",
  nrow(input_mat),
  "samples x",
  ncol(input_mat),
  "genes\n"
)


############################################################
# 5. CHECK GOOD SAMPLES AND GENES
############################################################

gsg <- goodSamplesGenes(
  input_mat,
  verbose = 3
)

if (!gsg$allOK) {
  
  input_mat <- input_mat[
    gsg$goodSamples,
    gsg$goodGenes
  ]
  
}


############################################################
# 6. SAMPLE CLUSTERING FOR OUTLIER DETECTION
############################################################

sampleTree <- hclust(
  dist(input_mat),
  method = "average"
)

pdf(
  file.path(
    output_dir,
    "Sample_clustering.pdf"
  ),
  width = 14,
  height = 8
)

plot(
  sampleTree,
  main = "Sample Clustering for Outlier Detection",
  sub = "",
  xlab = "",
  cex = 0.7
)

dev.off()


############################################################
# 7. CREATE SAMPLE METADATA
#
# Example sample:
# Cold_Tan23_1
#
# stress = Cold
# study = Tan23
# replicate = 1
############################################################

sample_info <- data.frame(
  Sample = rownames(input_mat),
  stringsAsFactors = FALSE
)

sample_split <- strsplit(
  sample_info$Sample,
  "_"
)

sample_info$stress <- sapply(
  sample_split,
  function(x) x[1]
)

sample_info$study <- sapply(
  sample_split,
  function(x) {
    if (length(x) >= 3) x[length(x) - 1] else NA
  }
)

sample_info$replicate <- sapply(
  sample_split,
  function(x) {
    if (length(x) >= 2) x[length(x)] else NA
  }
)

rownames(sample_info) <- sample_info$Sample

print(sample_info)


############################################################
# 8. CREATE BINARY STRESS TRAIT MATRIX
############################################################

stress_traits <- model.matrix(
  ~ 0 + stress,
  data = sample_info
)

# Remove "stress" prefix

colnames(stress_traits) <- sub(
  "^stress",
  "",
  colnames(stress_traits)
)

datTraits <- as.data.frame(
  stress_traits
)

# Ensure all columns are numeric

datTraits[] <- lapply(
  datTraits,
  as.numeric
)

rownames(datTraits) <- sample_info$Sample

# Match sample order

datTraits <- datTraits[
  match(
    rownames(input_mat),
    rownames(datTraits)
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(input_mat),
    rownames(datTraits)
  )
)

print(head(datTraits))


############################################################
# 9. SOFT-THRESHOLD POWER ANALYSIS
############################################################

powers <- c(
  1:10,
  seq(
    from = 12,
    to = 20,
    by = 2
  )
)

sft <- pickSoftThreshold(
  input_mat,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)


############################################################
# 10. PLOT SOFT-THRESHOLD RESULTS
############################################################

pdf(
  file.path(
    output_dir,
    "Soft_threshold.pdf"
  ),
  width = 12,
  height = 6
)

par(
  mfrow = c(1, 2)
)

cex1 <- 0.9

# Scale-free topology

plot(
  sft$fitIndices[, 1],
  -sign(
    sft$fitIndices[, 3]
  ) *
    sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale Independence"
)

text(
  sft$fitIndices[, 1],
  -sign(
    sft$fitIndices[, 3]
  ) *
    sft$fitIndices[, 2],
  labels = powers,
  cex = cex1
)

abline(
  h = 0.85,
  lty = 2
)

# Mean connectivity

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean Connectivity"
)

text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1
)

dev.off()


############################################################
# 11. SELECT SOFT POWER
#
# Check Soft_threshold.pdf before choosing this value
############################################################

picked_power <- 14


############################################################
# 12. RUN WGCNA
############################################################

temp_cor <- cor

cor <- WGCNA::cor

netwk <- blockwiseModules(
  
  input_mat,
  
  power = picked_power,
  
  networkType = "signed",
  
  TOMType = "signed",
  
  deepSplit = 2,
  
  pamRespectsDendro = FALSE,
  
  minModuleSize = 30,
  
  maxBlockSize = 4000,
  
  reassignThreshold = 0,
  
  mergeCutHeight = 0.25,
  
  saveTOMs = TRUE,
  
  saveTOMFileBase = file.path(
    output_dir,
    "TOM"
  ),
  
  numericLabels = TRUE,
  
  verbose = 3
)

# Restore normal cor function

cor <- temp_cor


############################################################
# 13. CONVERT MODULE LABELS TO COLORS
############################################################

mergedColors <- labels2colors(
  netwk$colors
)


############################################################
# 14. PLOT GENE DENDROGRAM
############################################################

pdf(
  file.path(
    output_dir,
    "Gene_dendrogram_modules.pdf"
  ),
  width = 16,
  height = 10
)

plotDendroAndColors(
  
  netwk$dendrograms[[1]],
  
  mergedColors[
    netwk$blockGenes[[1]]
  ],
  
  "Module colors",
  
  dendroLabels = FALSE,
  
  hang = 0.03,
  
  addGuide = TRUE,
  
  guideHang = 0.05
)

dev.off()


############################################################
# 15. CREATE GENE MODULE TABLE
############################################################

module_df <- data.frame(
  
  gene_id = colnames(input_mat),
  
  module = mergedColors,
  
  stringsAsFactors = FALSE
)

write.csv(
  module_df,
  file.path(
    output_dir,
    "Gene_modules.csv"
  ),
  row.names = FALSE
)


############################################################
# 16. COUNT GENES IN EACH MODULE
############################################################

module_sizes <- module_df %>%
  
  count(
    module,
    name = "n_genes"
  ) %>%
  
  arrange(
    desc(n_genes)
  )

print(module_sizes)

write.csv(
  module_sizes,
  file.path(
    output_dir,
    "Module_gene_counts.csv"
  ),
  row.names = FALSE
)


############################################################
# 17. CREATE CLEAN MODULE EIGENGENES
#
# IMPORTANT:
# Do NOT add treatment/sample names directly to MEs.
# MEs must remain completely numeric.
############################################################

MEs <- moduleEigengenes(
  
  input_mat,
  
  colors = mergedColors
  
)$eigengenes

MEs <- orderMEs(MEs)

# Match traits to eigengene samples

datTraits <- datTraits[
  match(
    rownames(MEs),
    rownames(datTraits)
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(MEs),
    rownames(datTraits)
  )
)


############################################################
# 18. MODULE EIGENGENE PLOT ACROSS SAMPLES
############################################################

module_order <- gsub(
  "^ME",
  "",
  colnames(MEs)
)

# Make a COPY for plotting

MEs_plot <- as.data.frame(MEs)

MEs_plot$treatment <- rownames(MEs_plot)

mME <- MEs_plot %>%
  
  pivot_longer(
    cols = -treatment,
    names_to = "module",
    values_to = "value"
  ) %>%
  
  mutate(
    
    module = gsub(
      "^ME",
      "",
      module
    ),
    
    module = factor(
      module,
      levels = module_order
    )
  )

p_ME <- mME %>%
  
  ggplot(
    aes(
      x = treatment,
      y = module,
      fill = value
    )
  ) +
  
  geom_tile() +
  
  theme_bw() +
  
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  ) +
  
  labs(
    title = "Module Eigengene Expression Across Samples",
    x = "Samples",
    y = "Modules",
    fill = "Eigengene"
  )

ggsave(
  file.path(
    output_dir,
    "Module_eigengene_samples.pdf"
  ),
  p_ME,
  width = 18,
  height = 10
)


############################################################
# 19. MODULE-STRESS CORRELATION
############################################################

moduleTraitCor <- WGCNA::cor(
  
  MEs,
  
  datTraits,
  
  use = "pairwise.complete.obs"
)

moduleTraitPvalue <- corPvalueStudent(
  
  moduleTraitCor,
  
  nSamples = nrow(input_mat)
)

write.csv(
  moduleTraitCor,
  file.path(
    output_dir,
    "Module_stress_correlations.csv"
  )
)

write.csv(
  moduleTraitPvalue,
  file.path(
    output_dir,
    "Module_stress_pvalues.csv"
  )
)


############################################################
# 20. PLOT ALL MODULE-STRESS RELATIONSHIPS
############################################################

textMatrix <- paste0(
  
  signif(
    moduleTraitCor,
    2
  ),
  
  "\n(",
  
  signif(
    moduleTraitPvalue,
    2
  ),
  
  ")"
)

dim(textMatrix) <- dim(moduleTraitCor)

pdf(
  file.path(
    output_dir,
    "All_Module_Stress_Heatmap.pdf"
  ),
  width = 14,
  height = 12
)

par(
  mar = c(10, 12, 4, 2)
)

labeledHeatmap(
  
  Matrix = moduleTraitCor,
  
  xLabels = colnames(datTraits),
  
  yLabels = rownames(moduleTraitCor),
  
  ySymbols = gsub(
    "^ME",
    "",
    rownames(moduleTraitCor)
  ),
  
  colorLabels = FALSE,
  
  colors = blueWhiteRed(50),
  
  textMatrix = textMatrix,
  
  setStdMargins = FALSE,
  
  cex.text = 0.6,
  
  zlim = c(-1, 1),
  
  main = "Module-Stress Relationships"
)

dev.off()


############################################################
# 21. SELECT TOP 10 STRESS-ASSOCIATED MODULES
#
# Modules are ranked by their strongest absolute correlation
# with ANY stress.
############################################################

module_max_cor <- apply(
  
  abs(moduleTraitCor),
  
  1,
  
  max,
  
  na.rm = TRUE
)

# Remove grey module if present

module_max_cor <- module_max_cor[
  !grepl(
    "^MEgrey$",
    names(module_max_cor)
  )
]

top10_ME_names <- names(
  
  sort(
    module_max_cor,
    decreasing = TRUE
  )[1:min(
    10,
    length(module_max_cor)
  )]
)

print(top10_ME_names)


############################################################
# 22. EXTRACT TOP 10 MODULES
############################################################

#moduleTraitCor_top10 <- moduleTraitCor[
  
#  top10_ME_names,
  
#  ,
  
#  drop = FALSE
#]

#moduleTraitPvalue_top10 <- moduleTraitPvalue[
  
#  top10_ME_names,
  
#  ,
  
#  drop = FALSE
#]
############################################################
# 21. SELECT TOP 10 MODULES BASED ON NUMBER OF GENES
############################################################

module_sizes_no_grey <- module_sizes %>%
  filter(module != "grey") %>%
  arrange(desc(n_genes))

top10_modules <- module_sizes_no_grey %>%
  slice_head(n = 10) %>%
  pull(module)

# Convert module colors to WGCNA eigengene names
top10_ME_names <- paste0("ME", top10_modules)

print("Top 10 modules based on gene number:")
print(module_sizes_no_grey %>% slice_head(n = 10))


############################################################
# 22. EXTRACT MODULE-STRESS CORRELATIONS
# FOR THE TOP 10 LARGEST MODULES
############################################################

moduleTraitCor_top10 <- moduleTraitCor[
  top10_ME_names,
  ,
  drop = FALSE
]

moduleTraitPvalue_top10 <- moduleTraitPvalue[
  top10_ME_names,
  ,
  drop = FALSE
]

############################################################
# 23. CREATE TOP 10 HEATMAP LABELS
############################################################

textMatrix_top10 <- paste0(
  
  signif(
    moduleTraitCor_top10,
    2
  ),
  
  "\n(",
  
  signif(
    moduleTraitPvalue_top10,
    2
  ),
  
  ")"
)

dim(textMatrix_top10) <- dim(
  moduleTraitCor_top10
)


############################################################
# 24. PLOT TOP 10 MODULE-STRESS HEATMAP
############################################################

pdf(
  file.path(
    output_dir,
    "Top10_Module_Stress_Heatmap.pdf"
  ),
  width = 14,
  height = 8
)

par(
  mar = c(10, 12, 4, 2)
)

labeledHeatmap(
  
  Matrix = moduleTraitCor_top10,
  
  xLabels = colnames(datTraits),
  
  yLabels = gsub(
    "^ME",
    "",
    rownames(moduleTraitCor_top10)
  ),
  
  ySymbols = gsub(
    "^ME",
    "",
    rownames(moduleTraitCor_top10)
  ),
  
  colorLabels = FALSE,
  
  colors = blueWhiteRed(50),
  
  textMatrix = textMatrix_top10,
  
  setStdMargins = FALSE,
  
  cex.text = 0.7,
  
  zlim = c(-1, 1),
  
  main = "Top 10 Stress-Associated WGCNA Modules"
)

dev.off()


############################################################
# 25. CALCULATE ADJACENCY MATRIX
#
# Used for intramodular connectivity.
############################################################

adjacency_matrix <- adjacency(
  
  input_mat,
  
  power = picked_power,
  
  type = "signed"
)


############################################################
# 26. CALCULATE INTRAMODULAR CONNECTIVITY
############################################################

connectivity <- intramodularConnectivity(
  
  adjacency_matrix,
  
  colors = mergedColors
)

connectivity$gene_id <- colnames(input_mat)

connectivity$module <- mergedColors


############################################################
# 27. CALCULATE MODULE MEMBERSHIP (kME)
############################################################

geneModuleMembership <- WGCNA::cor(
  
  input_mat,
  
  MEs,
  
  use = "pairwise.complete.obs"
)

# Get MM for the module to which each gene belongs

MM_own_module <- rep(
  NA,
  ncol(input_mat)
)

for (i in seq_len(ncol(input_mat))) {
  
  module_name <- mergedColors[i]
  
  ME_name <- paste0(
    "ME",
    module_name
  )
  
  if (ME_name %in% colnames(geneModuleMembership)) {
    
    MM_own_module[i] <-
      
      geneModuleMembership[
        i,
        ME_name
      ]
  }
}

connectivity$MM <- MM_own_module

connectivity$absMM <- abs(
  connectivity$MM
)


############################################################
# 28. SAVE COMPLETE CONNECTIVITY TABLE
############################################################

connectivity <- connectivity %>%
  
  select(
    gene_id,
    module,
    everything()
  )

write.csv(
  
  connectivity,
  
  file.path(
    output_dir,
    "All_gene_connectivity.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 29. IDENTIFY HUB GENES
#
# Definition:
# - Not grey
# - Top 10 genes by kWithin within each module
# - |MM| >= 0.8
############################################################

hub_genes <- connectivity %>%
  
  filter(
    module != "grey"
  ) %>%
  
  group_by(
    module
  ) %>%
  
  arrange(
    desc(kWithin),
    .by_group = TRUE
  ) %>%
  
  mutate(
    kWithin_rank = row_number()
  ) %>%
  
  filter(
    kWithin_rank <= 10,
    absMM >= 0.8
  ) %>%
  
  ungroup()

write.csv(
  
  hub_genes,
  
  file.path(
    output_dir,
    "Hub_genes.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 30. TOP 20 HUB CANDIDATES PER MODULE
############################################################

top20_hubs <- connectivity %>%
  
  filter(
    module != "grey"
  ) %>%
  
  group_by(
    module
  ) %>%
  
  arrange(
    desc(kWithin),
    .by_group = TRUE
  ) %>%
  
  slice_head(
    n = 20
  ) %>%
  
  ungroup()

write.csv(
  
  top20_hubs,
  
  file.path(
    output_dir,
    "Top20_hub_candidates_per_module.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 31. IDENTIFY EDGE / PERIPHERAL GENES
#
# Bottom 10% of intramodular connectivity
############################################################

edge_genes <- connectivity %>%
  
  filter(
    module != "grey"
  ) %>%
  
  group_by(
    module
  ) %>%
  
  mutate(
    
    kWithin_percentile = percent_rank(
      kWithin
    )
    
  ) %>%
  
  filter(
    kWithin_percentile <= 0.10
  ) %>%
  
  ungroup()

write.csv(
  
  edge_genes,
  
  file.path(
    output_dir,
    "Edge_peripheral_genes.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 32. CALCULATE GENE SIGNIFICANCE FOR ALL STRESSES
############################################################

geneTraitCor <- WGCNA::cor(
  
  input_mat,
  
  datTraits,
  
  use = "pairwise.complete.obs"
)

geneTraitPvalue <- corPvalueStudent(
  
  geneTraitCor,
  
  nSamples = nrow(input_mat)
)


############################################################
# 33. CREATE STRESS-SPECIFIC HUB GENES
#
# Criteria:
#
# |MM| >= 0.8
# |Gene-Stress Correlation| >= 0.5
############################################################

stress_hub_summary <- list()

for (stress_name in colnames(datTraits)) {
  
  stress_hubs <- connectivity
  
  stress_hubs$GS <-
    
    geneTraitCor[
      ,
      stress_name
    ]
  
  stress_hubs$GS_pvalue <-
    
    geneTraitPvalue[
      ,
      stress_name
    ]
  
  stress_hubs$absGS <- abs(
    stress_hubs$GS
  )
  
  stress_hubs <- stress_hubs %>%
    
    filter(
      
      module != "grey",
      
      absMM >= 0.8,
      
      absGS >= 0.5
      
    ) %>%
    
    arrange(
      
      desc(absGS),
      
      desc(kWithin)
      
    )
  
  stress_hub_summary[[stress_name]] <-
    
    stress_hubs
  
  write.csv(
    
    stress_hubs,
    
    file.path(
      
      output_dir,
      
      paste0(
        "Hub_genes_",
        stress_name,
        ".csv"
      )
      
    ),
    
    row.names = FALSE
  )
}


############################################################
# 34. SAVE COMPLETE R WORKSPACE
############################################################

save(
  
  data,
  
  input_mat,
  
  datTraits,
  
  sample_info,
  
  netwk,
  
  mergedColors,
  
  module_df,
  
  module_sizes,
  
  MEs,
  
  moduleTraitCor,
  
  moduleTraitPvalue,
  
  connectivity,
  
  hub_genes,
  
  edge_genes,
  
  geneTraitCor,
  
  geneTraitPvalue,
  
  file = file.path(
    output_dir,
    "WGCNA_complete_workspace.RData"
  )
)


############################################################
# 35. ANALYSIS COMPLETE
############################################################

cat("\n")
cat("============================================\n")
cat("WGCNA ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================\n")
cat("\n")

cat(
  "Results saved in:",
  output_dir,
  "\n"
)

cat(
  "Number of modules:",
  length(
    unique(
      mergedColors[
        mergedColors != "grey"
      ]
    )
  ),
  "\n"
)

cat(
  "Number of hub genes:",
  nrow(hub_genes),
  "\n"
)

cat(
  "Number of peripheral genes:",
  nrow(edge_genes),
  "\n"
)
############################################################
# 36. CREATE GENE CO-EXPRESSION NETWORK
############################################################

library(igraph)

# Select module
target_module <- "turquoise"

# Number of top connected genes
n_genes_network <- 25

# Edge threshold
edge_threshold <- 0.20


############################################################
# SELECT TOP CONNECTED GENES
############################################################

network_genes <- connectivity %>%
  
  filter(module == target_module) %>%
  
  arrange(desc(kWithin)) %>%
  
  slice_head(n = n_genes_network) %>%
  
  pull(gene_id)

cat(
  "Number of genes selected:",
  length(network_genes),
  "\n"
)


############################################################
# EXTRACT EXPRESSION
############################################################

network_expression <- input_mat[
  ,
  network_genes,
  drop = FALSE
]


############################################################
# CALCULATE ADJACENCY
############################################################

network_adjacency <- adjacency(
  
  network_expression,
  
  power = picked_power,
  
  type = "signed"
)

# Remove weak edges

network_adjacency[
  network_adjacency < edge_threshold
] <- 0

# Remove self-connections

diag(network_adjacency) <- 0


############################################################
# CREATE IGRAPH NETWORK
############################################################

g <- graph_from_adjacency_matrix(
  
  network_adjacency,
  
  mode = "undirected",
  
  weighted = TRUE,
  
  diag = FALSE
)


############################################################
# REMOVE ISOLATED GENES
############################################################

g <- delete_vertices(
  
  g,
  
  V(g)[degree(g) == 0]
)


############################################################
# ADD GENE CONNECTIVITY INFORMATION
############################################################

gene_info <- connectivity %>%
  
  filter(
    gene_id %in% V(g)$name
  )

V(g)$kWithin <- gene_info$kWithin[
  
  match(
    V(g)$name,
    gene_info$gene_id
  )
]


############################################################
# SAVE NETWORK EDGES
############################################################

edge_list <- as_data_frame(
  
  g,
  
  what = "edges"
)

write.csv(
  
  edge_list,
  
  file.path(
    output_dir,
    paste0(
      target_module,
      "_network_edges.csv"
    )
  ),
  
  row.names = FALSE
)


############################################################
# SAVE NETWORK NODES
############################################################

node_list <- data.frame(
  
  gene_id = V(g)$name,
  
  kWithin = V(g)$kWithin,
  
  degree = degree(g),
  
  stringsAsFactors = FALSE
)

node_list <- node_list %>%
  
  arrange(
    desc(kWithin)
  )

write.csv(
  
  node_list,
  
  file.path(
    output_dir,
    paste0(
      target_module,
      "_network_nodes.csv"
    )
  ),
  
  row.names = FALSE
)


############################################################
# PLOT NETWORK
############################################################

pdf(
  
  file.path(
    output_dir,
    paste0(
      target_module,
      "_gene_network.pdf"
    )
  ),
  
  width = 14,
  
  height = 14
)

set.seed(123)

plot(
  
  g,
  
  layout = layout_with_fr(g),
  
  vertex.size =
    4 +
    12 *
    (
      V(g)$kWithin /
        max(V(g)$kWithin)
    ),
  
  vertex.label.cex = 0.5,
  
  edge.width =
    E(g)$weight * 3,
  
  main = paste(
    "Gene Co-expression Network:",
    target_module
  )
)

dev.off()
