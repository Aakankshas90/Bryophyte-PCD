# ============================================================
# Annotate Marchantia polymorpha DESeq2 results
# Using MpTak_v6.1r2 annotation resources
#
# Inputs:
#   - DESeq2 result CSV files
#   - MpTak_v6.1r2 GFF
#   - MpTak_v6.1r2 KEGG
#   - MpTak_v6.1r2 KOG
#   - MpTak_v6.1r2 InterPro
#
# Output:
#   results/Marchantia_polymorpha/DESeq2/annotated/
#
# One row = one DEG/gene in one comparison.
# Multiple functional hits are collapsed with "; "
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  stop(
    "Package 'rtracklayer' is required.\n",
    "Install with:\n",
    "BiocManager::install('rtracklayer')"
  )
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Package 'data.table' is required.\n",
    "Install with:\n",
    "install.packages('data.table')"
  )
}

library(rtracklayer)
library(data.table)


# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

gff_file <- "references/Marchantia_polymorpha/MpTak_v6.1r2.gff"

kegg_file <- "references/Marchantia_polymorpha/MpTak_v6.1r2.KEGG.tsv"

kog_file <- "references/Marchantia_polymorpha/MpTak_v6.1r2.kog.summary.tsv"

interpro_file <- "references/Marchantia_polymorpha/MpTak_v6.1r2.interpro.tsv"


deg_dir <- "results/Marchantia_polymorpha/DESeq2"

output_dir <- file.path(
  deg_dir,
  "annotated"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

collapse_unique <- function(x) {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      nzchar(trimws(x))
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  x <- unique(x)
  
  paste(
    x,
    collapse = "; "
  )
}


# Convert transcript/protein ID to gene ID
#
# Examples:
# Mp1g00070.1 -> Mp1g00070
# Mp1g00070.2 -> Mp1g00070
#
# This is appropriate because the DEG results are gene-level,
# while KEGG/KOG/InterPro contain transcript/protein IDs.

to_gene_id <- function(x) {
  
  x <- as.character(x)
  
  sub(
    "\\.[0-9]+$",
    "",
    x
  )
}


# Decode URL-encoded GFF attributes such as:
# %2C -> ,
# %3B -> ;
#
# Protect against malformed values.

decode_gff <- function(x) {
  
  x <- as.character(x)
  
  out <- x
  
  ok <- !is.na(x)
  
  out[ok] <- vapply(
    x[ok],
    function(z) {
      tryCatch(
        utils::URLdecode(z),
        error = function(e) z
      )
    },
    character(1)
  )
  
  out
}


# ------------------------------------------------------------
# 3. Read GFF3
# ------------------------------------------------------------

cat("\n============================================\n")
cat("READING MARCHANTIA v6.1r2 GFF3\n")
cat("============================================\n")

gff <- readGFF(gff_file)

cat(
  "Total GFF records:",
  nrow(gff),
  "\n"
)

gene_gff <- gff[
  gff$type == "gene",
]

cat(
  "Gene records:",
  nrow(gene_gff),
  "\n"
)


# ------------------------------------------------------------
# 4. Extract gene-level annotation
# ------------------------------------------------------------

gff_fields <- c(
  "ID",
  "Name",
  "MapolyID",
  "symbol",
  "description",
  "product",
  "locus_type",
  "Note",
  "note",
  "Dbxref"
)

gff_fields <- intersect(
  gff_fields,
  colnames(gene_gff)
)

gene_annotation <- as.data.table(
  gene_gff[
    ,
    gff_fields,
    drop = FALSE
  ]
)


# Rename GFF ID to gene_id

setnames(
  gene_annotation,
  "ID",
  "gene_id"
)


# ------------------------------------------------------------
# 5. Clean GFF annotation columns
# ------------------------------------------------------------

for (col in colnames(gene_annotation)) {
  
  if (col != "gene_id") {
    
    if (is.list(gene_annotation[[col]])) {
      
      gene_annotation[
        ,
        (col) := vapply(
          get(col),
          collapse_unique,
          character(1)
        )
      ]
      
    }
    
    gene_annotation[
      ,
      (col) := decode_gff(
        get(col)
      )
    ]
  }
}


# Check duplicate gene IDs

if (anyDuplicated(gene_annotation$gene_id)) {
  
  duplicated_ids <- unique(
    gene_annotation$gene_id[
      duplicated(
        gene_annotation$gene_id
      )
    ]
  )
  
  stop(
    "Duplicate gene IDs found in GFF3:\n",
    paste(
      head(duplicated_ids, 20),
      collapse = ", "
    )
  )
}


cat(
  "Gene annotation table:",
  nrow(gene_annotation),
  "genes\n"
)

cat(
  "GFF annotation fields:\n"
)

print(
  colnames(gene_annotation)
)


# ------------------------------------------------------------
# 6. Read KEGG annotation
# ------------------------------------------------------------

cat("\n============================================\n")
cat("READING KEGG ANNOTATION\n")
cat("============================================\n")

kegg <- fread(
  kegg_file,
  sep = "\t",
  header = FALSE,
  quote = ""
)

cat(
  "KEGG rows:",
  nrow(kegg),
  "\n"
)

cat(
  "KEGG columns:",
  ncol(kegg),
  "\n"
)

print(
  head(kegg)
)


# The observed format is:
#
# Mp1g00070.1   KEGG   K10683   BRCA1-associated RING domain protein 1
#
# Therefore:

if (ncol(kegg) < 4) {
  
  stop(
    "Unexpected KEGG file format: fewer than 4 columns."
  )
}

kegg <- kegg[
  ,
  1:4
]

colnames(kegg) <- c(
  "transcript_id",
  "database",
  "KEGG_ID",
  "KEGG_description"
)


# Convert transcript ID -> gene ID

kegg[
  ,
  gene_id := to_gene_id(transcript_id)
]


# Collapse multiple KEGG hits per gene

kegg_gene <- kegg[
  ,
  .(
    KEGG_ID = collapse_unique(KEGG_ID),
    KEGG_description = collapse_unique(
      KEGG_description
    )
  ),
  by = gene_id
]


cat(
  "Genes with KEGG annotation:",
  nrow(kegg_gene),
  "\n"
)


# ------------------------------------------------------------
# 7. Read KOG annotation
# ------------------------------------------------------------

cat("\n============================================\n")
cat("READING KOG ANNOTATION\n")
cat("============================================\n")

# The KOG file has a header beginning with "#".
# Do NOT use comment.char="#" because that removes the header.
# Instead, read the file normally and remove the leading "#"
# from the first column name.

kog <- fread(
  kog_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  check.names = FALSE
)

# Remove leading "#" from the first column name

colnames(kog)[1] <- sub(
  "^#\\s*",
  "",
  colnames(kog)[1]
)

cat(
  "KOG rows:",
  nrow(kog),
  "\n"
)

cat(
  "KOG columns:\n"
)

print(
  colnames(kog)
)

print(
  head(kog)
)


# Check expected columns

required_kog <- c(
  "query_id",
  "accession",
  "definition",
  "category_codes",
  "categories"
)

missing_kog <- setdiff(
  required_kog,
  colnames(kog)
)

if (length(missing_kog) > 0) {
  
  stop(
    "KOG file is missing expected columns: ",
    paste(
      missing_kog,
      collapse = ", "
    )
  )
}


# Convert transcript/protein ID -> gene ID

kog[
  ,
  gene_id := to_gene_id(query_id)
]


# Collapse multiple transcript/KOG hits per gene

kog_gene <- kog[
  ,
  .(
    KOG_ID = collapse_unique(
      accession
    ),
    
    KOG_description = collapse_unique(
      definition
    ),
    
    KOG_category_code = collapse_unique(
      category_codes
    ),
    
    KOG_category = collapse_unique(
      categories
    )
  ),
  by = gene_id
]


cat(
  "Genes with KOG annotation:",
  nrow(kog_gene),
  "\n"
)
# ------------------------------------------------------------
# 8. Read InterPro annotation
# ------------------------------------------------------------

cat("\n============================================\n")
cat("READING INTERPRO ANNOTATION\n")
cat("============================================\n")

interpro <- fread(
  interpro_file,
  sep = "\t",
  header = FALSE,
  quote = ""
)

cat(
  "InterPro rows:",
  nrow(interpro),
  "\n"
)

cat(
  "InterPro columns:",
  ncol(interpro),
  "\n"
)

print(
  head(interpro)
)


# Your observed file has:
#
# Mp2g16790.1   Coils   Coil    Coil
#
# Therefore:
#
# column 1 = transcript/protein ID
# column 2 = annotation/database
# column 3 = accession
# column 4 = description


if (ncol(interpro) < 4) {
  
  stop(
    "Unexpected InterPro file format: fewer than 4 columns."
  )
}

interpro <- interpro[
  ,
  1:4
]

colnames(interpro) <- c(
  "transcript_id",
  "InterPro_database",
  "InterPro_ID",
  "InterPro_description"
)


interpro[
  ,
  gene_id := to_gene_id(
    transcript_id
  )
]


# ------------------------------------------------------------
# 9. Collapse InterPro hits per gene
# ------------------------------------------------------------

interpro_gene <- interpro[
  ,
  .(
    InterPro_database =
      collapse_unique(
        InterPro_database
      ),
    
    InterPro_ID =
      collapse_unique(
        InterPro_ID
      ),
    
    InterPro_description =
      collapse_unique(
        InterPro_description
      )
  ),
  by = gene_id
]


cat(
  "Genes with InterPro annotation:",
  nrow(interpro_gene),
  "\n"
)


# ------------------------------------------------------------
# 10. Create unified functional annotation table
# ------------------------------------------------------------

cat("\n============================================\n")
cat("MERGING FUNCTIONAL ANNOTATIONS\n")
cat("============================================\n")


functional_annotation <- merge(
  kegg_gene,
  kog_gene,
  by = "gene_id",
  all = TRUE
)

functional_annotation <- merge(
  functional_annotation,
  interpro_gene,
  by = "gene_id",
  all = TRUE
)


cat(
  "Genes in functional annotation table:",
  nrow(functional_annotation),
  "\n"
)


# ------------------------------------------------------------
# 11. Merge GFF + functional annotation
# ------------------------------------------------------------

annotation <- merge(
  gene_annotation,
  functional_annotation,
  by = "gene_id",
  all = TRUE
)


cat(
  "Total annotated gene records:",
  nrow(annotation),
  "\n"
)


# ------------------------------------------------------------
# 12. Save complete annotation table
# ------------------------------------------------------------

annotation_file <- file.path(
  output_dir,
  "Marchantia_v6.1r2_gene_annotation.tsv"
)

fwrite(
  annotation,
  annotation_file,
  sep = "\t",
  na = ""
)

cat(
  "\nSaved complete annotation table:\n",
  annotation_file,
  "\n"
)


# ------------------------------------------------------------
# 13. Find all DEG files
# ------------------------------------------------------------

deg_files <- list.files(
  deg_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)


# Don't accidentally process files already inside annotated/

deg_files <- deg_files[
  !grepl(
    "/annotated/",
    deg_files
  )
]


if (length(deg_files) == 0) {
  
  stop(
    "No CSV DEG files found in:\n",
    deg_dir
  )
}


cat("\n============================================\n")
cat("DEG FILES FOUND\n")
cat("============================================\n")

print(
  basename(deg_files)
)


# ------------------------------------------------------------
# 14. Annotate each DEG file
# ------------------------------------------------------------

summary_list <- vector(
  "list",
  length(deg_files)
)


all_annotated <- vector(
  "list",
  length(deg_files)
)


for (i in seq_along(deg_files)) {
  
  deg_file <- deg_files[i]
  
  cat("\n--------------------------------------------\n")
  cat(
    "Processing:",
    basename(deg_file),
    "\n"
  )
  cat("--------------------------------------------\n")
  
  
  # --------------------------------------------------------
  # Read DEG file
  # --------------------------------------------------------
  
  deg <- fread(
    deg_file
  )
  
  
  if (!"gene_id" %in% colnames(deg)) {
    
    stop(
      "File does not contain 'gene_id': ",
      deg_file
    )
  }
  
  
  cat(
    "DEG rows:",
    nrow(deg),
    "\n"
  )
  
  
  # --------------------------------------------------------
  # Check duplicate DEG IDs
  # --------------------------------------------------------
  
  if (anyDuplicated(deg$gene_id)) {
    
    warning(
      "Duplicate gene IDs detected in ",
      basename(deg_file),
      ". ",
      "Annotation will preserve all DEG rows."
    )
  }
  
  
  # --------------------------------------------------------
  # Merge annotation
  # --------------------------------------------------------
  
  annotated <- merge(
    deg,
    annotation,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  
  # --------------------------------------------------------
  # Add comparison name
  # --------------------------------------------------------
  
  comparison <- tools::file_path_sans_ext(
    basename(deg_file)
  )
  
  
  annotated[
    ,
    comparison := comparison
  ]
  
  
  # Put comparison after gene_id
  
  setcolorder(
    annotated,
    c(
      "gene_id",
      "comparison",
      setdiff(
        colnames(annotated),
        c(
          "gene_id",
          "comparison"
        )
      )
    )
  )
  
  
  # --------------------------------------------------------
  # Annotation coverage
  # --------------------------------------------------------
  
  n_total <- nrow(annotated)
  
  n_gff <- sum(
    !is.na(
      annotated$MapolyID
    )
  )
  
  n_kegg <- sum(
    !is.na(
      annotated$KEGG_ID
    )
  )
  
  n_kog <- sum(
    !is.na(
      annotated$KOG_ID
    )
  )
  
  n_interpro <- sum(
    !is.na(
      annotated$InterPro_ID
    )
  )
  
  
  # --------------------------------------------------------
  # Save annotated DEG file
  # --------------------------------------------------------
  
  output_file <- file.path(
    output_dir,
    paste0(
      comparison,
      "_annotated.tsv"
    )
  )
  
  
  fwrite(
    annotated,
    output_file,
    sep = "\t",
    na = ""
  )
  
  
  cat(
    "Saved:",
    output_file,
    "\n"
  )
  
  
  # --------------------------------------------------------
  # Store summary
  # --------------------------------------------------------
  
  summary_list[[i]] <- data.table(
    comparison = comparison,
    DEG_rows = n_total,
    GFF_annotated = n_gff,
    GFF_percent = round(
      100 * n_gff / n_total,
      2
    ),
    KEGG_annotated = n_kegg,
    KEGG_percent = round(
      100 * n_kegg / n_total,
      2
    ),
    KOG_annotated = n_kog,
    KOG_percent = round(
      100 * n_kog / n_total,
      2
    ),
    InterPro_annotated = n_interpro,
    InterPro_percent = round(
      100 * n_interpro / n_total,
      2
    )
  )
  
  
  all_annotated[[i]] <- annotated
}


# ------------------------------------------------------------
# 15. Combine all annotated DEG tables
# ------------------------------------------------------------

cat("\n============================================\n")
cat("COMBINING ALL ANNOTATED DEG FILES\n")
cat("============================================\n")


combined_deg <- rbindlist(
  all_annotated,
  fill = TRUE,
  use.names = TRUE
)


combined_file <- file.path(
  output_dir,
  "ALL_Marchantia_DEGs_annotated.tsv"
)


fwrite(
  combined_deg,
  combined_file,
  sep = "\t",
  na = ""
)


cat(
  "Combined table saved:\n",
  combined_file,
  "\n"
)

cat(
  "Total rows:",
  nrow(combined_deg),
  "\n"
)


# ------------------------------------------------------------
# 16. Save annotation coverage summary
# ------------------------------------------------------------

annotation_summary <- rbindlist(
  summary_list,
  fill = TRUE
)


summary_file <- file.path(
  output_dir,
  "annotation_coverage_summary.tsv"
)


fwrite(
  annotation_summary,
  summary_file,
  sep = "\t",
  na = ""
)


# ------------------------------------------------------------
# 17. Overall annotation summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("ANNOTATION COVERAGE SUMMARY\n")
cat("============================================\n")

print(
  annotation_summary
)


# ------------------------------------------------------------
# 18. Final message
# ------------------------------------------------------------

cat("\n============================================\n")
cat("ANNOTATION COMPLETE\n")
cat("============================================\n")

cat(
  "Annotated DEG files:",
  length(deg_files),
  "\n"
)

cat(
  "Output directory:\n",
  output_dir,
  "\n"
)

cat(
  "\nMain combined file:\n",
  combined_file,
  "\n"
)

cat(
  "\nCoverage summary:\n",
  summary_file,
  "\n"
)

cat("\nDone.\n")

