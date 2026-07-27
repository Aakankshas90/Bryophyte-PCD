#!/usr/bin/env Rscript

############################################################
# Create tx2gene table from GFF3/GTF annotation
#
# Project: Bryophyte-PCD
#
# Input:
#   GFF3 or GTF annotation file
#
# Output:
#   tx2gene.csv
#   tx2gene_summary.txt
#
############################################################

suppressPackageStartupMessages({

    library(rtracklayer)
    library(dplyr)

})


############################################################
# USER SETTINGS
############################################################

# Change only these two lines for each species

gff_file <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/references/Arabidopsis_thaliana/annotation.gff3"

output_dir <- "/Users/aakanksha/Desktop/github/Bryophyte-PCD/references/Arabidopsis_thaliana"


############################################################

cat("\n=========================================\n")
cat("MAKE tx2gene TABLE\n")
cat("Started:", Sys.time(), "\n")
cat("=========================================\n\n")


if(!file.exists(gff_file)){
    stop(
        "Annotation file not found:\n",
        gff_file
    )
}


dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


############################################################
# Import annotation
############################################################

cat("Reading annotation...\n")

annotation <- import(gff_file)

cat("Loaded",
    length(annotation),
    "features\n\n")


############################################################
# Keep transcript features
############################################################

cat("Searching for transcript features...\n")


transcripts <- annotation[
    annotation$type %in%
        c(
            "mRNA",
            "transcript",
            "RNA",
            "mrna"
        )
]


if(length(transcripts)==0){

    stop(
        "\nNo transcript features detected.\n",
        "Check annotation format.\n"
    )

}


cat(
    "Found",
    length(transcripts),
    "transcripts\n\n"
)



############################################################
# Convert to dataframe
############################################################

tx <- as.data.frame(transcripts)


############################################################
# Detect transcript identifier
############################################################

find_column <- function(columns, candidates){

    match <- candidates[
        candidates %in% columns
    ]

    if(length(match)==0){
        return(NULL)
    }

    return(match[1])

}



tx_id_col <- find_column(
    colnames(tx),
    c(
        "transcript_id",
        "transcript",
        "ID",
        "Name"
    )
)


gene_id_col <- find_column(
    colnames(tx),
    c(
        "gene_id",
        "gene",
        "Parent",
        "geneID",
        "gene_name"
    )
)



if(is.null(tx_id_col)){

    stop(
        "\nCannot identify transcript ID column.\n",
        "Available columns:\n",
        paste(colnames(tx), collapse=", ")
    )

}


if(is.null(gene_id_col)){

    stop(
        "\nCannot identify gene ID column.\n",
        "Available columns:\n",
        paste(colnames(tx), collapse=", ")
    )

}



cat("Transcript ID column:",
    tx_id_col,
    "\n")

cat("Gene ID column:",
    gene_id_col,
    "\n\n")



############################################################
# Create tx2gene
############################################################

tx2gene <- tx %>%

    dplyr::select(
        transcript = all_of(tx_id_col),
        gene = all_of(gene_id_col)
    ) %>%

    distinct()



############################################################
# Clean IDs
############################################################

clean_ids <- function(x){

    x <- gsub(
        "^(gene:|transcript:|rna-)",
        "",
        x
    )

    x <- gsub(
        ";.*$",
        "",
        x
    )

    return(x)

}



tx2gene$transcript <-
    clean_ids(tx2gene$transcript)


tx2gene$gene <-
    clean_ids(tx2gene$gene)



############################################################
# Remove missing values
############################################################

tx2gene <- tx2gene %>%

    filter(
        !is.na(transcript),
        !is.na(gene),
        transcript!="",
        gene!=""
    )



############################################################
# Summary
############################################################

n_transcripts <-
    nrow(tx2gene)

n_genes <-
    length(unique(tx2gene$gene))


cat("-----------------------------------------\n")
cat("SUMMARY\n")
cat("-----------------------------------------\n")

cat(
    "Transcripts:",
    n_transcripts,
    "\n"
)

cat(
    "Genes:",
    n_genes,
    "\n"
)



############################################################
# Save tx2gene
############################################################

outfile <-
    file.path(
        output_dir,
        "tx2gene.csv"
    )


write.csv(
    tx2gene,
    outfile,
    row.names = FALSE,
    quote = FALSE
)



############################################################
# Save summary
############################################################

summary_file <-
    file.path(
        output_dir,
        "tx2gene_summary.txt"
    )


sink(summary_file)

cat("tx2gene generation summary\n\n")

cat(
    "Annotation:\n",
    gff_file,
    "\n\n"
)

cat(
    "Transcripts:",
    n_transcripts,
    "\n"
)

cat(
    "Genes:",
    n_genes,
    "\n\n"
)

cat(
    "Generated:",
    Sys.time(),
    "\n"
)

sink()



cat("\nSaved:\n")
cat(outfile,"\n")
cat(summary_file,"\n")

cat("\nFinished successfully.\n")