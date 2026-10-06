# Comparative transcriptomic analysis of stress responses and programmed cell death-associated genes in *Marchantia polymorpha*

This repository contains the scripts, metadata, and reference files for the comparative reanalysis of publicly available RNA-seq datasets from *Marchantia polymorpha* exposed to diverse abiotic and biotic stresses.

The study integrates transcriptomic data from three independent studies and evaluates differential gene expression across nine stress-versus-control comparisons. Particular focus was placed on a curated set of 14 *M. polymorpha* programmed cell death (PCD)-associated genes, recurrent differential expression across stress conditions, and stress-associated gene co-expression modules identified using weighted gene co-expression network analysis (WGCNA).

## Main analyses

* RNA-seq quantification and transcript-level summarisation using Salmon and tximport
* Differential expression analysis using DESeq2 and apeglm
* Principal component analysis and sample-level assessment
* Identification and comparison of DEGs across nine stress conditions
* Analysis of the 34 genes shared across all nine comparisons, including expression direction and magnitude
* Analysis of 14 curated *M. polymorpha* PCD-associated genes
* WGCNA to identify stress-associated co-expression modules
* Functional interpretation using Gene Ontology enrichment
* Genome annotation and functional annotation of *M. polymorpha* genes

## Repository structure

```text
metadata/      Sample sheets and gene annotation files
references/    Reference genomes, and annotations
scripts/       R and shell scripts for individual analyses
```

## Reproducibility

Sample metadata and analysis scripts are retained in the repository to document the computational workflow.

## Data sources

All transcriptomic datasets analysed in this study were obtained from publicly available studies. The corresponding study-specific metadata and sample information are provided in `metadata/`.
