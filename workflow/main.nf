#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "data/test_fastq/*.fq"

process FASTQC {

    container 'biocontainers/fastqc:v0.11.9_cv8'

    publishDir "results/fastqc", mode: 'copy'

    input:
    path reads

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${reads}
    """
}

workflow {

    Channel
        .fromPath(params.reads)
        .set { read_files }

    FASTQC(read_files)
}