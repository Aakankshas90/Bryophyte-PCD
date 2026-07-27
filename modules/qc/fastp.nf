/*
===========================================================
FASTP MODULE
Author: Aakanksha S
Project: Bryophyte-PCD
===========================================================
*/

process FASTP {

    tag "${sample_id}"

    label 'process_medium'

    publishDir "${params.outdir}/fastp", mode: 'copy'

    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input:

    tuple val(sample_id), path(read1), path(read2)

    output:

    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fastq.gz"),
          path("${sample_id}_R2.trimmed.fastq.gz"),
          emit: reads

    path "*.html", emit: html

    path "*.json", emit: json

    script:

    """
    fastp \
        --thread ${task.cpus} \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_R1.trimmed.fastq.gz \
        --out2 ${sample_id}_R2.trimmed.fastq.gz \
        --html ${sample_id}.fastp.html \
        --json ${sample_id}.fastp.json
    """
}