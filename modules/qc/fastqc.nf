/*
===========================================================
FASTQC MODULE
Author: Aakanksha
Project: Bryophyte-PCD
===========================================================
*/

process FASTQC {

    tag "${sample_id}"

    label 'process_medium'

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:

    tuple val(sample_id), path(reads)

    output:

    tuple val(sample_id), path(reads), emit: reads

    path "*.html", emit: html

    path "*.zip", emit: zip

    script:

    """
    fastqc \
        --threads ${task.cpus} \
        ${reads}
    """
}