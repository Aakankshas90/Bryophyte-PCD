/*
===========================================================
STAR GENOME INDEX
===========================================================
*/

process STAR_INDEX {

    tag "STAR_INDEX"

    label 'process_high'

    publishDir "${params.outdir}/star_index", mode: 'copy'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_1'

    input:

    path genome_fasta

    path genome_gtf

    output:

    path "STAR_index", emit: index

    script:

    """
    mkdir STAR_index

    STAR \
        --runMode genomeGenerate \
        --runThreadN ${task.cpus} \
        --genomeDir STAR_index \
        --genomeFastaFiles ${genome_fasta} \
        --sjdbGTFfile ${genome_gtf}
    """
}