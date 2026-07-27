/*
===========================================================
FEATURECOUNTS
===========================================================
*/

process FEATURECOUNTS {

    tag "${sample_id}"

    label 'process_medium'

    publishDir "${params.outdir}/counts", mode: 'copy'

    container 'quay.io/biocontainers/subread:2.1.1--he4a0461_0'

    input:

    tuple val(sample_id),
          val(species),
          path(bam),
          path(bai)

    path annotation_gtf

    output:

    tuple val(sample_id),
          val(species),
          path("${sample_id}.featureCounts.txt"),
          emit: counts

    path("${sample_id}.featureCounts.summary"), emit: summary

    script:

    """
    featureCounts \
        -T ${task.cpus} \
        -p \
        -a ${annotation_gtf} \
        -o ${sample_id}.featureCounts.txt \
        ${bam}
    """
}