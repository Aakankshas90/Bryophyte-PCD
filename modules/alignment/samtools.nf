/*
===========================================================
SAMTOOLS
Sort, Index and Alignment Statistics
===========================================================
*/

process SAMTOOLS {

    tag "${sample_id}"

    label 'process_high'

    publishDir "${params.outdir}/bam", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    input:

    tuple val(sample_id),
          val(species),
          path(bam)

    output:

    tuple val(sample_id),
          val(species),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai"),
          emit: bam

    path("${sample_id}.flagstat.txt"), emit: flagstat

    path("${sample_id}.idxstats.txt"), emit: idxstats

    script:

    """
    samtools sort \
        -@ ${task.cpus} \
        -o ${sample_id}.sorted.bam \
        ${bam}

    samtools index \
        ${sample_id}.sorted.bam

    samtools flagstat \
        ${sample_id}.sorted.bam \
        > ${sample_id}.flagstat.txt

    samtools idxstats \
        ${sample_id}.sorted.bam \
        > ${sample_id}.idxstats.txt
    """
}