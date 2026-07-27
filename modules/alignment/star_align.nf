/*
===========================================================
STAR ALIGNMENT
===========================================================
*/

process STAR_ALIGN {

    tag "${sample_id}"

    label 'process_high'

    publishDir "${params.outdir}/alignment", mode: 'copy'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_1'

    input:

    tuple val(sample_id),
          val(species),
          path(read1),
          path(read2)

    path star_index

    output:

    tuple val(sample_id),
          val(species),
          path("${sample_id}.sorted.bam"),
          emit: bam

    path "${sample_id}.Log.final.out", emit: log

    script:

    """
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${star_index} \
        --readFilesIn ${read1} ${read2} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample_id}.

    mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.sorted.bam
    """
}