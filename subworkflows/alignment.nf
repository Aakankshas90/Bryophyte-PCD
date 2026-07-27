/*
===========================================================
ALIGNMENT SUBWORKFLOW
===========================================================
*/

include { STAR_INDEX } from '../modules/alignment/star_index'
include { STAR_ALIGN } from '../modules/alignment/star_align'
include { SAMTOOLS } from '../modules/alignment/samtools'

workflow ALIGNMENT {

    take:

    reads_ch

    genome_fasta

    annotation_gtf

    main:

    star_index = STAR_INDEX(
        genome_fasta,
        annotation_gtf
    )

    aligned = STAR_ALIGN(
        reads_ch,
        star_index.index
    )

    bam = SAMTOOLS(
        aligned.bam
    )

    emit:

    bam_files = bam.bam

    flagstat = bam.flagstat

    idxstats = bam.idxstats

}