/*
===========================================================
QC SUBWORKFLOW
===========================================================
*/

include { FASTQC } from '../modules/qc/fastqc'
include { FASTP  } from '../modules/qc/fastp'

workflow QC {

    take:

    reads_ch

    main:

    raw_qc = FASTQC(reads_ch)

    trimmed = FASTP(reads_ch)

    emit:

    trimmed_reads = trimmed.reads

    fastqc_html = raw_qc.html

    fastqc_zip = raw_qc.zip

    fastp_html = trimmed.html

    fastp_json = trimmed.json

}