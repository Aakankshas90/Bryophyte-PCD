#!/usr/bin/env bash

###############################################################################
# Script Name : 02_fastp_trim.sh
# Description : Adapter trimming and quality filtering using fastp
# Supports    : Single-end and paired-end FASTQ files
# Project     : Bryophyte-PCD
###############################################################################

set -euo pipefail


########################################
############ USER SETTINGS #############
########################################

PROJECT_DIR="/Users/aakanksha/Desktop/github/Bryophyte-PCD"

DATA_DIR="$PROJECT_DIR/data/Arabidopsis"

RAW_DIR="$DATA_DIR/raw"

TRIM_DIR="/Volumes/AakankshaHd/PCD/Arabidopsis/trimmed"

REPORT_DIR="$DATA_DIR/reports"

LOG_DIR="$DATA_DIR/logs"

DONE_DIR="$LOG_DIR/fastp_completed"

THREADS=8

########################################


mkdir -p "$TRIM_DIR"
mkdir -p "$REPORT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$DONE_DIR"


LOGFILE="$LOG_DIR/02_fastp_trim.log"



########################################
############ FUNCTIONS #################
########################################


timestamp(){

    date "+%Y-%m-%d %H:%M:%S"

}



log(){

    echo "$(timestamp) : $1" | tee -a "$LOGFILE"

}



########################################
############ CHECK FASTP ###############
########################################


if ! command -v fastp >/dev/null 2>&1
then

    log "ERROR: fastp not found in PATH."

    exit 1

fi



########################################
############ FIND FILES ################
########################################


mapfile -t FILES < <(find "$RAW_DIR" -maxdepth 1 -name "*.fastq.gz" | sort)


if [[ ${#FILES[@]} -eq 0 ]]
then

    log "ERROR: No FASTQ files found in:"
    log "$RAW_DIR"

    exit 1

fi



########################################
############ COUNTERS #################
########################################


TOTAL=${#FILES[@]}

COUNT=0
SUCCESS=0
SKIPPED=0
FAILED=0



########################################
############ START ####################
########################################


log "===================================================="
log "FASTP TRIMMING STARTED"
log "Started : $(timestamp)"
log "Input   : $RAW_DIR"
log "Output  : $TRIM_DIR"
log "===================================================="



########################################
############ MAIN LOOP #################
########################################


for FILE in "${FILES[@]}"
do


    BASENAME=$(basename "$FILE")



    # Skip R2 because it is processed with R1

    if [[ "$BASENAME" =~ (_2|_R2)\.fastq\.gz$ ]]
    then
        continue
    fi



    COUNT=$((COUNT+1))



    ####################################
    # Detect paired-end or single-end
    ####################################


    if [[ "$BASENAME" =~ _1\.fastq\.gz$ ]]
    then

        SAMPLE=${BASENAME%_1.fastq.gz}

        R1="$FILE"

        R2="$RAW_DIR/${SAMPLE}_2.fastq.gz"

        MODE="PE"



    elif [[ "$BASENAME" =~ _R1\.fastq\.gz$ ]]
    then

        SAMPLE=${BASENAME%_R1.fastq.gz}

        R1="$FILE"

        R2="$RAW_DIR/${SAMPLE}_R2.fastq.gz"

        MODE="PE"



    else

        SAMPLE=${BASENAME%.fastq.gz}

        R1="$FILE"

        R2=""

        MODE="SE"


    fi



    DONE="$DONE_DIR/${SAMPLE}.done"



    OUT1="$TRIM_DIR/${SAMPLE}_1.trim.fastq.gz"

    OUT2="$TRIM_DIR/${SAMPLE}_2.trim.fastq.gz"

    OUTSE="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"



    HTML="$REPORT_DIR/${SAMPLE}.fastp.html"

    JSON="$REPORT_DIR/${SAMPLE}.fastp.json"



    log "--------------------------------------------"

    log "Sample $COUNT / $TOTAL"

    log "$SAMPLE ($MODE)"



    ####################################
    # Resume support
    ####################################


    if [[ -f "$DONE" ]]
    then

        log "Already completed. Skipping."

        SKIPPED=$((SKIPPED+1))

        continue

    fi



    START=$(date +%s)



    ####################################
    # Paired-end trimming
    ####################################


    if [[ "$MODE" == "PE" ]]
    then


        if [[ ! -f "$R2" ]]
        then

            log "ERROR: Missing R2 file:"
            log "$R2"

            FAILED=$((FAILED+1))

            continue

        fi



        log "Running paired-end fastp"



        fastp \
            -i "$R1" \
            -I "$R2" \
            -o "$OUT1" \
            -O "$OUT2" \
            --detect_adapter_for_pe \
            --thread "$THREADS" \
            --html "$HTML" \
            --json "$JSON"



    ####################################
    # Single-end trimming
    ####################################


    else


        log "Running single-end fastp"



        fastp \
            -i "$R1" \
            -o "$OUTSE" \
            --thread "$THREADS" \
            --html "$HTML" \
            --json "$JSON"



    fi



    END=$(date +%s)

    RUNTIME=$((END-START))



    ####################################
    # Completion marker
    ####################################


    touch "$DONE"


    SUCCESS=$((SUCCESS+1))


    log "Completed: $SAMPLE"

    log "Runtime: ${RUNTIME} seconds"



done



########################################
############ SUMMARY ###################
########################################


log "===================================================="

log "FASTP SUMMARY"

log "Processed : $COUNT"

log "Successful: $SUCCESS"

log "Skipped   : $SKIPPED"

log "Failed    : $FAILED"

log "Finished  : $(timestamp)"

log "===================================================="