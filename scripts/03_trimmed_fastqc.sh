#!/usr/bin/env bash

###############################################################################
# Script Name : 03_trimmed_fastqc.sh
# Description : Run FastQC on fastp-trimmed FASTQ files
# Supports    : Single-end and paired-end reads
###############################################################################

set -uo pipefail


########################################
############ USER SETTINGS #############
########################################

PROJECT_DIR="/Users/aakanksha/Desktop/github/Bryophyte-PCD"

DATA_DIR="$PROJECT_DIR/data/Arabidopsis"

THREADS=8

########################################


RAW_DIR="/Volumes/AakankshaHd/PCD/Marchantia/trimmed"
OUT_DIR="$DATA_DIR/fastqc_trimmed"

LOG_DIR="$DATA_DIR/logs"
DONE_DIR="$LOG_DIR/trimmed_fastqc_completed"


mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$DONE_DIR"


LOGFILE="$LOG_DIR/03_trimmed_fastqc.log"



########################################
############ COLORS ####################
########################################

RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
YELLOW=$(tput setaf 3)
BLUE=$(tput setaf 4)
CYAN=$(tput setaf 6)
RESET=$(tput sgr0)


########################################

timestamp() {
    date "+%Y-%m-%d %H:%M:%S"
}


log() {

    echo -e "$1"

    echo "$(timestamp) : $(echo -e "$1" | sed 's/\x1b\[[0-9;]*m//g')" >> "$LOGFILE"

}



########################################
########### CHECK FASTQC ###############
########################################


if ! command -v fastqc >/dev/null 2>&1
then
    log "${RED}FastQC not found in PATH.${RESET}"
    exit 1
fi



########################################
######## FIND FASTQ FILES ##############
########################################


FILES=("$RAW_DIR"/*.fastq.gz)


if [[ ! -e "${FILES[0]}" ]]
then
    log "${RED}No FASTQ files found in:${RESET}"
    log "$RAW_DIR"
    exit 1
fi



TOTAL=${#FILES[@]}


SUCCESS=0
SKIPPED=0
FAILED=0
COUNT=0


PIPE_START=$(date +%s)



log "${BLUE}====================================================${RESET}"
log "${BLUE}RAW FASTQC STARTED${RESET}"
log "${BLUE}Started : $(timestamp)${RESET}"
log "${BLUE}====================================================${RESET}"



########################################
############ PROCESS FILES ############
########################################


for FILE in "${FILES[@]}"
do

    BASENAME=$(basename "$FILE")


    # Skip R2 files because they are handled with R1
    if [[ "$BASENAME" =~ (_2|_R2)\.fastq\.gz$ ]]
    then
        continue
    fi



    COUNT=$((COUNT+1))


    ####################################
    # Detect paired-end
    ####################################


    if [[ "$BASENAME" =~ (_1|_R1)\.fastq\.gz$ ]]
    then

        if [[ "$BASENAME" =~ _1\.fastq\.gz$ ]]
        then
            SAMPLE=${BASENAME%_1.fastq.gz}
            R1="$FILE"
            R2="$RAW_DIR/${SAMPLE}_2.fastq.gz"

        else
            SAMPLE=${BASENAME%_R1.fastq.gz}
            R1="$FILE"
            R2="$RAW_DIR/${SAMPLE}_R2.fastq.gz"
        fi


        MODE="PE"


    else

        SAMPLE=${BASENAME%.fastq.gz}

        R1="$FILE"
        R2=""

        MODE="SE"

    fi



    DONE="$DONE_DIR/${SAMPLE}.done"



    log ""
    log "${CYAN}--------------------------------------------${RESET}"
    log "${CYAN}Sample ${COUNT}${RESET}"
    log "${CYAN}$SAMPLE ($MODE)${RESET}"
    log "${CYAN}--------------------------------------------${RESET}"



    if [[ -f "$DONE" ]]
    then
        log "${YELLOW}Already completed. Skipping.${RESET}"
        SKIPPED=$((SKIPPED+1))
        continue
    fi



    START=$(date +%s)



    ####################################
    # Run FastQC
    ####################################


    if [[ "$MODE" == "PE" ]]
    then

        if [[ ! -f "$R2" ]]
        then
            log "${RED}Missing mate:${RESET} $R2"
            FAILED=$((FAILED+1))
            continue
        fi


        fastqc \
            -t "$THREADS" \
            -o "$OUT_DIR" \
            "$R1" \
            "$R2" >> "$LOGFILE" 2>&1


        OUT1="$OUT_DIR/${SAMPLE}_1_fastqc.zip"
        OUT2="$OUT_DIR/${SAMPLE}_2_fastqc.zip"



        if [[ -f "$OUT1" && -f "$OUT2" ]]
        then
            touch "$DONE"
            SUCCESS=$((SUCCESS+1))
            log "${GREEN}Completed successfully.${RESET}"

        else
            FAILED=$((FAILED+1))
            log "${RED}FastQC output missing.${RESET}"
        fi



    else


        fastqc \
            -t "$THREADS" \
            -o "$OUT_DIR" \
            "$R1" >> "$LOGFILE" 2>&1


        OUT="$OUT_DIR/${SAMPLE}_fastqc.zip"


        if [[ -f "$OUT" ]]
        then
            touch "$DONE"
            SUCCESS=$((SUCCESS+1))
            log "${GREEN}Completed successfully.${RESET}"

        else
            FAILED=$((FAILED+1))
            log "${RED}FastQC output missing.${RESET}"
        fi


    fi



    END=$(date +%s)

    log "Time : $((END-START)) seconds"


done



########################################
############ SUMMARY ###################
########################################


PIPE_END=$(date +%s)


log ""
log "${BLUE}====================================================${RESET}"
log "${BLUE}SUMMARY${RESET}"
log "${BLUE}====================================================${RESET}"

log "Total samples : $COUNT"
log "Successful    : $SUCCESS"
log "Skipped       : $SKIPPED"
log "Failed        : $FAILED"
log "Runtime       : $((PIPE_END-PIPE_START)) seconds"

log "${BLUE}Finished : $(timestamp)${RESET}"