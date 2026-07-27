#!/usr/bin/env bash

###############################################################################
# Script Name : 04_multiqc.sh
# Description : Aggregate FastQC and fastp reports using MultiQC
###############################################################################

set -uo pipefail


########################################
############ USER SETTINGS #############
########################################

PROJECT_DIR="/Users/aakanksha/Desktop/github/Bryophyte-PCD"

DATA_DIR="$PROJECT_DIR/data/Marchantia_polymorpha"

########################################


RAW_FASTQC="$DATA_DIR/fastqc_raw"
TRIM_FASTQC="$DATA_DIR/fastqc_trimmed"
FASTP_REPORTS="$DATA_DIR/reports"

OUT_DIR="$DATA_DIR/multiqc"

LOG_DIR="$DATA_DIR/logs"

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"


LOGFILE="$LOG_DIR/04_multiqc.log"



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
########### CHECK MULTIQC ##############
########################################


if ! command -v multiqc >/dev/null 2>&1
then

    log "${RED}MultiQC not found in PATH.${RESET}"

    exit 1

fi



########################################
############ RESUME CHECK ##############
########################################


if [[ -f "$OUT_DIR/multiqc_report.html" ]]
then

    log "${YELLOW}MultiQC report already exists.${RESET}"
    log "${YELLOW}Skipping.${RESET}"

    exit 0

fi



########################################
############ RUN MULTIQC ##############
########################################


log "${BLUE}====================================================${RESET}"
log "${BLUE}MULTIQC STARTED${RESET}"
log "${BLUE}Started : $(timestamp)${RESET}"
log "${BLUE}====================================================${RESET}"



START=$(date +%s)



multiqc \
    "$RAW_FASTQC" \
    "$TRIM_FASTQC" \
    "$FASTP_REPORTS" \
    -o "$OUT_DIR" \
    >> "$LOGFILE" 2>&1



STATUS=$?



if [[ "$STATUS" -eq 0 && -f "$OUT_DIR/multiqc_report.html" ]]
then

    END=$(date +%s)

    log "${GREEN}MultiQC completed successfully.${RESET}"
    log "Time : $((END-START)) seconds"

else

    log "${RED}MultiQC failed.${RESET}"

    exit 1

fi



########################################
############ SUMMARY ###################
########################################


log ""
log "${BLUE}====================================================${RESET}"
log "${BLUE}OUTPUT${RESET}"
log "${BLUE}====================================================${RESET}"

log "Report:"
log "$OUT_DIR/multiqc_report.html"

log "Finished : $(timestamp)"
