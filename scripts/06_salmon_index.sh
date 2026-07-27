#!/usr/bin/env bash

###############################################################################
# Script Name : 06_salmon_index.sh
# Description : Build a decoy-aware Salmon index
###############################################################################

set -uo pipefail

#######################################
########### USER SETTINGS #############
#######################################

PROJECT_DIR="/Users/aakanksha/Desktop/github/Bryophyte-PCD"

REF_DIR="$PROJECT_DIR/references/Arabidopsis_thaliana"

THREADS=8

#######################################

GENOME="$REF_DIR/genome.fa"
TRANSCRIPTS="$REF_DIR/transcripts.fa"

INDEX_DIR="$REF_DIR/salmon_index"

TMP_DIR="$REF_DIR/tmp"

LOG_DIR="$REF_DIR/logs"

mkdir -p "$INDEX_DIR"
mkdir -p "$TMP_DIR"
mkdir -p "$LOG_DIR"

LOGFILE="$LOG_DIR/06_salmon_index.log"

#######################################
############ COLORS ###################
#######################################

RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
YELLOW=$(tput setaf 3)
BLUE=$(tput setaf 4)
CYAN=$(tput setaf 6)
RESET=$(tput sgr0)

#######################################

timestamp(){

date "+%Y-%m-%d %H:%M:%S"

}

log(){

echo -e "$1"

echo "$(timestamp) : $(echo -e "$1" | sed 's/\x1b\[[0-9;]*m//g')" >> "$LOGFILE"

}

#######################################
######## DEPENDENCY CHECK #############
#######################################

for prog in salmon grep sed awk
do
    if ! command -v "$prog" >/dev/null
    then
        log "${RED}$prog not installed.${RESET}"
        exit 1
    fi
done

#######################################

log "${BLUE}====================================================${RESET}"
log "${BLUE}SALMON DECOY-AWARE INDEX${RESET}"
log "${BLUE}Started : $(timestamp)${RESET}"
log "${BLUE}====================================================${RESET}"

#######################################
########## RESUME CHECK ###############
#######################################

if [[ -f "$INDEX_DIR/index.ssi" ]]
then

    log "${GREEN}Index already exists.${RESET}"
    log "${GREEN}Skipping indexing.${RESET}"

    exit 0

fi

#######################################
######## FILE CHECK ###################
#######################################

for f in "$GENOME" "$TRANSCRIPTS"
do
    if [[ ! -f "$f" ]]
    then

        log "${RED}Missing:${RESET} $f"

        exit 1

    fi
done

#######################################
######## CREATE DECOY LIST ############
#######################################

log ""
log "${CYAN}Creating decoy list...${RESET}"

grep "^>" "$GENOME" \
| cut -d' ' -f1 \
| sed 's/>//' \
> "$TMP_DIR/decoys.txt"

#######################################
######## CONCATENATE FASTA ############
#######################################

log ""
log "${CYAN}Preparing gentrome...${RESET}"

cat "$TRANSCRIPTS" "$GENOME" > "$TMP_DIR/gentrome.fa"

#######################################
######## BUILD INDEX ##################
#######################################

START=$(date +%s)

log ""
log "${CYAN}Running Salmon index...${RESET}"

salmon index \
-t "$TMP_DIR/gentrome.fa" \
-d "$TMP_DIR/decoys.txt" \
-p "$THREADS" \
-i "$INDEX_DIR" \
2>&1 | tee -a "$LOGFILE"

STATUS=$?

END=$(date +%s)

#######################################
########### CHECK #####################
#######################################

if [[ $STATUS -ne 0 ]]
then

    log "${RED}Salmon indexing failed.${RESET}"

    exit 1

fi

if [[ ! -f "$INDEX_DIR/index.ssi" ]]
then

    log "${RED}Index incomplete.${RESET}"

    exit 1

fi

#######################################
############ SUMMARY ##################
#######################################

TIME=$((END-START))

log ""
log "${GREEN}Index built successfully.${RESET}"

log "Location : $INDEX_DIR"

log "Runtime  : ${TIME} seconds"

log "Finished : $(timestamp)"

log "${BLUE}====================================================${RESET}"