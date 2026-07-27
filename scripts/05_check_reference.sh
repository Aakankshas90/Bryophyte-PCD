#!/usr/bin/env bash

###############################################################################
# Script Name : 05_check_reference.sh
# Description : Verify genome, transcriptome and annotation before Salmon
###############################################################################

set -uo pipefail

########################################
############ USER SETTINGS #############
########################################

PROJECT_DIR="/Users/aakanksha/Desktop/github/Bryophyte-PCD"

REF_DIR="$PROJECT_DIR/references/Marchantia_polymorpha"

########################################

GENOME="$REF_DIR/genome.fa"
TRANSCRIPTS="$REF_DIR/transcripts.fa"
ANNOTATION="$REF_DIR/annotation.gff3"

LOG_DIR="$REF_DIR/logs"

mkdir -p "$LOG_DIR"

LOGFILE="$LOG_DIR/05_reference_check.log"

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
########### DEPENDENCY CHECK ###########
########################################

for prog in seqkit gffread awk grep
do
    if ! command -v "$prog" >/dev/null 2>&1
    then
        log "${RED}$prog not found.${RESET}"
        exit 1
    fi
done

########################################
############ FILE CHECK ################
########################################

log "${BLUE}====================================================${RESET}"
log "${BLUE}REFERENCE VALIDATION${RESET}"
log "${BLUE}Started : $(timestamp)${RESET}"
log "${BLUE}====================================================${RESET}"

for f in "$GENOME" "$TRANSCRIPTS" "$ANNOTATION"
do
    if [[ ! -f "$f" ]]
    then
        log "${RED}Missing:${RESET} $f"
        exit 1
    fi
done

log "${GREEN}All reference files found.${RESET}"

########################################
########### GENOME STATS ###############
########################################

log ""
log "${CYAN}Genome statistics${RESET}"

seqkit stats "$GENOME" | tee -a "$LOGFILE"

########################################
######## TRANSCRIPT STATS ##############
########################################

log ""
log "${CYAN}Transcript statistics${RESET}"

seqkit stats "$TRANSCRIPTS" | tee -a "$LOGFILE"

########################################
########### GFF SUMMARY ################
########################################

log ""
log "${CYAN}Annotation summary${RESET}"

GENES=$(grep -cv "^#" "$ANNOTATION")
MRNA=$(grep -c $'\tmRNA\t' "$ANNOTATION" || true)
EXONS=$(grep -c $'\texon\t' "$ANNOTATION" || true)
CDS=$(grep -c $'\tCDS\t' "$ANNOTATION" || true)

log "Total annotation entries : $GENES"
log "mRNA features            : $MRNA"
log "Exons                    : $EXONS"
log "CDS                      : $CDS"

########################################
###### CHROMOSOME NAME CHECK ###########
########################################

log ""
log "${CYAN}Checking chromosome names...${RESET}"

grep "^>" "$GENOME" \
| sed 's/>//' \
| awk '{print $1}' \
| sort -u > /tmp/genome_names.txt

grep -v "^#" "$ANNOTATION" \
| cut -f1 \
| sort -u > /tmp/gff_names.txt

MISSING=$(comm -23 /tmp/gff_names.txt /tmp/genome_names.txt | wc -l)

if [[ "$MISSING" -eq 0 ]]
then
    log "${GREEN}Genome and annotation chromosome names match.${RESET}"
else
    log "${RED}$MISSING chromosome names present in GFF but absent in genome.${RESET}"
    comm -23 /tmp/gff_names.txt /tmp/genome_names.txt | tee -a "$LOGFILE"
fi

########################################
###### TRANSCRIPT ID CHECK #############
########################################

log ""
log "${CYAN}Checking transcript IDs...${RESET}"

TRANS_FASTA=$(grep "^>" "$TRANSCRIPTS" | wc -l)

TRANS_GFF=$(grep $'\tmRNA\t' "$ANNOTATION" | wc -l)

log "Transcript FASTA entries : $TRANS_FASTA"
log "mRNA entries in GFF      : $TRANS_GFF"

########################################
########### GFF PARSE TEST #############
########################################

log ""
log "${CYAN}Testing annotation integrity...${RESET}"

if gffread "$ANNOTATION" -T -o /dev/null >/dev/null 2>&1
then
    log "${GREEN}GFF successfully parsed.${RESET}"
else
    log "${RED}GFF parsing failed.${RESET}"
    exit 1
fi

########################################
############ CLEANUP ###################
########################################

rm -f /tmp/genome_names.txt
rm -f /tmp/gff_names.txt

########################################
############ SUMMARY ###################
########################################

log ""
log "${BLUE}====================================================${RESET}"
log "${GREEN}Reference check completed.${RESET}"
log "Finished : $(timestamp)"
log "${BLUE}====================================================${RESET}"
