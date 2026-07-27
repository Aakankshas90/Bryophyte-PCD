#!/usr/bin/env bash

###############################################################################
# Script Name : 00_validate_fastq.sh
# Description : Validate raw FASTQ files before RNA-seq processing
# Supports    : Single-end and paired-end reads
###############################################################################

set -uo pipefail


########################################
############ USER SETTINGS #############
########################################

PROJECT_DIR="/Users/aakanksha/Desktop/github/Bryophyte-PCD"

DATA_DIR="$PROJECT_DIR/data/Arabidopsis/Time-point-exp"

########################################


RAW_DIR="$DATA_DIR/raw"
LOG_DIR="$DATA_DIR/logs"

mkdir -p "$LOG_DIR"

LOGFILE="$LOG_DIR/00_validate_fastq.log"


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
########### CHECK FILES ################
########################################


FILES=("$RAW_DIR"/*.fastq.gz)


if [[ ! -e "${FILES[0]}" ]]
then
    log "${RED}No FASTQ files found in:${RESET}"
    log "$RAW_DIR"
    exit 1
fi



log "${BLUE}====================================================${RESET}"
log "${BLUE}FASTQ VALIDATION STARTED${RESET}"
log "${BLUE}Started : $(timestamp)${RESET}"
log "${BLUE}====================================================${RESET}"


TOTAL=0
PASSED=0
FAILED=0



########################################
######## FASTQ VALIDATION FUNCTION #####
########################################


validate_fastq () {

    FILE="$1"

    NAME=$(basename "$FILE")


    log "${CYAN}Checking $NAME${RESET}"


    ####################################
    # gzip test
    ####################################


    if gzip -t "$FILE" 2>/dev/null
    then
        log "${GREEN}gzip integrity OK${RESET}"
    else
        log "${RED}gzip integrity FAILED${RESET}"
        return 1
    fi



    ####################################
    # FASTQ format check
    ####################################


    awk '

    {

        line = NR % 4


        if(line==1 && substr($0,1,1)!="@")
            exit 1


        if(line==3 && substr($0,1,1)!="+")
            exit 2

    }


    END {

        if(NR % 4 != 0)
            exit 3

    }


    ' <(gzip -dc "$FILE")



    STATUS=$?


    case $STATUS in

        0)
            log "${GREEN}FASTQ format OK${RESET}"
            ;;

        1)
            log "${RED}FASTQ header error${RESET}"
            return 1
            ;;

        2)
            log "${RED}Missing '+' separator${RESET}"
            return 1
            ;;

        3)
            log "${RED}FASTQ line count invalid${RESET}"
            return 1
            ;;

        *)
            log "${RED}Unknown FASTQ error${RESET}"
            return 1
            ;;

    esac


    return 0

}



########################################
########### MAIN LOOP ##################
########################################


for FILE in "${FILES[@]}"
do


    BASENAME=$(basename "$FILE")


    # skip R2, processed with R1
    if [[ "$BASENAME" =~ (_2|_R2)\.fastq\.gz$ ]]
    then
        continue
    fi



    TOTAL=$((TOTAL+1))


    ####################################
    # Paired-end detection
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



    log ""
    log "${BLUE}--------------------------------------------${RESET}"
    log "${BLUE}$SAMPLE ($MODE)${RESET}"
    log "${BLUE}--------------------------------------------${RESET}"



    SAMPLE_OK=true



    ####################################
    # paired-end mate check
    ####################################


    if [[ "$MODE" == "PE" ]]
    then

        if [[ ! -f "$R2" ]]
        then
            log "${RED}Missing mate:${RESET} $R2"
            SAMPLE_OK=false
        fi

    fi



    ####################################
    # validate files
    ####################################


    if ! validate_fastq "$R1"
    then
        SAMPLE_OK=false
    fi



    if [[ "$MODE" == "PE" ]]
    then

        if ! validate_fastq "$R2"
        then
            SAMPLE_OK=false
        fi

    fi



    if $SAMPLE_OK
    then

        PASSED=$((PASSED+1))

        log "${GREEN}Sample PASSED${RESET}"

    else

        FAILED=$((FAILED+1))

        log "${RED}Sample FAILED${RESET}"

    fi



done



########################################
############ SUMMARY ###################
########################################


log ""
log "${BLUE}====================================================${RESET}"
log "${BLUE}SUMMARY${RESET}"
log "${BLUE}====================================================${RESET}"

log "Samples checked : $TOTAL"
log "Passed          : $PASSED"
log "Failed          : $FAILED"

log "${BLUE}Finished : $(timestamp)${RESET}"
