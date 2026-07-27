#!/usr/bin/env bash

###############################################################################
# Script Name : 07_salmon_quant.sh
#
# Description :
#   Transcript abundance quantification using Salmon
#
# Supports:
#   Single-end and paired-end FASTQ files
#
###############################################################################

set -euo pipefail


########################################
############ USER SETTINGS #############
########################################

THREADS=8

PROJECT_DIR="/Users/aakanksha/Desktop/github/Bryophyte-PCD"

TRIMMED_DIR="/Volumes/AakankshaHd/PCD/Marchantia/trimmed"

INDEX_DIR="${PROJECT_DIR}/references/Marchantia_polymorpha/salmon_index"

OUTPUT_DIR="/Volumes/AakankshaHd/PCD/Marchantia/trimmed/salmon"

LOG_DIR="${PROJECT_DIR}/data/Marchantia/logs"



########################################


mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"


LOGFILE="${LOG_DIR}/07_salmon_quant.log"


exec > >(tee -a "$LOGFILE")
exec 2>&1



echo "===================================================="
echo "SALMON QUANTIFICATION STARTED"
echo "Started : $(date)"
echo "===================================================="



########################################
############ CHECKS ####################
########################################


if ! command -v salmon >/dev/null 2>&1
then
    echo "ERROR: Salmon not found"
    exit 1
fi


if [ ! -d "$INDEX_DIR" ]
then
    echo "ERROR: Salmon index missing:"
    echo "$INDEX_DIR"
    exit 1
fi



########################################
############ FIND FASTQ ################
########################################


FILES=("$TRIMMED_DIR"/*.fastq.gz)


if [[ ! -e "${FILES[0]}" ]]
then
    echo "ERROR: No FASTQ files found"
    exit 1
fi



TOTAL=0
COMPLETED=0
SKIPPED=0
FAILED=0



########################################
############ MAIN LOOP #################
########################################


for FILE in "${FILES[@]}"
do

    BASENAME=$(basename "$FILE")



    ####################################
    # Skip R2 files
    ####################################

    if [[ "$BASENAME" =~ (_2|_R2)\.trim\.fastq\.gz$ ]]
    then
        continue
    fi



    ####################################
    # Detect paired-end
    ####################################


    MODE=""

    if [[ "$BASENAME" =~ _1\.trim\.fastq\.gz$ ]]
    then

        SAMPLE=${BASENAME%_1.trim.fastq.gz}

        R1="$FILE"
        R2="${TRIMMED_DIR}/${SAMPLE}_2.trim.fastq.gz"

        MODE="PE"



    elif [[ "$BASENAME" =~ _R1\.trim\.fastq\.gz$ ]]
    then

        SAMPLE=${BASENAME%_R1.trim.fastq.gz}

        R1="$FILE"
        R2="${TRIMMED_DIR}/${SAMPLE}_R2.trim.fastq.gz"

        MODE="PE"



    else

        SAMPLE=${BASENAME%.trim.fastq.gz}

        R1="$FILE"
        MODE="SE"

    fi



    OUTDIR="${OUTPUT_DIR}/${SAMPLE}"

    TOTAL=$((TOTAL+1))



    echo
    echo "-----------------------------------------"
    echo "Sample : $SAMPLE"
    echo "Mode   : $MODE"
    echo "-----------------------------------------"



    ####################################
    # Skip completed
    ####################################


    if [ -f "${OUTDIR}/quant.sf" ]
    then

        echo "Already completed. Skipping."

        SKIPPED=$((SKIPPED+1))

        continue

    fi



    mkdir -p "$OUTDIR"



    ####################################
    # Run Salmon
    ####################################


    if [ "$MODE" == "PE" ]
    then

        if [ ! -f "$R2" ]
        then
            echo "Missing R2:"
            echo "$R2"

            FAILED=$((FAILED+1))

            continue
        fi


        salmon quant \
            -i "$INDEX_DIR" \
            -l A \
            -1 "$R1" \
            -2 "$R2" \
            -p "$THREADS" \
            --validateMappings \
            --gcBias \
            -o "$OUTDIR"


    else


        salmon quant \
            -i "$INDEX_DIR" \
            -l A \
            -r "$R1" \
            -p "$THREADS" \
            --validateMappings \
            --gcBias \
            -o "$OUTDIR"


    fi



    if [ -f "${OUTDIR}/quant.sf" ]
    then

        echo "SUCCESS: $SAMPLE"

        COMPLETED=$((COMPLETED+1))

    else

        echo "FAILED: $SAMPLE"

        FAILED=$((FAILED+1))

    fi



done



echo
echo "===================================================="
echo "SALMON QUANTIFICATION COMPLETE"
echo "Finished : $(date)"
echo "===================================================="

echo "Total samples : $TOTAL"
echo "Completed     : $COMPLETED"
echo "Skipped       : $SKIPPED"
echo "Failed        : $FAILED"

echo "===================================================="