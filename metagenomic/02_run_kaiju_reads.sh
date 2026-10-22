#!/bin/bash
# ==============================================================================
# Script Name: 02_run_kaiju_reads.sh
# Author: Sofia Fonseca Yepez ann automated writing with AI
# Description: Automated and iterative pipeline for paired-end shotgun metagenomes:
#                   Classifier: Kaiju (greedy mode with sensitive parameters for eukaryotes).
#                   Database: nr_euk
# ==============================================================================

# --- 1. CONFIGURATION & DIRECTORIES ---

INPUT_DIR="/"
OUTPUT_DIR="/"
KAJU_DB="/home/sfonseca/kaiju_db_nr_euk.fmi"
NODES_DMP="/home/sfonseca/nodes.dmp"
THREADS=8

# --- Fin de la configuración ---

mkdir -p "${OUTPUT_DIR}"

# --- Lista de identificadores (SRR) ---

muestras=(
    "SRR11786911"
    "SRR12425620"
    "SRR12425871"
    "SRR12426445"
    "SRR12426459"
    "SRR12426461"
    "SRR12426463"
    "SRR1509790"
    "SRR1509794"
    "SRR1509796"
    "SRR1509797"
    "SRR1509799"
    "SRR1509800"
)

echo "PROCESSING SAMPLE ${#muestras[@]} WITH Kaiju..."
echo "Modo: Greedy with score 70 "
echo "Data base: ${KAJU_DB}"
echo "Output files: ${OUTPUT_DIR}"
echo "----------------------------------------"

# ---Process each sample ---

for muestra in "${muestras[@]}"; do
    echo "Process sample: ${muestra}"

    R1="${INPUT_DIR}/${muestra}_R1_paired.fq.gz"
    R2="${INPUT_DIR}/${muestra}_R2_paired.fq.gz"

    if [[ ! -f "$R1" ]] || [[ ! -f "$R2" ]]; then
        echo "  ¡WARNING! Files ${muestra} not found. Searching:"
        echo "    ${R1}"
        echo "    ${R2}"
        echo "  Omitting this sample."
        continue
    fi

    OUTPUT_FILE="${OUTPUT_DIR}/${muestra}_kaiju.out"

    
    echo "  Comando: kaiju -t ${NODES_DMP} -f ${KAJU_DB} -i ${R1} -j ${R2} -o ${OUTPUT_FILE} -z ${THREADS} -v -s 70 -E 0.01"

    kaiju -t "${NODES_DMP}" \
          -f "${KAJU_DB}" \
          -i "${R1}" \
          -j "${R2}" \
          -o "${OUTPUT_FILE}" \
          -z "${THREADS}" \
          -v \
          -s 70 \
          -E 0.01

    if [ $? -eq 0 ]; then
        echo "  Sample ${muestra} processed successfully. File: ${OUTPUT_FILE}"
    else
        echo "  Error processing the sample ${muestra}. Review the files and parameters."
    fi
    echo "----------------------------------------"

done

echo "Processing completed."
