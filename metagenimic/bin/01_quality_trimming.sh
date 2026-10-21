#!/bin/bash
# ==============================================================================
# Script Name: 01_quality_trimming.sh
# Author: Sofia Fonseca Yepez ann automated writing with AI Gemini
# Description: Automated and iterative pipeline for paired-end shotgun metagenomes:
#              1. Pre-trimming Quality Control (FastQC + MultiQC)
#              2. Conservative Adapter & Quality Trimming (Trimmomatic v0.39)
#              3. Post-trimming Quality Control (FastQC + MultiQC)
# ==============================================================================

# ------------------------------
# 1. CONFIGURATION & DIRECTORIES
# ------------------------------
THREADS=4
RAW_DIR="."           # Directory containing raw *_1.fastq.gz and *_2.fastq.gz
ADAPTERS="/"          # Path to Illumina TruSeq3-PE adapters

# Output directories
OUTPUT_DIR="trimmo_kaiju_serial"
FASTQC_RAW_DIR="${OUTPUT_DIR}/01_fastqc_raw"
MULTIQC_RAW_DIR="${OUTPUT_DIR}/01_multiqc_raw"
TRIMMED_DIR="${OUTPUT_DIR}/02_trimmed_reads"
FASTQC_TRIM_DIR="${OUTPUT_DIR}/03_fastqc_trimmed"
MULTIQC_TRIM_DIR="${OUTPUT_DIR}/03_multiqc_trimmed"
LOG_DIR="${OUTPUT_DIR}/logs_pipeline_$(date +%Y%m%d_%H%M%S)"
SUMMARY_DIR="${OUTPUT_DIR}/resumen"

# Create directory structure
mkdir -p "$FASTQC_RAW_DIR" "$MULTIQC_RAW_DIR" "$TRIMMED_DIR" \
         "$FASTQC_TRIM_DIR" "$MULTIQC_TRIM_DIR" "$LOG_DIR" "$SUMMARY_DIR"

# Complete list of the 13 ETNP-OMZ metagenomic samples (PRJNA632347 & PRJNA254808)
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

# Summary status file
RESUMEN="${SUMMARY_DIR}/estado_muestras_pipeline.txt"
echo "SAMPLE | PRE_FASTQC | TRIMMOMATIC | POST_FASTQC | R1_PAIRED_SIZE | R2_PAIRED_SIZE | TIME" > "$RESUMEN"
echo "-------|------------|-------------|-------------|----------------|----------------|------" >> "$RESUMEN"

echo "============================================================"
echo "STARTING INTEGRATED QC & TRIMMOMATIC PIPELINE (SERIAL)"
echo "============================================================"
echo "✅ FastQC:      $(which fastqc)"
echo "✅ Trimmomatic: $(which trimmomatic)"
echo "✅ MultiQC:     $(which multiqc)"
echo "✅ Threads:     $THREADS"
echo "✅ Adapters:    $ADAPTERS"
echo "✅ Output Dir:  $OUTPUT_DIR"
echo "============================================================"

total_muestras=${#muestras[@]}
actual=0
exitosas=0
fallidas=0

# ==============================================================================
# 2. ITERATIVE PROCESSING PER SAMPLE (FASTQC PRE -> TRIMMOMATIC -> FASTQC POST)
# ==============================================================================
for muestra in "${muestras[@]}"; do
    actual=$((actual + 1))
    start_time=$(date +%s)

    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  PROCESSING SAMPLE $actual of $total_muestras: $muestra"
    echo "╚════════════════════════════════════════════════════════════╝"

    R1_RAW="${RAW_DIR}/${muestra}_1.fastq.gz"
    R2_RAW="${RAW_DIR}/${muestra}_2.fastq.gz"

    # Verify input files exist
    if [[ ! -f "$R1_RAW" ]] || [[ ! -f "$R2_RAW" ]]; then
        echo "   ❌ ERROR: Raw FASTQ files not found for $muestra ($R1_RAW, $R2_RAW)"
        echo "$muestra | MISSING_RAW | SKIPPED | SKIPPED | - | - | -" >> "$RESUMEN"
        fallidas=$((fallidas + 1))
        continue
    fi

    size1=$(stat -c%s "$R1_RAW" 2>/dev/null | numfmt --to=iec)
    size2=$(stat -c%s "$R2_RAW" 2>/dev/null | numfmt --to=iec)
    echo "   📁 Raw Input: R1 ($size1) | R2 ($size2)"

    # ----------------------------------------------------------
    # STEP 2A: PRE-TRIMMING FASTQC
    # ----------------------------------------------------------
    echo "   🔍 [1/3] Running Pre-Trimming FastQC for $muestra..."
    fastqc -t "$THREADS" -o "$FASTQC_RAW_DIR" "$R1_RAW" "$R2_RAW" \
        > "${LOG_DIR}/${muestra}_fastqc_raw.log" 2>&1

    if [ $? -eq 0 ]; then
        pre_qc_status="OK"
        echo "      ✅ Pre-Trimming FastQC completed."
    else
        pre_qc_status="FAIL"
        echo "      ⚠️ Warning: Pre-Trimming FastQC encountered an issue. Check log."
    fi

    # ----------------------------------------------------------
    # STEP 2B: TRIMMOMATIC ADAPTER & QUALITY FILTERING
    # ----------------------------------------------------------
    echo "   ✂️  [2/3] Running Trimmomatic PE for $muestra..."
    R1_paired="${TRIMMED_DIR}/${muestra}_R1_paired.fq.gz"
    R1_unpaired="${TRIMMED_DIR}/${muestra}_R1_unpaired.fq.gz"
    R2_paired="${TRIMMED_DIR}/${muestra}_R2_paired.fq.gz"
    R2_unpaired="${TRIMMED_DIR}/${muestra}_R2_unpaired.fq.gz"

    trimmomatic PE -threads "$THREADS" \
        "$R1_RAW" "$R2_RAW" \
        "$R1_paired" "$R1_unpaired" \
        "$R2_paired" "$R2_unpaired" \
        ILLUMINACLIP:${ADAPTERS}:2:30:7 \
        LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:50 \
        > "${LOG_DIR}/${muestra}_trimmomatic.log" 2>&1

    trim_exit=$?
    size_p1=$(stat -c%s "$R1_paired" 2>/dev/null || echo 0)
    size_p2=$(stat -c%s "$R2_paired" 2>/dev/null || echo 0)

    if [ $trim_exit -eq 0 ] && [ "$size_p1" -gt 1000 ] && [ "$size_p2" -gt 1000 ]; then
        trim_status="OK"
        echo "      ✅ Trimmomatic completed: R1_paired ($(numfmt --to=iec $size_p1)) | R2_paired ($(numfmt --to=iec $size_p2))"
    else
        trim_status="FAIL"
        echo "      ❌ ERROR: Trimmomatic failed or produced empty files for $muestra."
        echo "$muestra | $pre_qc_status | FAIL | SKIPPED | - | - | -" >> "$RESUMEN"
        fallidas=$((fallidas + 1))
        continue
    fi

    # ----------------------------------------------------------
    # STEP 2C: POST-TRIMMING FASTQC (on paired reads)
    # ----------------------------------------------------------
    echo "   🔍 [3/3] Running Post-Trimming FastQC for $muestra..."
    fastqc -t "$THREADS" -o "$FASTQC_TRIM_DIR" "$R1_paired" "$R2_paired" \
        > "${LOG_DIR}/${muestra}_fastqc_trimmed.log" 2>&1

    if [ $? -eq 0 ]; then
        post_qc_status="OK"
        echo "      ✅ Post-Trimming FastQC completed."
        exitosas=$((exitosas + 1))
    else
        post_qc_status="FAIL"
        echo "      ⚠️ Warning: Post-Trimming FastQC encountered an issue. Check log."
        fallidas=$((fallidas + 1))
    fi

    # Calculate sample processing duration
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    duration_min=$((duration / 60))
    duration_sec=$((duration % 60))

    echo "$muestra | $pre_qc_status | $trim_status | $post_qc_status | $(numfmt --to=iec $size_p1) | $(numfmt --to=iec $size_p2) | ${duration_min}m ${duration_sec}s" >> "$RESUMEN"

    echo "   💾 Current RAM usage:"
    free -h | grep "^Mem:" | awk '{print "      " $3 " used / " $2 " total"}'
    sleep 3
done

# ==============================================================================
# 3. GLOBAL AGGREGATION WITH MULTIQC (PRE- AND POST-TRIMMING)
# ==============================================================================
echo ""
echo "============================================================"
echo "📊 GENERATING CONSOLIDATED MULTIQC REPORTS"
echo "============================================================"

echo "   🔄 [1/2] Running MultiQC on Pre-Trimming FastQC outputs..."
multiqc "$FASTQC_RAW_DIR" \
    -o "$MULTIQC_RAW_DIR" \
    -n "multiqc_pre_trimmomatic_report.html" \
    > "${LOG_DIR}/multiqc_raw.log" 2>&1

if [ $? -eq 0 ]; then
    echo "      ✅ Pre-Trimming MultiQC report saved in: ${MULTIQC_RAW_DIR}/multiqc_pre_trimmomatic_report.html"
else
    echo "      ❌ Error running Pre-Trimming MultiQC. Check ${LOG_DIR}/multiqc_raw.log"
fi

echo "   🔄 [2/2] Running MultiQC on Post-Trimming FastQC & Trimmomatic logs..."
multiqc "$FASTQC_TRIM_DIR" "$LOG_DIR" \
    -o "$MULTIQC_TRIM_DIR" \
    -n "multiqc_post_trimmomatic_report.html" \
    > "${LOG_DIR}/multiqc_trimmed.log" 2>&1

if [ $? -eq 0 ]; then
    echo "      ✅ Post-Trimming MultiQC report saved in: ${MULTIQC_TRIM_DIR}/multiqc_post_trimmomatic_report.html"
else
    echo "      ❌ Error running Post-Trimming MultiQC. Check ${LOG_DIR}/multiqc_trimmed.log"
fi

# ==============================================================================
# 4. FINAL PIPELINE SUMMARY
# ==============================================================================
echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║                 PIPELINE EXECUTION COMPLETE                ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo "   Total samples: $total_muestras | ✅ Successful: $exitosas | ❌ Failed: $fallidas"
echo "   📋 Summary table: $RESUMEN"
echo ""
cat "$RESUMEN" | column -t -s '|'
