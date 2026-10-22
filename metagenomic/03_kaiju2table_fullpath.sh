#!/bin/bash
# ==============================================================================
# Script Name: 01_quality_trimming.sh
# Author: Sofia Fonseca Yepez ann automated writing with AI
# Description: Automated and iterative pipeline for paired-end shotgun metagenomes
#               Sequence alignment against the database
#              

# ==============================================================================

# ------------------------------
# 1. CONFIGURATION & DIRECTORIES
# ------------------------------

TAXDIR="/"           # Directories with nodes.dmp y names.dmp
SAMPLES_DIR="/"        # Directorio donde están los archivos .kaiju.out
OUTPUT_DIR="/"       # Directorio donde se guardarán los resultados

# Crear directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Archivos de taxonomía
NODES="$TAXDIR/nodes.dmp"
NAMES="$TAXDIR/names.dmp"

# Verificar que existan los archivos de taxonomía
if [[ ! -f "$NODES" || ! -f "$NAMES" ]]; then
    echo "Error: No se encuentran nodes.dmp o names.dmp en $TAXDIR"
    exit 1
fi

# Definir lista de muestras
SAMPLES=(
    "SRR11786911_kaiju.out"
    "SRR12425620_kaiju.out"
    "SRR12425871_kaiju.out"
    "SRR12426445_kaiju.out"
    "SRR12426459_kaiju.out"
    "SRR12426461_kaiju.out"
    "SRR12426463_kaiju.out"
    "SRR1509790_kaiju.out"
    "SRR1509794_kaiju.out"
    "SRR1509796_kaiju.out"
    "SRR1509797_kaiju.out"
    "SRR1509799_kaiju.out"
    "SRR1509800_kaiju.out"
)

# Iterar sobre cada muestra
for SAMPLE in "${SAMPLES[@]}"; do
    BASENAME="${SAMPLE%_kaiju.out}"          # Quita la extensión .kaiju.out
    INPUT="$SAMPLES_DIR/$SAMPLE"

    # Verificar que el archivo de entrada existe
    if [[ ! -f "$INPUT" ]]; then
        echo "Advertencia: No se encuentra $INPUT. Se omite."
        continue
    fi

    echo "Procesando $BASENAME ..."


    echo "  Generando archivo con ruta completa (-p)..."
    kaiju2table -t "$NODES" -n "$NAMES" -r species \
        -p -o "$OUTPUT_DIR/${BASENAME}_fullpath.tsv" \
        "$INPUT" 2>> "$OUTPUT_DIR/${BASENAME}_fullpath.log"

    # Verificar que se generó correctamente
    if [[ $? -eq 0 ]]; then
        echo "  ✅ Archivo con ruta completa generado: ${BASENAME}_fullpath.tsv"

      else
        echo "  ❌ Hubo un error al generar: ${BASENAME}_fullpath.tsv (Revisa el log)"
    fi


    # Pequeña pausa para no sobrecargar el sistema
    sleep 1

done

echo "¡Procesamiento de todas las muestras completado!"
echo "Los resultados están en: $OUTPUT_DIR"

# Mostrar resumen de archivos generados
echo ""
echo "Archivos generados:"
ls -lh "$OUTPUT_DIR" | grep -E "\.tsv$"

