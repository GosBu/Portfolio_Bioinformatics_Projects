#!/bin/bash

# Aktywacja środowiska Conda 
eval "$(conda shell.bash hook)"
conda activate SRA_tools

# Lista próbek 
SAMPLES=("SRR23609077" "SRR23609078" "SRR23609079" "SRR23609080" "SRR23609081" \
         "SRR23609082" "SRR23609083" "SRR23609084" "SRR23609085" "SRR23609086")

# Pętla do pobierania każdej próbki
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Pobieram próbkę: $SAMPLE"
    fastq-dump "$SAMPLE" --split-3 --skip-technical --gzip
done
