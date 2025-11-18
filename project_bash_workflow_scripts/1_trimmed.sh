#!/bin/bash

source activate SRA_tools

# Lista próbek
samples=("SRR25629154" "SRR25629153")


# Liczba rdzeni
THREADS=4

# Parametry Trimmomatic
PARAMS="LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50"

# Petla po wszystkich próbkach
for sample in "${samples[@]}"
do
echo "Przycinam próbkę: $sample"

trimmomatic PE -phred33 -threads $THREADS \
${sample}_1.fastq.gz ${sample}_2.fastq.gz \
${sample}_1_trimmed_paired.fastq.gz /dev/null \
${sample}_2_trimmed_paired.fastq.gz /dev/null \
$PARAMS
echo "Zakończono: $sample"
done
