#!/bin/bash

# Dalsze kroki po mapowaniu - sortowanie, nadanie indeksu

# Ścieżka do indeksu genomu referencyjnego
ref_genome="/mnt/c/Users/magor/25_PJATK_Genomika/genome_index/coli_index"

# Tablica do próbek 2 i 3 powtórzenia biologicznego
samples=('SRR25629154' 'SRR25629153')

for sample in "${samples[@]}"; do
    echo "Mapowanie próbki: $sample"

    # Mapowanie i konwersja do BAM
    bwa mem -t 4 "$ref_genome" "${sample}_1_trimmed_paired.fastq.gz" "${sample}_2_trimmed_paired.fastq.gz" |
        samtools view -bS > "${sample}.bam"

    # Sortowanie pliku BAM
    samtools sort "${sample}.bam" -o "${sample}.sorted.bam"

    # Indeksowanie pliku BAM
    samtools index "${sample}.sorted.bam"

    echo "Zakończono: $sample"
    echo "--------------------------"
done
