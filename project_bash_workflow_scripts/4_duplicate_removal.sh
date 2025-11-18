#!/bin/bash

# Oznaczenie i usunięcie duplikatów PCR z plików .bam (odczyty sparowane - PE)

# Ścieżka do indeksu genomu referencyjnego
ref_genome="/mnt/c/Users/magor/25_PJATK_Genomika/genome_index/coli_index"

# Tablica do próbek 2 i 3 powtórzenia biologicznego
samples=('Ecoli_rep2' 'Ecoli_rep3')

for sample in "${samples[@]}"; do
    echo "Mapowanie próbki: $sample"
    
    # Konwersja SAM -> BAM
    samtools view -bS "${sample}.sam" > "${sample}.bam"

    # Naprawa sparowanych końców odczytów
    samtools fixmate -m "${sample}.bam" "${sample}_fixmate.bam"

    # Sortowanie pliku
    samtools sort "${sample}_fixmate.bam" -o "${sample}_fixmate_sorted.bam"

    # Indeksowanie
    samtools index "${sample}_fixmate_sorted.bam"

    # Oznaczenie i usunięcie duplikatów
    samtools markdup -r "${sample}_fixmate_sorted.bam" "${sample}_fixmate_sorted_dedup.bam"

    echo "Zakończono: $sample"
    echo "--------------------------"
done

echo "Proces zakończony"
