#!/bin/bash

# Mapowanie przyciętych odczytów dla E. coli K-12 przy użyciu BWA MEM

# Ścieżka do indeksu genomu referencyjnego
ref_genome="/mnt/c/Users/magor/25_PJATK_Genomika/genome_index/coli_index"


# Powtórzenie 2: SRR25629154
bwa mem -t 4 "$ref_genome" "SRR25629154_1_trimmed_paired.fastq.gz" "SRR25629154_2_trimmed_paired.fastq.gz" > "Ecoli_rep2.sam"

# Powtórzenie 3: SRR25629153
bwa mem -t 4 "$ref_genome" "SRR25629153_1_trimmed_paired.fastq.gz" "SRR25629153_2_trimmed_paired.fastq.gz" > "Ecoli_rep3.sam"
