#!/bin/bash

# Ścieżka do genomu referencyjnego
genom_ref="/mnt/c/Users/magor/25_PJATK_Genomika/genome_index/GCF_000005845.2_ASM584v2_genomic.fna"

# Tablica dla prób
samples=('Ecoli_rep2' 'Ecoli_rep3')


# Przetworzenie 2 powtórzenia biologicznego
sample=${samples[0]}

# Szacowanie stopnia pokrycia odczytów 2 powtórzenia biologicznego
bcftools mpileup -O b -o raw_2.bcf -f $genom_ref -q 20 -Q 30 ${sample}_fixmate_sorted_dedup.bam

# Identyfikacja SNV (warianty pojedynczych nt)
bcftools call -m -v --ploidy 1 -o variants_2.vcf raw_2.bcf


# Przetworzenie 3 powtórzenia biologicznego
sample=${samples[1]}

# Szacowanie stopnia pokrycia odczytów 3 powtórzenia biologicznego
bcftools mpileup -O b -o raw_3.bcf -f $genom_ref -q 20 -Q 30 ${sample}_fixmate_sorted_dedup.bam

# Identyfikacja SNV (warianty pojedynczych nt)
bcftools call -m -v --ploidy 1 -o variants_3.vcf raw_3.bcf