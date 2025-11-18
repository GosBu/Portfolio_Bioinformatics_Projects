#!/bin/bash

# -------------------------------
# Skrypt do przetwarzania danych NGS
# -------------------------------

# Aktywacja środowiska SRA Tools
conda activate SRA_tools

# Pobranie plików fastq
fastq-dump SRR1536581 --split-3 --skip-technical --gzip

# Tworzenie katalogu na raporty FastQC
mkdir -p final_project_quality_reports

# Analiza jakości odczytów FastQC
fastqc -t 4 -o final_project_quality_reports/ *.fastq.gz

# -------------------------------
# Trimmomatic - przycinanie sekwencji
# -------------------------------
source activate trimed

samples=("SRR1536581")
THREADS=4
PARAMS="LEADING:8 TRAILING:8 SLIDINGWINDOW:5:20 MINLEN:40"

for sample in "${samples[@]}"; do
    echo "Przycinam próbkę: $sample"
    trimmomatic SE -phred33 -threads $THREADS \
        ${sample}.fastq.gz \
        ${sample}_trimmed.fastq.gz \
        $PARAMS
    echo "Zakończono: $sample"
done

# FastQC po przycięciu
fastqc -t 4 -o final_project_quality_reports/ *trimmed.fastq.gz

# -------------------------------
# Usuwanie adapterów Cutadapt
# -------------------------------
cutadapt -a AGATCGGAAGAGC -o SRR1536581_cutadapt.fastq.gz SRR1536581_trimmed.fastq.gz
fastqc -t 4 -o final_project_quality_reports/ *cutadapt.fastq.gz

# -------------------------------
# MultiQC raport zbiorczy
# -------------------------------
mkdir -p final_project_quality_reports/multiqc_raports
multiqc -o final_project_quality_reports/multiqc_raports -f -v final_project_quality_reports

# -------------------------------
# Mapowanie BWA MEM
# -------------------------------
conda activate mapping
mkdir -p genome_index

# Jeśli masz skompresowany genom referencyjny
# gunzip hg38.fa.gz
# BWA index, jeśli indeks nie istnieje
# bwa index hg38.fa

bwa mem -t 4 "/ścieżka/do/genome_index/Homo_sapiens_assembly38.fasta" \
    "SRR1536581_cutadapt.fastq.gz" > "Human_rep.sam"

# -------------------------------
# Sortowanie i indeksowanie BAM
# -------------------------------
samtools sort Human_rep.sam -o Human_rep.sorted.bam
samtools index Human_rep.sorted.bam
samtools flagstat Human_rep.sorted.bam > mapping_stats.txt

# -------------------------------
# Usuwanie duplikatów
# -------------------------------
samtools markdup -r Human_rep.sorted.bam Human_rep.sorted_dedup.bam
samtools index Human_rep.sorted_dedup.bam

# -------------------------------
# Indeksowanie genomu do IGV
# -------------------------------
samtools faidx "/ścieżka/do/genome_index/Homo_sapiens_assembly38.fasta"

# -------------------------------
# Identyfikacja wariantów BCFtools
# -------------------------------
conda activate BCF_tools
bcftools mpileup -O b -o raw_1.bcf -f "/ścieżka/do/genome_index/Homo_sapiens_assembly38.fasta" -q 20 -Q 30 Human_rep.sorted_dedup.bam
bcftools call -m -v --ploidy 1 -o final_project_bcf_variants_1.vcf raw_1.bcf

# Statystyki wariantów
bcftools stats final_project_bcf_variants_1.vcf > _final_project_bcf_variant_stats.txt
bcftools view -v snps final_project_bcf_variants_1.vcf | grep -v "^#" | wc -l
bcftools view -v indels final_project_bcf_variants_1.vcf | grep -v "^#" | wc -l

# -------------------------------
# Filtrowanie wariantów
# -------------------------------
bcftools filter -s LOWQUAL -e '%QUAL<20 || DP<10' final_project_bcf_variants_1.vcf -Ov -o filtered_variants.vcf
bcftools stats filtered_variants.vcf > variant_stats.txt
bcftools view -v snps filtered_variants.vcf | grep -v "^#" | wc -l
bcftools view -v indels filtered_variants.vcf | grep -v "^#" | wc -l

echo "Pipeline zakończony pomyślnie."
