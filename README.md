# Bioinformatics Portfolio / Portfolio Bioinformatyczne

---

[Polish version below]

---

## **Repository Purpose**

This portfolio showcases practical projects in bioinformatics, genetic data
analysis, and machine learning. Each project is based on real biological data
(DNA sequences, SNPs, genes, and mutations) and demonstrates skills in:

- Working with large-scale biological datasets (Spark, Cloud)
- Performing bioinformatics analyses (Pandas, Biopython, scikit-learn)
- Interpreting biological results through data-driven insights

---

## Projects

### 1. **Genetic Data Analysis in SQL**
`project_sql/`

Analysis of patient and genetic test data using SQL queries. Focused on
identifying tests with multiple patients and analysing the frequency of variants
across different test types (SNP Array, NGS, WES, WGS).

**Technologies:** SQL, Relational Databases  
**Contents:** `project_sql.md`

---

### 2. **NGS Analysis Pipeline in Snakemake**
`project_snakemake/`  

An automated workflow for processing NGS data (FASTQ → QC → Alignment → Variant
Calling) using Snakemake and Conda environments.

**Technologies:** Snakemake, Conda, BWA, Samtools, FastQC, MultiQC  
**Contents:** `Snakefile`, `config.yaml`, `envs/`, `data/`, `results/`, `README.md`,
`dag_full.png`, `report.html`

---

### 3. **NGS QC & Trimming Pipeline in Nextflow**
`project_nextflow/`

An automated workflow for quality control and trimming of NGS data
(FASTQ → QC → MultiQC → Trimming) using Nextflow and configurable parameters.
The workflow uses channels, params objects, and log.info for process tracking.

**Technologies:** Nextflow, FASTQC, MultiQC, Trimmomatic
**Contents:** `PRODATA_7.nf`, `timeline.html`, `report.html`, `trace.txt`, `results/`

---

### 4. **Genetic Variant Analysis using Apache Spark (1000 Genomes Project)**
`project_spark_snp_analysis/`  

A project utilising Apache Spark for large-scale genetic data analysis (VCF
files from the 1000 Genomes Project). The goal was to count SNPs, identify
variants in the **CFH** gene, and visualise the variant distribution.

**Technologies:** PySpark, Biopython, Pandas, Matplotlib  
**Environment:** Google Cloud VM (8 vCPU, 62 GB RAM)  
**Results:** SNP distribution histogram, CFH gene analysis
**Contents:** `SNP_chr1_CFH_Spark_GoogleCloud`, `SNP_chr1_CFH_Spark_GoogleCloud`,
`Zał.1_Klastry_Projekt_Liczba_wariantów_6468094`.
`Zał.2_Klastry_Projekt_Całkowita_liczba_SNP_chrom_1`,
`Zał.3_Klastry_Projekt_Liczba_SNP_w_genie_CFH`,
`Zał.5_Klastry_Projekt_ściągnięcie_archiwum`,
`Zał.6_Klastry_Projekt_Zrzut_ekranu_Przepustowość`,
`Zał.7_Klastry_Projekt_Zrzut_ekranu_Wykorzystanie_miejsca_na_dysku`,
`Zał.8_Klastry_Projekt_Zrzut_ekranu_Wykorzystanie_procesora`

---

### 5. **Basic DNA Sequence Analysis (Python)**
`project_python_dna_basic_analysis/`  

A Python script for analysing DNA sequences in FASTA format. Counts nucleotide
frequencies, calculates GC content, and generates bar plots of nucleotide
composition (A, C, G, T, GC).

**Technologies:** Python, Biopython, Pandas, Matplotlib
**Contents:** `basic_dna_analysis.py`, `basic_dna_analysis-code_description.txt`,
`clean_sequences.fasta`, `ncbi_reference_sequences.txt`, `sequences.txt`

---

### 6. **Genomic Workflow (Bash)**
`project_bash_workflow_scripts/`  

A collection of Bash scripts automating NGS data analysis for *E. coli K-12*,
from read trimming to variant detection.
Designed for reproducibility and modular execution.

**Workflow Steps:**
- Download reads (`fastq-dump`)
- Read trimming (`Trimmomatic`)
- Mapping (`BWA`)
- Sorting and indexing (`Samtools`)
- Duplicate removal (`Samtools markdup`)
- Variant calling (`bcftools`)

**Technologies:** Bash, BWA, Samtools, bcftools, Conda  
**Contents:** `NGS_reads_download.sh`, `1_trimmed.sh`, `2_mapping.sh`,
`3_sorting_&iand_indexing.sh`, `4_duplicate_removal.sh`, `5_variant_calling.sh`

---

### 7. **ML_mutation_classifier**
`project_ML_mutation_classifier/`

A project aimed at creating a classification model to predict whether a given
mutation in a DNA sequence is pathogenic. Focused on genetic sequences
associated with hereditary diseases.

**Goal:** Build a model to classify DNA mutations as pathogenic or benign.
**Technologies:** Python, scikit-learn, Pandas, NumPy, Matplotlib
**Results:** Model performance metrics (accuracy, precision, recall), visualisation
of mutation effects
**Contents:** `mutation_classifier_cross_validation.py`, `mutation_classifier.py`,
`mutation_classifier_interpretacja_wynikow.md`

---

### 8. **NGS Data Analysis**
`project_NGS_Human_Variants/`

A project demonstrating skills acquired in the course Bioinformatics in
Genomics and Transcriptomics. The goal is to perform a complete NGS data
analysis workflow, from raw read quality control to downstream analysis.

**Goal:** Showcase a full NGS workflow, including QC, processing, mapping, and
specific genomic or transcriptomic analyses.
**Technologies:** Bash, Python, R, FastQC, MultiQC, Trimmomatic, Cutadapt,
BWA/STAR, Samtools, IGV, Matplotlib, ggplot2
**Contents:** `NGS_final_project_dokumentacja.pdf`, `NGS_final_project_dokumentacja.Rmd`,
`NGS_final_project_full_code.txt`, `NGS_pipeline.sh`

---

### 9. **Bioinformatic FASTA Pipeline (Bash)**
`project_linux/`

A project demonstrating practical proficiency in Linux CLI and bioinformatics
tools. Confirms ability to automate tasks using Bash scripts and build
repeatable, secure analytical workflows for genomic data—essential for
IT/Data Science roles.

**Goal:** Set up working environment, manage genomic files (Insulin gene orthologs),
process data, format FASTA sequences, and run automated analysis via CLI.
**Technologies:** Bash, Linux CLI, AWK, sed, grep
**Contents:** `README.md`, `linux_download_repo.sh`,
`bioinformatyczny_pipeline_bash_obróbka_sekwencji_FASTA.txt`

---

### 10. **Data Analysis in R**
`project_R`

A project demonstrating practical proficiency in R programming and data
analysis. Confirms ability to perform a complete data workflow, from
import and cleaning to exploration, statistical analysis, visualisation, and
reporting—showcasing both technical accuracy and analytical insight.

**Technologies:** R, tidyverse, ggplot2, RMarkdown
**Contents:** `Analiza_Danych_Genomowych_R.pdf`, `Analiza_Danych_Genomowych_R.Rmd`

---

## Core Technical Competencies

**Programming & Scripting**  
![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-121011?logo=gnu-bash&logoColor=white)
![SQL](https://img.shields.io/badge/SQL-4479A1?logo=postgresql&logoColor=white)

**Bioinformatics**  
![Snakemake](https://img.shields.io/badge/Snakemake-3C8DBC?logo=snakemake&logoColor=white)
![BWA](https://img.shields.io/badge/BWA-0288D1)
![Samtools](https://img.shields.io/badge/Samtools-009688)
![bcftools](https://img.shields.io/badge/bcftools-4CAF50)
![FastQC](https://img.shields.io/badge/FastQC-1976D2)
![MultiQC](https://img.shields.io/badge/MultiQC-6A1B9A)

**Data Analysis & Visualization**  
![Pandas](https://img.shields.io/badge/Pandas-150458?logo=pandas&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557C?logo=plotly&logoColor=white)
![Seaborn](https://img.shields.io/badge/Seaborn-3776AB?logo=python&logoColor=white)
![scikit-learn](https://img.shields.io/badge/scikit--learn-F7931E?logo=scikit-learn&logoColor=white)

**Cloud & Big Data**  
![Apache Spark](https://img.shields.io/badge/Apache%20Spark-E25A1C?logo=apachespark&logoColor=white)
![Google Cloud](https://img.shields.io/badge/Google%20Cloud-4285F4?logo=googlecloud&logoColor=white)
![AWS](https://img.shields.io/badge/AWS-FF9900?logo=amazonaws&logoColor=white)

**Reproducibility & Environments**  
![Conda](https://img.shields.io/badge/Conda-44A833?logo=anaconda&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-2496ED?logo=docker&logoColor=white)
![Git](https://img.shields.io/badge/Git-F05032?logo=git&logoColor=white)

---

## About the Author

**Gosia** - MSc in Biotechnology (2016) and Postgraduate Studies in
Bioinformatics (PJATK, 2025). Passionate about connecting code, data, and the
question *“why?”* into one coherent story.

📍 Currently based in Europe, Poland  
📫 Contact: **[magorzata.bujak@gmail.com](mailto:magorzata.bujak@gmail.com)**

---

![Last Update](https://img.shields.io/badge/Updated-November%202025-blue)

---

# Portfolio Bioinformatyczne

---

## **Cel repozytorium**

To portfolio prezentuje praktyczne projekty z zakresu bioinformatyki, analizy
danych genetycznych i uczenia maszynowego. Każdy projekt opiera się na rzeczywistych danych biologicznych (DNA, SNP, geny, mutacje) i pokazuje
umiejętność:

- pracy z dużymi zbiorami danych (Spark, Cloud),
- analizy danych biologicznych (Pandas, Biopython, scikit-learn),
- interpretacji wyników w kontekście biologicznym.

---

## Projekty

### 1. **Analiza danych genetycznych w SQL**  
`project_sql/`  

Analiza danych pacjentów i wyników testów genetycznych przy użyciu zapytań SQL.  
Skupiono się na identyfikacji testów z wieloma pacjentami i analizie częstości
wariantów w różnych typach testów (SNP Array, NGS, WES, WGS).  

**Technologie:** SQL, relacyjne bazy danych  
**Plik:** `project_sql.md`

---

### 2. **Pipeline analizy NGS w Snakemake**  
`project_snakemake/`  

Zautomatyzowany workflow przetwarzania danych NGS (FASTQ → QC → Alignment →
Variant Calling) z wykorzystaniem Snakemake i środowisk Conda.  

**Technologie:** Snakemake, Conda, BWA, Samtools, FastQC, MultiQC  
**Zawartość:** `Snakefile`, `config.yaml`, `envs/`, `data/`, `results/`, `README.md`, `dag_full.png`, `report.html`

---

### 3. **Pipeline QC i trymowania NGS w Nextflow**
`project_nextflow/`

Zautomatyzowany workflow do kontroli jakości i trymowania danych NGS
(FASTQ → QC → MultiQC → Trimming) z wykorzystaniem Nextflow i parametrów konfiguracyjnych.
Workflow używa kanałów, obiektów params oraz log.info do monitorowania przebiegu procesów.

**Technologie:** Nextflow, FASTQC, MultiQC, Trimmomatic
**Zawartość:** `PRODATA_7.nf`, `config.yaml`, `timeline.html`, `report.html`, `trace.txt`

---

### 4. **Analiza wariantów genetycznych w Apache Spark (1000 Genomes Project)**  
`project_spark_snp_analysis/`  

Projekt wykorzystujący Apache Spark do analizy dużych danych genetycznych (plik
VCF z projektu 1000 Genomes). Celem było policzenie SNP, identyfikacja
wariantów w genie CFH i wizualizacja rozkładu wariantów.  

**Technologie:** PySpark, Biopython, Pandas, Matplotlib  
**Środowisko:** Google Cloud VM (8 vCPU, 62 GB RAM)  
**Wynik:** histogram rozkładu SNP, analiza genu CFH

---

### 5. **Podstawowa analiza sekwencji DNA (Python)**  
`dna_basic_analysis/`

Skrypt analizujący sekwencje DNA w formacie FASTA. Zlicza częstość nukleotydów,
oblicza zawartość GC, generuje wykresy procentowej zawartości A, C, G, T oraz
GC.  

**Technologie:** Python, Biopython, Pandas, Matplotlib  
**Plik:** `basic_dna_analysis.py`  
**Wynik:** `basic_dna_analysis-code_description.txt`, `clean_sequences.fasta`, `ncbi_reference_sequences.txt`, `sequences.txt`

---

### 6. **Genomic Workflow (Bash)**  
`bash_workflow_scripts/`  

Zestaw skryptów Bash automatyzujących analizę danych NGS dla *E. coli K-12* od
przycinania FASTQ po wykrywanie wariantów.  

**Etapy:**
- Pobieranie odczytów (`fastq-dump`)
- Trymowanie odczytów (`Trimmomatic`)
- Mapowanie (`BWA`)
- Sortowanie i indeksowanie (`Samtools`)
- Usuwanie duplikatów (`Samtools markdup`)
- Analiza wariantów (`bcftools`)

**Technologie:** Bash, BWA, Samtools, bcftools, Conda  
**Pliki:** `NGS_reads_download.sh`, `1_trimmed.sh`, `2_mapping.sh`, `3_sorting_&iand_indexing.sh`, `4_duplicate_removal.sh`, `5_variant_calling.sh`

---

## Kluczowe kompetencje techniczne

**Programowanie i Skrypty**  
![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-121011?logo=gnu-bash&logoColor=white)
![SQL](https://img.shields.io/badge/SQL-4479A1?logo=postgresql&logoColor=white)

**Bioinformatyka**  
![Snakemake](https://img.shields.io/badge/Snakemake-3C8DBC?logo=snakemake&logoColor=white)
![BWA](https://img.shields.io/badge/BWA-0288D1)
![Samtools](https://img.shields.io/badge/Samtools-009688)
![bcftools](https://img.shields.io/badge/bcftools-4CAF50)
![FastQC](https://img.shields.io/badge/FastQC-1976D2)
![MultiQC](https://img.shields.io/badge/MultiQC-6A1B9A)

**Analiza i Wizualizacja Danych**  
![Pandas](https://img.shields.io/badge/Pandas-150458?logo=pandas&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557C?logo=plotly&logoColor=white)
![Seaborn](https://img.shields.io/badge/Seaborn-3776AB?logo=python&logoColor=white)
![scikit-learn](https://img.shields.io/badge/scikit--learn-F7931E?logo=scikit-learn&logoColor=white)

**Chmury i Big Data**  
![Apache Spark](https://img.shields.io/badge/Apache%20Spark-E25A1C?logo=apachespark&logoColor=white)
![Google Cloud](https://img.shields.io/badge/Google%20Cloud-4285F4?logo=googlecloud&logoColor=white)
![AWS](https://img.shields.io/badge/AWS-FF9900?logo=amazonaws&logoColor=white)

**RPotwarzalność i Środowiska**  
![Conda](https://img.shields.io/badge/Conda-44A833?logo=anaconda&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-2496ED?logo=docker&logoColor=white)
![Git](https://img.shields.io/badge/Git-F05032?logo=git&logoColor=white)

---

## O autorce

**Gosia** - mgr inż. Biotechnologii (2016), absolwentka Studiów Podyplomowych
z Bioinformatyki (PJATK, 2025). Łącze dane, kod i pytanie *„dlaczego?”* w jedną
spójną historię.

📍 Obecnie: Europa, Polska  
📫 Kontakt: **[magorzata.bujak@gmail.com](mailto:magorzata.bujak@gmail.com)**

---

![Last Update](https://img.shields.io/badge/Updated-November%202025-blue)
