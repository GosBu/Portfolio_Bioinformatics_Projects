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
**File:** `project_sql`

---

### 2. **NGS Analysis Pipeline in Snakemake**
`project_snakemake/`  

An automated workflow for processing NGS data (FASTQ → QC → Alignment → Variant
Calling) using Snakemake and Conda environments.

**Technologies:** Snakemake, Conda, BWA, Samtools, FastQC, MultiQC  
**Contents:** `Snakefile`, `config.yaml`, `envs/`, `data/`, `results/`

---

### 3. **Genetic Variant Analysis using Apache Spark (1000 Genomes Project)**
`project_spark_snp_analysis/`  

A project utilising Apache Spark for large-scale genetic data analysis (VCF
files from the 1000 Genomes Project). The goal was to count SNPs, identify
variants in the **CFH** gene, and visualise the variant distribution.

**Technologies:** PySpark, Biopython, Pandas, Matplotlib  
**Environment:** Google Cloud VM (8 vCPU, 62 GB RAM)  
**Results:** SNP distribution histogram, CFH gene analysis

---

### 4. **Basic DNA Sequence Analysis (Python)**
`dna_basic_analysis/`  

A Python script for analysing DNA sequences in FASTA format. Counts nucleotide
frequencies, calculates GC content, and generates bar plots of nucleotide
composition (A, C, G, T, GC).

**Technologies:** Python, Biopython, Pandas, Matplotlib  
**File:** `basic_dna_analysis.py`  
**Output:** `wyniki.csv` + bar charts

---

### 5. **Genomic Workflow (Bash)**
`bash_workflow_scripts/`  

A collection of Bash scripts automating NGS data analysis for *E. coli K-12*,
from read trimming to variant detection.
Designed for reproducibility and modular execution.

**Workflow Steps:**
- Read trimming (`Trimmomatic`)
- Mapping (`BWA`)
- Sorting and indexing (`Samtools`)
- Duplicate removal (`Samtools markdup`)
- Variant calling (`bcftools`)

**Technologies:** Bash, BWA, Samtools, bcftools, Conda  
**Files:** `trymowanie.sh`, `mapowanie.sh`, `sortowanie_indeksowanie.sh`, `usuwanie_duplikatow.sh`, `analiza_wariantow.sh`

---

## Core Technical Competencies

**Programming:**  
![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-121011?logo=gnu-bash&logoColor=white)

**Bioinformatics:**  
![Snakemake](https://img.shields.io/badge/Snakemake-3C8DBC?logo=snakemake&logoColor=white)
![Biopython](https://img.shields.io/badge/Biopython-3776AB?logo=python&logoColor=white)
![Samtools](https://img.shields.io/badge/Samtools-009688?logo=data:image/svg+xml;base64,PHN2ZyBmaWxsPSIjZmZmIiB2aWV3Qm94PSIwIDAgMTIgMTIiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+PHJlY3QgeD0iMSIgeT0iMSIgd2lkdGg9IjEwIiBoZWlnaHQ9IjEwIiByeD0iMSIvPjwvc3ZnPg==)

**Data Analysis & Visualization:**  
![Pandas](https://img.shields.io/badge/Pandas-150458?logo=pandas&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557C?logo=plotly&logoColor=white)
![Seaborn](https://img.shields.io/badge/Seaborn-3776AB?logo=python&logoColor=white)

**Cloud / Big Data:**  
![Apache Spark](https://img.shields.io/badge/Apache%20Spark-E25A1C?logo=apachespark&logoColor=white)
![Google Cloud](https://img.shields.io/badge/Google%20Cloud-4285F4?logo=googlecloud&logoColor=white)
![AWS](https://img.shields.io/badge/AWS-FF9900?logo=amazonaws&logoColor=white)


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

To portfolio prezentuje praktyczne projekty z zakresu bioinformatyki, analizy danych genetycznych i uczenia maszynowego.  
Każdy projekt opiera się na rzeczywistych danych biologicznych (DNA, SNP, geny, mutacje) i pokazuje umiejętność:

- pracy z dużymi zbiorami danych (Spark, Cloud),
- analizy danych biologicznych (Pandas, Biopython, scikit-learn),
- interpretacji wyników w kontekście biologicznym.

---

## Projekty

### 1. **Analiza danych genetycznych w SQL**  
`project_sql/`  
Analiza danych pacjentów i wyników testów genetycznych przy użyciu zapytań SQL.  
Skupiono się na identyfikacji testów z wieloma pacjentami i analizie częstości wariantów w różnych typach testów (SNP Array, NGS, WES, WGS).  

**Technologie:** SQL, relacyjne bazy danych  
**Plik:** `pd4657-DATAB-projekt.sql`

---

### 2. **Pipeline analizy NGS w Snakemake**  
`project_snakemake/`  
Zautomatyzowany workflow przetwarzania danych NGS (FASTQ → QC → Alignment → Variant Calling) z wykorzystaniem Snakemake i środowisk Conda.  

**Technologie:** Snakemake, Conda, BWA, Samtools, FastQC, MultiQC  
**Zawartość:** `Snakefile`, `config.yaml`, `envs/`, `data/`, `results/`

---

### 3. **Analiza wariantów genetycznych w Spark (1000 Genomes Project)**  
`project_spark_snp_analysis/`  
Projekt wykorzystujący Apache Spark do analizy dużych danych genetycznych (plik VCF z projektu 1000 Genomes).  
Celem było policzenie SNP, identyfikacja wariantów w genie CFH i wizualizacja rozkładu wariantów.  

**Technologie:** PySpark, Biopython, Pandas, Matplotlib  
**Środowisko:** Google Cloud VM (8 vCPU, 62 GB RAM)  
**Wynik:** histogram rozkładu SNP, analiza genu CFH

---

### 4. **Podstawowa analiza sekwencji DNA (Python)**  
`dna_basic_analysis/`  
Skrypt analizujący sekwencje DNA w formacie FASTA.  
Zlicza częstość nukleotydów, oblicza zawartość GC, generuje wykresy procentowej zawartości A, C, G, T oraz GC.  

**Technologie:** Python, Biopython, Pandas, Matplotlib  
**Plik:** `basic_dna_analysis.py`  
**Wynik:** `wyniki.csv` + wykresy słupkowe  

---

### 5. **Genomic Workflow (Bash)**  
`bash_workflow_scripts/`  
Zestaw skryptów Bash automatyzujących analizę danych NGS dla *E. coli K-12* od przycinania FASTQ po wykrywanie wariantów.  

**Etapy:**
- Trymowanie odczytów (`Trimmomatic`)
- Mapowanie (`BWA`)
- Sortowanie i indeksowanie (`Samtools`)
- Usuwanie duplikatów (`Samtools markdup`)
- Analiza wariantów (`bcftools`)

**Technologie:** Bash, BWA, Samtools, bcftools, Conda  
**Pliki:** `trymowanie.sh`, `mapowanie.sh`, `sortowanie_indeksowanie.sh`, `usuwanie_duplikatow.sh`, `analiza_wariantow.sh`

---

## Kluczowe kompetencje techniczne

| Obszar | Technologie / Narzędzia |
|--------|--------------------------|
| **Programowanie** | Python (Pandas, Matplotlib, Biopython), Bash |
| **Big Data / Chmura** | Apache Spark, Google Cloud, AWS |
| **Analiza biologiczna** | Snakemake, FastQC, BWA, Samtools, bcftools, Conda |
| **Bazy danych** | SQL |
| **Wizualizacja danych** | Matplotlib, Seaborn |

---

## O autorce

**Gosia** - z wykształcenia mgr inż. Biotechnologii (2016), ukończyłam Studia Podyplomowe z Bioinformatyki (PJATK, 2025), aby zdobyć kluczowe umiejętności, które pomogą mi w przebranżowieniu się.  
Uwielbiam łączyć kod, dane i pytanie *„dlaczego?”* w jedną spójną historię.

📍 Obecnie: Europa, Polska  
📫 Kontakt: **[magorzata.bujak@gmail.com](mailto:magorzata.bujak@gmail.com)**

---

![Last Update](https://img.shields.io/badge/Updated-November%202025-blue)
