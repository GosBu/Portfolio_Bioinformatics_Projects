# Pipeline for SARS-CoV-2 NGS Data Analysis

Projekt Snakemake do analizy danych sekwencjonowania (NGS) wirusa SARS-CoV-2. Pipeline obejmuje jakość danych, przycięcie, mapowanie, pokrycie, warianty oraz raport MultiQC.

## Struktura katalogów

```

PRODATA\_projekt/
├── Snakefile
├── config.yaml
├── envs/
│   ├── fastqc-env.yaml
│   ├── trimmomatic-env.yaml
│   ├── bwa-env.yaml
│   ├── samtools-env.yaml
│   ├── bcftools-env.yaml
│   └── multiqc-env.yaml
├── data/
│   └── raw/
│       ├── sequence_NC_045512.2.fasta
│       └── <sample>\_<1|2>.fastq.gz
├── results/
│   ├── quality\_reports/
│   ├── index
│       └── genome\_index/
│   └── quality\_reports\_multiqc/
└── dag.pdf

````

## Wymagania

- Python 3.10
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Snakemake (instalacja poniżej)

## Instalacja środowiska

```bash
conda create -n snakemake python=3.10
conda activate snakemake
conda install -c conda-forge -c bioconda snakemake
````

## Konfiguracja

Zdefiniuj próbki i inne parametry w `config.yaml`:

```yaml
samples:
  - SRR23609077
  - SRR23609078
  - SRR23609079
  - SRR23609080
  - SRR23609081
  - SRR23609082
  - SRR23609083
  - SRR23609084
  - SRR23609085
  - SRR23609086
  
reads: ["1", "2"]
```

## Uruchomienie workflowu

### 1. Standardowo:

```bash
snakemake --use-conda --cores 4
```

### 2. Raport HTML z całego workflowu:

```bash
snakemake --report raport.html
```

### 3. Generacja grafu DAG:

```bash
snakemake --dag --forceall | dot -Tpdf > dag.pdf
```

## Etapy analizy

1. `fastqc` – analiza jakości sekwencji
2. `multiqc` – raport zbiorczy z jakości danych 
3. `trimmomatic` – przycięcie słabej jakości odczytów
4. `bwa_index_genome_ref` – indeksacja genomu referencyjnego
5. `bwa_mem_do_bam` – mapowanie FASTQ na genom (BAM)
6. `samtools_sort_index` – sortowanie i indeksacja BAM
7. `genome_coverage` – pokrycie genomu
8. `variant_calling` – wykrywanie wariantów (bcftools)
