# Analiza wariantów genetycznych

## 1. Opis danych i źródeł

**Dane:**
Projekt 1000 Genomes, plik w formacie VCF z wariantami genetycznymi dla chromosomu 1.

**Źródło:**
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

**Format pliku:**
VCF (Variant Call Format), kolumny: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.

**Cel analizy:**
- Policzyć liczbę wariantów SNP na chromosomie 1
- Zidentyfikować SNP w genie CFH (pozycje 196,621,007 – 196,716,634)
- Zwizualizować rozkład wariantów na chromosomie



## 2. Metody analizy i użyte narzędzia

**Chmura:**
Google Cloud Compute Engine

**Maszyna:**
Instancja Virtual Machine: C4, 8 vCPU, 62 GB pamięci; Linux 12

**Środowisko:**
- Apache Spark 4.0.0
- Python 3.11.2
- JAVA openjdk 17.0.16
- wirtualne środowisko o nazwie własnej: *venv_bio*

**Biblioteki:**
- PySpark - do efektywnego przetwarzania dużych plików (np. VCF); działa na RDD, które pozwalają analizować bardzo duże pliki w rozproszony sposób
- Pandas - do pracy z danymi tabelarycznymi (np. filtrowanie)
- Biopython - ułatwia pracę z plikami VCF oraz analizę sekwencji np.wariantów genetycznych (np. SNP)
- Matplotlib - do wizualizacji wyników


**Metody analizy:**

- Zainstalowanie potrzebnych bibliotek (`sudo apt-get`) i utworzenie oraz aktywacja witrualnego środowiska o nazwie *venv_bio*, aby zainstalować pozostałe biblioteki (wszystkie poza java) (`pip install`)
- Pobranie danych (`wget`)
- Rozpakowanie pobranego archiwum (`unzip`)
- Aktywacja Spark (SparkContext)
- Wczytaj pliku VCF dla chromosomu 1 do RDD
- Usunięcie nagłówków (`.filter()` i funkcja `lambda`)
- Policzenie wszystkich wariantów w pliku (`.count()`)
- Podziel linie na kolumny (`.map()` i `lambda`)
- Wyciągnięcie interesujących kolumn i zapisanie każdej do krotki (`parsed.map()` i `lambda`)
- Filtracja wariantów SNP (REF i ALT długości 1) (`.filter`)
- Policzenie wszystkich SNP na chromosomie 1 (`.count()`)
- Wyznaczenie SNP w genie CFH na chromosomie 1
- Wizualizacja rozkładu SNP (histogram) (biblioteka matplot)



## 3. Wyniki i ich interpretacja

- (Zał. 1) Całkowita liczba wariantów SNP: **6 468 094**
- (Zał. 2) Całkowita liczba SNP na chromosomie 1: **6 196 151**
- (Zał. 3) Liczba SNP w genie CFH: **2 777** - to normalna liczba wariantów dla tego genu w populacji

- (Zał. 4) Histogram pokazuje rozmieszczenie SNP na chromosome 1 - s rozmieszczenie w histogramie może wskazywać na obszary bardziej podatne na mutacje.
- 


## 4. Wnioski

Gen CFH jest jednym z genów związanych z chorobami oczu (np. AMD). Analiza SNP pozwala ocenić częstość wariantów i ich potencjalne znaczenie.


## Wykresy i tabele

Załączniki:
- Zał. 1: Zrzuty ekranu: *Zał.1_Klastry_Projekt_Liczba_wariantów_6468094*
- Zał. 2: Zrzuty ekranu: *Zał.2_Klastry_Projekt_Całkowita_liczba_SNP_chrom_1*
- Zał. 3: Zrzuty ekranu: *Zał.3_Klastry_Projekt_Liczba_SNP_w_genie_CFH*
- Zał. 4: Histogram: *snp_distribution_chr1_sample.png*
- Zał. 5: Zrzuty ekranu: *Zał.5_Klastry_Projekt_ściągnięcie_archiwum*
- Zał. 6: Zrzuty ekranu: *Zał.6_Klastry_Projekt_Zrzut_ekranu_Przepustowość*
- Zał. 7: Zrzut ekranu: *Zał.7_Klastry_Projekt_Zrzut_ekranu_Wykorzystanie_miejsca_na_dysku*
- Zał. 8: Zrzut ekranu: *Zał.8_Klastry_Projekt_Zrzut_ekranu_Wykorzystanie_procesora*
