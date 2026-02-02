# Indian Population Structural Variant Pipeline

A modular, reproducible Python pipeline that analyzes **structural variants (SVs)** from the **1000 Genomes Project Phase 3** (GRCh38 assembly) and compares their allele frequencies in **South Asian (Indian)** populations vs global ones, focusing on genes involved in nutrition and metabolism.

## Why This Project?

I am passionate about **nutrigenomics** — how genetic variation shapes individual responses to diet, particularly in metabolic diseases like type 2 diabetes (which disproportionately affects Indian populations).  

Structural variants (deletions, duplications, inversions) often have larger functional impacts than SNPs, yet they are understudied in South Asian cohorts. This pipeline explores whether SVs in key metabolic genes show population-specific patterns that could inform **personalized nutrition** strategies (e.g., low-carb recommendations for individuals with certain variants).

This work aligns with research at labs like **CCMB Hyderabad** (genetic variation in Indian populations, health disparities) and builds on my earlier projects on TCF7L2 variants and ketogenic diet gene expression.

## Biological Question

**Do structural variants in nutrition- and metabolism-related genes show different allele frequencies in South Asian (Indian) populations compared to other global groups, and what might this imply for diet-related health outcomes (e.g., diabetes risk or response to low-carb diets)?**

## Dataset

- **Source**: 1000 Genomes Project Phase 3 integrated structural variant callset (lifted to GRCh38/hg38)
- **File**: `ALL.wgs.integrated_sv_map_v2_GRCh38.20130502.svs.genotypes.vcf.gz` (~18 MB compressed)
- **URL**: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/
- **Samples**: 2,504 individuals (including South Asian super-population = SAS)
- **Content**: ~68,700 structural variants with genotypes and population-specific allele frequencies (SAS_AF, EUR_AF, AFR_AF, etc.)

## Approach & Methods

1. **Download & Setup**  
   - Automated HTTP download + gzip decompression to `output/`

2. **Parsing**  
   - Line-by-line streaming read (memory-efficient for large VCF)
   - Extract key INFO tags: SVTYPE, SVLEN, END, SAS_AF, EUR_AF, AFR_AF
   - Overlap check: variant position (POS–END) overlaps predefined gene regions (±100 kb)

3. **Filtering**  
   - Keep only DEL and DUP types
   - Match against 9 nutrition/metabolism genes (TCF7L2, GCK, PCK1, G6PC, PDK4, FASN, PPARGC1A, HMGCS2, LCT)

4. **Output & Visualization**  
   - CSV table of matching variants
   - Text summary with overview
   - Two plots:
     - Population frequency bar plot (SAS vs EUR vs AFR per gene)
     - SV size distribution histogram (absolute length per gene)

## Results

- **Variants found**: 35 structural variants (mostly deletions) overlapping gene regions after widening search windows to ±100 kb.
- **Key genes with hits**:
  - LCT (lactose tolerance): highest number and frequencies
  - TCF7L2 (diabetes risk): several deletions
  - PDK4 (ketogenic adaptation): deletions
  - FASN, HMGCS2, others: fewer but present
- **Population patterns**: Some SVs show higher allele frequency in SAS (South Asian/Indian) compared to EUR/AFR — potential relevance for diet response or metabolic health disparities in India.
- **Plots**:
  - `sv_population_frequencies.png`: Visual comparison of allele frequencies across populations per gene.
  - `sv_size_distribution.png`: Distribution of SV lengths (skipped in some runs due to too few unique values).

## Errors Encountered & How They Were Solved

- **cyvcf2 installation failed on Windows**  
  → Switched to pure Python line-by-line parsing (no C dependencies, more transparent).

- **Gene names not present in INFO tags**  
  → Added coordinate-based overlap check using hg38 gene regions (from Ensembl/NCBI).

- **KDE crash in seaborn histplot ("dataset input should have multiple elements")**  
  → Disabled KDE (`kde=False`) when few unique SVLEN values; added check for minimum data points.

- **KeyError: 'SAS_AF_num' during sorting**  
  → Removed unused temp column; now convert `SAS_AF` to numeric on-the-fly for sorting.

- **Indentation & import errors during development**  
  → Fixed with consistent 4-space indentation, proper imports, and PyCharm reformat (Ctrl+Alt+L).

These errors taught me real-world debugging, resilience, and the importance of modular, testable code.

## Skills Demonstrated

- Modular Python architecture (config, utils, parser, visualization separation)
- Efficient streaming file parsing of large genomic VCFs
- Genomic interval overlap logic (coordinate-based filtering)
- Population genetics data handling (allele frequencies, super-populations)
- Data visualization (matplotlib/seaborn static plots)
- Error handling, logging, reproducible workflows (requirements.txt)
- Public data integration (FTP download from EBI/1000 Genomes)
- Debugging under deadline pressure

## How to Run

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   
2. Run the pipeline:
   ```bash
      python main.py

## Outputs appear in output/:
- sv_results.csv: Table of matching SVs (sorted by SAS_AF)
- sv_summary.txt: Human-readable overview 
- sv_population_frequencies.png: Bar plot of frequencies per gene/population 
- sv_size_distribution.png: Histogram of SV lengths per gene

## Basic workflow:
Download VCF → Unzip → Parse line-by-line → Filter by gene regions & SV type → Extract population frequencies → Save table + plots.

## Future Prospects
- Add exact gene annotation (Ensembl VEP API)
- Integrate GWAS Catalog or UK Biobank for functional impact 
- Build interactive dashboard (Streamlit/Plotly Dash)
- Extend to tandem repeats or epigenomic data in Indian cohorts 
- Compare with clinical datasets (e.g., diabetes cohorts)

***This project sparked my interest in population-specific genomics and its link to nutrition — I hope it contributes to better understanding of health disparities in India.***

Feedback welcome!
Gayathri Jonnalagadda — February 2026
