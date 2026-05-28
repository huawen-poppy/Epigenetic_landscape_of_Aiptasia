# Aiptasia Nanopore DRS Analysis

Analysis code for Oxford Nanopore direct RNA sequencing (DRS) data from *Aiptasia diaphana* under symbiotic and aposymbiotic conditions.

## Repository overview

This repository contains the scripts used for:

- preprocessing Nanopore DRS data
- reference-guided transcriptome reconstruction and quantification
- differential gene and transcript expression analysis
- differential isoform usage and isoform switch analysis
- poly(A) tail-length analysis
- transcript-resolved epitranscriptomic analysis
- Gene Ontology enrichment analysis

## Project structure

```bash
scripts/
├── 01_preprocessing/
├── 02_isoquant/
├── 03_deg_det/
├── 04_dtu_isoswitch/
├── 05_polya/
├── 06_epitranscriptome/
└── 07_go_enrichment/
```

Each folder contains the code for the corresponding analysis step.

## Script modules

### `01_preprocessing`
Preprocessing of Nanopore DRS reads, including input preparation, alignment-related steps, and intermediate file generation.

### `02_isoquant`
Reference-guided transcriptome reconstruction and quantification using IsoQuant.

### `03_deg_det`
Gene-level and transcript-level differential expression analyses.

### `04_dtu_isoswitch`
Differential transcript usage, isoform switching, and switch consequence analyses.

### `05_polya`
Transcript-level poly(A) tail-length estimation and downstream comparative analyses.

### `06_epitranscriptome`
Transcript-resolved epitranscriptomic analysis, including candidate RNA modification analyses.

### `07_go_enrichment`
Gene Ontology enrichment analyses for selected gene or transcript sets.

## Analysis workflow

The analyses were developed for *Aiptasia diaphana* Nanopore DRS data from two biological conditions:

- symbiotic
- aposymbiotic

Typical workflow:

1. preprocessing  
2. IsoQuant-based transcriptome reconstruction and quantification  
3. DEG/DET analysis  
4. DTU / isoform switch analysis  
5. poly(A) tail-length analysis  
6. epitranscriptomic analysis  
7. GO enrichment  

## Main software

- IsoQuant
- edgeR
- IsoformSwitchAnalyzeR
- eggNOG-mapper
- TandemMod
- Python
- R

## Notes

- Outputs from one step are typically used as inputs for downstream analyses.
- File paths, software environments, and computational settings may need to be adapted before reuse.
- This repository is intended to document the analysis workflow used in the associated manuscript.

## Contact
- **Constance**
- **huawen.zhong@kaust.edu.sa**
