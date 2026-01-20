# Microbiome Analysis for Medeiros et al., 2024

[![DOI](https://img.shields.io/badge/DOI-10.3389%2Ffnins.2024.1309075-blue)](https://doi.org/10.3389/fnins.2024.1309075)

Analysis code for the paper: **"Slowing Alzheimer's disease progression through probiotic supplementation"** published in *Frontiers in Neuroscience*.

## Overview

This repository contains R code for processing QIIME2 metagenomic sequencing data to analyze bacterial phylum-level abundance differences in Alzheimer's disease patients. The analysis compares Bacteroidetes and Firmicutes relative abundance across baseline, control (ADC), and probiotic (ADP) treatment cohorts, investigating the gut-brain axis in neurodegeneration.

## Repository Structure

```
├── README.md                           # This file
├── analysis/
│   └── R-Processing-and-Figure.md      # Main analysis workflow
├── data/                               # Input data (not included)
│   └── README.md                       # Data requirements
└── figures/                            # Output figures
    └── README.md                       # Figure descriptions
```

## Requirements

### R Version
- R >= 4.1.0

### R Packages

Install all dependencies using `pacman`:

```r
install.packages("pacman")
library(pacman)
pacman::p_load(
  decontam, phyloseq, data.table, ggplot2, BiocManager,
  qiime2R, DESeq2, tidyverse, RColorBrewer, viridis,
  vegan, pheatmap, patchwork, ggpubr, lme4, nlme,
  microViz, rstatix, ANCOMBC, microbiome
)
```

**Note:** Some packages require Bioconductor. Install BiocManager first if needed:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("phyloseq", "DESeq2", "ANCOMBC", "microbiome"))
```

## Input Data

The analysis requires the following QIIME2 output files (not included in this repository):

| File | Description |
|------|-------------|
| `table.qza` | Feature table (ASV/OTU counts) |
| `rooted-tree.qza` | Phylogenetic tree |
| `taxonomy.qza` | Taxonomic classifications |
| `Alzheimers_metadata.txt` | Sample metadata with cohort assignments |

Place these files in the `data/` directory before running the analysis.

## Analysis Workflow

1. **Data Import**: Load QIIME2 artifacts using `qiime2R`
2. **Taxonomic Aggregation**: Aggregate features at phylum level
3. **Normalization**: Calculate relative abundances
4. **Subset**: Focus on Bacteroidetes and Firmicutes
5. **Statistical Testing**: Pairwise t-tests with Bonferroni correction
6. **Visualization**: Generate publication-ready box plots

## Key Findings

- Significant difference in Bacteroidetes abundance between Probiotic and Baseline groups (p = 0.015, adjusted p = 0.046)
- No significant differences in Firmicutes abundance across cohorts
- Data confirmed to be normally distributed (Shapiro-Wilk test, p > 0.05 for all groups)

## Output

The analysis generates:
- `phylum_comp.png` - Box plot of phylum relative abundance by cohort
- `phylum_comp.tiff` - High-resolution TIFF for publication

## Usage

1. Clone this repository
2. Place required data files in `data/`
3. Open `analysis/R-Processing-and-Figure.md` in RStudio
4. Run the code chunks sequentially

## Authors

**Code:** Todd Testerman

**Paper:** Destynie Medeiros, Kristina McMurry, Melissa Pfeiffer, Kayla Newsome, Todd Testerman, Joerg Graf, Adam C. Silver, and Paola Sacchetti

## Citation

If you use this code or find it helpful, please cite the original paper:

> Medeiros D, McMurry K, Pfeiffer M, Newsome K, Testerman T, Graf J, Silver AC and Sacchetti P (2024) Slowing Alzheimer's disease progression through probiotic supplementation. *Front. Neurosci.* 18:1309075. doi: [10.3389/fnins.2024.1309075](https://doi.org/10.3389/fnins.2024.1309075)

### BibTeX

```bibtex
@article{medeiros2024slowing,
  title={Slowing Alzheimer's disease progression through probiotic supplementation},
  author={Medeiros, Destynie and McMurry, Kristina and Pfeiffer, Melissa and Newsome, Kayla and Testerman, Todd and Graf, Joerg and Silver, Adam C and Sacchetti, Paola},
  journal={Frontiers in Neuroscience},
  volume={18},
  pages={1309075},
  year={2024},
  publisher={Frontiers},
  doi={10.3389/fnins.2024.1309075}
}
```

## License

This project is available for academic and research use.
