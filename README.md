# MARG-signature

This repository contains the R scripts supporting the main computational analyses reported in the manuscript:

**“Machine learning-based integration identifies a macrophage activation-related prognostic signature and therapeutic implications in lung adenocarcinoma.”**

## Overview

The scripts in this repository cover the main data-processing and computational analysis procedures used in the study, including:

- TCGA-LUAD and GEO data preprocessing
- Macrophage activation-related score calculation
- Differential expression analysis
- Prognostic signature construction and validation
- Single-cell RNA sequencing quality control
- Pseudotime trajectory analysis
- Immune response-related analyses
- Drug sensitivity prediction using CTRP and PRISM data

Only code directly relevant to the results reported in the manuscript is included. Raw datasets, intermediate files, unpublished exploratory analyses, and code related to ongoing studies are not provided.

## Repository structure

```text
MARG-signature/
├── 2-RawData/
│   ├── 1-TCGA/
│   │   ├── 1-PD.R
│   │   └── 2-TCGA-LUAD.R
│   ├── 2-GEO/
│   │   ├── A-GEO Datasets.R
│   │   ├── B-GeneID.R
│   │   ├── C-Combined Data.R
│   │   └── GSE11969/code/agilent_raw.R
│   └── 3-表型基因/
│       └── 表型处理.R
├── 3.1-Score_Construction/
│   └── code/ssgsea.R
├── 3.2-Score_Validation/
│   └── code/ssgsea.R
├── 4-Diffanalysis/
│   └── code/deseq2.R
├── 11-Drug/
│   └── code/
│       ├── 1-tide.R
│       ├── 2-tcia.R
│       ├── 3-IMV.R
│       └── 4-CTRP-PRISM.R
├── 13-SC_QC/
│   └── SC_QC.R
└── 14-Sc_track/
    └── code/
        ├── monocle.R
        └── orderCells.R
