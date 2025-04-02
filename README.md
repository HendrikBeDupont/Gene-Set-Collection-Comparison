# Gene Set Collection Comparison for Functional Enrichment Analysis in vEDS

This repository contains an R-based analysis pipeline designed to perform **functional bioinformatics analysis** on transcriptomic data, with a specific focus on comparing gene set collections. Using a dataset of fibroblasts derived from vEDS patients (GSE239914), the pipeline identifies **differentially expressed genes (DEGs)**, conducts **gene set enrichment analysis**, and evaluates how different collections influence biological interpretation.

In addition to the standard enrichment workflow, this pipeline introduces a systematic **comparison of five widely used gene set databases**—GO:BP, KEGG, Reactome, MSigDB (Hallmark), and WikiPathways—based on pathway count, vEDS relevance, and clustering interpretability. This comparison helps highlight the impact of database selection on biological insights, particularly in rare disease research.

# Overview

This project includes the following key steps:

✔ **Data Import and Preprocessing** of RNA-seq differential expression data  
✔ **Volcano Plot Generation** for DEG visualisation and thresholding  
✔ **Over-representation analysis (ORA)** across multiple gene set collections  
✔ **Pathway Relevance Scoring** using a custom vEDS-related keyword screen  
✔ **Functional Clustering Analysis** to evaluate redundancy and interpretability  
✔ **Protein-Protein Interaction (PPI) Network Analysis** using Cytoscape + stringApp  
✔ **Comparative Assessment** of biological coverage and clarity of enrichment output  

All analyses are performed in **R**, and integration with **Cytoscape** supports PPI visualisation and enrichment mapping.

# Dataset

- **GSE239914** – Fibroblast expression profiles from 18 vEDS patients and 36 healthy controls (available via GEO).

# Authors

Originally developed and maintained by [Hendrik Dupont](https://github.com/HendrikBeDupont), based on a pipeline by [Dr. Summer-Kutmon](https://github.com/mkutmon), with project-specific extensions for gene set comparison and rare disease pathway evaluation.

# Purpose

This project was carried out as part of a research assignment to answer the question:  
**"How do different gene set collections affect biological interpretation, and which is most suitable for studying vEDS?"**

# Keywords

vEDS, gene set enrichment, transcriptomics, pathway analysis, extracellular matrix
