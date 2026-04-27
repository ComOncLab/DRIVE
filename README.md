### DRIVE:  database of drug-induced transcriptomic response, alternative splicing, and splicing-derived neoantigens in cancer cell lines

Although drugs, alternative splicing, and neoantigens are intimately interconnected, existing resources fail to integrate these three dimensions. To bridge this gap, we constructed DRIVE (***D***rug-induced ***R***NA profile, **I**soform ***V***ariability and ***E***pitopes), a comprehensive resource integrating, drug induced transcriptomic responses, alternative splicing and neoantigen landscapes. We systematically retrieved drug-treated and matched control cancer cell line transcriptomic datasets from the GEO database, utilizing Large Language Models (LLMs) combined with rigorous manual inspection to ensure high-fidelity metadata. By implementing a standardized pipeline for raw data processing, we quantified drug-induced differential gene expression, AS alterations, and predicted splice-neoantigens across thousands of conditions. 

#### Contents

This repository contains code for the processing pipeline of raw data, reproducible data, code, and analysis reports for the main results in the paper, as well as code for building the DRIVE Shiny App.

1. Processing pipeline of raw data
   - [processing matched treatment and control samples](./scripts/raw_data_processing/run_gse): Extract FASTQ, QC, mapping, Transcript and gene expression quantification, HLA typing.
   - [translate into protein sequences](./scripts/raw_data_processing/run_jcast): predict peptides using [Jcast](https://github.com/ed-lau/jcast)
   - [predicting HLA binding](./scripts/raw_data_processing/run_binding): predicting HLA binding using NetMHCpan4.1
2. The data for analysis in manuscript and shiny app is available in Zenodo with DOI: 
3. The analysis report can be view online in https://comonclab.github.io/DRIVE/

