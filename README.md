Sex-Specific Transcriptomic and Behavioural Responses of Parhyale hawaiensis to Phenanthrene

Description
This repository contains the R scripts and analytical pipelines used to evaluate the molecular and behavioural toxicity of phenanthrene (PHE) in the marine amphipod Parhyale hawaiensis. The code integrates survival analysis, dose-response modelling for behavioural impairments (feeding suppression, mating latency), and sex-specific differential gene expression analysis.

Data Availability
The raw transcriptomic sequence read data analysed by these scripts reside in the NCBI Sequence Read Archive (SRA) under BioProject accession PRJNA1424690 (BioSamples SAMN55391704 to SAMN55391793).

System Requirements
The analysis requires a standard computer running R (version 4.4.2 or higher).
Primary packages required include:
DESeq2 (v1.42)
survival (for Kaplan-Meier analysis)
drc (for dose-response modelling)
clusterProfiler / topGO (for functional enrichment)
ggplot2 (for data visualisation)

Repository Structure
The repository divides the analysis into five sequential scripts:

01_Import_and_QC.R: Imports HTSeq-count matrices, filters low-count transcripts, and applies variance-stabilised transformations (VST) to generate Principal Component Analysis (PCA) plots.

02_DESeq_Modelling.R: Executes the Likelihood Ratio Test (LRT) using the multifactorial design (~ Sex + Treatment + Sex:Treatment). This script identifies transcripts exhibiting significant sex-specific responses to phenanthrene and calculates pairwise contrasts (e.g., Female Phenanthrene vs. Female DMSO).

03_Functional_Enrichment.R: Performs Gene Ontology (GO) enrichment analysis on the differentially expressed genes to identify significantly impacted biological processes (e.g., chitin metabolism, neurological function).

04_Figure_Generation.R: Contains the plotting functions to generate publication-ready graphics, including dual volcano plots, hierarchically clustered heatmaps, and boxplots for specific molecular markers (e.g., NPF receptor, chitinase).

05_Behavioural_Statistics.R: Processes the in vivo phenotypic data. It calculates the 96 h LC50 values, models feeding suppression EC50s using bias-reduced logistic regression, and performs log-rank tests for mating latency.

Usage
Users must execute the scripts sequentially. Ensure the working directory contains the raw count matrices and the associated metadata table (sample_metadata.csv) before running script 01.
