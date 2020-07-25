# SARS-CoV2 Longitudinal Systems-level analyses 
COVID-19 study from acute to recovered patients: single-level analyses and omics integration.

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Repo description](#repo-description)
* [Single-level analyses](#single-level-analyses)
* [Omics integration](#omics-integration)
* [Figures](#figures)
* [Setup](#setup)

## General info
This project used multiple omics data:
- Plasma protein expression (Olink - NPX values)
- Cell abundance (CyTOF - Grid cell abundance)
	
## Dependencies
Project is created with:
* RStudio version: 3.6.0

## Repo description
- ```MOFA_Helsinki.R``` used to train MOFA model and some functions to uncover sources of variation 
- ```MEM_clean.R``` used for mixed-effect modelling 
- ```cell_corr.R``` uses cell abundance to build spearman correlation plots
- ```Single_features.R``` used to plot trends for all datasets
- ```Sample_tpt.R``` used to plot general timeline of sample aquisition for Olink and CyTOF

## Single-level analyses
### Grid
- Grid is a supervised learning algorithm based on t-SNE implementation. It uses the manual classification of cells in CyTOF samples and then the automatic classification of cells in new samples through the use of machine learning techniques based on these manual classifications.
- To run Grid, install it using:
```
$ pip install cellgrid
```
### Mixed-Effect modelling
- Partial-bayesian mixed-effect modelling to account for covariates 
- ```MEM_clean.R``` extracts confounding variables for downstream use 

## Omics integration
### MOFA
- Multi-Omics Factor Analysis (MOFA) was used in this study in order to deconvolute the main sources of variation in the differents sets of data mentioned above. For more information, read their published Methods paper [Argelaguet et al. (2018)](https://www.embopress.org/doi/10.15252/msb.20178124). 
- MOFA is publicly accessible here: https://github.com/bioFAM/MOFA 

## Figures
### Spearman Correlation 
- ```cell_corr.R``` uses cell abundance dataframe built from Grid that is sub-setted by grouped days for correlation
- Option to re-order one matrix in accordance to the other for comparison purposes

## Setup
To run this project, install it locally using devtools:

```
$ install.packages('devtools')
$ library(devtools)
$ install_github('rodriluc/SARS-CoV2_study')
$ library(SARS-CoV2_study)
```
