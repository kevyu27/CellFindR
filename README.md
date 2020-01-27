# CellFindR
scRNAseq tools for cell clustering

## Introduction
CellFindR is an algorithm that iteratively applies Louvain clustering to single cell datasets in an automated fashion to identify biologically meaningful cell subpopulations.  It utilizes many of the functions within Seurat and is implemented as a wraparound on that package. CellFindR's outputs include matricies with gene expression information for each subpopulation as well as automated generation of violin plots and feature plots for the top differentially expressed genes.  It also includes additional data visualization functionalities.  CellfindR was developed by Kevin Yu in the Tward lab at UCSF.

## How to use:
CellFindR uses wrapper functions on Seurat and other R packages. Please install the required packages. 
To use:
run all the CellfindR functions included in the 'CellFindR.R file'. 
open template 'cellfindR_function_scripts.R' to walk through how to run your data through the pipeline. 

### Version 2.0.0 1/26/2020
updated to work on Seurat 3.1


### Version 1.0.0 08/20/2019
published onto github! 

### Reference:

If you find CellfindR useful, please reference:

Development of the Mouse and Human Cochlea at Single Cell Resolution
Kevin Shengyang Yu, Stacey M. Frumm, Jason S. Park, Katharine Lee, Daniel M. Wong, Lauren Byrnes, Sarah M. Knox, Julie B. Sneddon, Aaron D. Tward

doi: https://doi.org/10.1101/739680
BioRxiv https://www.biorxiv.org/content/10.1101/739680v2!
 

**any inquiries or debugging questions ** please contact

 *kyu@ucsf.edu* or
 *aaron.tward@ucsf.edu*
 
 
