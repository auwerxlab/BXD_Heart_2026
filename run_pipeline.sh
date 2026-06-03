#!/bin/bash


# Run script from the project main directory by executing "bash run_pipeline.sh"


#########################
## Create conda environment
#########################

ENV_NAME="BXD_Heart_R_4_4_3"
R_VERSION="4.3.3"

# Make conda available in this non-interactive shell
source "$(conda info --base)/etc/profile.d/conda.sh"

# Create the environment only if it does not exist
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  conda create -y -n "$ENV_NAME" -c conda-forge "r-base=$R_VERSION"
fi

# Activate it for the rest of the script
conda activate "$ENV_NAME"


#########################
## Phenotypes
#########################
# phenotypic data preparation and analysis
Rscript --vanilla ./Scripts/Phenotypes/01_pheno_data_preparation_and_analysis.R


#########################
## Transcriptomics
#########################
# raw fq data QC, STAR-mapping and QC, counts table extraction
Rscript --vanilla ./Scripts/Transcriptomics/01_RNAseq_processing.R
# differential expression analysis (DEA), followed by GSEA
Rscript --vanilla ./Scripts/Transcriptomics/02_RNAseq_DEA_GSEA.R
# single-cell deconvolution
Rscript --vanilla ./Scripts/Transcriptomics/03_RNAseq_single_cell_deconvolution.R


#########################
## Proteomics
#########################
# data preparation
Rscript --vanilla ./Scripts/Proteomics/01_proteomics_data_prep.R
# differential abundance analysis (DAA), followed by GSEA
Rscript --vanilla ./Scripts/Proteomics/02_proteomics_DAA_GSEA.R


#########################
## Lipidomics
#########################
# data preparation, lipids annotation, data aggregation
Rscript --vanilla ./Scripts/Lipidomics/01_lipidomics_data_prep.R
# differential abundance analysis (DAA), followed by GSEA
Rscript --vanilla ./Scripts/Lipidomics/02_lipidomics_DAA_GSEA.R
# WGCNA applied to lipidomic data to compute co-abundance networks
Rscript --vanilla ./Scripts/Lipidomics/03_lipidomics_WGCNA.R


#########################
## Multi-omics
#########################
# Within- and across-omics features correlation
Rscript --vanilla ./Scripts/Multiomics/01_omics_correlation.R
# GSEA using as ranking the correlation of selected gene with all other 
# LV-expressed genes
Rscript --vanilla ./Scripts/Multiomics/02_gene_gene_corr_GSEA.R
# Multivariate (Mahalanobis) distance for strain outliers detection, within and 
# across lipidomics and phenotypic layers
Rscript --vanilla ./Scripts/Multiomics/03_mahalanobis_and_z_score.R


#########################
## QTL mapping
#########################
# BXD genetic markers and genetic maps preparation
Rscript --vanilla ./Scripts/QTL_and_misc/01_QTL_genetic_markers_preparation.R
# Phenotypic data preparation for QTL mapping
Rscript --vanilla ./Scripts/QTL_and_misc/02_QTL_preparation_phenotypes.R
# Transcriptomics data preparation for QTL mapping
Rscript --vanilla ./Scripts/QTL_and_misc/03_QTL_preparation_RNAseq.R
# Proteomics data preparation for QTL mapping
Rscript --vanilla ./Scripts/QTL_and_misc/04_QTL_preparation_proteomic.R
# Lipidomics data preparation for QTL mapping
Rscript --vanilla ./Scripts/QTL_and_misc/05_QTL_preparation_lipidomic.R
# QTL mapping
Rscript --vanilla ./Scripts/QTL_and_misc/06_QTL_mapping.R
# QTL results extraction and formatting
Rscript --vanilla ./Scripts/QTL_and_misc/07_QTL_results_extraction.R


#########################
## Other systems genetics 
#########################
# Genes physical and QTL colocalization analysis within lQTL peak's CI
Rscript --vanilla ./Scripts/QTL_and_misc/08_colocalization_analysis_lipids.R
# Features (phenotypes, lipids) observed variance decomposition within and 
# across diets
Rscript --vanilla ./Scripts/QTL_and_misc/09_variance_decomposition.R
# QTL pleiotropy analysis at lQTL peak loci, with genes displaying colocalized 
# (e/p)QTL peaks
Rscript --vanilla ./Scripts/QTL_and_misc/10_pleiotropy.R
# Forward and reverse mediation - partial scan only with lipid/gene features 
# with overlapping QTL peaks and full scan for selected lipids testing all 
# possible gene mediators (all LV-expressed genes)
Rscript --vanilla ./Scripts/QTL_and_misc/11_mediation.R






