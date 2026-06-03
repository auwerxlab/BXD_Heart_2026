

library(plyr)
library(data.table)


source("./Scripts/QTL_and_misc/QTL_fun.R")


dirQTL_input <- "./Data/QTL_mapping/qtl_processed_inputs"
if(!dir.exists(dirQTL_input)){
  dir.create(dirQTL_input, recursive = T)
}


################################################################################
## QTL preparation
################################################################################

################################################################################
## Transcriptomic data for eQTL. mRNA expression.
################################################################################


metadata               <- readRDS("./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")
mRNA_counts            <- readRDS("./Data/RNAseq_processing/mRNA_counts/geneLevel/STAR_RNAseq_count_matrix.RDS")
metadata$strainDiet    <- paste0(metadata$strain, "_", metadata$diet)
metadata               <- metadata[metadata$sample_name %in% colnames(mRNA_counts), ]
metadata               <- metadata[match(colnames(mRNA_counts), metadata$sample_name), ]
rownames(metadata)     <- metadata$sample_name
colnames(metadata)     <- gsub("\\.", "_", gsub("\\[", "", gsub("\\]\\.|\\/", "_", colnames(metadata))))
colnames(metadata)     <- make.names(colnames(metadata))
mRNA_counts            <- as.matrix(mRNA_counts)
stopifnot(identical(rownames(metadata), colnames(mRNA_counts)))
y                      <- edgeR::DGEList(counts = mRNA_counts, samples = meta)
keep.exprs             <- edgeR::filterByExpr(y, group = "strain_diet")
table(keep.exprs)
y                      <- y[keep.exprs, , keep.lib.sizes = FALSE]
y                      <- edgeR::calcNormFactors(y, method = "TMM")
mRNA_counts.TMM        <- edgeR::cpm(y)
dim(mRNA_counts.TMM)

overwrite  <- F
tt <- lapply(c("HFD", "CD", "all"), function(diet){
  # diet <- "HFD"
  
  saveDir <- paste0(dirQTL_input, "/RNAseq", "/", ifelse(diet == "all", "CD_HFD", diet))
  if(!dir.exists(saveDir)){
    dir.create(saveDir, recursive = T)
  }
  if(file.exists(paste0(saveDir, "/BXD.json")) & !overwrite){
    cat("QTL files already created, skipped...\n")
    return()
  } else{
    cat("Creating QTL files...\n")
  }
  
  
  df.counts <- t(mRNA_counts.TMM)
  if(diet != "all"){
    df.counts <- df.counts[grepl(diet, rownames(df.counts)), ]
  }
  
  tmp       <- meta[meta$correctSampleName %in% rownames(df.counts), ]
  geno      <- readRDS(paste0(dirQTL_input, "/BXD_geno.RDS"))
  geno.diet <- geno[, "marker", drop = F]
  for(i in 1:nrow(tmp)){
    # print(i)
    # i <- 1
    geno.diet <- cbind(geno.diet, geno[, tmp$correctStrain[i], drop = F])
    colnames(geno.diet)[ncol(geno.diet)] <- paste0(tmp$correctSampleName[i])
  }
  crossinfo                 <- meta[meta$correctSampleName %in% rownames(df.counts), ]
  crossinfo                 <- dplyr::select(crossinfo, c("correctSampleName", "correctStrain", "strainConditionReplicate", "diet"))
  colnames(crossinfo)       <- c("id", "strain", "strainReplicate", "diet")
  crossinfo$cross_direction <- "BxD"
  
  samplescovar              <- as.data.frame(dplyr::select(crossinfo, c("id", "strain", "diet", "strainReplicate")))
  samplescovar$sex          <- "m"
  rownames(samplescovar)    <- samplescovar$id
  crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -c("strain", "diet", "strainReplicate")))
  rownames(crossinfo)       <- crossinfo$id
  
  
  # transform to normal distribution
  norm                             <- transform2normal_byCol_INT(df.counts)
  transf.type.df                   <- norm$transfType.df
  # remove phenotypes with less than 15 numeric measurements
  transf.df                        <- norm$transf.df
  tmpLogi                          <- apply(transf.df, 2, function(kk) (sum(!is.na(kk)) + sum(!is.infinite(kk) & !is.na(kk))) > 15)
  # table(tmpLogi)
  transf.df                        <- transf.df[, tmpLogi]
  rownames(transf.df)              <- rownames(df.counts)
  df.counts                        <- transf.df
  
  phenocovar                  <- data.frame(phenoID = colnames(df.counts), stringsAsFactors = F)
  rownames(phenocovar)        <- phenocovar$phenoID
  df.counts                   <- df.counts[, rownames(phenocovar)]
  df.counts[is.na(df.counts)] <- NA
  
  file.copy(from = paste0(dirQTL_input, "/BXD_gmap.csv"), to = saveDir)
  file.copy(from = paste0(dirQTL_input, "/BXD_pmap.csv"), to = saveDir)
  na_strings <- c("-", "NA")      # the default (see ?qtl2::write_control_file)
  
  qtl2::write_control_file(paste0(saveDir, "/BXD.json"),
                           description     = "BXD mouse data from GeneNetwork",
                           crosstype       = "risib",
                           geno_file       = "BXD_geno.csv",
                           geno_transposed = TRUE,
                           geno_codes      = list(B = 1, D = 2),
                           xchr            = "X",
                           sex_covar       = "sex",
                           sex_codes       = c(f = "female", m = "male"),
                           gmap_file       = "BXD_gmap.csv",
                           pmap_file       = "BXD_pmap.csv",
                           pheno_file      = "BXD_pheno.csv",
                           covar_file      = "BXD_covar.csv",
                           phenocovar_file = "BXD_phenocovar.csv",
                           crossinfo_file  = "BXD_crossinfo.csv",
                           crossinfo_codes = c("BxD" = 0),
                           alleles         = c("B", "D"),
                           na.strings      = na_strings,
                           overwrite       = TRUE)
  cat("\n\n")
  NA
})


################################################################################
## (3) Transcriptomic data with cell fractions for cf-eQTL (or cf-tQTL)
##     https://www.nature.com/articles/s41588-021-00909-9
################################################################################

# load new sc-deconvolution results - PCA projections
df.scdec.pca.coord.list <- lapply(c("CD", "HFD", "CD_HFD"), function(dd){
  # dd <- "CD"
  
  files.pca           <- list.files("./Data/singleCell_deconvolution/cells_proportion_pca/", full.names = T)
  df.scdec.pca.coord  <- files.pca[(grepl("__all_genes__", files.pca) & grepl("__inferred_cell_size", files.pca) & 
                                      grepl("__music2__", files.pca) & grepl("\\/pca_projections__", files.pca) &
                                      grepl(paste0("bulkRNAseq_", dd, "___scRNAseq"), files.pca))]
  df.scdec.pca.coord  <- readRDS(df.scdec.pca.coord)
  df.scdec.pca.expVar <- files.pca[(grepl("__all_genes__", files.pca) & grepl("__inferred_cell_size", files.pca) & 
                                      grepl("__music2__", files.pca) & grepl("\\/pca_explained_var__", files.pca) &
                                      grepl(paste0("bulkRNAseq_", dd, "___scRNAseq"), files.pca))]
  df.scdec.pca.expVar <- readRDS(df.scdec.pca.expVar)
  pcs.select          <- names(df.scdec.pca.expVar)[1:which(cumsum(df.scdec.pca.expVar) > 90)[1]]
  df.scdec.pca.coord  <- df.scdec.pca.coord[, pcs.select]
  df.scdec.pca.coord  <- as.data.frame(df.scdec.pca.coord)
  colnames(df.scdec.pca.coord) <- paste0("cf__", colnames(df.scdec.pca.coord))
  summary(df.scdec.pca.coord)
  round(apply(df.scdec.pca.coord, 2, mean), 4)
  round(apply(df.scdec.pca.coord, 2, sd), 4)
  # print(dim(df.scdec.pca.coord))
  df.scdec.pca.coord
  
})
names(df.scdec.pca.coord.list) <- c("CD", "HFD", "CD_HFD")


dir_source <- paste0(dirQTL_input, "/RNAseq")
dir_dest   <- paste0(dirQTL_input, "/RNAseq_cf_pca")

# copy data folder from previous step
if(dir.exists(dir_dest)){
  system(paste0("rm -r ", dir_dest))
}
cmd <- paste0("cp -r ", dir_source, "/ ", dir_dest)
system(cmd)

tt <- lapply(dir(dir_dest, full.names = T), function(x){
  # x <- dir(dir_dest, full.names = T)[1]
  print(x)
  
  diet.x    <- gsub(".*\\/", "", x)
  covarFile <- list.files(x, pattern = "_covar.csv", full.names = T)
  covar.df  <- qtl2::fread_csv(covarFile)
  rowsOrder <- rownames(covar.df)
  
  covar.df <- merge(covar.df, df.scdec.pca.coord.list[[diet.x]], by = "row.names", all.x = T, all.y = F)
  
  colnames(covar.df)[colnames(covar.df) == "Row.names"] <- "id"
  covar.df <- covar.df[match(rowsOrder, covar.df$id), ]
  stopifnot(all(rowsOrder == covar.df$id))
  
  print(dim(df.scdec.pca.coord.list[[diet.x]]))
  print(dim(covar.df))
  print(head(covar.df, 2))
  print("")
  
  qtl2convert::write2csv(covar.df, paste0(x, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
  
})


