
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
## Proteomic data aggregated from Coon for pQTL
## protein-level intensities
################################################################################



###########
## prepare for QTL mapping
###########

prot.df         <- readRDS("./Data/input_data/proteomics/coon_formatted/raw_intensity_table_samples_excluded.RDS")
samples.meta.df <- readRDS("./Data/input_data/proteomics/coon_formatted/sample_metadata.RDS")
prot.meta.df    <- readRDS("./Data/input_data/proteomics/coon_formatted/proteins_metadata.RDS")

prot.df.aggr <- reshape2::melt(as.matrix(prot.df))
prot.df.aggr <- merge(prot.df.aggr, samples.meta.df, by.x = "Var2", by.y = "sample", all.x = T, all.y = F)
prot.df.aggr <- data.table(prot.df.aggr)[, list(value = mean(value, na.rm = T)),
                                         by = c("strain", "diet", "Var1")]
prot.df.aggr$strainDiet   <- paste0(prot.df.aggr$strain, "_", prot.df.aggr$diet)
prot.df.aggr.meta         <- unique(dplyr::select(prot.df.aggr, c("strain", "diet", "strainDiet")))
prot.df.aggr.meta$sample  <- prot.df.aggr.meta$strainDiet
prot.df.aggr.meta         <- type.convert(as.data.frame(prot.df.aggr.meta), as.is = T)
prot.df.aggr              <- dcast(prot.df.aggr, Var1~strainDiet, value.var = "value")
prot.df.aggr              <- as.data.frame(prot.df.aggr)
prot.df.aggr[is.na(prot.df.aggr)] <- NA
rownames(prot.df.aggr)    <- prot.df.aggr$Var1
prot.df.aggr              <- dplyr::select(prot.df.aggr, -"Var1")

# since all samples of a strains_diet condition were matched to a unique batch, we can add batch meta to the aggregated data metadata
prot.df.aggr.meta           <- merge(prot.df.aggr.meta, unique(dplyr::select(samples.meta.df, c("batch", "strain_diet", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2"))), by.x = "strainDiet", by.y = "strain_diet", all = F)
rownames(samples.meta.df)   <- samples.meta.df$sample
rownames(prot.df.aggr.meta) <- prot.df.aggr.meta$strainDiet

samples.meta.df   <- samples.meta.df[match(colnames(prot.df), rownames(samples.meta.df)), ]
prot.df.aggr.meta <- prot.df.aggr.meta[match(colnames(prot.df.aggr), rownames(prot.df.aggr.meta)), ]

identical(samples.meta.df$sample, colnames(prot.df))
all(samples.meta.df$sample %in% colnames(prot.df))

identical(colnames(prot.df), rownames(samples.meta.df))
all(samples.meta.df$sample %in% colnames(prot.df))

identical(prot.df.aggr.meta$strainDiet, colnames(prot.df.aggr))
all(prot.df.aggr.meta$strainDiet %in% colnames(prot.df.aggr))

overwrite  <- F
tt         <- lapply(c("HFD", "CD", "all"), function(diet){
  # diet <- "all"
  
  tt <- lapply(c("individuals", "average"), function(aggr){
    # aggr <- "individuals"
    # aggr <- "average"
    print(c(diet, aggr))
    
    saveDir <- paste0(dirQTL_input, "/", aggr, "_coonProteomic", "/", ifelse(diet == "all", "CD_HFD", diet))
    if(!dir.exists(saveDir)){
      dir.create(saveDir, recursive = T)
    }
    if(file.exists(paste0(saveDir, "/BXD.json")) & !overwrite){
      cat("QTL files already created, skipped...\n")
      return()
    } else{
      cat("Creating QTL files...\n")
    }
    
    if(aggr == "individuals"){
      df.intensity <- t(prot.df)
    } else{
      df.intensity <- t(prot.df.aggr)
    }
    
    
    if(diet != "all"){
      df.intensity <- df.intensity[grepl(diet, rownames(df.intensity)), ]
    }
    df.intensity <- df.intensity[, !is.na(colnames(df.intensity))]
    
    if(aggr == "individuals"){
      tmp       <- samples.meta.df[samples.meta.df$sample %in% rownames(df.intensity), ]
    } else{
      tmp       <- prot.df.aggr.meta[prot.df.aggr.meta$sample %in% rownames(df.intensity), ]
    }
    
    geno      <- readRDS(paste0(dirQTL_input, "/BXD_geno.RDS"))
    geno.diet <- geno[, "marker", drop = F]
    for(i in 1:nrow(tmp)){
      # print(i)
      geno.diet <- cbind(geno.diet, geno[, as.character(tmp$strain[i]), drop = F])
      colnames(geno.diet)[ncol(geno.diet)] <- paste0(tmp$sample[i])
    }
    
    if(aggr == "individuals"){
      crossinfo                 <- samples.meta.df[samples.meta.df$sample %in% rownames(df.intensity), ]
      crossinfo                 <- dplyr::select(crossinfo, c("strain", "diet", "dietReplicate", "organ", "organRegion", "sampleUniqueID", "strain_diet", "sample", "sample_name", "batch", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2"))
      colnames(crossinfo)       <- c("strain", "diet", "dietReplicate", "organ", "organRegion", "sampleUniqueID", "strain_diet", "id", "sample_name", "batch", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2")
      crossinfo$cross_direction <- "BxD"
      
      samplescovar              <- as.data.frame(dplyr::select(crossinfo, c("id", "sample_name", "strain", "diet", "dietReplicate", "organ", "organRegion", "sampleUniqueID", "strain_diet", "batch", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2")))
      samplescovar$sex          <- "m"
      rownames(samplescovar)    <- samplescovar$id
      crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -c("strain", "sample_name", "diet", "dietReplicate", "organ", "organRegion", "sampleUniqueID", "strain_diet", "batch", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2")))
      rownames(crossinfo)       <- crossinfo$id
    } else{
      crossinfo                 <- prot.df.aggr.meta[prot.df.aggr.meta$sample %in% rownames(df.intensity), ]
      crossinfo                 <- dplyr::select(crossinfo, c("strain", "diet", "strainDiet", "sample", "batch", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2"))
      colnames(crossinfo)       <- c("strain", "diet", "strainDiet", "id", "batch", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2")
      crossinfo$cross_direction <- "BxD"
      
      samplescovar              <- as.data.frame(dplyr::select(crossinfo, c("id", "strain", "diet", "strainDiet", "batch", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2")))
      samplescovar$sex          <- "m"
      rownames(samplescovar)    <- samplescovar$id
      crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -c("strain", "diet", "strainDiet", "batch", "batch_median_raw", "batch_median_log2", "batch_avg_raw", "batch_avg_log2")))
      rownames(crossinfo)       <- crossinfo$id
    }
    
    
    saveRDS(df.intensity, paste0(saveDir, "/pheno_original_distribution.RDS"))
    
    
    # transform to normal distribution
    norm                             <- transform2normal_byCol_INT(df.intensity)
    transf.type.df                   <- norm$transfType.df
    
    # remove phenotypes with less than 15 strains with numeric measurements
    transf.df                        <- norm$transf.df
    strainsVec                       <- gsub("\\..*", "", rownames(df.intensity))
    tmpLogi                          <- apply(transf.df, 2, function(kk){
      # tmpLogi <- (sum(!is.na(kk)) + sum(!is.infinite(kk) & !is.na(kk))) > 15
      tmpLogi <- !is.na(kk) & (!is.infinite(kk) & !is.na(kk))
      length(unique(strainsVec[tmpLogi])) > 15
    })
    # table(tmpLogi)
    transf.df                        <- transf.df[, tmpLogi]
    rownames(transf.df)              <- rownames(df.intensity)
    df.intensity                     <- transf.df
    
    phenocovar                        <- data.frame(phenoID = colnames(df.intensity), stringsAsFactors = F)
    rownames(phenocovar)              <- phenocovar$phenoID
    df.intensity                      <- df.intensity[, rownames(phenocovar)]
    df.intensity[is.na(df.intensity)] <- NA
    
    samplescovar <- samplescovar[match(rownames(df.intensity), rownames(samplescovar)), ]
    crossinfo    <- crossinfo[match(rownames(df.intensity), rownames(crossinfo)), ]
    
    if(!identical(rownames(phenocovar), colnames(df.intensity))){
      stop("ERROR: non-matching dimnames 1")
    }
    if(!identical(rownames(samplescovar), rownames(df.intensity))){
      stop("ERROR: non-matching dimnames 2")
    }
    if(!identical(rownames(crossinfo), rownames(df.intensity))){
      stop("ERROR: non-matching dimnames 3")
    }
    if(!identical(colnames(geno.diet)[-1], rownames(df.intensity))){
      stop("ERROR: non-matching dimnames 4")
    }
    
    
    # identical(rownames(phenocovar), colnames(df.intensity))
    
    saveRDS(norm, paste0(saveDir, "/pheno_normal_distribution.RDS"))
    # saveRDS(transf.type.df, paste0(saveDir, "/pheno_normal_transformation_type.RDS"))
    qtl2convert::write2csv(geno.diet, paste0(saveDir, "/BXD_geno.csv"), overwrite = TRUE, comment = "Genotypes for BXD data")
    qtl2convert::write2csv(samplescovar, paste0(saveDir, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
    qtl2convert::write2csv(crossinfo, paste0(saveDir, "/BXD_crossinfo.csv"), overwrite = TRUE, comment = paste0("Cross info for BXD data\n#", "(all lines formed from cross between female B and male D)"))
    qtl2convert::write2csv(phenocovar, paste0(saveDir, "/BXD_phenocovar.csv"), overwrite = TRUE, comment = "Phenotype covariates (metadata) for BXD phenotype data")
    qtl2convert::write2csv(cbind(id = rownames(df.intensity), df.intensity), paste0(saveDir, "/BXD_pheno.csv"), overwrite = TRUE, comment = "BXD phenotype data")
    
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
  NA
})




################################################################################
## (5c) Average Coon proteomic data with cell fractions for cf-pQTL (or cf-pQTL)
##     https://www.nature.com/articles/s41588-021-00909-9
################################################################################

# load new sc-deconvolution results - PCA projections
df.scdec.pca.coord.list <- lapply(c("CD", "HFD", "CD_HFD"), function(dd){
  # dd <- "CD"
  
  files.pca           <- list.files("./Data/singleCell_deconvolution/cells_proportion_pca/", full.names = T)
  df.scdec.pca.coord  <- files.pca[(grepl("__all_genes__", files.pca) & grepl("__inferred_cell_size", files.pca) & 
                                      grepl("__music2__", files.pca) & grepl("\\/avg_pca_projections__", files.pca) &
                                      grepl(paste0("bulkRNAseq_", dd, "___scRNAseq"), files.pca))]
  df.scdec.pca.coord  <- readRDS(df.scdec.pca.coord)
  df.scdec.pca.expVar <- files.pca[(grepl("__all_genes__", files.pca) & grepl("__inferred_cell_size", files.pca) & 
                                      grepl("__music2__", files.pca) & grepl("\\/avg_pca_explained_var__", files.pca) &
                                      grepl(paste0("bulkRNAseq_", dd, "___scRNAseq"), files.pca))]
  df.scdec.pca.expVar <- readRDS(df.scdec.pca.expVar)
  pcs.select          <- names(df.scdec.pca.expVar)[1:which(cumsum(df.scdec.pca.expVar) > 90)[1]]
  df.scdec.pca.coord  <- df.scdec.pca.coord[, pcs.select]
  df.scdec.pca.coord  <- as.data.frame(df.scdec.pca.coord)
  colnames(df.scdec.pca.coord) <- paste0("cf__", colnames(df.scdec.pca.coord))
  summary(df.scdec.pca.coord)
  round(apply(df.scdec.pca.coord, 2, mean), 4)
  round(apply(df.scdec.pca.coord, 2, sd), 4)
  df.scdec.pca.coord
})
names(df.scdec.pca.coord.list) <- c("CD", "HFD", "CD_HFD")




dir_source <- paste0(dirQTL_input, "/average_coonProteomic")
dir_dest   <- paste0(dirQTL_input, "/average_coonProteomic_cf_pca")

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
  
  # rownames(covar.df) <- covar.df$Row.names
  colnames(covar.df)[colnames(covar.df) == "Row.names"] <- "id"
  covar.df           <- covar.df[match(rowsOrder, covar.df$id), ]
  stopifnot(all(rowsOrder == covar.df$id))
  
  print(dim(df.scdec.pca.coord.list[[diet.x]]))
  print(dim(covar.df))
  print(head(covar.df, 2))
  print("")
  
  qtl2convert::write2csv(covar.df, paste0(x, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
})



################################################################################
## (5d) Individuals Coon proteomic data with cell fractions for cf-eQTL (or cf-tQTL)
##     https://www.nature.com/articles/s41588-021-00909-9
## Cell fractions were estimated from pooled RNAseq samples and therefore not 
## available for all samples. We assume similar tissue composition for all samples
## in the same condition by using the same cell fraction estimated from pooled
## RNAseq data
################################################################################


# load new sc-deconvolution results - PCA projections
df.scdec.pca.coord.list <- lapply(c("CD", "HFD", "CD_HFD"), function(dd){
  # dd <- "CD"
  
  files.pca           <- list.files("./Data/singleCell_deconvolution/cells_proportion_pca/", full.names = T)
  df.scdec.pca.coord  <- files.pca[(grepl("__all_genes__", files.pca) & grepl("__inferred_cell_size", files.pca) & 
                                      grepl("__music2__", files.pca) & grepl("\\/avg_pca_projections__", files.pca) &
                                      grepl(paste0("bulkRNAseq_", dd, "___scRNAseq"), files.pca))]
  df.scdec.pca.coord  <- readRDS(df.scdec.pca.coord)
  df.scdec.pca.expVar <- files.pca[(grepl("__all_genes__", files.pca) & grepl("__inferred_cell_size", files.pca) & 
                                      grepl("__music2__", files.pca) & grepl("\\/avg_pca_explained_var__", files.pca) &
                                      grepl(paste0("bulkRNAseq_", dd, "___scRNAseq"), files.pca))]
  df.scdec.pca.expVar <- readRDS(df.scdec.pca.expVar)
  pcs.select          <- names(df.scdec.pca.expVar)[1:which(cumsum(df.scdec.pca.expVar) > 90)[1]]
  df.scdec.pca.coord  <- df.scdec.pca.coord[, pcs.select]
  df.scdec.pca.coord  <- as.data.frame(df.scdec.pca.coord)
  colnames(df.scdec.pca.coord) <- paste0("cf__", colnames(df.scdec.pca.coord))
  summary(df.scdec.pca.coord)
  round(apply(df.scdec.pca.coord, 2, mean), 4)
  round(apply(df.scdec.pca.coord, 2, sd), 4)
  df.scdec.pca.coord
  
})
names(df.scdec.pca.coord.list) <- c("CD", "HFD", "CD_HFD")



dir_source <- paste0(dirQTL_input, "/individuals_coonProteomic")
dir_dest   <- paste0(dirQTL_input, "/individuals_coonProteomic_cf_pca")

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
  covar.df  <- cbind("id" = rownames(covar.df), covar.df)
  rowsOrder <- rownames(covar.df)
  colsOrder <- colnames(covar.df)
  
  covar.df           <- merge(covar.df, df.scdec.pca.coord.list[[diet.x]], by.x = "strain_diet", by.y = "row.names", all.x = T, all.y = F)
  covar.df           <- dplyr::select(covar.df, c(colsOrder, colnames(covar.df)[!(colnames(covar.df) %in% colsOrder)]))
  covar.df           <- covar.df[match(rowsOrder, covar.df$id), ]
  rownames(covar.df) <- covar.df$id
  stopifnot(all(rowsOrder == covar.df$id))
  
  print(dim(df.scdec.pca.coord.list[[diet.x]]))
  print(dim(covar.df))
  print(head(covar.df, 2))
  print("")
  
  qtl2convert::write2csv(covar.df, paste0(x, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
  
})


