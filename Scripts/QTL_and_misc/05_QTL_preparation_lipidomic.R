

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
## Newly normalized lipidomic data and derived lipidomic features
################################################################################

lipid.derived.list  <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")

# get metadata table
meta            <- lipid.derived.list$MetaData
meta$SampleName[grepl("DBA|C57", meta$SampleName)] <- gsub("^BXD", "", meta$SampleName[grepl("DBA|C57", meta$SampleName)])
rownames(meta)  <- meta$SampleName

# retain tables to perform QTL mapping on
lipid.derived.list <- lipid.derived.list[grepl("derived_|Filt_MassNorm_BatchNorm_Data_NArm", names(lipid.derived.list))]
rownames(lipid.derived.list$Filt_MassNorm_BatchNorm_Data_NArm) <- lipid.derived.list$Filt_MassNorm_BatchNorm_Data_NArm$Identifier

# retain only sample columns
lipid.derived.list <- lapply(lipid.derived.list, function(x){
  cols_select <- colnames(x)[grepl("BXD|DBA|C57", colnames(x))]
  dplyr::select(x, dplyr::all_of(cols_select))
})

# split between individual and avg data
lipid.derived.list <- split(lipid.derived.list, ifelse(grepl("derived_avgAggr", names(lipid.derived.list)), "average", "individuals"))

# reorder columns to match them
# rename DBA and C57 columns
lipid.derived.list <- lapply(lipid.derived.list, function(x){
  # x <- lipid.derived.list[[1]]
  
  # check if always the same nb. of columns is available
  stopifnot(length(unique(unlist(lapply(x, ncol)))) == 1)
  
  # reorder
  cols_order <- colnames(x[[1]])
  x          <- lapply(x, function(y) dplyr::select(y, dplyr::all_of(cols_order)))
  stopifnot(all(sapply(lapply(x, colnames), function(z) identical(z, colnames(x[[1]])))))
  
  # rename DBA and C57 columns
  x <- lapply(x, function(y){
    colnames(y) <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", colnames(y)))
    y
  })
  
  x
})


# retain only metadata entries for which data are available
meta$SampleName <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", meta$SampleName))
meta$strainDiet <- paste0(meta$Strain, "_", meta$Diet)
meta            <- meta[meta$strainDiet %in% colnames(lipid.derived.list$average[[1]]) & meta$SampleName %in% colnames(lipid.derived.list$individuals[[1]]), ]

stopifnot(all(colnames(lipid.derived.list$average[[1]]) %in% meta$strainDiet))
stopifnot(all(colnames(lipid.derived.list$individuals[[1]]) %in% meta$SampleName))

overwrite  <- F
na_strings <- c("-", "NA")
tt         <- lapply(c("HFD", "CD", "all"), function(diet){
  # diet <- "HFD"
  
  tt <- lapply(c("average", "individuals"), function(zz){
    # zz <- "average"
    
    lipid.df.list <- lipid.derived.list[[zz]]
    tt            <- lapply(names(lipid.df.list), function(yy){
      # yy <- names(lipid.df.list)[1]
      
      cat(diet, " - ", zz, " - ", yy, "\n")
      
      saveDir <- paste0(dirQTL_input, "/", zz, "_lipidomic_", yy, "/", ifelse(diet == "all", "CD_HFD", diet))
      if(!dir.exists(saveDir)){
        dir.create(saveDir, recursive = T)
      }
      if(file.exists(paste0(saveDir, "/BXD.json")) & !overwrite){
        cat("QTL files already created, skipped...\n")
        return()
      } else{
        cat("Creating QTL files...\n")
      }
      
      df.expr <- as.matrix(t(lipid.df.list[[yy]]))
      if(diet != "all"){
        df.expr <- df.expr[grepl(diet, rownames(df.expr)), ]
      }
      
      if(zz == "average"){

        tmp.meta           <- unique(dplyr::select(meta, c("Strain", "Diet", "strainDiet")))
        colnames(tmp.meta) <- c("strain", "diet", "strainDiet")
        rownames(tmp.meta) <- tmp.meta$strainDiet
        tmp.meta           <- tmp.meta[tmp.meta$strainDiet %in% rownames(df.expr), ]
        stopifnot(nrow(tmp.meta) == nrow(df.expr))
        geno               <- readRDS(paste0(dirQTL_input, "/BXD_geno.RDS"))
        geno.diet          <- geno[, "marker", drop = F]
        for(i in 1:nrow(tmp.meta)){
          # print(i)
          # i <- 6
          geno.diet <- cbind(geno.diet, geno[, tmp.meta$strain[i], drop = F])
          colnames(geno.diet)[ncol(geno.diet)] <- paste0(tmp.meta$strainDiet[i])
        }
        crossinfo                 <- tmp.meta[tmp.meta$strainDiet %in% rownames(df.expr), ]
        crossinfo                 <- dplyr::select(crossinfo, c("strainDiet", "strain", "diet"))
        colnames(crossinfo)       <- c("id", "strain", "diet")
        crossinfo$cross_direction <- "BxD"
        
        samplescovar              <- as.data.frame(dplyr::select(crossinfo, c("id", "strain", "diet")))
        samplescovar$sex          <- "m"
        rownames(samplescovar)    <- samplescovar$id
        crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -c("strain", "diet")))
        rownames(crossinfo)       <- crossinfo$id
        
      } else if(zz == "individuals"){
        
        tmp.meta           <- unique(meta[meta$SampleName %in% rownames(df.expr), ])
        rownames(tmp.meta) <- tmp.meta$SampleName
        stopifnot(nrow(tmp.meta) == nrow(df.expr))
        geno               <- readRDS(paste0(dirQTL_input, "/BXD_geno.RDS"))
        geno.diet          <- geno[, "marker", drop = F]
        for(i in 1:nrow(tmp.meta)){
          # print(i)
          # i <- 1
          geno.diet <- cbind(geno.diet, geno[, tmp.meta$Strain[i], drop = F])
          colnames(geno.diet)[ncol(geno.diet)] <- paste0(tmp.meta$SampleName[i])
        }
        crossinfo                 <- meta[meta$SampleName %in% rownames(df.expr), ]
        crossinfo                 <- dplyr::select(crossinfo, c("SampleName", "Strain", "Batch", "Diet"))
        colnames(crossinfo)       <- c("id", "strain", "batch", "diet")
        crossinfo$cross_direction <- "BxD"
        
        samplescovar              <- as.data.frame(dplyr::select(crossinfo, c("id", "strain", "diet", "batch")))
        samplescovar$sex          <- "m"
        rownames(samplescovar)    <- samplescovar$id
        crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -c("strain", "diet", "batch")))
        rownames(crossinfo)       <- crossinfo$id
      }
      
      
      saveRDS(df.expr, paste0(saveDir, "/pheno_original_distribution.RDS"))
      
      # transform to normal distribution
      norm                             <- transform2normal_byCol_INT(df.expr)
      transf.type.df                   <- norm$transfType.df
      # remove phenotypes with less than 15 strains with phenotypic measurements
      transf.df                        <- norm$transf.df
      strain_vec                       <- gsub("_.*", "", rownames(df.expr))
      tmpLogi                          <- apply(transf.df, 2, function(kk) length(unique(strain_vec[!is.na(kk) & !is.infinite(kk) & !is.na(kk)])) > 15)
      transf.df                        <- transf.df[, tmpLogi]
      rownames(transf.df)              <- rownames(df.expr)
      df.expr                          <- transf.df
      
      phenocovar                  <- data.frame(phenoID = colnames(df.expr), stringsAsFactors = F)
      rownames(phenocovar)        <- phenocovar$phenoID
      df.expr                     <- df.expr[colnames(geno.diet)[-1], rownames(phenocovar)]
      df.expr[is.na(df.expr)]     <- NA
      
      
      stopifnot(identical(colnames(geno.diet)[-1], rownames(df.expr)))
      stopifnot(identical(colnames(geno.diet)[-1], rownames(samplescovar)))
      stopifnot(identical(colnames(geno.diet)[-1], rownames(crossinfo)))
      stopifnot(identical(rownames(phenocovar), colnames(df.expr)))
      
      saveRDS(norm, paste0(saveDir, "/pheno_normal_distribution.RDS"))
      saveRDS(transf.type.df, paste0(saveDir, "/pheno_normal_transformation_type.RDS"))
      qtl2convert::write2csv(geno.diet, paste0(saveDir, "/BXD_geno.csv"), overwrite = TRUE, comment = "Genotypes for BXD data")
      qtl2convert::write2csv(samplescovar, paste0(saveDir, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
      qtl2convert::write2csv(crossinfo, paste0(saveDir, "/BXD_crossinfo.csv"), overwrite = TRUE, comment = paste0("Cross info for BXD data\n#", "(all lines formed from cross between female B and male D)"))
      qtl2convert::write2csv(phenocovar, paste0(saveDir, "/BXD_phenocovar.csv"), overwrite = TRUE, comment = "Phenotype covariates (metadata) for BXD phenotype data")
      qtl2convert::write2csv(cbind(id = rownames(df.expr), df.expr), paste0(saveDir, "/BXD_pheno.csv"), overwrite = TRUE, comment = "BXD phenotype data")
      
      file.copy(from = paste0(dirQTL_input, "/BXD_gmap.csv"), to = saveDir)
      file.copy(from = paste0(dirQTL_input, "/BXD_pmap.csv"), to = saveDir)
      
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
     
      NA
    })
    NA
  })
  NA
})


