
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
## Phenotypic data - averaged phenotypes
################################################################################

overwrite     <- F
pheno.mean.df <- readRDS("./Data/phenome/phenoData_transformation/pheno_formatted_mean_allValues.RDS")
pheno.mean.df <- as.data.frame(pheno.mean.df)

tt <- lapply(c("HFD", "CD", "all"), function(diet){
  # diet <- "HFD"
  
  saveDir <- paste0(dirQTL_input, "/average_pheno_mean/", ifelse(diet == "all", "CD_HFD", diet))
  if(!dir.exists(saveDir)){
    dir.create(saveDir, recursive = T)
  }
  if(file.exists(paste0(saveDir, "/BXD.json")) & !overwrite){
    cat("QTL files already created, skipped...\n")
    return()
  } else{
    cat("Creating QTL files...\n")
  }
  
  if(diet == "all"){
    pheno.mean.diet  <- reshape2::melt(pheno.mean.df, id.vars = c("strain", "diet"))
  } else{
    pheno.mean.diet  <- reshape2::melt(pheno.mean.df[pheno.mean.df$diet == diet, ], id.vars = c("strain", "diet"))
  }
  pheno.mean.diet$sample    <- paste0(pheno.mean.diet$strain, "_", pheno.mean.diet$diet)
  
  tmp       <- unique(dplyr::select(pheno.mean.diet, dplyr::all_of(c("strain", "diet"))))
  geno      <- readRDS(paste0(dirQTL_input, "/BXD_geno.RDS"))
  geno.diet <- geno[, "marker", drop = F]
  for(i in 1:nrow(tmp)){
    geno.diet <- cbind(geno.diet, geno[, tmp$strain[i], drop = F])
    colnames(geno.diet)[ncol(geno.diet)] <- paste0(tmp$strain[i], "_", tmp$diet[i])
  }
  crossinfo                 <- cbind(id = paste0(tmp$strain, "_", tmp$diet), tmp)
  colnames(crossinfo)       <- c("id", "strain", "diet")
  crossinfo$cross_direction <- "BxD"
  
  samplescovar              <- as.data.frame(dplyr::select(crossinfo, dplyr::all_of(c("id", "strain", "diet"))))
  samplescovar$sex          <- "m"
  rownames(samplescovar)    <- samplescovar$id
  crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -dplyr::all_of(c("strain", "diet"))))
  rownames(crossinfo)       <- crossinfo$id
  
  phenocovar                       <- type.convert(as.data.frame(unique(dplyr::select(pheno.mean.diet, dplyr::all_of(c("variable"))))), as.is = T)
  colnames(phenocovar)             <- c("phenoID")
  rownames(phenocovar)             <- phenocovar$phenoID
  pheno.mean.diet           <- reshape2::dcast(pheno.mean.diet, sample~variable, value.var = "value")
  
  tmp.df <- pheno.mean.diet
  rownames(tmp.df) <- tmp.df$sample
  tmp.df <- tmp.df[, -1, drop = F]
  saveRDS(tmp.df, paste0(saveDir, "/pheno_original_distribution.RDS"))
  
  # transform to normal distribution
  norm <- transform2normal_byCol_INT(pheno.mean.diet[, -1, drop = F])
  
  transf.type.df                   <- norm$transfType.df
  # remove phenotypes with less than 15 numeric measurements
  transf.df                        <- norm$transf.df
  
  tmpLogi                          <- apply(transf.df, 2, function(kk) (sum(!is.na(kk)) + sum(!is.infinite(kk) & !is.na(kk))) > 15)
  transf.df                        <- transf.df[, tmpLogi, drop = F]
  pheno.mean.diet           <- cbind(pheno.mean.diet[, 1, drop = F], transf.df)
  rownames(pheno.mean.diet) <- pheno.mean.diet$sample
  phenocovar                       <- phenocovar[phenocovar$phenoID %in% colnames(pheno.mean.diet), , drop = F]
  pheno.mean.diet           <- pheno.mean.diet[, rownames(phenocovar), drop = F]
  pheno.mean.diet[is.na(pheno.mean.diet)] <- NA
  
  saveRDS(norm, paste0(saveDir, "/pheno_normal_distribution.RDS"))
  saveRDS(transf.type.df, paste0(saveDir, "/pheno_normal_transformation_type.RDS"))
  qtl2convert::write2csv(geno.diet, paste0(saveDir, "/BXD_geno.csv"), overwrite = TRUE, comment = "Genotypes for BXD data")
  qtl2convert::write2csv(samplescovar, paste0(saveDir, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
  qtl2convert::write2csv(crossinfo, paste0(saveDir, "/BXD_crossinfo.csv"), overwrite = TRUE, comment = paste0("Cross info for BXD data\n#", "(all lines formed from cross between female B and male D)"))
  qtl2convert::write2csv(phenocovar, paste0(saveDir, "/BXD_phenocovar.csv"), overwrite = TRUE, comment = "Phenotype covariates (metadata) for BXD phenotype data")
  qtl2convert::write2csv(cbind(id = rownames(pheno.mean.diet), pheno.mean.diet), paste0(saveDir, "/BXD_pheno.csv"), overwrite = TRUE, comment = "BXD phenotype data")
  
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
  NA
})



################################################################################
## Phenotypic data - individuals phenotypes
################################################################################


overwrite   <- F
pheno.indiv <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_withMeta.RDS")
pheno.indiv <- dplyr::select(pheno.indiv, dplyr::all_of(c("strain", "diet", "dietReplicate", "strainReplicate", colnames(pheno.indiv)[6:ncol(pheno.indiv)])))

lapply(c("HFD", "CD", "all"), function(diet){
  # diet <- "all"
  
  saveDir <- paste0(dirQTL_input, "/individuals_pheno", "/", ifelse(diet == "all", "CD_HFD", diet))
  if(!dir.exists(saveDir)){
    dir.create(saveDir, recursive = T)
  }
  if(file.exists(paste0(saveDir, "/BXD.json")) & !overwrite){
    cat("QTL files already created, skipped...\n")
    return()
  } else{
    cat("Creating QTL files...\n")
  }
  
  if(diet == "all"){
    pheno.indiv.diet  <- reshape2::melt(pheno.indiv, id.vars = c("strain", "diet", "dietReplicate", "strainReplicate"))
  } else{
    pheno.indiv.diet  <- reshape2::melt(pheno.indiv[pheno.indiv$diet == diet, ], id.vars = c("strain", "diet", "dietReplicate", "strainReplicate"))
  }
  
  pheno.indiv.diet$sample    <- paste0(pheno.indiv.diet$strain, "_", pheno.indiv.diet$diet, "_", pheno.indiv.diet$strainReplicate, "_", pheno.indiv.diet$dietReplicate)
  
  tmp       <- unique(dplyr::select(pheno.indiv.diet, dplyr::all_of(c("strain", "diet", "dietReplicate", "strainReplicate"))))
  geno      <- readRDS(paste0(dirQTL_input, "/BXD_geno.RDS"))
  geno.diet <- geno[, "marker", drop = F]
  for(i in 1:nrow(tmp)){
    # print(i)
    geno.diet <- cbind(geno.diet, geno[, tmp$strain[i], drop = F])
    colnames(geno.diet)[ncol(geno.diet)] <- paste0(tmp$strain[i], "_", tmp$diet[i], "_", tmp$strainReplicate[i], "_", tmp$dietReplicate[i])
  }
  crossinfo                 <- cbind(id = paste0(tmp$strain, "_", tmp$diet, "_", tmp$strainReplicate, "_", tmp$dietReplicate), tmp)
  colnames(crossinfo)       <- c("id", "strain", "diet", "dietReplicate", "strainReplicate")
  crossinfo$cross_direction <- "BxD"
  crossinfo$sample_name     <- paste0(crossinfo$strain, "_", crossinfo$strainReplicate, "_", crossinfo$diet, crossinfo$dietReplicate)
  
  samplescovar              <- as.data.frame(dplyr::select(crossinfo, dplyr::all_of(c("id", "sample_name", "strain", "diet", "dietReplicate", "strainReplicate"))))
  samplescovar$sex          <- "m"
  rownames(samplescovar)    <- samplescovar$id
  crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -dplyr::all_of(c("strain", "sample_name", "diet", "dietReplicate", "strainReplicate"))))
  rownames(crossinfo)       <- crossinfo$id
  
  phenocovar                 <- type.convert(as.data.frame(unique(dplyr::select(pheno.indiv.diet, dplyr::all_of(c("variable"))))), as.is = T)
  colnames(phenocovar)       <- c("phenoID")
  rownames(phenocovar)       <- phenocovar$phenoID
  pheno.indiv.diet           <- type.convert(as.data.frame(dcast(data.table(pheno.indiv.diet), sample~variable, value.var = "value")), as.is = T)
  
  df.tmp           <- pheno.indiv.diet
  rownames(df.tmp) <- df.tmp$sample
  df.tmp           <- df.tmp[, -1, drop = F]
  saveRDS(df.tmp, paste0(saveDir, "/pheno_original_distribution.RDS"))
  
  # transform to normal distribution
  norm <- transform2normal_byCol_INT(pheno.indiv.diet[, -1, drop = F])
  
  transf.type.df                   <- norm$transfType.df
  # remove phenotypes with less than 15 numeric measurements
  transf.df                        <- norm$transf.df
  tmpLogi                          <- apply(transf.df, 2, function(kk) (sum(!is.na(kk)) + sum(!is.infinite(kk) & !is.na(kk))) > 15)
  transf.df                        <- transf.df[, tmpLogi, drop = F]
  pheno.indiv.diet                 <- cbind(pheno.indiv.diet[, 1, drop = F], transf.df)
  rownames(pheno.indiv.diet)       <- pheno.indiv.diet$sample
  phenocovar                       <- phenocovar[phenocovar$phenoID %in% colnames(pheno.indiv.diet), , drop = F]
  pheno.indiv.diet                 <- pheno.indiv.diet[, rownames(phenocovar), drop = F]
  pheno.indiv.diet[is.na(pheno.indiv.diet)] <- NA
  
  saveRDS(norm, paste0(saveDir, "/pheno_normal_distribution.RDS"))
  saveRDS(transf.type.df, paste0(saveDir, "/pheno_normal_transformation_type.RDS"))
  qtl2convert::write2csv(geno.diet, paste0(saveDir, "/BXD_geno.csv"), overwrite = TRUE, comment = "Genotypes for BXD data")
  qtl2convert::write2csv(samplescovar, paste0(saveDir, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
  qtl2convert::write2csv(crossinfo, paste0(saveDir, "/BXD_crossinfo.csv"), overwrite = TRUE, comment = paste0("Cross info for BXD data\n#", "(all lines formed from cross between female B and male D)"))
  qtl2convert::write2csv(phenocovar, paste0(saveDir, "/BXD_phenocovar.csv"), overwrite = TRUE, comment = "Phenotype covariates (metadata) for BXD phenotype data")
  qtl2convert::write2csv(cbind(id = rownames(pheno.indiv.diet), pheno.indiv.diet), paste0(saveDir, "/BXD_pheno.csv"), overwrite = TRUE, comment = "BXD phenotype data")
  
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
})



################################################################################
## Fasted study phenotypic data - averaged phenotypes
################################################################################

# keep strains overlapping with fed strain
strains.keep <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_withMeta.RDS")
table(strains.keep$strain == strains.keep$strain)
strains.keep <- stringr::str_sort(unique(strains.keep$strain), numeric = T)

overwrite <- F
tt        <- lapply(c("allStrains", "commonStrains"), function(ss){
  # ss <- "commonStrains"
  
  pheno.mean.df <- readRDS("./Data/phenome/fasted_phenoData_transformation/fasted_pheno_formatted_mean_allValues.RDS")
  pheno.mean.df <- as.data.frame(pheno.mean.df)
  
  
  if(ss == "commonStrains"){
    strains.noOverlap   <- strains.keep[!(strains.keep %in% unique(pheno.mean.df$strain))]
    strains.noIntersect <- unique(pheno.mean.df$strain)[!(unique(pheno.mean.df$strain) %in% strains.keep)]
    if(length(strains.noOverlap) > 0){
      cat("The following strains to keep are not found in the input data (N =", length(strains.noOverlap), "):", paste(strains.noOverlap, collapse = "; "), "\n")
    }
    if(length(strains.noIntersect) > 0){
      cat("The following strains are removed from the input data (N =", length(strains.noIntersect), "):",  paste(strains.noIntersect, collapse = "; "), "\n")
    }
    pheno.mean.df <- pheno.mean.df[pheno.mean.df$strain %in% strains.keep, ]
  }
  
  
  tt <- lapply(c("HFD", "CD", "all"), function(diet){
    # diet <- "HFD"
    
    saveDir <- paste0(dirQTL_input, "/average_fasted_pheno_mean_", ss, "/", ifelse(diet == "all", "CD_HFD", diet))
    if(!dir.exists(saveDir)){
      dir.create(saveDir, recursive = T)
    }
    if(file.exists(paste0(saveDir, "/BXD.json")) & !overwrite){
      cat("QTL files already created, skipped...\n")
      return()
    } else{
      cat("Creating QTL files...\n")
    }
    
    if(diet == "all"){
      pheno.mean.diet  <- reshape2::melt(pheno.mean.df, id.vars = c("strain", "diet"))
    } else{
      pheno.mean.diet  <- reshape2::melt(pheno.mean.df[pheno.mean.df$diet == diet, ], id.vars = c("strain", "diet"))
    }
    # unique(pheno.mean.diet$diet)
    
    pheno.mean.diet$sample    <- paste0(pheno.mean.diet$strain, "_", pheno.mean.diet$diet)
    tmp                       <- unique(dplyr::select(pheno.mean.diet, dplyr::all_of(c("strain", "diet"))))
    
    if(ss == "commonStrains"){
      geno  <- readRDS(paste0(dirQTL_input, "/BXD_geno.RDS"))
    } else{
      geno        <- readRDS(paste0(dirQTL_input, "/BXD_geno_all_BXD_strains.RDS"))
      cols.remove <- c("Chr", "cM", "Mb", "Mb_mm9", "Mb_mm10")
      cols.remove <- cols.remove[cols.remove %in% colnames(geno)]
      geno        <- dplyr::select(geno, -dplyr::all_of(cols.remove))
      colnames(geno)[colnames(geno) == "Locus"] <- "marker"
    }
    stopifnot(all(tmp$strain %in% colnames(geno)))
    
    geno.diet <- geno[, "marker", drop = F]
    for(i in 1:nrow(tmp)){
      geno.diet <- cbind(geno.diet, geno[, tmp$strain[i], drop = F])
      colnames(geno.diet)[ncol(geno.diet)] <- paste0(tmp$strain[i], "_", tmp$diet[i])
    }
    crossinfo                 <- cbind(id = paste0(tmp$strain, "_", tmp$diet), tmp)
    colnames(crossinfo)       <- c("id", "strain", "diet")
    crossinfo$cross_direction <- "BxD"
    
    samplescovar              <- as.data.frame(dplyr::select(crossinfo, dplyr::all_of(c("id", "strain", "diet"))))
    samplescovar$sex          <- "m"
    rownames(samplescovar)    <- samplescovar$id
    crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -dplyr::all_of(c("strain", "diet"))))
    rownames(crossinfo)       <- crossinfo$id
    
    phenocovar                <- type.convert(as.data.frame(unique(dplyr::select(pheno.mean.diet, dplyr::all_of(c("variable"))))), as.is = T)
    colnames(phenocovar)      <- c("phenoID")
    rownames(phenocovar)      <- phenocovar$phenoID
    pheno.mean.diet           <- reshape2::dcast(pheno.mean.diet, sample~variable, value.var = "value")
    
    tmp.df           <- pheno.mean.diet
    rownames(tmp.df) <- tmp.df$sample
    tmp.df           <- tmp.df[, -1]
    saveRDS(tmp.df, paste0(saveDir, "/pheno_original_distribution.RDS"))
    
    # transform to normal distribution
    norm <- transform2normal_byCol_INT(pheno.mean.diet[, -1, drop = F])
    
    transf.type.df                   <- norm$transfType.df
    # remove phenotypes with less than 15 numeric measurements
    transf.df                        <- norm$transf.df
    
    tmpLogi                          <- apply(transf.df, 2, function(kk) (sum(!is.na(kk)) + sum(!is.infinite(kk) & !is.na(kk))) > 15)
    transf.df                        <- transf.df[, tmpLogi, drop = F]
    pheno.mean.diet           <- cbind(pheno.mean.diet[, 1, drop = F], transf.df)
    rownames(pheno.mean.diet) <- pheno.mean.diet$sample
    phenocovar                       <- phenocovar[phenocovar$phenoID %in% colnames(pheno.mean.diet), , drop = F]
    pheno.mean.diet           <- pheno.mean.diet[, rownames(phenocovar), drop = F]
    pheno.mean.diet[is.na(pheno.mean.diet)] <- NA
    
    saveRDS(norm, paste0(saveDir, "/pheno_normal_distribution.RDS"))
    saveRDS(transf.type.df, paste0(saveDir, "/pheno_normal_transformation_type.RDS"))
    qtl2convert::write2csv(geno.diet, paste0(saveDir, "/BXD_geno.csv"), overwrite = TRUE, comment = "Genotypes for BXD data")
    qtl2convert::write2csv(samplescovar, paste0(saveDir, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
    qtl2convert::write2csv(crossinfo, paste0(saveDir, "/BXD_crossinfo.csv"), overwrite = TRUE, comment = paste0("Cross info for BXD data\n#", "(all lines formed from cross between female B and male D)"))
    qtl2convert::write2csv(phenocovar, paste0(saveDir, "/BXD_phenocovar.csv"), overwrite = TRUE, comment = "Phenotype covariates (metadata) for BXD phenotype data")
    qtl2convert::write2csv(cbind(id = rownames(pheno.mean.diet), pheno.mean.diet), paste0(saveDir, "/BXD_pheno.csv"), overwrite = TRUE, comment = "BXD phenotype data")
    
    
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
    NA
  })
})



################################################################################
## Fasted study phenotypic data - individuals phenotypes
################################################################################



# keep strains overlapping with fed strain
strains.keep <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_withMeta.RDS")
strains.keep <- strains.keep[!strains.keep$excludeFromAnalysis, ]
table(strains.keep$strain == strains.keep$strain)
strains.keep <- stringr::str_sort(unique(strains.keep$strain), numeric = T)

overwrite <- F
tt        <- lapply(c("allStrains", "commonStrains"), function(ss){
  # ss <- "allStrains"
  
  pheno.indiv <- read.table("./Data/input_data/BXD_fasted_MPD/Auwerx1_submit.csv", sep = ",", header = T, stringsAsFactors = F)
  
  
  
  if(ss == "commonStrains"){
    strains.noOverlap   <- strains.keep[!(strains.keep %in% unique(pheno.indiv$strain))]
    strains.noIntersect <- unique(pheno.indiv$strain)[!(unique(pheno.indiv$strain) %in% strains.keep)]
    if(length(strains.noOverlap) > 0){
      cat("The following strains to keep are not found in the input data (N =", length(strains.noOverlap), "):", paste(strains.noOverlap, collapse = "; "), "\n")
    }
    if(length(strains.noIntersect) > 0){
      cat("The following strains are removed from the input data (N =", length(strains.noIntersect), "):",  paste(strains.noIntersect, collapse = "; "), "\n")
    }
    pheno.indiv <- pheno.indiv[pheno.indiv$strain %in% strains.keep, ]
  }
  
  lapply(c("HFD", "CD", "all"), function(diet){
    # diet <- "all"
    
    saveDir <- paste0(dirQTL_input, "/individuals_fasted_pheno_", ss, "/", ifelse(diet == "all", "CD_HFD", diet))
    if(!dir.exists(saveDir)){
      dir.create(saveDir, recursive = T)
    }
    if(file.exists(paste0(saveDir, "/BXD.json")) & !overwrite){
      cat("QTL files already created, skipped...\n")
      return()
    } else{
      cat("Creating QTL files...\n")
    }
    
    if(diet == "all"){
      pheno.indiv.diet  <- reshape2::melt(pheno.indiv, id.vars = c("strain", "diet", "id"))
    } else{
      pheno.indiv.diet  <- reshape2::melt(pheno.indiv[pheno.indiv$diet == diet, ], id.vars = c("strain", "diet", "id"))
    }
    
    pheno.indiv.diet$sample    <- paste0(pheno.indiv.diet$strain, "_", pheno.indiv.diet$diet, "_", gsub("[^[:alnum:]]", "", pheno.indiv.diet$id))
    
    if(diet == "all"){
      stopifnot(length(unique(pheno.indiv.diet$sample)) == nrow(pheno.indiv))
    } else{
      stopifnot(length(unique(pheno.indiv.diet$sample)) == nrow(pheno.indiv[pheno.indiv$diet == diet, ]))
    }
    
    
    tmp <- unique(dplyr::select(pheno.indiv.diet, dplyr::all_of(c("strain", "diet", "id"))))
    
    if(ss == "commonStrains"){
      geno  <- readRDS(paste0(dirQTL_input, "/BXD_geno.RDS"))
    } else{
      # t1 <- readRDS(paste0(dirQTL_input, "/BXD_geno_all_BXD_strains.RDS"))
      # t2 <- readRDS(paste0(dirQTL_input, "/BXD_geno_all_BXD_strains_new_geno.RDS"))
      geno        <- readRDS(paste0(dirQTL_input, "/BXD_geno_all_BXD_strains.RDS"))
      cols.remove <- c("Chr", "cM", "Mb", "Mb_mm9", "Mb_mm10")
      cols.remove <- cols.remove[cols.remove %in% colnames(geno)]
      geno        <- dplyr::select(geno, -dplyr::all_of(cols.remove))
      colnames(geno)[colnames(geno) == "Locus"] <- "marker"
    }
    stopifnot(all(tmp$strain %in% colnames(geno)))
    
    geno.diet <- geno[, "marker", drop = F]
    for(i in 1:nrow(tmp)){
      # print(i)
      geno.diet <- cbind(geno.diet, geno[, tmp$strain[i], drop = F])
      colnames(geno.diet)[ncol(geno.diet)] <- paste0(tmp$strain[i], "_", tmp$diet[i], "_", gsub("[^[:alnum:]]", "", tmp$id[i]))
    }
    crossinfo                 <- cbind(id = paste0(tmp$strain, "_", tmp$diet, "_", gsub("[^[:alnum:]]", "", tmp$id)), dplyr::select(tmp, -dplyr::all_of(c("id"))))
    colnames(crossinfo)       <- c("id", "strain", "diet")
    crossinfo$cross_direction <- "BxD"
    crossinfo$sample_name     <- crossinfo$id
    
    samplescovar              <- as.data.frame(dplyr::select(crossinfo, dplyr::all_of(c("id", "sample_name", "strain", "diet"))))
    samplescovar$sex          <- "m"
    rownames(samplescovar)    <- samplescovar$id
    crossinfo                 <- as.data.frame(dplyr::select(crossinfo, -dplyr::all_of(c("strain", "sample_name", "diet"))))
    rownames(crossinfo)       <- crossinfo$id
    
    phenocovar                 <- type.convert(as.data.frame(unique(dplyr::select(pheno.indiv.diet, dplyr::all_of(c("variable"))))), as.is = T)
    colnames(phenocovar)       <- c("phenoID")
    rownames(phenocovar)       <- phenocovar$phenoID
    pheno.indiv.diet           <- type.convert(as.data.frame(dcast(data.table(pheno.indiv.diet), sample~variable, value.var = "value")), as.is = T)
    
    df.tmp           <- pheno.indiv.diet
    rownames(df.tmp) <- df.tmp$sample
    df.tmp           <- df.tmp[, -1, drop = F]
    saveRDS(df.tmp, paste0(saveDir, "/pheno_original_distribution.RDS"))
    
    # transform to normal distribution
    norm <- transform2normal_byCol_INT(pheno.indiv.diet[, -1, drop = F])
    
    transf.type.df                   <- norm$transfType.df
    # remove phenotypes with less than 15 numeric measurements
    transf.df                        <- norm$transf.df
    tmpLogi                          <- apply(transf.df, 2, function(kk) (sum(!is.na(kk)) + sum(!is.infinite(kk) & !is.na(kk))) > 15)
    transf.df                        <- transf.df[, tmpLogi, drop = F]
    pheno.indiv.diet                 <- cbind(pheno.indiv.diet[, 1, drop = F], transf.df)
    rownames(pheno.indiv.diet)       <- pheno.indiv.diet$sample
    phenocovar                       <- phenocovar[phenocovar$phenoID %in% colnames(pheno.indiv.diet), , drop = F]
    pheno.indiv.diet                 <- pheno.indiv.diet[, rownames(phenocovar), drop = F]
    pheno.indiv.diet                 <- pheno.indiv.diet[colnames(geno.diet)[-1], ]
    pheno.indiv.diet[is.na(pheno.indiv.diet)] <- NA
    
    stopifnot(identical(colnames(geno.diet)[-1], rownames(samplescovar)))
    stopifnot(identical(colnames(geno.diet)[-1], rownames(pheno.indiv.diet)))
    stopifnot(identical(colnames(geno.diet)[-1], rownames(crossinfo)))
    stopifnot(identical(colnames(pheno.indiv.diet), rownames(phenocovar)))
    
    saveRDS(norm, paste0(saveDir, "/pheno_normal_distribution.RDS"))
    saveRDS(transf.type.df, paste0(saveDir, "/pheno_normal_transformation_type.RDS"))
    qtl2convert::write2csv(geno.diet, paste0(saveDir, "/BXD_geno.csv"), overwrite = TRUE, comment = "Genotypes for BXD data")
    qtl2convert::write2csv(samplescovar, paste0(saveDir, "/BXD_covar.csv"), overwrite = TRUE, comment = paste0("Sample covariates"))
    qtl2convert::write2csv(crossinfo, paste0(saveDir, "/BXD_crossinfo.csv"), overwrite = TRUE, comment = paste0("Cross info for BXD data\n#", "(all lines formed from cross between female B and male D)"))
    qtl2convert::write2csv(phenocovar, paste0(saveDir, "/BXD_phenocovar.csv"), overwrite = TRUE, comment = "Phenotype covariates (metadata) for BXD phenotype data")
    qtl2convert::write2csv(cbind(id = rownames(pheno.indiv.diet), pheno.indiv.diet), paste0(saveDir, "/BXD_pheno.csv"), overwrite = TRUE, comment = "BXD phenotype data")
    
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
    NA
  })
})





