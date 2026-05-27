library(data.table)


pheno_ind_files <- c("./Data/input_data/phenotypic_data/pheno_formatted_withMeta.RDS",
                     "./Data/phenome/phenoData_transformation/pheno_derived_traits_1.RDS")
pheno.ind.df.list <- lapply(pheno_ind_files, function(y){
  # y <- pheno_ind_files[1]
  # print(y)
  df.y            <- readRDS(y)
  rownames(df.y)  <- paste0(df.y$strain, "_", df.y$diet, "_", df.y$strainReplicate, "_", df.y$dietReplicate)
  df.y            <- df.y[!df.y$excludeFromAnalysis, ]
  df.y            <- df.y[, -c(1:19), drop = F]
  df.y            <- cbind("sample_id" = rownames(df.y), df.y)
  df.y
})
pheno.ind.df <- Reduce(function(...) merge(..., all = T, by = c("sample_id")), pheno.ind.df.list)
rownames(pheno.ind.df) <- pheno.ind.df$sample_id
pheno.ind.df <- dplyr::select(pheno.ind.df, -dplyr::all_of("sample_id"))
pheno.ind.df[is.na(pheno.ind.df)] <- NA


pheno_ind_fasted_files     <- c("./Data/input_data/BXD_fasted_MPD/Auwerx1_submit.csv",
                                "./Data/phenome/fasted_phenoData_transformation/fasted_pheno_derived_traits_1.RDS")
pheno.ind.fasted.df.list   <- lapply(pheno_ind_fasted_files, function(y){
  # y <- pheno_ind_fasted_files[2]
  # print(y)
  if(grepl("csv$", y)){
    df.y <- read.table(y, sep = ",", header = T, stringsAsFactors = F)
  } else{
    df.y <- readRDS(y)
  }
  # keep only commoin strains
  df.y            <- df.y[df.y$strain %in% unique(gsub("_.*", "", rownames(pheno.ind.df))), ]
  stopifnot(all(unique(gsub("_.*", "", rownames(pheno.ind.df))) %in% unique(df.y$strain)))
  rownames(df.y)  <- paste0(df.y$strain, "_", df.y$diet, "_", df.y$id)
  cols_rm         <- c("id", "strain", "sex", "diet")
  df.y            <- dplyr::select(df.y, -dplyr::all_of(cols_rm[cols_rm %in% colnames(df.y)]))
  df.y            <- cbind("sample_id" = rownames(df.y), df.y)
  df.y
})
pheno.ind.fasted.df <- Reduce(function(...) merge(..., all = T, by = c("sample_id")), pheno.ind.fasted.df.list)
rownames(pheno.ind.fasted.df) <- pheno.ind.fasted.df$sample_id
pheno.ind.fasted.df <- dplyr::select(pheno.ind.fasted.df, -dplyr::all_of("sample_id"))
pheno.ind.fasted.df[is.na(pheno.ind.fasted.df)] <- NA
colnames(pheno.ind.fasted.df) <- paste0("fasted__", colnames(pheno.ind.fasted.df))




lipid.derived.list         <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")
lipid.df                   <- lipid.derived.list$Filt_MassNorm_BatchNorm_Data_NArm
rownames(lipid.df)         <- lipid.df$Identifier
lipid.df                   <- lipid.df[!lipid.df$isUnknown, ]
lipid.df                   <- lipid.df[, -c(1:11)]
colnames(lipid.df)         <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", colnames(lipid.df)))
lipid.df                   <- as.data.frame(t(lipid.df))
lipid.df                   <- log2(lipid.df + 1)
lipid.df                   <- lipid.df[, colSums(is.na(lipid.df)) < 0.5 * nrow(lipid.df)]




features_list.all <- list("pheno"        = pheno.ind.df,
                          "pheno_fasted" = pheno.ind.fasted.df,
                          "protein"      = prot.df,
                          "lipid"        = lipid.df)



kinmat_file  <- list.files("./Data/QTL_mapping/qtl_outputs/average_lipidomic", pattern = "kinship_overall", full.names = T)
kinmat_file  <- kinmat_file[grepl("CD", kinmat_file)]
kinmat       <- readRDS(readRDS)
diag(kinmat) <- 1 # make sure diag is 1 (numerical approx can make it slightly deviated from 1)


saveDir <- "./Data/heritability/lme4qtl_var_decomp"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}


nbCores <- 50
tt      <- lapply(names(features_list.all), function(x){
  # x <- names(features_list.all)[2]
  
  cat("Decomposing variability for all features of input type:", x, "...\n")
  
  
  df.x        <- features_list.all[[x]]
  var.df.list <- pbmcapply::pbmclapply(mc.cores = nbCores, X = colnames(df.x), FUN = function(y){
    # y <- colnames(df.x)[1]
    # y <- "HexCer[NS] d18:0_22:1_RT_12.391"
    
    df.y              <- dplyr::select(df.x, dplyr::all_of(y))
    df.y$strain       <- factor(gsub("_.*", "", rownames(df.y)))
    df.y$diet         <- factor(ifelse(grepl("CD", rownames(df.y)), "CD", "HFD"))
    df.y$strain_diet  <- factor(paste0(df.y$strain, "_", df.y$diet))
    colnames(df.y)[1] <- "feature"
    df.y              <- df.y[!is.na(df.y$feature), ]
    head(df.y)
    
    kin.mat           <- kinmat_list$CD
    rownames(kin.mat) <- colnames(kin.mat) <- gsub("_.*", "", rownames(kin.mat))
    
    kin.mat_CD_HFD           <- as.matrix(Matrix::kronecker(Matrix::Diagonal(2), kin.mat))
    rownames(kin.mat_CD_HFD) <- colnames(kin.mat_CD_HFD) <- c(paste0(rownames(kin.mat), "_CD"), paste0(rownames(kin.mat), "_HFD"))
    kin.mat_CD_HFD           <- kin.mat_CD_HFD[order(rownames(kin.mat_CD_HFD)), order(rownames(kin.mat_CD_HFD))]
    
    # kin.mat_CD_HFD[1:5, 1:5]
    
    kin.mat        <- kin.mat[rownames(kin.mat) %in% df.y$strain, colnames(kin.mat) %in% df.y$strain]
    kin.mat_CD_HFD <- kin.mat_CD_HFD[rownames(kin.mat_CD_HFD) %in% df.y$strain_diet, colnames(kin.mat_CD_HFD) %in% df.y$strain_diet]
    
    fit.list          <- list("diet_CD"     = lme4qtl::relmatLmer(feature ~ 1 + (1|strain), df.y[df.y$diet == "CD" & df.y$strain %in% rownames(kin.mat), ], relmat = list(strain = kin.mat)),
                              "diet_HFD"    = lme4qtl::relmatLmer(feature ~ 1 + (1|strain), df.y[df.y$diet == "HFD" & df.y$strain %in% rownames(kin.mat), ], relmat = list(strain = kin.mat)),
                              "diet_CD_HFD" = lme4qtl::relmatLmer(feature ~ 1 + diet + (1|strain) + (1|strain_diet), df.y[df.y$strain %in% rownames(kin.mat), ], relmat = list(strain = kin.mat, strain_diet = kin.mat_CD_HFD)))
    
    
    
    
    # exctact perc. of explained variance
    var.df.y <- lapply(names(fit.list), function(kk){
      # kk <- names(fit.list)[1]
      
      fit            <- fit.list[[kk]]
      vc             <- as.data.frame(lme4::VarCorr(fit))
      vc$effect_type <- "rndeff"
      
      # compute fixed term variance (same as done by insight::get_variance, insight::get_variance_fixed and same as suggested by chatGPT)
      beta  <- lme4::fixef(fit)
      if(!all(names(beta) == "(Intercept)")){
        X            <- model.matrix(terms(fit), data = df.y)
        # https://people.math.ethz.ch/~maechler/MEMo-pages/lMMwR.pdf?utm_source=chatgpt.com
        fitted_fixed <- as.vector(X %*% beta)
        var_fixed    <- var(fitted_fixed)
        fixed_row    <- data.frame(grp = "diet", var1 = NA, var2 = NA, vcov = var_fixed, sdcor = sd(fitted_fixed), effect_type = "fixeff", stringsAsFactors = FALSE)
        vc           <- rbind(vc, fixed_row)
      }
      vc$H2        <- vc$vcov / sum(vc$vcov) 
      vc$grp       <- gsub("strain_diet", "strain:diet", vc$grp)
      rownames(vc) <- paste0(vc$grp, "__", vc$effect_type)
      
      
      out           <- cbind("model_type" = kk, as.data.frame(t(vc[, "H2", drop = F])))
      colnames(out) <- tolower(colnames(out))
      colnames(out)[colnames(out) == "residual__rndeff"] <- "residual"
      out
      
    })
    var.df.y              <- as.data.frame(rbindlist(var.df.y, use.names = T, fill = T))
    var.df.y$feature_id   <- y
    var.df.y$feature_type <- x
    
    cols_order <- c("feature_type", "feature_id", "model_type", "diet__fixeff", "strain__rndeff", "strain:diet__rndeff", "residual")
    cols_order <- cols_order[cols_order %in% colnames(var.df.y)] 
    stopifnot(all(colnames(var.df.y) %in% cols_order))
    
    var.df.y     <- dplyr::select(var.df.y, dplyr::all_of(cols_order))
    var.df.y
  })
  var.df <- as.data.frame(rbindlist(var.df.list))
  
  saveRDS(var.df, paste0(saveDir, "/lme4qtl_var_decomp___", x, ".RDS"))
  NA
})







