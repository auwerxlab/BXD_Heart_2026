
library(data.table)


source("./Scripts/Lipidomics/lipidomics_fun.R")

################################################################################
## Compute mahalanobis distances to measure how far a sample is from a 
## multivariate distribution taking into account the covariance between variables.
################################################################################



pheno.all <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_list_all_fed_fasted.RDS")


# Remove FS as highly covaring/correlated, makes matrix not inversible (collinearity)
# echocardiography:        26 weeks old
# fasted__dist_treadmill1: 23 weeks old
# fasted__dist_treadmill2: 25 weeks old
# fasted__dist_wheel:      23 weeks old
pheno.keep   <- c("A_duration_[ms]", "BW_Week27_echo_[g]", "DTE_[ms]", "EF_[%]", "IVCT_[ms]", "IVRT_[ms]",
                  "IVS_diastole_[mm]", "IVS_systole_[mm]", "LV_mass_echo_corrected/Tibia_length_[mg/mm]",
                  "LVID_diastole_[mm]", "LVID_systole_[mm]", "LVPW_diastole_[mm]", "LVPW_systole_[mm]",
                  "MV_A_Vel_[mm/s]", "MV_E_Vel_[mm/s]", "MV_E/A_correct",
                  "fasted__dist_treadmill1", "fasted__dist_treadmill2", "fasted__dist_wheel")
tmp.df.pheno <- dplyr::select(pheno.all$avg_fed_and_fasted, dplyr::all_of(pheno.keep))

lipid.derived.list  <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")
meta.lipids.df      <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_metadata.RDS")
lipid.df            <- lipid.derived.list$derived_avgAggr_data
lipid.df            <- lipid.df[!lipid.df$isUnknown, ]
if(all(!(c("Lipid.Grandparent.Class", "Lipid.Parent.Class") %in% colnames(lipid.df)))){
  lipid.df  <- merge(lipid.df, meta.lipids.df, by = "Identifier", all = F)
}
meta_cols             <- colnames(lipid.df)[!grepl("bxd|dba|c57", colnames(lipid.df), ignore.case = T)]
rownames(lipid.df)    <- lipid.df$Identifier
lipid.df              <- type.convert(reshape2::melt(lipid.df, id.vars = meta_cols), as.is = T)
tot_pool              <- data.table(lipid.df)[, list(tot_pool  = sum(value, na.rm = T)), by = c("variable")]
tot_class             <- data.table(lipid.df)[, list(tot_class = sum(value, na.rm = T)), by = c("variable", "Lipid.Parent.Class")]
lipid_perc            <- merge(tot_class, tot_pool, by = "variable")
lipid_perc$perc_class <- (lipid_perc$tot_class / lipid_perc$tot_pool) * 100
# Merge classes representing < 1% of total lipid pool. Those are dangerous as a small difference here is in proportion huge (e.g. a change 0.1% > 0.2% --> 0.2/0.1 = 2, 30% > 40% --> 40/30 = 1.3).
# The differences for lowly abundant lipids are exaggerated, especially when transforming compositional data with CLR or ILR (involves logs)
# Minor lipids are also often highly variable (multiplicative variability: 0.01% → 0.03% = +200%) and I observed that they are strongly correlated (e.g. lyso-lipid classes) --> extremely tight cluster for major lipids and very wide clusters for rare lipids.
classes_merge         <- unique(lipid_perc$Lipid.Parent.Class[lipid_perc$perc_class < 1])
lipid_perc$class      <- ifelse(lipid_perc$Lipid.Parent.Class %in% classes_merge, "rest", lipid_perc$Lipid.Parent.Class)
lipid_perc.aggr       <- data.table(lipid_perc)[, list(tot_perc_class = sum(perc_class)), by = c("variable", "class")]
lipid_perc            <- as.data.frame(dcast(lipid_perc.aggr, variable~class, value.var = "tot_perc_class"))
rownames(lipid_perc)  <- lipid_perc$variable
lipid_perc            <- lipid_perc[, -1]
rownames(lipid_perc)  <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", rownames(lipid_perc)))
lipid_perc            <- lipid_perc[rownames(tmp.df.pheno), ]
colnames(lipid_perc)  <- paste0("parent_class__", colnames(lipid_perc), "__perc_tot_pool")


lipid_sum            <- merge(tot_class, tot_pool, by = "variable")
lipid_sum$class      <- ifelse(lipid_sum$Lipid.Parent.Class %in% classes_merge, "rest", lipid_sum$Lipid.Parent.Class)
lipid_sum.aggr       <- data.table(lipid_sum)[, list(tot_sum_class = sum(tot_class)), by = c("variable", "class")]
lipid_sum            <- as.data.frame(dcast(lipid_sum.aggr, variable~class, value.var = "tot_sum_class"))
rownames(lipid_sum)  <- lipid_sum$variable
lipid_sum            <- lipid_sum[, -1]
rownames(lipid_sum)  <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", rownames(lipid_sum)))
lipid_sum            <- lipid_sum[rownames(tmp.df.pheno), ]
lipid_sum            <- log2(lipid_sum + 1)

lipid.avg.df            <- lipid.derived.list$derived_avgAggr_data[!lipid.derived.list$derived_avgAggr_data$isUnknown, ]
lipid.avg.df            <- dplyr::select(lipid.avg.df, dplyr::all_of(colnames(lipid.avg.df)[grepl("BXD|C57|DBA", colnames(lipid.avg.df), ignore.case = T)]))
colnames(lipid.avg.df)  <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", colnames(lipid.avg.df)))
lipid.avg.df            <- as.data.frame(t(lipid.avg.df))
lipid.avg.df[lipid.avg.df == Inf]  <- NA
lipid.avg.df[lipid.avg.df == -Inf] <- NA
lipid.avg.df            <- log2(lipid.avg.df + 1)
lipid.avg.df            <- lipid.avg.df[rownames(pheno.avg.df), ]
remove(lipid.derived.list)



stopifnot(identical(rownames(tmp.df.pheno), rownames(lipid_sum)))
stopifnot(identical(rownames(tmp.df.pheno), rownames(lipid.avg.df)))
stopifnot(identical(rownames(tmp.df.pheno), rownames(lipid_perc)))

mat_list <- list("phenotypes"                           = tmp.df.pheno,
                 "lipid_pclasses_sum"                   = lipid_sum,
                 "lipid_pclasses_perc"                  = compositions::ilr(lipid_perc),
                 "lipid_species"                        = lipid.avg.df[, apply(lipid.avg.df, 2, function(x) all(!is.na(x)))],
                 "phenotypes__and__lipid_pclasses_sum"  = cbind(tmp.df.pheno, lipid_sum),
                 "phenotypes__and__lipid_pclasses_perc" = cbind(tmp.df.pheno, compositions::ilr(lipid_perc)),
                 "phenotypes__and__lipid_species"       = cbind(tmp.df.pheno, lipid.avg.df[, apply(lipid.avg.df, 2, function(x) all(!is.na(x)))]))

mahal.dist.list <- lapply(names(mat_list), function(x){
  # x <- names(mat_list)[1]
  
  tryCatch({
    
    mat_all  <- mat_list[[x]]
    mat_cd   <- mat_all[grepl("CD", rownames(mat_all)), ]
    pca_cd   <- prcomp(mat_cd, center = TRUE, scale. = TRUE)
    expl     <- cumsum(pca_cd$sdev^2) / sum(pca_cd$sdev^2)
    nb_pc    <- which(expl >= 0.9)[1]  # number of PCs to retain
    # Project both datasets onto first k PCs
    mat_cd   <- predict(pca_cd, newdata = mat_all[grepl("CD", rownames(mat_all)), ])[, 1:nb_pc]
    mat_hfd  <- predict(pca_cd, newdata = mat_all[grepl("HFD", rownames(mat_all)), ])[, 1:nb_pc]
    mat_all  <- rbind(mat_cd, mat_hfd)
    
    cov_mat  <- cov(mat_cd)
    med_vec  <- apply(mat_cd, 2, median)
    
    md2     <- mahalanobis(mat_all, center = med_vec, cov = cov_mat)  # squared distances
    md      <- sqrt(md2)
    p       <- ncol(mat_all)
    pvals   <- pchisq(md2, df = p, lower.tail = FALSE) # https://www.r-bloggers.com/2021/08/how-to-calculate-mahalanobis-distance-in-r/
    alpha   <- 0.05
    
    stopifnot(identical(names(md2), names(md)))
    stopifnot(identical(names(md2), names(pvals)))
    
    out <- data.frame(sample_id = names(md2),
                      md2       = unname(md2),
                      md        = unname(md),
                      pval      = unname(pvals),
                      stringsAsFactors = F)
    out$input_data_type <- x
    out$baseline_type   <- "med_CD_baseline_pca"
    out$cov_mat_type    <- "cov_CD"
    out
  }, error = function(e){
    cat("\n\n")
    print(e)
    cat("\n\n")
    NULL
  })
})
names(mahal.dist.list) <- names(mat_list)
mahal.dist.list        <- lapply(mahal.dist.list, function(x) rbindlist(unlist(x, recursive = F)))

saveDir <- "./Data/lipidomics/pheno_mahalanobis_distance"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}
saveRDS(mahal.dist.list, paste0(saveDir, "/pheno_lipidomic_mahalanobis_distance_list.RDS"))








################################################################################
## Phenotypic z-scores
################################################################################



pheno.all <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_list_all_fed_fasted.RDS")
df.pheno  <- pheno.all$avg_fed_and_fasted


lipid.derived.list  <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")

df.lipid            <- lipid.derived.list$derived_avgAggr_data[!lipid.derived.list$derived_avgAggr_data$isUnknown, ]
df.lipid            <- dplyr::select(df.lipid, dplyr::all_of(colnames(df.lipid)[grepl("BXD|C57|DBA", colnames(df.lipid), ignore.case = T)]))
colnames(df.lipid)  <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", colnames(df.lipid)))
df.lipid            <- as.data.frame(t(df.lipid))
df.lipid[df.lipid == Inf]  <- NA
df.lipid[df.lipid == -Inf] <- NA
df.lipid            <- log2(df.lipid + 1)
df.lipid            <- df.lipid[rownames(pheno.avg.df), ]

df.lipid.class            <- as.data.frame(t(lipid.derived.list$derived_avgAggr_lipid_class_sum))
df.lipid.class[df.lipid.class == Inf]  <- NA
df.lipid.class[df.lipid.class == -Inf] <- NA
rownames(df.lipid.class)  <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", rownames(df.lipid.class)))
df.lipid.class            <- df.lipid.class[rownames(df.lipid), ]
df.lipid.class            <- df.lipid.class[, grepl("derived_lipid_class_sum", colnames(df.lipid.class))]


df.pheno.zscore       <- scale(df.pheno)
df.lipid.zscore       <- scale(df.lipid)
df.lipid.class.zscore <- scale(df.lipid.class)
df.lipid.zscore       <- df.lipid.zscore[rownames(df.pheno.zscore), ]
df.lipid.class.zscore <- df.lipid.class.zscore[rownames(df.pheno.zscore), ]
stopifnot(all(rownames(df.lipid.zscore) %in% rownames(df.pheno.zscore)))
stopifnot(all(rownames(df.pheno.zscore) %in% rownames(df.lipid.zscore)))
stopifnot(all(rownames(df.pheno.zscore) %in% rownames(df.lipid.class.zscore)))

df.all.zscore <- as.data.frame(cbind(df.pheno.zscore, df.lipid.zscore, df.lipid.class.zscore))
df.all.zscore <- cbind("strain" = gsub("_.*", "", rownames(df.all.zscore)), "diet" = gsub(".*_", "", rownames(df.all.zscore)), df.all.zscore)


saveDir <- "./Data/lipidomics/pheno_zscore"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}
saveRDS(df.all.zscore, paste0(saveDir, "/pheno_lipidomic_zscores.RDS"))











