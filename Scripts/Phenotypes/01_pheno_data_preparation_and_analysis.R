# clear working space
rm(list=ls())
# # clear all plots
# dev.off()
# clear console
cat("\014")


library(data.table)

 
################################################################################
## saver data as RDS
################################################################################

saveDir <- "./Data/input_data/phenotypic_data"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

pheno.ind.df <- fread("./Data/pheno_formatted.txt", data.table = F, stringsAsFactors = F)
saveRDS(pheno.ind.df, paste0(saveDir, "/pheno_formatted_withMeta.RDS"))

################################################################################
## Download fasted study data
## Dataset "Auwerx1" on Mouse Phenome Database (MPD)
# https://pubmed.ncbi.nlm.nih.gov/27284200/
################################################################################

dirSave <- "./Data/input_data/BXD_fasted_MPD"
if(!dir.exists(dirSave)){
  dir.create(dirSave, recursive = T)
}
if(!file.exists(paste0(dirSave, "/Auwerx1_submit.csv"))){
  system(paste0("wget -O ", dirSave, "/Auwerx1_submit.csv https://phenome.jax.org/getcurationfile?id=135&mode=download"))
}


################################################################################
## aggregate by strain (average-aggregation)
################################################################################

dirSave <- "./Data/phenome/phenoData_transformation"
if(!dir.exists(dirSave)){
  dir.create(dirSave, recursive = T)
}

pheno.ind.df <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_withMeta.RDS")

pheno_columns  <- colnames(pheno.ind.df)
pheno_columns  <- pheno_columns[-c(1:5)]
pheno.avg      <- data.table(pheno.ind.df)[, lapply(.SD, function(ff) mean(ff, na.rm = T)), by = c("strain", "diet"), .SDcols = pheno_columns]
pheno.avg[is.na(pheno.avg)] <- NA

saveRDS(pheno.avg, paste0(dirSave, "/pheno_formatted_mean_allValues.csv"))


################################################################################
## compute mean for fasted study
################################################################################


pheno.indiv.fasted <- read.table("./Data/input_data/BXD_fasted_MPD/Auwerx1_submit.csv", sep = ",", header = T, stringsAsFactors = F)

pheno_columns         <- colnames(pheno.indiv.fasted)
pheno_columns         <- pheno_columns[-c(1:4)]
pheno.avg.fasted      <- data.table(pheno.indiv.fasted)[, lapply(.SD, function(ff) mean(ff, na.rm = T)), by = c("strain", "diet"), .SDcols = pheno_columns]

pheno.avg.fasted[is.na(pheno.avg.fasted)]       <- NA


dirSave <- "./Data/phenome/fasted_phenoData_transformation"
if(!dir.exists(dirSave)){
  dir.create(dirSave, recursive = T)
}
write.table(pheno.avg.fasted, paste0(dirSave, "/fasted_pheno_formatted_mean_allValues.csv"), sep = ",", row.names = F)
saveRDS(pheno.avg.fasted, paste0(dirSave, "/fasted_pheno_formatted_mean_allValues.RDS"))


################################################################################
## compute global (independent of strain) phenotypes diet fold change (of 
## averaged pheno) and diet effect size
## Note on LFC computation - diet contrast all strains:
##
## Now:
## log2((s1_HFD + s2_HFD) / 2) - log2((s1_CD + s2_CD) / 2) = 
## log2(((s1_HFD + s2_HFD) / 2) / ((s1_CD + s2_CD) / 2))
################################################################################

pheno.individuals <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_withMeta.RDS")
pheno.individuals <- reshape2::melt(pheno.individuals, id.vars = colnames(pheno.individuals)[1:5])
pheno.individuals <- dplyr::select(pheno.individuals, c("strain", "diet", "variable", "value"))

computeLFC <- function(value_num, diet_code){
  log2(mean(value_num[diet_code == "HFD" & !is.na(value_num)]) / mean(value_num[diet_code == "CD" & !is.na(value_num)]))
}

computeFC <- function(value_num, diet_code){
  LFC <- mean(value_num[diet_code == "HFD" & !is.na(value_num)]) / mean(value_num[diet_code == "CD" & !is.na(value_num)])
}

computePval <- function(value_num, diet_code){
  tryCatch({
    stats::t.test(x = value_num[diet_code == "CD" & !is.na(value_num)], y = value_num[diet_code == "HFD" & !is.na(value_num)])$p.value
  }, error = function(e){
    as.numeric(NA)
  })
}

overall.pheno.DEA              <- data.table(pheno.individuals)[, list(LFC  = computeLFC(value, diet),
                                                                       FC   = computeFC(value, diet),
                                                                       pVal = computePval(value, diet)),
                                                                by = c("variable")]
colnames(overall.pheno.DEA)[1] <- "phenotype"


saveDir <- paste0(dataPath, "/phenome/diet_and_genetic_effects")
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}
saveRDS(overall.pheno.DEA, paste0(saveDir, "/pheno_FC_diet_allStrains.RDS"))
write.table(overall.pheno.DEA, paste0(saveDir, "/pheno_FC_diet_allStrains.csv"), sep = ",", row.names = F)




################################################################################
## prepare object with all phenotypic data
################################################################################



pheno_avg_files     <- list("fed"    = c("./Data/phenome/phenoData_transformation/pheno_formatted_mean_allValues.RDS"),
                            "fasted" = c("./Data/phenome/fasted_phenoData_transformation/fasted_pheno_formatted_mean_allValues.RDS"))
pheno.avg.df.list <- lapply(names(pheno_avg_files), function(x){
  # x <- names(pheno_avg_files)[2]
  lapply(pheno_avg_files[[x]], function(y){
    # y <- pheno_avg_files[[x]][1]
    # print(y)
    df.y            <- readRDS(y)
    strain_col      <- colnames(df.y)[grepl("strain", colnames(df.y), ignore.case = T)]
    diet_col        <- colnames(df.y)[grepl("^diet", colnames(df.y), ignore.case = T)]
    df.y$strainDiet <- paste0(df.y[[strain_col]], "_", df.y[[diet_col]])
    df.y            <- dplyr::select(df.y, -dplyr::all_of(c(strain_col, diet_col)))
    if(x == "fasted"){
      colnames(df.y)[colnames(df.y) != "strainDiet"] <- paste0("fasted__", colnames(df.y)[colnames(df.y) != "strainDiet"])
    }
    df.y
  })
})
pheno.avg.df <- Reduce(function(...) merge(..., all = T, by = c("strainDiet")), unlist(pheno.avg.df.list, recursive = F))
pheno.avg.df <- pheno.avg.df[pheno.avg.df$strainDiet %in% pheno.avg.df.list[[1]][[1]]$strainDiet, ]
pheno.avg.df <- as.data.frame(pheno.avg.df)
rownames(pheno.avg.df) <- pheno.avg.df$strainDiet
pheno.avg.df <- dplyr::select(pheno.avg.df, -dplyr::all_of("strainDiet"))
pheno.avg.df[is.na(pheno.avg.df)] <- NA
remove(pheno_avg_files, pheno.avg.df.list)


#################################
## individuals phenotypic data
#################################

pheno.ind.df.fed            <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_withMeta.RDS")
# pheno.ind.df.fed            <- readRDS("./pheno_formatted_withMeta.RDS")
pheno.ind.df.fed$strainDiet <- paste0(pheno.ind.df.fed$strain, "_", pheno.ind.df.fed$diet)
pheno.ind.df.fed            <- pheno.ind.df.fed[pheno.ind.df.fed$strainDiet %in% rownames(pheno.avg.df), ]
rownames(pheno.ind.df.fed)  <- paste0(pheno.ind.df.fed$strain, "_", pheno.ind.df.fed$diet, "_", pheno.ind.df.fed$strainReplicate, "_", pheno.ind.df.fed$dietReplicate)
pheno.ind.df.fed            <- dplyr::select(pheno.ind.df.fed, -dplyr::all_of(c("strain", "diet", "strainReplicate", "dietReplicate", "sex")))
pheno.ind.df.fed            <- unique(as.data.frame(pheno.ind.df.fed))
pheno.ind.df.fed[is.na(pheno.ind.df.fed)] <- NA


pheno.ind.df.fasted            <- read.table("./Data/input_data/BXD_fasted_MPD/Auwerx1_submit.csv", sep = ",", header = T, stringsAsFactors = F)
pheno.ind.df.fasted$strainDiet <- paste0(pheno.ind.df.fasted$strain, "_", pheno.ind.df.fasted$diet)
pheno.ind.df.fasted            <- pheno.ind.df.fasted[pheno.ind.df.fasted$strainDiet %in% rownames(pheno.avg.df), ]
rownames(pheno.ind.df.fasted)  <- paste0(pheno.ind.df.fasted$strain, "_", pheno.ind.df.fasted$diet, "_", pheno.ind.df.fasted$id)
pheno.ind.df.fasted            <- dplyr::select(pheno.ind.df.fasted, -dplyr::all_of(c("strain", "sex", "id", "diet")))
pheno.ind.df.fasted            <- unique(as.data.frame(pheno.ind.df.fasted))
pheno.ind.df.fasted[is.na(pheno.ind.df.fasted)] <- NA
colnames(pheno.ind.df.fasted)  <- paste0("fasted__", colnames(pheno.ind.df.fasted))

pheno.all <- list(avg_fed_and_fasted = pheno.avg.df,
                  ind_fed            = pheno.ind.df.fed,
                  ind_fasted         = pheno.ind.df.fasted)

saveRDS(pheno.all, "./Data/input_data/phenotypic_data/pheno_formatted_list_all_fed_fasted.RDS")
# saveRDS(pheno.all, "./pheno_formatted_list_all_fed_fasted.RDS")





################################################################################
## linear mixed models and models comparison
################################################################################


library(lme4)
library(rlang)

pheno.all <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_list_all_fed_fasted.RDS")

saveDir <- "./Plots/pheno_exploration/pheno_fc_new"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}

# restrict to just some phenotypes
cols_order <- c("LV_mass_echo_corrected/Tibia_length_[mg/mm]", "BW_Week27_echo_[g]",
                "LVID_diastole_[mm]", "LVID_systole_[mm]", "LVPW_diastole_[mm]", "LVPW_systole_[mm]", "IVS_diastole_[mm]", "IVS_systole_[mm]", 
                "EF_[%]", "FS_[%]",
                "IVRT_[ms]", "IVCT_[ms]", "DTE_[ms]", "A_duration_[ms]", "MV_A_Vel_[mm/s]", "MV_E_Vel_[mm/s]", "MV_E/A_correct")

tmp.df <- type.convert(reshape2::melt(as.matrix(pheno.all$ind_fed)), as.is = T)
tmp.df <- tmp.df[!is.na(tmp.df$value), ]
tmp.df <- tmp.df[!(tmp.df$Var2 %in% "BW_Week27_echo_[g]"), ]
tmp.df <- tmp.df[tmp.df$Var2 %in% cols_order, ]
cols_order[!c(cols_order %in% tmp.df$Var2)]

tmp.df$diet   <- ifelse(grepl("CD", tmp.df$Var1), "CD", "HFD")
tmp.df$strain <- gsub("_.*", "", tmp.df$Var1)

# model diet fixed effects and strain random effects, to see if some traits are affected by the strain component
nbCores             <- 60
set.seed(123)
stats.full_model.df <- parallel::mclapply(mc.cores = min(c(nbCores, length(unique(tmp.df$Var2)))), X = unique(tmp.df$Var2), FUN = function(x){
  # x <- unique(tmp.df$Var2)[1]
  
  df.x             <- tmp.df[tmp.df$Var2 == x, ]
  df.x             <- df.x[!is.na(df.x$value), ]
  
  # same random slope for each random term
  set.seed(123)
  fit0 <- lmerTest::lmer(value ~ diet + (1 | strain), data = df.x, REML = F)
  sum0 <- summary(fit0)
  # diet effect may differ by strain - use also a random slope
  set.seed(123)
  fit1 <- lmerTest::lmer(value ~ diet + (1 + diet | strain), data = df.x, REML = F)
  sum1 <- summary(fit1)
  
  
  # test benefit of adding random term
  rt0 <- lmerTest::ranova(fit0)
  rt1 <- lmerTest::ranova(fit1)
  
  # compare models
  aov  <- anova(fit0, fit1)
  # better for LMMs
  comp <- pbkrtest::PBmodcomp(largeModel = fit1, smallModel = fit0, nsim = 1000, seed = 123) # more reliable p-value for the random-slope term, use a parametric bootstrap
  
  data.frame(trait_id    = x,
             model_0     = "value ~ diet + (1 | strain)",
             model_1     = "value ~ diet + (1 + diet | strain)",
             coef_intercept_0 = sum0$coefficients["(Intercept)", "Estimate"],
             pval_intercept_0 = sum0$coefficients["(Intercept)", "Pr(>|t|)"],
             coef_diet_0      = sum0$coefficients["dietHFD", "Estimate"],
             pval_diet_0      = sum0$coefficients["dietHFD", "Pr(>|t|)"],
             coef_intercept_1 = sum1$coefficients["(Intercept)", "Estimate"],
             pval_intercept_1 = sum1$coefficients["(Intercept)", "Pr(>|t|)"],
             coef_diet_1      = sum1$coefficients["dietHFD", "Estimate"],
             pval_diet_1      = sum1$coefficients["dietHFD", "Pr(>|t|)"],
             pval_raov0       = rt0$`Pr(>Chisq)`[2],
             raov0_term       = "(1 | strain)",
             pval_raov1       = rt1$`Pr(>Chisq)`[2],
             raov1_term       = "(1 + diet | strain)",
             pval_aov         = aov$`Pr(>Chisq)`[2],
             pval_bootstrap   = comp$test["PBtest", "p.value"],
             stringsAsFactors = F)
  
})
stats.full_model.df <- do.call(rbind, stats.full_model.df)


stats.df <- data.table(tmp.df)[, list(log2FC = log2(mean(value[diet == "HFD"])) - log2(mean(value[diet == "CD"])),
                                      pval   = ifelse(sum(diet == "HFD") >= 2 & sum(diet == "CD") >= 2,
                                                      t.test(x = value[diet == "HFD"], y = value[diet == "CD"])$p.value,
                                                      as.numeric(NA)),
                                      n_CD   = sum(diet == "CD"),
                                      n_HFD  = sum(diet == "HFD")), by = c("strain", "Var2")]
stats.df$pval.adj <- p.adjust(stats.df$pval, "BH")
stats.df$star     <- gtools::stars.pval(stats.df$pval)
stats.df$star.adj <- gtools::stars.pval(stats.df$pval.adj)


tmp.df           <- as.data.frame(dcast(stats.df, strain~Var2, value.var = "log2FC"))
rownames(tmp.df) <- tmp.df[, 1]
tmp.df           <- tmp.df[, -1]
hc               <- hclust(stats::dist(tmp.df), method = "ward.D2")

qn <- quantile(stats.df$log2FC, c(0.02, 0.98))
stats.df$log2FC[stats.df$log2FC < qn[1]] <- qn[1]
stats.df$log2FC[stats.df$log2FC > qn[2]] <- qn[2]

vars_order <- c("IVCT_[ms]", "MV_E/A_correct", "MV_A_Vel_[mm/s]", "MV_E_Vel_[mm/s]", "DTE_[ms]", "IVRT_[ms]", "A_duration_[ms]",
                "EF_[%]", "FS_[%]", "LVID_diastole_[mm]", "LVID_systole_[mm]", "IVS_diastole_[mm]", "IVS_systole_[mm]",
                "LVPW_diastole_[mm]", "LVPW_systole_[mm]", "BW_Week27_echo_[g]", "LV_mass_echo_corrected/Tibia_length_[mg/mm]")
stats.df$Var2   <- factor(stats.df$Var2, levels = rev(vars_order[vars_order %in% unique(stats.df$Var2)]))
stats.df$strain <- factor(stats.df$strain, levels = hc$labels[hc$order])
stats.df$log2FC[is.na(stats.df$pval)] <- NA



pl <- ggplot(stats.df, aes(x = strain, y = Var2, fill = log2FC)) +
  geom_tile() +
  geom_text(aes(label = star)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(low = "blue", high = "red", na.value = "#000000") +
  coord_equal() +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "#2f2f2f"))
# pl
ggsave(paste0(saveDir, "/strain_diet_FC_heatmap__nominal_pVal.pdf"), plot = pl, width = 10, height = 5, useDingbats = F, limitsize = F)

pl <- ggplot(stats.df, aes(x = strain, y = Var2, fill = log2FC)) +
  geom_tile() +
  geom_text(aes(label = star.adj)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(low = "blue", high = "red", na.value = "#000000") +
  coord_equal() +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "#2f2f2f"))
# pl
ggsave(paste0(saveDir, "/strain_diet_FC_heatmap__adjusted_pVal.pdf"), plot = pl, width = 10, height = 5, useDingbats = F, limitsize = F)

stats.full_model.df$trait_id            <- factor(stats.full_model.df$trait_id, levels = rev(vars_order[vars_order %in% unique(stats.full_model.df$trait_id)]))
stats.full_model.df$pval_bootstrap_adj  <- p.adjust(stats.full_model.df$pval_bootstrap, "BH")
stats.full_model.df$pval_raov0_adj      <- p.adjust(stats.full_model.df$pval_raov0, "BH")
stats.full_model.df$star_bootstrap      <- gtools::stars.pval(stats.full_model.df$pval_bootstrap)
stats.full_model.df$star_bootstrap_adj  <- gtools::stars.pval(stats.full_model.df$pval_bootstrap_adj)
stats.full_model.df$star_raov0_term     <- gtools::stars.pval(stats.full_model.df$pval_raov0)
stats.full_model.df$star_raov0_term_adj <- gtools::stars.pval(stats.full_model.df$pval_raov0_adj)

pl.1 <- ggplot(stats.full_model.df, aes(x = 1, y = trait_id, fill = -log10(stats.full_model.df$pval_bootstrap))) +
  geom_tile() +
  geom_text(aes(label = star_bootstrap)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = pals::magma(200)) +
  coord_equal()
# pl.1
ggsave(paste0(saveDir, "/strain_diet_strain_lmm_parametric_bootstrap__nominal_pVal.pdf"), plot = pl.1, width = 10, height = 5, useDingbats = F, limitsize = F)

pl.1.adj <- ggplot(stats.full_model.df, aes(x = 1, y = trait_id, fill = -log10(stats.full_model.df$pval_bootstrap_adj))) +
  geom_tile() +
  geom_text(aes(label = star_bootstrap_adj)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = pals::magma(200)) +
  coord_equal()
# pl.1.adj
ggsave(paste0(saveDir, "/strain_diet_strain_lmm_parametric_bootstrap__adjusted_pVal.pdf"), plot = pl.1.adj, width = 10, height = 5, useDingbats = F, limitsize = F)


pl.2 <- ggplot(stats.full_model.df, aes(x = 1, y = trait_id, fill = -log10(stats.full_model.df$pval_raov0))) +
  geom_tile() +
  geom_text(aes(label = star_raov0_term)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = pals::magma(200)) +
  coord_equal()
# pl.2
ggsave(paste0(saveDir, "/strain_strain_lmm_strain_reff_ranova__nominal_pVal.pdf"), plot = pl.2, width = 10, height = 5, useDingbats = F, limitsize = F)


pl.2.adj <- ggplot(stats.full_model.df, aes(x = 1, y = trait_id, fill = -log10(stats.full_model.df$pval_raov0_adj))) +
  geom_tile() +
  geom_text(aes(label = star_raov0_term_adj)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = pals::magma(200)) +
  coord_equal()
# pl.2.adj
ggsave(paste0(saveDir, "/strain_strain_lmm_strain_reff_ranova__adjusted_pVal.pdf"), plot = pl.2.adj, width = 10, height = 5, useDingbats = F, limitsize = F)



################################################################################
## phenotypes pairwise correlation and linear models testing phenotype*diet 
## interaction coeff
################################################################################



pheno.all <- readRDS("./Data/input_data/phenotypic_data/pheno_formatted_list_all_fed_fasted.RDS")

# restrict to just some phenotypes
cols_order <- c("LV_mass_echo_corrected/Tibia_length_[mg/mm]", "BW_Week27_echo_[g]",
                "LVID_diastole_[mm]", "LVID_systole_[mm]", "LVPW_diastole_[mm]", "LVPW_systole_[mm]", "IVS_diastole_[mm]", "IVS_systole_[mm]", 
                "EF_[%]", "FS_[%]",
                "IVRT_[ms]", "IVCT_[ms]", "DTE_[ms]", "A_duration_[ms]", "MV_A_Vel_[mm/s]", "MV_E_Vel_[mm/s]", "MV_E/A_correct")

saveDir <- "./Plots/pheno_exploration/pairwise_corr_heatmaps"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}

tt <- lapply(c("avg_data", "ind_data"), function(zz){
  # zz <- "avg_data"
  # zz <- "ind_data"
  
  
  diets_run <- c("CD", "HFD", "CD_HFD")
  mat_list  <- lapply(diets_run, function(x){
    # x <- "HFD"
    # x <- "CD_HFD"
    
    
    if(zz == "avg_data"){
      tmp.df   <- pheno.all$avg_fed_and_fasted
    } else{
      tmp.df   <- pheno.all$ind_fed
    }
    tmp.df <- dplyr::select(tmp.df, dplyr::all_of(cols_order[cols_order %in% colnames(tmp.df)]))
    
    if(x != "CD_HFD"){
      tmp.df <- tmp.df[grepl(x, rownames(tmp.df)), ]
    }
    
    
  
    ###########
    ## test difference in CD and HFD slopes
    ###########
    
    if(x == "CD_HFD"){
      
      # Note:
      # Use z-scores to avoid out-of scale coefficients
      slopes_list_interaction <- lapply(colnames(tmp.df), function(k1){
        # k1 <- colnames(tmp.df)[1]
        out <- lapply(colnames(tmp.df), function(k2){
          # k2 <- colnames(tmp.df)[1]


          lm_df      <- data.frame(y = scale(tmp.df[[k1]]), x = scale(tmp.df[[k2]]))
          lm_df$diet <- factor(ifelse(grepl("CD", rownames(tmp.df)), "CD", "HFD"), levels = c("CD", "HFD"))
          lm_fit     <- lm(y~x*diet, lm_df)
          lm_summary <- summary(lm_fit)
          coef_mat   <- lm_summary$coefficients
          data.frame(trait_1                 = k1,
                     trait_2                 = k2,
                     r2                      = lm_summary$r.squared,
                     r2_adj                  = lm_summary$adj.r.squared,
                     intercept               = coef_mat["(Intercept)", "Estimate"],
                     coef                    = coef_mat["x", "Estimate"],
                     dietHFD_coef            = coef_mat["dietHFD", "Estimate"],
                     dietHFD_inter_coef      = coef_mat["x:dietHFD", "Estimate"],
                     intercept_pval          = coef_mat["(Intercept)", "Pr(>|t|)"],
                     coef_pval               = coef_mat["x", "Pr(>|t|)"],
                     dietHFD_coef_pval       = coef_mat["dietHFD", "Pr(>|t|)"],
                     dietHFD_inter_coef_pval = coef_mat["x:dietHFD", "Pr(>|t|)"])
          
        })
        do.call(rbind, out)
      })
      slopes_inter_df <- do.call(rbind, slopes_list_interaction)
      
      coef_mat_inter                <- reshape2::dcast(slopes_inter_df, trait_1~trait_2, value.var = "dietHFD_inter_coef")
      rownames(coef_mat_inter)      <- coef_mat_inter[, 1]
      coef_mat_inter                <- as.matrix(coef_mat_inter[, -1])
      diag(coef_mat_inter)          <- NA
      coef_pval_mat_inter           <- reshape2::dcast(slopes_inter_df, trait_1~trait_2, value.var = "dietHFD_inter_coef_pval")
      rownames(coef_pval_mat_inter) <- coef_pval_mat_inter[, 1]
      coef_pval_mat_inter           <- as.matrix(coef_pval_mat_inter[, -1])
      diag(coef_pval_mat_inter)     <- NA
      coef_pval_adj_inter           <- matrix(NA, nrow = nrow(coef_pval_mat_inter), ncol = ncol(coef_pval_mat_inter), dimnames = dimnames(coef_pval_mat_inter))
      lt_idx                        <- lower.tri(coef_pval_mat_inter, diag = FALSE)
      lt_idx                        <- which(lt_idx & !is.na(coef_pval_mat_inter))
      coef_pval_adj_inter[lt_idx]   <- p.adjust(coef_pval_mat_inter[lt_idx], method = "BH")
      coef_pval_adj_inter[upper.tri(coef_pval_adj_inter)] <- t(coef_pval_adj_inter)[upper.tri(coef_pval_adj_inter)]
      coef_star_mat_inter           <- gtools::stars.pval(coef_pval_adj_inter)
      rownames(coef_star_mat_inter) <- colnames(coef_star_mat_inter) <- rownames(coef_pval_adj_inter)
      
    } else{
      
      slopes_inter_df     <- NULL
      coef_mat_inter      <- NULL
      coef_pval_mat_inter <- NULL
      coef_pval_adj_inter <- NULL
      coef_star_mat_inter <- NULL
      
    }
    
    
    ###########
    ## pearson correlation
    ###########
    
    corr_obj           <- WGCNA::corAndPvalue(tmp.df, method = "pearson")
    cor_mat            <- corr_obj$cor
    pval_mat           <- corr_obj$p
    pval_adj           <- matrix(NA, nrow = nrow(pval_mat), ncol = ncol(pval_mat), dimnames = dimnames(pval_mat))
    lt_idx             <- lower.tri(pval_mat, diag = FALSE)
    lt_idx             <- which(lt_idx & !is.na(pval_mat))
    pval_adj[lt_idx]   <- p.adjust(pval_mat[lt_idx], method = "BH")
    pval_adj[upper.tri(pval_adj)] <- t(pval_adj)[upper.tri(pval_adj)]
    star_mat           <- gtools::stars.pval(pval_adj)
    rownames(star_mat) <- colnames(star_mat) <- rownames(pval_adj)
    
    
    
    corr_obj.strain    <- WGCNA::corAndPvalue(t(tmp.df), method = "pearson")
    cor_mat.strain     <- corr_obj.strain$cor
    pval_mat.strain    <- corr_obj.strain$p
    pval_adj.strain    <- matrix(NA, nrow = nrow(pval_mat.strain), ncol = ncol(pval_mat.strain), dimnames = dimnames(pval_mat.strain))
    lt_idx             <- lower.tri(pval_mat.strain, diag = FALSE)
    lt_idx             <- which(lt_idx & !is.na(pval_mat.strain))
    pval_adj.strain[lt_idx]   <- p.adjust(pval_mat.strain[lt_idx], method = "BH")
    pval_adj.strain[upper.tri(pval_adj.strain)] <- t(pval_adj.strain)[upper.tri(pval_adj.strain)]
    star_mat.strain           <- gtools::stars.pval(pval_adj.strain)
    rownames(star_mat.strain) <- colnames(star_mat.strain) <- rownames(pval_adj.strain)
    
    
    ###########
    ## plot
    ###########
    
    
    if(x == "CD_HFD"){
      
      coef_mat_inter_plot       <- coef_mat_inter
      diag(coef_mat_inter_plot) <- NA
      qn                        <- quantile(coef_mat_inter_plot, c(0.02, 0.98), na.rm = T)
      coef_mat_inter_plot[coef_mat_inter_plot < qn[1] & !is.na(coef_mat_inter_plot)] <- qn[1]
      coef_mat_inter_plot[coef_mat_inter_plot > qn[2] & !is.na(coef_mat_inter_plot)] <- qn[2]
      
      
      paletteLength <- 200
      color_vec     <- rev(pals::brewer.rdbu(paletteLength))
      breaks_vec    <- c(seq(min(coef_mat_inter_plot, na.rm = T), 0, length.out = ceiling(paletteLength / 2) + 1),
                         seq(max(coef_mat_inter_plot, na.rm = T) / paletteLength, max(coef_mat_inter_plot, na.rm = T), length.out = floor(paletteLength / 2)))
      
      pheno_order          <- cols_order[cols_order %in% rownames(coef_mat_inter_plot)]
      coef_mat_inter_plot  <- coef_mat_inter_plot[pheno_order, pheno_order]
      coef_mat_inter       <- coef_mat_inter[pheno_order, pheno_order]
      coef_pval_mat_inter  <- coef_pval_mat_inter[pheno_order, pheno_order]
      coef_pval_adj_inter  <- coef_pval_adj_inter[pheno_order, pheno_order]
      coef_star_mat_inter  <- coef_star_mat_inter[pheno_order, pheno_order]
      
      pdf(paste0(saveDir, "/lm_interaction_coeff_heatmap__", x, "__", zz, "__bounded_colorscale.pdf"), width = 14, height = 8, useDingbats = F)
      print(pheatmap::pheatmap(coef_mat_inter_plot,
                               display_numbers   = coef_star_mat_inter,
                               border_color      = NA,
                               clustering_method = "ward.D2",
                               cluster_rows      = F,
                               cluster_cols      = F,
                               cellheight        = 15,
                               cellwidth         = 15,
                               number_color      = "#000000",
                               angle_col         = 45,
                               color             = color_vec,
                               breaks            = breaks_vec,
                               na_col            = "#FFFFFF",
                               main              = paste0("Diet interaction coefficient")))
      dev.off()
    }
    
    
    cor_mat_plot       <- cor_mat
    diag(cor_mat_plot) <- NA
    qn                 <- quantile(cor_mat_plot, c(0.02, 0.98), na.rm = T)
    cor_mat_plot[cor_mat_plot < qn[1] & !is.na(cor_mat_plot)] <- qn[1]
    cor_mat_plot[cor_mat_plot > qn[2] & !is.na(cor_mat_plot)] <- qn[2]
    
    cor_mat_plot.strain       <- cor_mat.strain
    diag(cor_mat_plot.strain) <- NA
    qn                        <- quantile(cor_mat_plot.strain, c(0.02, 0.98), na.rm = T)
    cor_mat_plot.strain[cor_mat_plot.strain < qn[1] & !is.na(cor_mat_plot.strain)] <- qn[1]
    cor_mat_plot.strain[cor_mat_plot.strain > qn[2] & !is.na(cor_mat_plot.strain)] <- qn[2]
    
    
    paletteLength <- 200
    color_vec     <- rev(pals::brewer.rdbu(paletteLength))
    breaks_vec    <- c(seq(min(cor_mat_plot, na.rm = T), 0, length.out = ceiling(paletteLength / 2) + 1),
                       seq(max(cor_mat_plot, na.rm = T) / paletteLength, max(cor_mat_plot, na.rm = T), length.out = floor(paletteLength / 2)))
    
    pheno_order   <- cols_order[cols_order %in% rownames(cor_mat_plot)]
    cor_mat_plot  <- cor_mat_plot[pheno_order, pheno_order]
    cor_mat       <- cor_mat[pheno_order, pheno_order]
    pval_mat      <- pval_mat[pheno_order, pheno_order]
    pval_adj      <- pval_adj[pheno_order, pheno_order]
    star_mat      <- star_mat[pheno_order, pheno_order]
    
    pdf(paste0(saveDir, "/pearson_corr_heatmap__", x, "__", zz, "__bounded_colorscale.pdf"), width = 14, height = 8, useDingbats = F)
    print(pheatmap::pheatmap(cor_mat_plot,
                             display_numbers   = star_mat,
                             border_color      = NA,
                             clustering_method = "ward.D2",
                             cluster_rows      = F,
                             cluster_cols      = F,
                             cellheight        = 15,
                             cellwidth         = 15,
                             number_color      = "#000000",
                             angle_col         = 45,
                             color             = color_vec,
                             breaks            = breaks_vec,
                             na_col            = "#FFFFFF",
                             main              = paste0("Correlation in ", x)))
    dev.off()
    
    color_vec     <- rev(pals::magma(paletteLength))
    color_vec     <- pals::brewer.reds(paletteLength)
    
    pdf(paste0(saveDir, "/pearson_strain_corr_heatmap__", x, "__", zz, "__bounded_colorscale.pdf"), width = 14, height = 8, useDingbats = F)
    print(pheatmap::pheatmap(cor_mat_plot.strain,
                             display_numbers   = star_mat.strain,
                             border_color      = NA,
                             clustering_method = "ward.D2",
                             cluster_rows      = F,
                             cluster_cols      = F,
                             cellheight        = 15,
                             cellwidth         = 15,
                             number_color      = "#000000",
                             angle_col         = 45,
                             color             = color_vec,
                             # breaks            = breaks_vec,
                             na_col            = "#FFFFFF",
                             main              = paste0("Strain correlation in ", x)))
    dev.off()
    
    list(cor_mat          = cor_mat,
         cor_pval_mat     = pval_mat,
         cor_adj_pval_mat = pval_adj,
         star_adj_mat     = star_mat,
         
         lm_inter_coef_summary_df   = slopes_inter_df,
         lm_inter_coef_mat          = coef_mat_inter,
         lm_inter_coef_pval_mat     = coef_pval_mat_inter,
         lm_inter_coef_adj_pval_mat = coef_pval_adj_inter,
         lm_inter_coef_star_adj_mat = coef_star_mat_inter)
    
  })
  names(mat_list) <- diets_run
  
  
  star_mat <- mat_list$CD$star_adj_mat
  star_mat[upper.tri(star_mat)] <- mat_list$HFD$star_adj_mat[upper.tri(mat_list$HFD$star_adj_mat)]
  
  cor_mat <- mat_list$CD$cor_mat
  cor_mat[upper.tri(cor_mat)] <- mat_list$HFD$cor_mat[upper.tri(mat_list$HFD$cor_mat)]
  diag(cor_mat) <- NA
  
  paletteLength <- 200
  color_vec     <- rev(pals::brewer.rdbu(paletteLength))
  breaks_vec    <- c(seq(min(cor_mat, na.rm = T), 0, length.out = ceiling(paletteLength / 2) + 1),
                     seq(max(cor_mat, na.rm = T) / paletteLength, max(cor_mat, na.rm = T), length.out = floor(paletteLength / 2)))
  
  
  pdf(paste0(saveDir, "/pearson_corr_heatmap__CD_lower_diag_HFD_upper_diag__", zz, "__bounded_colorscale.pdf"), width = 14, height = 8, useDingbats = F)
  print(pheatmap::pheatmap(cor_mat,
                           display_numbers   = star_mat,
                           border_color      = NA,
                           clustering_method = "ward.D2",
                           cluster_rows      = F,
                           cluster_cols      = F,
                           cellheight        = 15,
                           cellwidth         = 15,
                           number_color      = "#000000",
                           angle_col         = 45,
                           color             = color_vec,
                           breaks            = breaks_vec,
                           na_col            = "#FFFFFF",
                           main              = "Pearson correlation - CD and HFD"))
  dev.off()
  
  
  mat_list
})












