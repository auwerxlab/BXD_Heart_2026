

library(data.table)
library(limma)


source("./Scripts/lipidomics_fun.R")


################################################################################
## load and format lipidomics raw data
################################################################################

sheetNames        <- openxlsx::getSheetNames("./Data/BXDMice_heart_lipidomicsData.xlsx")
lipid.list        <- lapply(sheetNames, function(x){
  print(x)
  openxlsx::read.xlsx("./Data/input_data/lipidomics/BDXMice_heart_lipidomicsData.xlsx", sheet = x)
})
names(lipid.list) <- sheetNames

saveDir <- "./Data/input_data/lipidomics"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
saveRDS(lipid.list, "./Data/input_data/lipidomics/BXD_heart_lipidomics_formatted.RDS")

################################################################################
## explore lipidomics raw data before normalization
################################################################################

lipid.list      <- readRDS("./Data/input_data/lipidomics/BXD_heart_lipidomics_formatted.RDS")
samples2exclude <- c("BXD89_CD_2")

tt <- lapply(c("rawData_all_samples", "rawData_outliers_rm"), function(kk){
  # x <- "all_samples"
  
  if(kk == "rawData_all_samples"){
    meta        <- lipid.list$MetaData
    raw.lipid   <- lipid.list$Raw_Data
    
  } else{
    meta        <- lipid.list$MetaData
    raw.lipid   <- lipid.list$Raw_Data
    stopifnot(all(samples2exclude %in% meta$SampleName))
    meta        <- meta[!(meta$SampleName %in% samples2exclude), ]
    raw.lipid   <- dplyr::select(raw.lipid, -samples2exclude) 
  }
  diagn.plots <- diagnosticPlots(raw.lipid, meta)
  
  saveFolder <- paste0("./Plots/lipidomics/QC_noUnknown/", kk)
  if(!dir.exists(saveFolder)){
    dir.create(saveFolder, recursive = T)
  }
  
  ggsave(paste0(saveFolder, "/", kk, "__duplicatedLipids_retentionTime.pdf"),                      plot = diagn.plots$pl.scatter.duplicatedLipids.retTime, width = 28, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__duplicatedLipids_intensityRange.pdf"),                     plot = diagn.plots$pl.scatter.duplLipids.range, width = 28, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__duplicatedLipids_intensityRange_ordered.pdf"),             plot = diagn.plots$pl.scatter.duplLipids.range.ordered, width = 28, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_nbSamples_byBatch.pdf"),                           plot = diagn.plots$pl.nbSamples.byBatch, width = 3, height = 1.5, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_samplesIntensity_byPolarization.pdf"),             plot = diagn.plots$pl.boxlot.byPolarization, width = 7, height = 3.5, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_samplesIntensity_overall.pdf"),                    plot = diagn.plots$pl.boxlot.overall, width = 7, height = 2.6, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__dotPlot_cLength_vs_Saturation_byPolarization.pdf"),        plot = diagn.plots$pl.dotplot.cLength.vs.saturation.byPolarization, width = 17, height = 4, useDingbats = F)
  # TODO
  ggsave(paste0(saveFolder, "/", kk, "__dotPlot_cLength_vs_Saturation_overall.pdf"),               plot = diagn.plots$pl.dotplot.cLength.vs.saturation.overall, width = 10, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__dotPlot_cLength_vs_Saturation_by_lipidClass.pdf"),         plot = diagn.plots$pl.dotplot.cLength.vs.saturation.byClass, width = 17, height = 30, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__histogram_saturationRatio_byPolarization.pdf"),            plot = diagn.plots$pl.histo.saturationRatio.byPolarization, width = 6, height = 2.4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__histogram_saturationRatio_overall.pdf"),                   plot = diagn.plots$pl.histo.saturationRatio.overall, width = 3.8, height = 2.4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_nbLipids_byLipidClass_byPolarization.pdf"),        plot = diagn.plots$pl.histo.nbLipidByClass.byPolarization, width = 5.6, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_nbLipids_byLipidClass_overall.pdf"),               plot = diagn.plots$pl.histo.nbLipidByClass.overall, width = 5.6, height = 2.6, useDingbats = F)
  
  tt <- lapply(names(diagn.plots$pl.density.list.byPolarization), function(x){
    # print(x)
    ggsave(paste0(saveFolder, "/", kk, "__densityPlot_colorBy_", x, "_byPolarization.pdf"), plot = diagn.plots$pl.density.list.byPolarization[[x]], width = 5.2, height = 2, useDingbats = F)
    NA
  })
  tt <- lapply(names(diagn.plots$pl.density.list.overall), function(x){
    # print(x)
    ggsave(paste0(saveFolder, "/", kk, "__densityPlot_colorBy_", x, "_overall.pdf"), plot = diagn.plots$pl.density.list.overall[[x]], width = 3, height = 2, useDingbats = F)
    NA
  })
  tt <- lapply(names(diagn.plots$pl.density.pkAnn.list.byPolarization), function(x){
    # print(x)
    ggsave(paste0(saveFolder, "/", kk, "__densityPlot_colorBy_", x, "_annotated_byPolarization.pdf"), plot = diagn.plots$pl.density.pkAnn.list.byPolarization[[x]], width = 5.2, height = 2, useDingbats = F)
    NA
  })
  tt <- lapply(names(diagn.plots$pl.density.pkAnn.list.overall), function(x){
    # print(x)
    ggsave(paste0(saveFolder, "/", kk, "__densityPlot_colorBy_", x, "_annotated_overall.pdf"), plot = diagn.plots$pl.density.pkAnn.list.overall[[x]], width = 3, height = 2, useDingbats = F)
    NA
  })
  
  tt <- lapply(names(diagn.plots$pl.pca.list.samples), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.pca.list.samples[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__pca_plot_samples_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.pca.list.samples[[x]][[y]], width = ifelse(y == "by_polarity", 6, 3.7), height = 2.4, useDingbats = F)
    })
    NA
  })
  tt <- lapply(names(diagn.plots$pl.mds.list.samples), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.mds.list.samples[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__mds_plot_samples_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.mds.list.samples[[x]][[y]], width = ifelse(y == "by_polarity", 6, 3.7), height = 2.4, useDingbats = F)
    })
    NA
  })
  tt <- lapply(names(diagn.plots$pl.umap.list.samples), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.umap.list.samples[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__umap_plot_samples_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.umap.list.samples[[x]][[y]], width = ifelse(y == "by_polarity", 6, 3.7), height = 2.4, useDingbats = F)
    })
    NA
  })
  
  tt <- lapply(names(diagn.plots$pl.pca.list.lipids), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.pca.list.lipids[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__pca_plot_lipids_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.pca.list.lipids[[x]][[y]], width = ifelse(y == "by_polarity", 6, 4.1), height = 2.4, useDingbats = F)
    })
    NA
  })
  tt <- lapply(names(diagn.plots$pl.mds.list.lipids), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.mds.list.lipids[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__mds_plot_lipids_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.mds.list.lipids[[x]][[y]], width = ifelse(y == "by_polarity", 6, 4.1), height = 2.4, useDingbats = F)
    })
    NA
  })
  tt <- lapply(names(diagn.plots$pl.umap.list.lipids), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.umap.list.lipids[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__umap_plot_lipids_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.umap.list.lipids[[x]][[y]], width = ifelse(y == "by_polarity", 6, 4.1), height = 2.4, useDingbats = F)
    })
    NA
  })
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_sample_vs_blank_control_diff.pdf"), plot = diagn.plots$pl.sample.blank.diff, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__scatterPlot_median_sample_vs_blank_control_diff.pdf"), plot = diagn.plots$pl.median.sample.blank.diff, width = 10, height = 3, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_sample_vs_pooled_control_diff.pdf"), plot = diagn.plots$pl.sample.pooled.diff, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_median_sample_vs_pooled_control_diff.pdf"), plot = diagn.plots$pl.median.sample.pooled.diff, width = 10, height = 3, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__scatterPlot_RT_vs_mz_noUnknowns.pdf"), plot = diagn.plots$pl.rt.vs.mz, width = 10, height = 9, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__scatterPlot_RT_vs_mz_withUnknowns.pdf"), plot = diagn.plots$pl.rt.vs.mz.withUnknown, width = 10, height = 9, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__densityPlot_pooled_controls.pdf"), plot = diagn.plots$pl.pooled_intensity.distrib.all, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__densityPlot_pooled_controls_noUnknowns.pdf"), plot = diagn.plots$pl.pooled_intensity.distrib.noUnknown, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__densityPlot_pooled_controls_unknowns.pdf"), plot = diagn.plots$pl.pooled_intensity.distrib.unknown, width = 16, height = 12, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__densityPlot_blank_controls.pdf"), plot = diagn.plots$pl.blank_intensity.distrib, width = 9, height = 5, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_sample_mass_variation_batches.pdf"), plot = diagn.plots$pl.weight_variation.batch, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_sample_mass_variation_overall.pdf"), plot = diagn.plots$pl.weight_variation.all, width = 25, height = 4, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__pca_plot_with_pooled_controls.pdf"), plot = diagn.plots$pl.pca.with_pooled_controls, width = 3.7, height = 2.6, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__pca_plot_with_pooled_controls_color_by_batch.pdf"), plot = diagn.plots$pl.pca.with_pooled_controls.batches, width = 6.65, height = 5, useDingbats = F)
  
})



################################################################################
## Data normalization and processing
################################################################################

lipid.list <- readRDS("./Data/input_data/lipidomics/BXD_heart_lipidomics_formatted.RDS")
lipid.list <- lipid.list[names(lipid.list) %in% c("MetaData", "Raw_Data")]
lipids.df  <- lipid.list$Raw_Data
meta.df    <- lipid.list$MetaData

#############################
## (1) exclude outliers
#############################

samples2exclude     <- c("BXD89_CD_2")
lipids.df           <- dplyr::select(lipids.df, -dplyr::all_of(samples2exclude))
message_1           <- paste0(length(samples2exclude), " outlier sample excluded...")
lipid.list$Raw_Data <- lipids.df

cols_meta   <- colnames(lipids.df)[!grepl("BXD|PWT|Blank|DBA|C57", colnames(lipids.df))]
cols_sample <- colnames(lipids.df)[!(colnames(lipids.df) %in% cols_meta) & !grepl("PWT|Blank", colnames(lipids.df))]
cols_pooled <- colnames(lipids.df)[!(colnames(lipids.df) %in% cols_meta) & grepl("PWT", colnames(lipids.df))]
cols_blank  <- colnames(lipids.df)[!(colnames(lipids.df) %in% cols_meta) & grepl("Blank", colnames(lipids.df))]


#############################
## (2) sample weight normalization
#############################


lipids.df.massNorm <- lipids.df
for(i in c(cols_sample, cols_pooled)){
  # i <- cols_sample[1]
  
  sample_mass             <- meta.df$Weight[meta.df$SampleName == i]
  lipids.df.massNorm[[i]] <- lipids.df.massNorm[[i]] / sample_mass
  
}
lipids.df.massNorm[is.na(lipids.df.massNorm)] <- NA
lipid.list$MassNorm_Data                      <- lipids.df.massNorm

remove(i, sample_mass)


#############################
## (4.1) set to NA all entries with exactly the same intensity for the same lipid 
##       for different retention times within a sample. Based on the raw data
## (4.2) (a) filter entries with less than a two-fold difference between the 
##           sample and batch-blank. Based on the raw data
##       (b) filter entries with less than a two-fold difference between the 
##           sample and blank avg across batches --> Useful for batch 3 that 
##           doesn't have a blank control. Based on the raw data
#############################

tmp.df            <- type.convert(reshape2::melt(lipids.df, id.vars = cols_meta, variable.name = "sample"), as.is = T)
tmp.df            <- tmp.df[tmp.df$sample %in% cols_sample, ]
tmp.df$unique_id  <- paste0(gsub("_RT_.*", "", tmp.df$Identifier), "__", tmp.df$sample, "__", tmp.df$value)
tmp.df$unique_id2 <- paste0(tmp.df$Identifier, "__", tmp.df$sample)
tmp.df$dupl       <- tmp.df$unique_id %in% tmp.df$unique_id[duplicated(tmp.df$unique_id)] & !grepl("PWT|Blank", tmp.df$sample)
entries.rm        <- tmp.df$unique_id2[tmp.df$dupl]
message_2         <- paste0(sum(tmp.df$dupl), " duplicated entries (",
                            sum(tmp.df$dupl & !tmp.df$isUnknown), " / ", round((sum(tmp.df$dupl & !tmp.df$isUnknown) / sum(!tmp.df$isUnknown)) * 100, 2), "% known; ",
                            sum(tmp.df$dupl & tmp.df$isUnknown), " / ", round((sum(tmp.df$dupl & tmp.df$isUnknown) / sum(tmp.df$isUnknown)) * 100, 2), "% unknown species entries",
                            ") set to NA. This affects ", 
                            length(unique(gsub("_RT_.*", "", tmp.df$Identifier)[tmp.df$dupl & !tmp.df$isUnknown])), " known and ",
                            length(unique(gsub("_RT_.*", "", tmp.df$Identifier)[tmp.df$dupl & tmp.df$isUnknown])), " unknown lipid species...")

##########
##########

blank_avg.df   <- data.frame(Identifier        = lipids.df$Identifier,
                             blank_avg_overall = apply(lipids.df[, grepl("Blank", colnames(lipids.df))], 1, function(x) mean(x, na.rm = T)), stringsAsFactors = F)
blank.df       <- lipids.df[, grepl("Identifier|Blank", colnames(lipids.df))]
blank.df       <- type.convert(reshape2::melt(blank.df, id.vars = c("Identifier"), value.name = "blank", variable.name = "blank_sample"), as.is = T)
blank.df$batch <- plyr::mapvalues(blank.df$blank_sample, from = meta.df$SampleName, to = meta.df$Batch, warn_missing = F)
tmp.df$batch   <- plyr::mapvalues(tmp.df$sample, from = meta.df$SampleName, to = meta.df$Batch, warn_missing = F)
tmp.df         <- merge(tmp.df, blank_avg.df, by = "Identifier", all = T)
tmp.df         <- merge(tmp.df, blank.df, by = c("Identifier", "batch"), all = T)

tmp.df$sample_blank_ratio     <- tmp.df$value / tmp.df$blank
tmp.df$sample_avg_blank_ratio <- tmp.df$value / tmp.df$blank_avg_overall

tmp.df$blank_rm <- (tmp.df$sample_blank_ratio <= 2 & !is.na(tmp.df$blank) & !is.na(tmp.df$value)) | (tmp.df$sample_avg_blank_ratio <= 2 & !is.na(tmp.df$sample_avg_blank_ratio) & !is.na(tmp.df$value))
entries.rm      <- unique(c(entries.rm, tmp.df$unique_id2[tmp.df$blank_rm]))

message_3    <- paste0(sum(tmp.df$blank_rm & !tmp.df$dupl), " intensities lower than blank (",
                       sum(tmp.df$blank_rm & !tmp.df$dupl & !tmp.df$isUnknown), " / ", round((sum(tmp.df$blank_rm & !tmp.df$dupl & !tmp.df$isUnknown) / sum(!tmp.df$isUnknown)) * 100, 2), "% known; ",
                       sum(tmp.df$blank_rm & !tmp.df$dupl & tmp.df$isUnknown), " / ", round((sum(tmp.df$blank_rm & !tmp.df$dupl & tmp.df$isUnknown) / sum(tmp.df$isUnknown)) * 100, 2), "% unknown species entries",
                       ") set to NA. This affects ", 
                       length(unique(gsub("_RT_.*", "", tmp.df$Identifier)[tmp.df$blank_rm & !tmp.df$dupl & !tmp.df$isUnknown])), " known and ",
                       length(unique(gsub("_RT_.*", "", tmp.df$Identifier)[tmp.df$blank_rm & !tmp.df$dupl & tmp.df$isUnknown])), " unknown lipid species...")

##########
##########

tables_filter <- names(lipid.list)[names(lipid.list) != "MetaData"]
filtered.list <- parallel::mclapply(mc.cores = length(tables_filter), X = tables_filter, FUN = function(x){
# filtered.list <- lapply(tables_filter, FUN = function(x){
  # x <- tables_filter[2]
  # print(x)
  
  df.x           <- lipid.list[[x]]
  df.x           <- type.convert(reshape2::melt(df.x, id.vars = cols_meta, variable.name = "sample"), as.is = T)
  df.x$unique_id <- paste0(df.x$Identifier, "__", df.x$sample)
  # stopifnot(all(entries.rm %in% df.x$unique_id))
  entries.rm[!(entries.rm %in% df.x$unique_id)]
  df.x$value[df.x$unique_id %in% entries.rm] <- NA
  df.x           <- reshape2::dcast(df.x, as.formula(paste0(paste(cols_meta, collapse = "+"), "~sample")), value.var = "value")
  df.x
})
names(filtered.list) <- paste0("Filt_", gsub("Raw_", "", tables_filter))
lipid.list           <- c(lipid.list, filtered.list)


remove(tmp.df, blank_avg.df, blank.df, entries.rm, tables_filter, filtered.list)


#############################
## (5) simple batch correction on different kind of filtered data tables
#############################

tables_batch_corr   <- names(lipid.list)[grepl("Filt_", names(lipid.list))]
batch_corr.list     <- lapply(tables_batch_corr, function(x){
  # x <- tables_batch_corr[1]
  # x <- tables_batch_corr[3]
  # print(x)
  
  df.x                <- lipid.list[[x]]
  colnames(df.x)      <- gsub("batch10_", "", colnames(df.x))
  overall_mean        <- unname(apply(dplyr::select(df.x, dplyr::all_of(cols_sample)), 1, function(x) mean(x, na.rm = T)))
  df.x.batchNorm      <- df.x[, !grepl("Blank|PWT", colnames(df.x))]
  normFactors.df      <- dplyr::select(df.x.batchNorm, dplyr::all_of(colnames(df.x.batchNorm)[!grepl("BXD|C57|DBA|PWT|Blank", colnames(df.x.batchNorm))]))
  for(i in unique(meta.df$Batch)){
    # i <- unique(meta.df$Batch)[1]
    
    samples_batch <- meta.df$SampleName[meta.df$Batch == i & !meta.df$is_control_sample]
    samples_batch <- samples_batch[samples_batch %in% colnames(df.x.batchNorm)]
    batch_mean    <- apply(dplyr::select(df.x, dplyr::all_of(samples_batch)), 1, function(x) ifelse(all(is.na(x)), NA, mean(x, na.rm = T)))
    # dplyr::select(df.x, dplyr::all_of(samples_batch))[169, ]
    norm_factors  <- overall_mean / batch_mean
    
    df.x.batchNorm[, samples_batch]       <- df.x.batchNorm[, samples_batch] * norm_factors
    normFactors.df[[paste0("batch_", i)]] <- norm_factors
  }
  df.x.batchNorm[is.na(df.x.batchNorm)] <- NA
  normFactors.df[is.na(normFactors.df)] <- NA
  
  out        <- list(df.x.batchNorm, normFactors.df)
  names(out) <- c(gsub("_Data", "_BatchNorm_Data", x),
                  paste0(x, "_BatchNorm_Factors"))
  out
})
batch_corr.list <- unlist(batch_corr.list, recursive = F)

lipid.list <- c(lipid.list, batch_corr.list)

orer_idx   <- c(which(!grepl("_BatchNorm_Factors", names(lipid.list))),
                which(grepl("_BatchNorm_Factors", names(lipid.list))))
lipid.list <- lipid.list[orer_idx]


remove(tables_batch_corr, batch_corr.list, orer_idx)


#############################
## (6) print data filtering stats
#############################

message_all <- paste0("Data filtering statistics:\n", paste(paste0("\t--> ", c(message_1, message_2, message_3)), collapse = "\n"), "\n\n")
cat(message_all)

#############################
## (7) save processed data
#############################

saveDir <- paste0("./Data/lipidomics/data_processing")
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
saveRDS(lipid.list, paste0(saveDir, "/BXD_heart_lipidomics_filtered_normalized.RDS"))
cat(message_all, file = paste0(saveDir, "/data_filtering_log.txt"), append = F)


################################################################################
## diagnostic plots after filtering and normalization
################################################################################

lipid.list   <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_filtered_normalized.RDS")
entries_plot <- names(lipid.list)[!(names(lipid.list) %in% c("MetaData", "Raw_Data"))]
entries_plot <- entries_plot[!grepl("_BatchNorm_Factors", entries_plot)]

nbCores <- 20
tt      <- parallel::mclapply(mc.cores = min(c(nbCores, length(entries_plot))), X = entries_plot, FUN = function(kk){
  # kk <- entries_plot[1]

  cat("\n\n***************************************************************\nCreating diagnostic plots for", kk, "...\n**************************************************************\n")
  
  diagn.plots <- diagnosticPlots(lipid.list[[kk]], lipid.list$MetaData[lipid.list$MetaData$SampleName %in% colnames(lipid.list[[kk]]), ])
  
  saveFolder <- paste0("./Plots/lipidomics/QC_noUnknown/", kk)
  if(!dir.exists(saveFolder)){
    dir.create(saveFolder, recursive = T)
  }
  
  
  ggsave(paste0(saveFolder, "/", kk, "__duplicatedLipids_retentionTime.pdf"),                      plot = diagn.plots$pl.scatter.duplicatedLipids.retTime, width = 28, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__duplicatedLipids_intensityRange.pdf"),                     plot = diagn.plots$pl.scatter.duplLipids.range, width = 28, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__duplicatedLipids_intensityRange_ordered.pdf"),             plot = diagn.plots$pl.scatter.duplLipids.range.ordered, width = 28, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_nbSamples_byBatch.pdf"),                           plot = diagn.plots$pl.nbSamples.byBatch, width = 3, height = 1.5, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_samplesIntensity_byPolarization.pdf"),             plot = diagn.plots$pl.boxlot.byPolarization, width = 7, height = 3.5, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_samplesIntensity_overall.pdf"),                    plot = diagn.plots$pl.boxlot.overall, width = 7, height = 2.6, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__dotPlot_cLength_vs_Saturation_byPolarization.pdf"),        plot = diagn.plots$pl.dotplot.cLength.vs.saturation.byPolarization, width = 17, height = 4, useDingbats = F)
  # TODO
  ggsave(paste0(saveFolder, "/", kk, "__dotPlot_cLength_vs_Saturation_overall.pdf"),               plot = diagn.plots$pl.dotplot.cLength.vs.saturation.overall, width = 10, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__dotPlot_cLength_vs_Saturation_by_lipidClass.pdf"),         plot = diagn.plots$pl.dotplot.cLength.vs.saturation.byClass, width = 17, height = 30, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__histogram_saturationRatio_byPolarization.pdf"),            plot = diagn.plots$pl.histo.saturationRatio.byPolarization, width = 6, height = 2.4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__histogram_saturationRatio_overall.pdf"),                   plot = diagn.plots$pl.histo.saturationRatio.overall, width = 3.8, height = 2.4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_nbLipids_byLipidClass_byPolarization.pdf"),        plot = diagn.plots$pl.histo.nbLipidByClass.byPolarization, width = 5.6, height = 4, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_nbLipids_byLipidClass_overall.pdf"),               plot = diagn.plots$pl.histo.nbLipidByClass.overall, width = 5.6, height = 2.6, useDingbats = F)
  
  tt <- lapply(names(diagn.plots$pl.density.list.byPolarization), function(x){
    # print(x)
    ggsave(paste0(saveFolder, "/", kk, "__densityPlot_colorBy_", x, "_byPolarization.pdf"), plot = diagn.plots$pl.density.list.byPolarization[[x]], width = 5.2, height = 2, useDingbats = F)
    NA
  })
  tt <- lapply(names(diagn.plots$pl.density.list.overall), function(x){
    # print(x)
    ggsave(paste0(saveFolder, "/", kk, "__densityPlot_colorBy_", x, "_overall.pdf"), plot = diagn.plots$pl.density.list.overall[[x]], width = 3, height = 2, useDingbats = F)
    NA
  })
  tt <- lapply(names(diagn.plots$pl.density.pkAnn.list.byPolarization), function(x){
    # print(x)
    ggsave(paste0(saveFolder, "/", kk, "__densityPlot_colorBy_", x, "_annotated_byPolarization.pdf"), plot = diagn.plots$pl.density.pkAnn.list.byPolarization[[x]], width = 5.2, height = 2, useDingbats = F)
    NA
  })
  tt <- lapply(names(diagn.plots$pl.density.pkAnn.list.overall), function(x){
    # print(x)
    ggsave(paste0(saveFolder, "/", kk, "__densityPlot_colorBy_", x, "_annotated_overall.pdf"), plot = diagn.plots$pl.density.pkAnn.list.overall[[x]], width = 3, height = 2, useDingbats = F)
    NA
  })
  
  tt <- lapply(names(diagn.plots$pl.pca.list.samples), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.pca.list.samples[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__pca_plot_samples_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.pca.list.samples[[x]][[y]], width = ifelse(y == "by_polarity", 6, 3.7), height = 2.4, useDingbats = F)
    })
    NA
  })
  tt <- lapply(names(diagn.plots$pl.mds.list.samples), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.mds.list.samples[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__mds_plot_samples_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.mds.list.samples[[x]][[y]], width = ifelse(y == "by_polarity", 6, 3.7), height = 2.4, useDingbats = F)
    })
    NA
  })
  tt <- lapply(names(diagn.plots$pl.umap.list.samples), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.umap.list.samples[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__umap_plot_samples_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.umap.list.samples[[x]][[y]], width = ifelse(y == "by_polarity", 6, 3.7), height = 2.4, useDingbats = F)
    })
    NA
  })
  
  tt <- lapply(names(diagn.plots$pl.pca.list.lipids), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.pca.list.lipids[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__pca_plot_lipids_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.pca.list.lipids[[x]][[y]], width = ifelse(y == "by_polarity", 6, 4.1), height = 2.4, useDingbats = F)
    })
    NA
  })
  tt <- lapply(names(diagn.plots$pl.mds.list.lipids), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.mds.list.lipids[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__mds_plot_lipids_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.mds.list.lipids[[x]][[y]], width = ifelse(y == "by_polarity", 6, 4.1), height = 2.4, useDingbats = F)
    })
    NA
  })
  tt <- lapply(names(diagn.plots$pl.umap.list.lipids), function(x){
    # print(x)
    tt <- lapply(names(diagn.plots$pl.umap.list.lipids[[x]]), function(y){
      ggsave(paste0(saveFolder, "/", kk, "__umap_plot_lipids_colorBy_", x, "_", ifelse(y == "by_polarity", "byPolarization", "overall"), ".pdf"), 
             plot = diagn.plots$pl.umap.list.lipids[[x]][[y]], width = ifelse(y == "by_polarity", 6, 4.1), height = 2.4, useDingbats = F)
    })
    NA
  })
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_sample_vs_blank_control_diff.pdf"), plot = diagn.plots$pl.sample.blank.diff, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__scatterPlot_median_sample_vs_blank_control_diff.pdf"), plot = diagn.plots$pl.median.sample.blank.diff, width = 10, height = 3, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_sample_vs_pooled_control_diff.pdf"), plot = diagn.plots$pl.sample.pooled.diff, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__boxplot_median_sample_vs_pooled_control_diff.pdf"), plot = diagn.plots$pl.median.sample.pooled.diff, width = 10, height = 3, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__scatterPlot_RT_vs_mz_noUnknowns.pdf"), plot = diagn.plots$pl.rt.vs.mz, width = 10, height = 9, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__scatterPlot_RT_vs_mz_withUnknowns.pdf"), plot = diagn.plots$pl.rt.vs.mz.withUnknown, width = 10, height = 9, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__densityPlot_pooled_controls.pdf"), plot = diagn.plots$pl.pooled_intensity.distrib.all, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__densityPlot_pooled_controls_noUnknowns.pdf"), plot = diagn.plots$pl.pooled_intensity.distrib.noUnknown, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__densityPlot_pooled_controls_unknowns.pdf"), plot = diagn.plots$pl.pooled_intensity.distrib.unknown, width = 16, height = 12, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__densityPlot_blank_controls.pdf"), plot = diagn.plots$pl.blank_intensity.distrib, width = 9, height = 5, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_sample_mass_variation_batches.pdf"), plot = diagn.plots$pl.weight_variation.batch, width = 16, height = 12, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__barPlot_sample_mass_variation_overall.pdf"), plot = diagn.plots$pl.weight_variation.all, width = 25, height = 4, useDingbats = F)
  
  ggsave(paste0(saveFolder, "/", kk, "__pca_plot_with_pooled_controls.pdf"), plot = diagn.plots$pl.pca.with_pooled_controls, width = 3.7, height = 2.6, useDingbats = F)
  ggsave(paste0(saveFolder, "/", kk, "__pca_plot_with_pooled_controls_color_by_batch.pdf"), plot = diagn.plots$pl.pca.with_pooled_controls.batches, width = 6.65, height = 5, useDingbats = F)
  
  
})


################################################################################
## extract lipid species database info using LipidSigR
## Note: LipidSigR is great but most of the analyses are not very flexible as
##       can be achieved with custom scripts
## https://lipidsig.bioinfomics.org/lipidsigr/index.html
################################################################################

lipid.derived.list <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")
lipids_vec         <- lipid.derived.list$Raw_Data$Identifier[!lipid.derived.list$Raw_Data$isUnknown]
lipids_vec         <- unique(gsub("_RT_.*", "", lipids_vec))
lipids_vec_cp      <- lipids_vec
lipids_vec         <- unique(gsub("[NS]", "", lipids_vec, fixed = T)) # Class letter meaning (NS vs NDS): The bracketed subclass NS historically means Non-hydroxy fatty acid + Sphingosine base (sphingosine = d18:1). If your base is d18:0 (dihydrosphingosine / sphinganine), the correct subclass is NDS (non-hydroxy fatty acid + dihydrosphingosine), not NS. Therefore Cer[NS] d18:0_16:0 is internally inconsistent: NS implies a d18:1 base but the name shows d18:0.
lipids_vec         <- unique(gsub("SP d", "Sphingosine ", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("^AC ", "CAR ", lipids_vec))
lipids_vec         <- unique(gsub("Plasmanyl-PC ", "PC(", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("Plasmanyl-PE ", "PE(", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("Plasmenyl-PC ", "PC(", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("Plasmenyl-PE ", "PE(", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("Plasmenyl-PE ", "PE(", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("Alkanyl-DG ", "DG(", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("Alkanyl-DG ", "DG(", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("Alkenyl-DG ", "DG(", lipids_vec, fixed = T))
lipids_vec         <- unique(gsub("Alkenyl-TG ", "TG(", lipids_vec, fixed = T))
lipids_vec[grepl("\\[OH\\] OH", lipids_vec)] <- paste0(gsub("\\[OH\\] OH-", " ", lipids_vec[grepl("\\[OH\\]", lipids_vec)]), ";O")
lipids_vec[grepl("Cer ", lipids_vec)]  <- gsub("_", "/", lipids_vec[grepl("Cer ", lipids_vec)]) # Ceramides always use "/"
lipids_vec[grepl("CerP ", lipids_vec)] <- gsub("_", "/", lipids_vec[grepl("CerP ", lipids_vec)]) # Ceramides always use "/"
lipids_vec[grepl("\\(", lipids_vec) & !grepl("\\)", lipids_vec)] <- paste0(lipids_vec[grepl("\\(", lipids_vec) & !grepl("\\)", lipids_vec)], ")")
names(lipids_vec)  <- lipids_vec_cp

parsed_lipids       <- rgoslin::parseLipidNames(lipids_vec)
table(parsed_lipids$Grammar == "NOT_PARSEABLE")
parsed_lipids       <- parsed_lipids[parsed_lipids$Grammar != "NOT_PARSEABLE", ]
lipid.df            <- lipid.derived.list$Raw_Data
lipid.df            <- cbind("feature" = lipid.df$Identifier, lipid.df)
lipid.df            <- lipid.df[, grepl("feature|BXD|DBA|C57", colnames(lipid.df), ignore.case = F) & !grepl("PWT|Blank", colnames(lipid.df), ignore.case = F)]
# average isomers (otherwise LipidSigR doesn't work)
lipid.df$feature    <- gsub("_RT_.*", "", lipid.df$feature)
lipid.df            <- lipid.df[lipid.df$feature %in% names(lipids_vec), ]
lipid.df$feature    <- lipids_vec[lipid.df$feature]
lipid.df.aggr       <- data.table(lipid.df)[, lapply(.SD, function(zz) mean(zz, na.rm = T)), by = c("feature")]
lipid.df            <- lipid.df.aggr
lipid.df            <- lipid.df[lipid.df$feature %in% parsed_lipids$Original.Name, ]
lipid.df            <- as.data.frame(lipid.df[match(parsed_lipids$Original.Name, lipid.df$feature), ])
meta.df             <- lipid.derived.list$MetaData
meta.df             <- meta.df[meta.df$SampleName %in% colnames(lipid.df), ]
meta.df             <- meta.df[match(colnames(lipid.df)[-1], meta.df$SampleName), ]
stopifnot(identical(colnames(lipid.df)[-1], meta.df$SampleName))
meta.df$sample_name <- meta.df$SampleName
meta.df$label_name  <- meta.df$SampleName
# meta.df$group       <- paste0(meta.df$Strain, "_", meta.df$Diet)
meta.df$group       <- meta.df$Diet
meta.df$pair        <- NA
meta.df             <- dplyr::select(meta.df, dplyr::all_of(c("sample_name", "label_name", "group", "pair")))
se                  <- LipidSigR::as_summarized_experiment(lipid.df, parsed_lipids, meta.df, se_type = "de_two", paired_sample = FALSE)
processed_se        <- LipidSigR::data_process(se, exclude_missing_pct = 100, replace_na_method = "mean", normalization = "none", transform = "none")
lipids_meta         <- as.data.frame(processed_se@elementMetadata)


# save another metadata table to match what done downstream
lipid_clean      <- data.frame(original_id = names(lipids_vec), clean_id = unname(lipids_vec), stringsAsFactors = F)
lipids.meta.df.2 <- merge(lipids_meta, dplyr::select(parsed_lipids, colnames(parsed_lipids)[!(colnames(parsed_lipids) %in% colnames(lipids_meta))]), by.x = "feature", by.y = "Original.Name")
lipids.meta.df.2 <- merge(lipid_clean, lipids.meta.df.2, by.x = "clean_id", by.y = "feature")
lipids.meta.df.2 <- lipids.meta.df.2[, c(which(grepl("id|normalized", colnames(lipids.meta.df.2), ignore.case = T)), which(!grepl("id|normalized", colnames(lipids.meta.df.2), ignore.case = T)))]


parsed_lipids <- cbind(feature_id = plyr::mapvalues(parsed_lipids$Original.Name, from = unname(lipids_vec), to = names(lipids_vec)), parsed_lipids)
lipids_meta   <- cbind(feature_id = plyr::mapvalues(lipids_meta$feature, from = unname(lipids_vec), to = names(lipids_vec)), lipids_meta)


lipids.meta.df <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_metadata.RDS")
lipids.meta.df <- cbind(feature_id = gsub("_RT_.*", "", lipids.meta.df$Identifier), lipids.meta.df)
lipids.meta.df <- lipids.meta.df[!lipids.meta.df$isUnknown, ]
lipids.meta.df$feature_id[!(lipids.meta.df$feature_id %in% lipids_meta$feature_id)]
lipids.meta.df <- merge(lipids.meta.df, lipids_meta, by = "feature_id")

saveRDS(parsed_lipids, "./Data/lipidomics/data_processing/goslin_lipid_meta.RDS")
saveRDS(lipids_meta, "./Data/lipidomics/data_processing/LipidSigR_lipid_meta.RDS")
saveRDS(lipids.meta.df, "./Data/lipidomics/data_processing/BXD_heart_lipidomics_metadata__with_LipidSigR_lipid_meta.RDS")
saveRDS(lipids.meta.df.2, "./Data/lipidomics/data_processing/BXD_heart_goslin_LipidSigR_lipids_meta.RDS")



################################################################################
## Save lipid species metadata
################################################################################

lipid.list   <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_filtered_normalized.RDS")
lipids.df    <- lipid.list$Filt_MassNorm_BatchNorm_Data

# GP = Glycerophospholipids (Phospholipids)
# GL = Glycerolipids
# SP = Sphingolipids
# ST = Sterol Lipids
# FA = Fatty Acyls
lipid_grandparents_map <- c("PC"     = "GP",
                            "PC"     = "GP",
                            "PE"     = "GP",
                            "PG"     = "GP",
                            "PI"     = "GP",
                            "PS"     = "GP",
                            "CL"     = "GP",
                            "MLCL"   = "GP",
                            "LysoPC" = "GP",
                            "LysoPE" = "GP",
                            "LysoPG" = "GP",
                            "LysoPI" = "GP",
                            
                            "DG"     = "GL",
                            "TG"     = "GL",
                            
                            "Cer"    = "SP",
                            "HexCer" = "SP",
                            "SM"     = "SP",
                            "SP"     = "SP",
                            
                            "CE"     = "ST",
                            
                            "AC"     = "FA")


chain_counts <- c("AC"             = 1,
                  "LysoPC"         = 1, 
                  "LysoPE"         = 1,
                  "LysoPG"         = 1, 
                  "LysoPI"         = 1,
                  "DG"             = 2, 
                  "Alkanyl-DG"     = 2,
                  "Alkenyl-DG"     = 2,
                  "PC"             = 2,
                  "PC[OH]"         = 2,
                  "PE"             = 2,
                  "PE-NMe"         = 2,
                  "PE-NMe2"        = 2,
                  "PG"             = 2,
                  "PI"             = 2,
                  "PI[OH]"         = 2,
                  "Plasmanyl-PC"   = 2,
                  "Plasmenyl-PC"   = 2,
                  "Plasmanyl-PE"   = 2,
                  "Plasmenyl-PE"   = 2,
                  "PS"             = 2,
                  "TG"             = 3,
                  "Alkenyl-TG"     = 3,
                  "CL"             = 4,
                  "MLCL"           = 3,
                  "CE"             = 1,
                  "Cer[NS]"        = 2,
                  "CerP"           = 2,
                  "HexCer[NS]"     = 2,
                  "SM"             = 1,
                  "SP"             = 1)
stopifnot(all(unique(lipids.df$Lipid.Class[!is.na(lipids.df$Lipid.Class)]) %in% names(chain_counts)))



lipids.meta.df <- dplyr::select(lipids.df, dplyr::all_of(c("Identifier", "Retention.Time..min.", "Quant.Ion", 
                                                           "Polarity", "Area..max.", "Identification", "Lipid.Class", 
                                                           "isUnknown", "cLength", "nbDoubleBonds", "Features.Found")))

rownames(lipids.meta.df)               <- lipids.meta.df$Identifier
lipids.meta.df$Lipid.Parent.Class      <- gsub("CerP", "Cer", gsub("Alkanyl-|Alkenyl-|Plasmanyl-|Plasmenyl-|-NMe2|-NMe|\\[.*", "", lipids.meta.df$Lipid.Class))
stopifnot(all(unique(lipids.meta.df$Lipid.Parent.Class[!is.na(lipids.meta.df$Lipid.Parent.Class)]) %in% names(lipid_grandparents_map)))
lipids.meta.df$Lipid.Grandparent.Class <- unname(lipid_grandparents_map[lipids.meta.df$Lipid.Parent.Class])
lipids.meta.df$sat_ratio               <- lipids.meta.df$cLength / lipids.meta.df$nbDoubleBonds
lipids.meta.df$nb_chains               <- as.integer(unname(chain_counts[lipids.meta.df$Lipid.Class]))
stopifnot(sum(is.na(lipids.meta.df$nb_chains) & !lipids.meta.df$isUnknown) == 0)
nb_chains_empirical                   <- unlist(lapply(lipids.meta.df$Identifier, function(x){
  elements.x <- unlist(strsplit(gsub("_RT.*", "", x), "_"))
  elements.x <- elements.x[grepl("[0-9]", elements.x)]
  length(elements.x)
}))
stopifnot(all(nb_chains_empirical[!lipids.meta.df$isUnknown] <= lipids.meta.df$nb_chains[!lipids.meta.df$isUnknown]))
stopifnot(all(nb_chains_empirical[!lipids.meta.df$isUnknown & nb_chains_empirical != 1] == lipids.meta.df$nb_chains[!lipids.meta.df$isUnknown & nb_chains_empirical != 1]))

saveRDS(lipids.meta.df, "./Data/lipidomics/data_processing/BXD_heart_lipidomics_metadata.RDS")


################################################################################
## Derive lipidomic features
################################################################################

lipid.list     <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_filtered_normalized.RDS")
lipids.df      <- lipid.list$Filt_MassNorm_BatchNorm_Data
sample_cols    <- colnames(lipids.df)[grepl("^BXD|^DBA|^C57", colnames(lipids.df))]
lipids.meta.df <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_metadata.RDS")


#############################
## Add lipid parent and grandparent class, lipid saturation ratio and lipid nb. of chain (by class)
#############################
tmp.df                <- lipids.df
rownames(tmp.df)      <- tmp.df$Identifier
tmp.df                <- dplyr::select(tmp.df, -c("Retention.Time..min.", "Quant.Ion", "Area..max.", "Features.Found", "Polarity", "Identification"))
tmp.df                <- tmp.df[, !grepl("PWT|Blank|RSD", colnames(tmp.df))]
colnames(tmp.df)[grepl("DBA|C57", colnames(tmp.df))] <- gsub("^BXD", "", colnames(tmp.df)[grepl("DBA|C57", colnames(tmp.df))])
tmp.df                <- cbind(Lipid.Parent.Class = gsub("CerP", "Cer", gsub("Alkanyl-|Alkenyl-|Plasmanyl-|Plasmenyl-|-NMe2|-NMe|\\[.*", "", tmp.df$Lipid.Class)), tmp.df)
stopifnot(all(unique(tmp.df$Lipid.Parent.Class[!is.na(tmp.df$Lipid.Parent.Class)]) %in% names(lipid_grandparents_map)))
tmp.df                <- cbind(Lipid.Grandparent.Class = unname(lipid_grandparents_map[tmp.df$Lipid.Parent.Class]), tmp.df)
tmp.df                <- cbind(sat_ratio = tmp.df$cLength / tmp.df$nbDoubleBonds, tmp.df)
tmp.df                <- cbind(nb_chains = as.integer(unname(chain_counts[tmp.df$Lipid.Class])), tmp.df)
# tmp.df                <- cbind(nb_chains_empirical = nb_chains_empirical, tmp.df)
stopifnot(sum(is.na(tmp.df$nb_chains) & !tmp.df$isUnknown) == 0)
nb_chains_empirical   <- unlist(lapply(tmp.df$Identifier, function(x){
  # x <- tmp.df$Identifier[40]
  elements.x <- unlist(strsplit(gsub("_RT.*", "", x), "_"))
  elements.x <- elements.x[grepl("[0-9]", elements.x)]
  length(elements.x)
}))
stopifnot(all(nb_chains_empirical[!tmp.df$isUnknown] <= tmp.df$nb_chains[!tmp.df$isUnknown]))
stopifnot(all(nb_chains_empirical[!tmp.df$isUnknown & nb_chains_empirical != 1] == tmp.df$nb_chains[!tmp.df$isUnknown & nb_chains_empirical != 1]))


sample_cols   <- colnames(tmp.df)[grepl("^BXD|^DBA|^C57", colnames(tmp.df))]
meta_cols.tmp <- colnames(tmp.df)[!grepl("^BXD|^DBA|^C57|PWT|Blank", colnames(tmp.df))]
meta_cols     <- c("Identifier", "Lipid.Grandparent.Class", "Lipid.Parent.Class", "Lipid.Class",
                   "isUnknown", "cLength", "nbDoubleBonds", "sat_ratio", "nb_chains")
stopifnot(all(meta_cols.tmp %in% meta_cols))

#############################
## (1) strain average-aggregated lipid intensities
#############################
tmp.df.melted               <- type.convert(reshape2::melt(tmp.df,
                                                           id.vars       = meta_cols,
                                                           value.name    = "intensity",
                                                           variable.name = "sample"), as.is = T)
tmp.df.melted               <- tmp.df.melted[grepl("BXD|C57|DBA", tmp.df.melted$sample), ]
tmp.df.melted$strainDiet    <- gsub("CD.*$", "CD", gsub("HFD.*$", "HFD", tmp.df.melted$sample))
lipids.df.aggr              <- as.data.frame(data.table(tmp.df.melted)[, list(avg_int = mean(intensity, na.rm = T)), by = c(meta_cols, "strainDiet")])
lipids.df.aggr              <- reshape2::dcast(lipids.df.aggr, as.formula(paste0(paste(meta_cols, collapse = "+"), "~strainDiet")), value.var = "avg_int")
lipids.df.aggr[is.na(lipids.df.aggr)] <- NA
rownames(lipids.df.aggr)    <- lipids.df.aggr$Identifier

#############################
## (2) total sum / average by class types
#############################

class_meta.df           <- unique(dplyr::select(tmp.df, dplyr::all_of(c("Lipid.Class", "Lipid.Parent.Class", "Lipid.Grandparent.Class"))))
rownames(class_meta.df) <- NULL
class_meta.df           <- class_meta.df[!is.na(class_meta.df$Lipid.Class), ]
for(i in colnames(class_meta.df)){
  class_meta.df[[paste0(i, ".nbLipids")]] <- unlist(lapply(class_meta.df[[i]], function(x) sum(!is.na(tmp.df[[i]]) & tmp.df[[i]] == x & !tmp.df$isUnknown, na.rm = T)))
}
print(knitr::kable(class_meta.df))

class_cols              <- c("class" = "Lipid.Class", "parent_class" = "Lipid.Parent.Class", "grandparent_class" = "Lipid.Grandparent.Class")
lipid_class_aggr_all.df <- lapply(names(class_cols), function(x){
  # x <- names(class_cols)[1]
  
  class_col           <- class_cols[[x]]
  aggr.list.x         <- lapply(c("avg", "sum"), function(y){
    # y <- "avg"
    if(y == "avg"){
      fun_aggr <- function(k){ mean(k, na.rm = T) }
    } else if(y == "sum"){
      fun_aggr <- function(k){ sum(k, na.rm = T) }
    }
    aggr.df.y           <- as.data.frame(data.table(tmp.df[!is.na(tmp.df[[class_col]]), ])[, lapply(.SD, fun_aggr), by = c(class_col), .SDcols = sample_cols])
    rownames(aggr.df.y) <- paste0(x, "_", aggr.df.y[[class_col]], "_", ifelse(y == "avg", "avg", "tot_sum"))
    aggr.df.y[, -1]
  })
  stopifnot(all(sapply(lapply(aggr.list.x, colnames), function(x) identical(x, colnames(aggr.list.x[[1]])))))
  aggr.df.x <- do.call(rbind, aggr.list.x)
  aggr.df.x
})
stopifnot(all(sapply(lapply(lipid_class_aggr_all.df, colnames), function(x) identical(x, colnames(lipid_class_aggr_all.df[[1]])))))
lipid_class_aggr_all.df <- do.call(rbind, lipid_class_aggr_all.df)
lipid_class_aggr_all.df <- as.data.frame(lipid_class_aggr_all.df)

#############################
## Concatenate derived features
#############################


lipid.list.derived   <- lipid.list[names(lipid.list) %in% c("MetaData", "Raw_Data", "Filt_MassNorm_BatchNorm_Data")]
lipid.list.derived$Filt_MassNorm_BatchNorm_Data_NArm <- lipids.df

derived_features.all <- list(derived_lipid_class_sum         = lipid_class_aggr_all.df[grepl("_sum$", rownames(lipid_class_aggr_all.df)), ],
                             derived_lipid_class_avg         = lipid_class_aggr_all.df[grepl("_avg$", rownames(lipid_class_aggr_all.df)), ])

stopifnot(all(sapply(lapply(derived_features.all, ncol), function(x) identical(x, ncol(derived_features.all[[1]])))))
derived_features.all <- lapply(derived_features.all, function(x) dplyr::select(x, dplyr::all_of(colnames(derived_features.all[[1]]))))
stopifnot(all(sapply(lapply(derived_features.all, colnames), function(x) identical(x, colnames(derived_features.all[[1]])))))
derived_feat_names   <- names(derived_features.all)
derived_features.all <- lapply(derived_feat_names, function(x){
  df.x           <- derived_features.all[[x]]
  rownames(df.x) <- paste0(x, "___", rownames(df.x))
  df.x
})
names(derived_features.all) <- derived_feat_names

lipid.list.derived <- c(lipid.list.derived, derived_features.all)

#############################
## Strain-averaged aggregated features
#############################

aggr.tmp.list             <- list(derived_avgAggr_data      = lipids.df.aggr)
derived_features.all.aggr <- derived_features.all[!grepl("PCA", names(derived_features.all))]
derived_features.all.aggr <- lapply(derived_features.all.aggr, function(x){
  # x <- derived_features.all[[1]]
  
  tmp.df.aggr.x            <- t(x)
  tmp.df.aggr.x            <- type.convert(as.data.frame(cbind(strainDiet = gsub("CD.*$", "CD", gsub("HFD.*$", "HFD", rownames(tmp.df.aggr.x))), tmp.df.aggr.x)), as.is = T)
  tmp.df.aggr.x$strainDiet <- as.character(tmp.df.aggr.x$strainDiet)
  df.aggr.x                <- data.table(tmp.df.aggr.x)[, lapply(.SD, function(zz) mean(zz[!is.infinite(zz)], na.rm = T)), by = c("strainDiet")]
  df.aggr.x                <- as.data.frame(df.aggr.x, stringsAsFactors = F)
  rownames(df.aggr.x)      <- as.character(df.aggr.x$strainDiet)
  df.aggr.x                <- dplyr::select(df.aggr.x, -dplyr::all_of("strainDiet"))
  df.aggr.x                <- as.data.frame(t(df.aggr.x))
  df.aggr.x
})
names(derived_features.all.aggr) <- gsub("derived_", "derived_avgAggr_", names(derived_features.all.aggr))
derived_features.all.aggr        <- c(aggr.tmp.list, derived_features.all.aggr)

lipid.list.derived <- c(lipid.list.derived, derived_features.all.aggr)

saveRDS(lipid.list.derived, "./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")


