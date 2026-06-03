library(limma)
library(ggplot2)
library(stringr)
library(data.table)



################################################################################
## Differential abundance analysis - DEA
################################################################################


prot.df         <- readRDS("./Data/input_data/proteomics/coon_formatted/raw_intensity_table_samples_excluded.RDS")
samples.meta.df <- readRDS("./Data/input_data/proteomics/coon_formatted/sample_metadata.RDS")
prot.meta.df    <- readRDS("./Data/input_data/proteomics/coon_formatted/proteins_metadata.RDS")

prot.df         <- prot.df[, colnames(prot.df) %in% samples.meta.df$sample]

# aggregate protein metadata by:
# (1) uniprot gene ID to have a correct table to use for merging in the DEA analysis on protein-level aggregated data
# (2) ensembl gene ID to have a correct table to use for merging in the DEA analysis on gene-level aggregated data
aggr_fun <- function(x){
  if(all(is.na(x))){
    return(as.character(NA))
  }
  paste(sort(unique(x[!is.na(x)])), collapse = ";")
}
prot.meta.df.aggr.byProt <- data.table(prot.meta.df[!is.na(prot.meta.df$uniprot_id), ])[, lapply(.SD, aggr_fun), by = c("uniprot_id")]
prot.meta.df.aggr.byProt <- type.convert(as.data.frame(prot.meta.df.aggr.byProt), as.is = T)
prot.meta.df.aggr.byGene <- data.table(prot.meta.df[!is.na(prot.meta.df$ensembl_gene_id), ])[, lapply(.SD, aggr_fun), by = c("ensembl_gene_id")]
prot.meta.df.aggr.byGene <- type.convert(as.data.frame(prot.meta.df.aggr.byGene), as.is = T)

# calculate completeness, i.e. in how many samples a protein is detected
tmp.df <- reshape2::melt(as.matrix(prot.df))
tmp.df <- merge(tmp.df, samples.meta.df, by.x = "Var2", by.y = "sample")
nb_samples      <- ncol(prot.df)
nb_samples_CD   <- sum(grepl("CD", colnames(prot.df)))
nb_samples_HFD  <- sum(grepl("HFD", colnames(prot.df)))
nb_strains      <- length(unique(samples.meta.df$strain))
completeness.df <- data.table(tmp.df)[, list(completeness_perc_samples      = (sum(!is.na(value)) / nb_samples) * 100,
                                             completeness_perc_samples_CD   = (sum(!is.na(value[diet == "CD"])) / nb_samples_CD) * 100,
                                             completeness_perc_samples_HFD  = (sum(!is.na(value[diet == "HFD"])) / nb_samples_HFD) * 100,
                                             completeness_perc_strains      = (length(unique(strain[!is.na(value)])) / nb_strains) * 100,
                                             completeness_perc_strains_CD   = (length(unique(strain[!is.na(value) & diet == "CD"])) / nb_strains) * 100,
                                             completeness_perc_strains_HFD  = (length(unique(strain[!is.na(value) & diet == "HFD"])) / nb_strains) * 100),
                                      by = c("Var1")]
completeness.df <- type.convert(as.data.frame(completeness.df), as.is = T)
completeness.df <- completeness.df[order(completeness.df$completeness_perc_samples, decreasing = T), ]
colnames(completeness.df)[colnames(completeness.df) == "Var1"] <- "prot_id"
remove(tmp.df)


# log2-transform to "normalize" data
prot.df.log2 <- log2(prot.df)

base::range(prot.df, na.rm = T)
base::range(prot.df.log2, na.rm = T)

# perform median-normalization to partially remove batch effects
sizefactor               <- matrixStats::colMedians(as.matrix(prot.df.log2), na.rm = TRUE)
prot.df.log2.median.norm <- sweep(prot.df.log2, 2, sizefactor)

# keep only proteins expressed in at least half of the samples (using completeness)
completeness_threshold_perc <- 50
tmpLogi.keep                <- (completeness.df$completeness_perc_samples >= completeness_threshold_perc & 
                                  completeness.df$completeness_perc_strains_CD >= completeness_threshold_perc & 
                                  completeness.df$completeness_perc_strains_HFD >= completeness_threshold_perc)
proteins.keep               <- completeness.df$prot_id[tmpLogi.keep]
prot.df.log2.filt           <- prot.df.log2[rownames(prot.df.log2) %in% proteins.keep, ]


# aggregate by gene ensembl id, to also do the analysis by gene (useful for GSEA to avoid having duplicated ranking for genes with multiple protein isoforms)
aggr_fun <- function(x){
  if(all(is.na(x))){
    return(as.numeric(NA))
  }
  mean(x, na.rm = T)
}
tmp.df                 <- merge(unique(dplyr::select(prot.meta.df[!is.na(prot.meta.df$ensembl_gene_id), ], c("uniprot_id", "ensembl_gene_id"))), 
                                prot.df.log2.filt, by.x = "uniprot_id", by.y = "row.names")
tmp.df                 <- dplyr::select(tmp.df, -c("uniprot_id"))
prot.df.log2.filt.aggr <- data.table(tmp.df)[, lapply(.SD, aggr_fun), by = c("ensembl_gene_id")]
prot.df.log2.filt.aggr <- type.convert(as.data.frame(prot.df.log2.filt.aggr), as.is = T)
rownames(prot.df.log2.filt.aggr) <- prot.df.log2.filt.aggr$ensembl_gene_id
prot.df.log2.filt.aggr <- dplyr::select(prot.df.log2.filt.aggr, -c("ensembl_gene_id"))
remove(tmp.df)

# calculate completeness for gene-level aggregated daza,i.e. in how many samples a gene is detected
tmp.df <- reshape2::melt(as.matrix(prot.df.log2.filt.aggr))
tmp.df <- merge(tmp.df, samples.meta.df, by.x = "Var2", by.y = "sample")
completeness.df.gene <- data.table(tmp.df)[, list(completeness_perc_samples      = (sum(!is.na(value)) / nb_samples) * 100,
                                                  completeness_perc_samples_CD   = (sum(!is.na(value[diet == "CD"])) / nb_samples_CD) * 100,
                                                  completeness_perc_samples_HFD  = (sum(!is.na(value[diet == "HFD"])) / nb_samples_HFD) * 100,
                                                  completeness_perc_strains      = (length(unique(strain[!is.na(value)])) / nb_strains) * 100,
                                                  completeness_perc_strains_CD   = (length(unique(strain[!is.na(value) & diet == "CD"])) / nb_strains) * 100,
                                                  completeness_perc_strains_HFD  = (length(unique(strain[!is.na(value) & diet == "HFD"])) / nb_strains) * 100),
                                           by = c("Var1")]
completeness.df.gene <- type.convert(as.data.frame(completeness.df.gene), as.is = T)
completeness.df.gene <- completeness.df.gene[order(completeness.df.gene$completeness_perc_samples, decreasing = T), ]
colnames(completeness.df.gene)[colnames(completeness.df.gene) == "Var1"] <- "ensembl_gene_id"
remove(tmp.df)

# match dimnames
rownames(samples.meta.df) <- samples.meta.df$sample
samples.meta.df           <- samples.meta.df[match(colnames(prot.df.log2), rownames(samples.meta.df)), ]

# check dimnames
stopifnot(identical(rownames(samples.meta.df), colnames(prot.df.log2)))
stopifnot(identical(rownames(samples.meta.df), colnames(prot.df.log2.filt)))
stopifnot(identical(rownames(samples.meta.df), colnames(prot.df.log2.filt.aggr)))

###############
## Create several data-inspection plots
###############

saveDir <- "./Plots/proteomics/coon_QC/data_completeness"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
completeness_threshold_perc <- 50
dfPlot      <- reshape2::melt(completeness.df, id.vars = "prot_id")
dfPlot$what <- gsub("completeness_perc_", "", dfPlot$variable)
map_vec     <- c("samples" = "all samples (CD & HFD)", "samples_CD" = "CD samples", "samples_HFD" = "HFD samples", "strains" = "strains (CD & HFD)", "strains_CD" = "strains (CD)", "strains_HFD" = "strains (HFD)")
dfPlot$what <- map_vec[dfPlot$what]
color_pal   <- c("all samples (CD & HFD)" = "#CF9FFF", "CD samples" = "#FF7518", "HFD samples" = "#FF3131", "strains (CD & HFD)" = "#0437F2", "strains (CD)" = "#32CD32", "strains (HFD)" = "#008080")
dfPlot$what <- factor(dfPlot$what, levels = names(color_pal))
pl <- ggplot(dfPlot, aes(x = value, fill = what)) +
  geom_vline(xintercept = completeness_threshold_perc, color = "black", size = 0.2, linetype = "dashed") +
  geom_histogram(bins = 100, alpha = 0.6) +
  scale_fill_manual(values = color_pal) +
  xlab("Protein expressed in samples [%]") +
  theme_classic() +
  labs(fill = "Completeness\nacross:")
# pl
ggsave(paste0(saveDir, "/completeness_histogram.pdf"), plot = pl, width = 7, height = 3, useDingbats = F, limitsize = F)

completeness_threshold_perc <- 50
dfPlot <- reshape2::melt(completeness.df, id.vars = "prot_id")
dfPlot <- lapply(unique(dfPlot$variable), function(k){
  out <- lapply(0:100, function(x){
    data.frame(perc_thresh          = x,
               nb_prot_above_thresh = sum(dfPlot$value >= x & dfPlot$variable == k),
               what                 = k,
               stringsAsFactors = F)
  })
  do.call(rbind, out)
})
dfPlot      <- do.call(rbind, dfPlot)
dfPlot$what <- gsub("completeness_perc_", "", dfPlot$what)
map_vec     <- c("samples" = "all samples (CD & HFD)", "samples_CD" = "CD samples", "samples_HFD" = "HFD samples", "strains" = "strains (CD & HFD)", "strains_CD" = "strains (CD)", "strains_HFD" = "strains (HFD)")
dfPlot$what <- map_vec[dfPlot$what]
hlines.df   <- lapply(unique(dfPlot$what), function(x){
  data.frame(what        = x,
             y_intercept = dfPlot$nb_prot_above_thresh[dfPlot$perc_thresh == completeness_threshold_perc & dfPlot$what == x],
             stringsAsFactors = F)
})
hlines.df      <- do.call(rbind, hlines.df)
color_pal      <- c("all samples (CD & HFD)" = "#CF9FFF", "CD samples" = "#FF7518", "HFD samples" = "#FF3131", "strains (CD & HFD)" = "#0437F2", "strains (CD)" = "#32CD32", "strains (HFD)" = "#008080")
dfPlot$what    <- factor(dfPlot$what, levels = names(color_pal))
hlines.df$what <- factor(hlines.df$what, levels = names(color_pal))

# pals::pal.bands(color_pal)
pl <- ggplot(dfPlot, aes(x = perc_thresh, y = nb_prot_above_thresh, color = what)) +
  geom_vline(xintercept = completeness_threshold_perc, color = "black", size = 0.2, linetype = "dashed") +
  geom_hline(data = hlines.df, aes(yintercept = y_intercept, color = what), size = 0.2, linetype = "dashed") +
  geom_line(size = 0.1) +
  geom_point(size = 0.1) +
  scale_color_manual(values = color_pal) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  xlab("Completeness threshold") +
  ylab("Nb. proteins above threshold") +
  theme_classic() +
  labs(color = "Completeness\nacross:")
pl
ggsave(paste0(saveDir, "/completeness_thresholding_nb_proteins.pdf"), plot = pl, width = 7, height = 5, useDingbats = F, limitsize = F)


saveDir <- "./Plots/proteomics/coon_QC/data_distribution"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

tt <- lapply(c("raw_data", "log2_data"), function(x){
  # x <- "raw_data"
  
  if(x == "raw_data"){
    dfPlot <- prot.df
  } else{
    dfPlot <- prot.df.log2
  }
  
  stopifnot(!is.unsorted(rev(unique(completeness.df$completeness_perc_samples))))
  dfPlot       <- type.convert(reshape2::melt(as.matrix(dfPlot)), as.is = T)
  dfPlot       <- dfPlot[!is.na(dfPlot$value), ]
  dfPlot       <- dfPlot[dfPlot$Var1 %in% completeness.df$prot_id[1:100], ]
  dfPlot$value <- log2(dfPlot$value + 1)
  pl <- ggplot(dfPlot, aes(x = value, color = Var1, group = Var1)) +
    geom_density(size = 0.2) +
    theme_classic() +
    theme(legend.position = "none")
  # pl
  ggsave(paste0(saveDir, "/distribution_top_100_complete_proteins__", x, ".pdf"), plot = pl, width = 4, height = 3.2, useDingbats = F, limitsize = F)
  
  
  pl <- ggplot(dfPlot, aes(x = value, color = Var2, group = Var2)) +
    geom_density(size = 0.1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    theme_classic() +
    theme(legend.position = "none")
  # pl
  ggsave(paste0(saveDir, "/distribution_samples__", x, ".pdf"), plot = pl, width = 4, height = 3.2, useDingbats = F, limitsize = F)
  
  NA
})



saveDir <- "./Plots/proteomics/coon_QC/PCA"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

tmpLogi.keep           <- completeness.df$completeness_perc_samples >= 50 & completeness.df$completeness_perc_strains_CD >= 50 & completeness.df$completeness_perc_strains_HFD >= 50
pca                    <- FactoMineR::PCA(t(prot.df[rownames(prot.df) %in% completeness.df$prot_id[tmpLogi.keep], ]), graph = F)
dfPlot                 <- cbind(pca$ind$coord, samples.meta.df[match(rownames(pca$ind$coord), samples.meta.df$sample), ])
pl <- ggplot(dfPlot, aes(x = Dim.1, y = Dim.2, color = batch, shape = diet)) +
  geom_point(size = 0.9) +
  xlab(paste0("Dim.1 (", round(pca$eig[1, "percentage of variance"], 2), "%)")) +
  ylab(paste0("Dim.2 (", round(pca$eig[2, "percentage of variance"], 2), "%)")) +
  theme_classic() +
  theme(legend.spacing.y = unit(0.1, 'cm'),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(0.1, "cm"),
        legend.margin    = margin(0.02,0,0,0, unit="cm")) +
  labs(color = "Batch", shape = "Diet")
# pl
ggsave(paste0(saveDir, "/pca_raw_data.pdf"), plot = pl, width = 4, height = 3.2, useDingbats = F, limitsize = F)

pca                    <- FactoMineR::PCA(t(prot.df.log2.filt), graph = F)
dfPlot                 <- cbind(pca$ind$coord, samples.meta.df[match(rownames(pca$ind$coord), samples.meta.df$sample), ])
pl <- ggplot(dfPlot, aes(x = Dim.1, y = Dim.2, color = batch, shape = diet)) +
  geom_point(size = 0.9) +
  xlab(paste0("Dim.1 (", round(pca$eig[1, "percentage of variance"], 2), "%)")) +
  ylab(paste0("Dim.2 (", round(pca$eig[2, "percentage of variance"], 2), "%)")) +
  theme_classic() +
  theme(legend.spacing.y = unit(0.1, 'cm'),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(0.1, "cm"),
        legend.margin    = margin(0.02,0,0,0, unit="cm")) +
  labs(color = "Batch", shape = "Diet")
# pl
ggsave(paste0(saveDir, "/pca_log2_data.pdf"), plot = pl, width = 4, height = 3.2, useDingbats = F, limitsize = F)



saveDir <- "./Plots/proteomics/coon_QC/median_inspection"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
dfPlot.1      <- type.convert(reshape2::melt(as.matrix(prot.df.log2)), as.is = T)
dfPlot.1      <- dfPlot.1[!is.na(dfPlot.1$value), ]
dfPlot.1$what <- "Log2 data"

dfPlot.2      <- type.convert(reshape2::melt(as.matrix(prot.df.log2)), as.is = T)
dfPlot.2      <- dfPlot.2[!is.na(dfPlot.2$value), ]
dfPlot.2$what <- "Log2 median-normalized data"

dfPlot <- rbind(dfPlot.1, dfPlot.2)

pl <- ggplot(dfPlot, aes(x = Var2, y = value)) +
  geom_boxplot(outlier.shape = NA, size = 0.2) +
  facet_wrap(~what, ncol = 1, scales = "free_y") +
  ylab("Intensity") +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 3.5),
        axis.title.x = element_blank())
# pl
ggsave(paste0(saveDir, "/samples_median__log2_data.pdf"), plot = pl, width = 16, height = 5, useDingbats = F, limitsize = F)

remove(pca, dfPlot, dfPlot.1, dfPlot.2, tt, saveDir, pl, completeness_threshold_perc)

###############
## Do DEA
##
## > Notes: 
## The batch effect (seems to be present from PCA) 
## can  not be distinguished from the strain effect because of the 
## bad experimental design, where all samples within one strain are in the same
## batch. But to partially overcome technical differences (batch differences in
## median intensity are very minimal), we account for the log2-median batch
## intensity (log2 because we also use log2 protein data) in the limma linear
## model. To account for the similarities between strains we use a limma block
## design.
##
## Gene-level analysis is useful for GSEA, to avoid duplicated rankings for
## gene(s) with multiple protein isoforms.
##
## We perform DEA on proteins expressed in at least 50% of samples and 50% of
## strains in CD and HFD.
##
## In limma, the data don't have to be normally distributed, the residuals yes.
##
## About duplicateCorrelation and blocking:
## https://support.bioconductor.org/p/125489/#125602
## https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/dupcor.html
###############

samples.meta.df$strain_diet_flat <- gsub("\\/", "_", samples.meta.df$strain_diet)
samples.meta.df$strain_flat      <- gsub("\\/", "_", samples.meta.df$strain)


contrast_types <- c("population_diet_effect", "strain_diet_effect")
tt             <- lapply(c("protein_level", "gene_level"), function(x){
  # x <- "protein_level"
  # x <- "gene_level"
  
  cat("\nPerforming DEA at", gsub("_", " ", x), "...\n")
  
  tt <- lapply(contrast_types, function(cc){
    # cc <- "population_diet_effect"
    cat("\tContrast type:", cc, "\n")
    
    
    if(x == "protein_level"){
      data.expr.use        <- prot.df.log2.filt
      completeness.use     <- completeness.df
      meta.features.use    <- prot.meta.df.aggr.byProt
      meta.mergeBy         <- "uniprot_id"
      completeness.mergeBy <- "prot_id"
    } else{
      data.expr.use        <- prot.df.log2.filt.aggr
      completeness.use     <- completeness.df.gene
      meta.features.use    <- prot.meta.df.aggr.byGene
      meta.mergeBy         <- "ensembl_gene_id"
      completeness.mergeBy <- "ensembl_gene_id"
    }
    completenes_col_map <- c("completeness_perc_samples"     = "perc_samples_detected",
                             "completeness_perc_samples_CD"  = "perc_CD_samples_detected",
                             "completeness_perc_samples_HFD" = "perc_HFD_samples_detected",
                             "completeness_perc_strains"     = "perc_strains_detected",
                             "completeness_perc_strains_CD"  = "perc_CD_strains_detected",
                             "completeness_perc_strains_HFD" = "perc_HFD_strains_detected")
    
    if(cc == "population_diet_effect"){
      
      # Which lipids respond to HFD on average?
      designFormula  <- as.formula(paste0("~0+diet+batch_median_log2")) # 0 means no intercept for the linear model
      mm             <- model.matrix(designFormula, data = samples.meta.df)
      colnames(mm)   <- gsub("diet", "", colnames(mm))
      contrasts_def  <- list("HFD_vs_CD_diet" = "(HFD)-(CD)")
      
    } else if(cc == "strain_diet_effect"){
      
      # How does strain X respond to diet?
      designFormula        <- as.formula(paste0("~0+strain_diet_flat+batch_median_log2"))
      mm                   <- model.matrix(designFormula, data = samples.meta.df)
      colnames(mm)         <- gsub("strain_diet_flat", "", colnames(mm))
      contrasts_def        <- lapply(unique(samples.meta.df$strain_flat), function(y) paste0("(", y, "_HFD)-(", y, "_CD)"))
      names(contrasts_def) <- paste0("HFD_vs_CD_strain_", unique(samples.meta.df$strain_flat))
      
    } 
    contrasts <- lapply(contrasts_def, function(k){
      # k <- contrasts_def[1]
      # k <- contrasts_def[grepl("BXD1$", names(contrasts_def))]
      eval(parse(text = paste0("out <- makeContrasts(", unname(k), ", levels = mm)")))
      colnames(out) <- names(k)
      out
    })
    
    dc                  <- limma::duplicateCorrelation(data.expr.use, design = mm, block = samples.meta.df$strain)
    fitModel            <- limma::lmFit(data.expr.use, design = mm, block = samples.meta.df$strain, correlation = dc$consensus.correlation)
    
    top_tables <- lapply(1:length(contrasts), function(k){
      # k <- 1
      fitModel.2.k            <- limma::eBayes(limma::contrasts.fit(fitModel, contrasts = contrasts[[k]]))
      table.k                 <- limma::topTable(fitModel.2.k, sort.by = "none", n = Inf)
      table.k[[meta.mergeBy]] <- rownames(table.k)
      cNames                  <- colnames(table.k)
      table.k                 <- merge(table.k, completeness.use, by.x = meta.mergeBy, by.y = completeness.mergeBy, all.x = T, all.y = F)
      for(i in names(completenes_col_map)){
        colnames(table.k)[colnames(table.k) == i] <- unname(completenes_col_map[[i]])
      }
      
      table.k                 <- merge(table.k, meta.features.use, by = meta.mergeBy, all.x = T, all.y = F)
      table.k                 <- dplyr::select(table.k, c(cNames, colnames(table.k)[!(colnames(table.k) %in% cNames)]))
      table.k$contrastID      <- contrasts_def[[k]]
      table.k$contrast        <- names(contrasts_def)[k]
      table.k                 <- table.k[order(table.k$adj.P.Val, table.k$P.Value), ]
      list(top_table     = table.k,
           fitted_bmodel = fitModel.2.k)
    })
    names(top_tables) <- names(contrasts)
    
    # to preserve compatibility with downstream scripts as before only global diet-level contrast was fitted
    if(cc == "population_diet_effect"){
      saveDir <- "./Data/proteomics/coon_DEA"
      if(!dir.exists(saveDir)){
        dir.create(saveDir, recursive = T)
      }
      topTable.df <- top_tables$HFD_vs_CD_diet$top_table
      saveRDS(topTable.df, paste0(saveDir, "/DEA_proteome_limma_table_all_contrasts__", x, ".RDS"))
      write.table(topTable.df, paste0(saveDir, "/DEA_proteome_limma_table_all_contrasts__", x, ".csv"), row.names = F, sep = ",")
    } else{
      saveDir <- "./Data/proteomics/coon_DEA/single_strains"
      if(!dir.exists(saveDir)){
        dir.create(saveDir, recursive = T)
      }
      topTable.df <- lapply(top_tables, function(w) w$top_table)
      saveRDS(topTable.df, paste0(saveDir, "/DEA_proteome_limma_table_all_contrasts__", x, "__contrast_type_", cc, ".RDS"))
    }
    
    
    ###############
    ## DEA fitting plots
    ###############
    
    # to preserve compatibility with downstream scripts as before only global diet-level contrast was fitted
    if(cc == "population_diet_effect"){
      saveDir <- "./Plots/proteomics/coon_DEA/MD_plots"
      if(!dir.exists(saveDir)){
        dir.create(saveDir, recursive = T)
      }
      pdf(paste0(saveDir, "/MD_plot__HFD_vs_CD__", x, ".pdf"), width = 5, height = 3.8, useDingbats = F)
      plotMD(top_tables$HFD_vs_CD_diet$fitted_bmodel, column = 1)
      dev.off()
    } else{
      saveDir <- paste0("./Plots/proteomics/coon_DEA/MD_plots/single_strains/", cc)
      if(!dir.exists(saveDir)){
        dir.create(saveDir, recursive = T)
      }
      tt <- lapply(names(top_tables), function(ww){
        # ww <- names(top_tables)[1]
        pdf(paste0(saveDir, "/MD_plot__", cc, "__", ww, ".pdf"), width = 5, height = 3.8, useDingbats = F)
        plotMD(top_tables[[ww]]$fitted_bmodel, column = 1)
        dev.off()
      })
      
    }
    
    
    # plotMD(fitModel, column = 3)
    
    NA
    
  })
})



################################################################################
## Genes set enrichment analysis - GSEA
################################################################################


top_table_files.prot        <- list.files("./Data/proteomics/coon_DEA/single_strains", full.names = T, pattern = "protein_level")
names(top_table_files.prot) <- gsub(".*contrast_type_|\\.RDS", "", top_table_files.prot, ignore.case = T)
top_table_files.prot        <- c("population_diet_effect" = "./Data/proteomics/coon_DEA/DEA_proteome_limma_table_all_contrasts__protein_level.RDS", top_table_files.prot)

top_table_files.gene        <- list.files("./Data/proteomics/coon_DEA/single_strains", full.names = T, pattern = "gene_level")
names(top_table_files.gene) <- gsub(".*contrast_type_|\\.RDS", "", top_table_files.gene, ignore.case = T)
top_table_files.gene        <- c("population_diet_effect" = "./Data/proteomics/coon_DEA/DEA_proteome_limma_table_all_contrasts__gene_level.RDS", top_table_files.gene)

top_table_files             <- list("protein_level" = top_table_files.prot,
                                    "gene_level"    = top_table_files.gene)
top_tables.list             <- lapply(top_table_files, function(x){
  # x <- top_table_files[[1]]
  out <- lapply(names(x), function(y){
    # y <- names(x)[1]
    if(y == "population_diet_effect"){
      list("HFD_vs_CD" = readRDS(x[[y]]))
    } else{
      # already a list
      readRDS(x[[y]])
    }
  })
  names(out) <- names(x)
  out
})

mitoSheets                   <- readxl::excel_sheets("./Data/input_data/mitochondrial_ressources_geneSets/Mouse.MitoCarta3.0.xls")
mouse.mito.genes             <- lapply(mitoSheets, function(x) readxl::read_excel("./Data/input_data/mitochondrial_ressources_geneSets/Mouse.MitoCarta3.0.xls", sheet = x))
names(mouse.mito.genes)      <- mitoSheets
mitocarta.pathways           <- mouse.mito.genes$`C MitoPathways`
colnames(mitocarta.pathways) <- c("gs_name", "gs_hierarchy", "gene_symbol")
mitocarta.pathways           <- tidyr::separate_rows(mitocarta.pathways, "gene_symbol", sep = ", ")
geneTable                    <- unique(dplyr::select(mouse.mito.genes$`B Mouse All Genes`, c("Symbol", "EnsemblGeneID")))
colnames(geneTable)          <- c("gene_symbol", "ensembl_gene_id")
geneTable                    <- tidyr::separate_rows(geneTable, c("ensembl_gene_id"), sep = "\\|")
geneTable                    <- unique(tidyr::separate_rows(geneTable, c("gene_symbol"), sep = ","))
mitocarta.pathways           <- tidyr::separate_rows(mitocarta.pathways, "gene_symbol", sep = ", ")
mitocarta.pathways           <- merge(mitocarta.pathways, geneTable, by = "gene_symbol")
mitocarta.pathways           <- unique(na.omit(mitocarta.pathways))
colnames(mitocarta.pathways) <- c("gene_symbol", "gs_name", "gs_hierarchy", "ensembl_gene")


msigdbr.geneSets <- list(BP           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP"),
                         CC           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC"),
                         MF           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "MF"),
                         REACTOME     = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"),
                         HPO          = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "HPO"),
                         HALLMARK     = msigdbr::msigdbr(species = "Mus musculus", category = "H"),
                         KEGG         = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG"),
                         WIKIPATHWAYS = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS"),
                         MITOCARTA    = mitocarta.pathways)
remove(geneTable, mitoSheets, mouse.mito.genes, mitocarta.pathways)

overwrite <- F
gsea.list <- lapply(c("logFC_ranking"), function(s){
  # s <- "logFC_ranking"
  
  lapply(names(top_tables.list), function(w){
    # w <- names(top_tables.list)[1]
    
    lapply(names(top_tables.list[[w]]), function(z){
      # z <- names(top_tables.list[[w]])[1]
      
      
      saveDir <- paste0("./Data/proteomics/coon_GSEA/", w, "/", s)
      if(!dir.exists(saveDir)){
        dir.create(saveDir, recursive = T)
      }
      out_file_rds <- paste0(saveDir, "/proteomics_GSEA_results___contrast_type_", z, ".RDS")
      if(file.exists(out_file_rds)){
        cat("\t\t--> Analysis already completed! Skipping...\n")
        return()
      }
      
      out.list <- lapply(names(top_tables.list[[w]][[z]]), function(k){
        # k <- names(top_tables.list[[w]][[z]])[1]
        
        cat("Running GSEA for:\n\t> contrast type:\t", z,
            "\n\t> contrast id:\t\t", k,
            "\n\t> DEA type:\t\t", w, 
            "\n\t> ranking type:\t\t", s, "...\n")
        
        tmp      <- top_tables.list[[w]][[z]][[k]]
        stopifnot(length(unique(tmp$contrastID)) == 1 & length(unique(tmp$contrast)) == 1)
        tmp      <- unique(dplyr::select(tmp, c("ensembl_gene_id", "uniprot_id", "logFC", "t", "adj.P.Val", "P.Value")))
        
        if(w == "gene_level"){
          
          tmp       <- tmp[!is.na(tmp$ensembl_gene_id) & !(tmp$ensembl_gene_id %in% c("", " ", "NA")), ]
          tmp       <- tmp[!is.na(tmp$P.Value) & !is.na(is.na(tmp$logFC)), ]
          # tt <- tmp[tmp$ensembl_gene_id %in% tmp$ensembl_gene_id[duplicated(tmp$ensembl_gene_id)], ]
          
          duplGenes <- unique(tmp$ensembl_gene_id[duplicated(tmp$ensembl_gene_id)])
          tmp       <- tmp[!(tmp$ensembl_gene_id %in% duplGenes), ]
          
        } else if(w == "protein_level"){
          
          tmp <- tmp[!is.na(tmp$uniprot_id) & !(tmp$uniprot_id %in% c("", " ", "NA")), ]
          tmp <- tidyr::separate_rows(tmp, "ensembl_gene_id", sep = "\\;")
          tmp <- tmp[!is.na(tmp$ensembl_gene_id) & !(tmp$ensembl_gene_id %in% c("", " ", "NA")), ]
          tmp <- tmp[!is.na(tmp$P.Value) & !is.na(is.na(tmp$logFC)), ]
          
          duplGenes <- unique(tmp$ensembl_gene_id[duplicated(tmp$ensembl_gene_id)])
          if(length(duplGenes) > 0){
            tmp.dupl  <- tmp[tmp$ensembl_gene_id %in% duplGenes, ]
            tmp.dupl  <- tmp.dupl[order(tmp.dupl$ensembl_gene_id, tmp.dupl$P.Value), ]
            # for each duplicated gene take:
            # (1a) if all not-significant and logFC discordant, discard the gene as the directionalit of non-significant genes is less interpretable
            # (1b) if not all the protein are significant, the significant ones
            # (2)  among the remaining ones, select the one with most significant change (optional, could first select the one without isoform tag or)
            tmp.dupl.filt <- plyr::ddply(tmp.dupl, c("ensembl_gene_id"), function(zz){
              # zz <- tmp.dupl[tmp.dupl$ensembl_gene_id == tmp.dupl$ensembl_gene_id[1], ]
              
              if(!all(zz$adj.P.Val >= 0.05)){
                zz <- zz[zz$adj.P.Val < 0.05, ]
              } else{
                # if all not-significant and logFC discordant, discard the gene as the directionality of non-significant genes is less interpretable
                if(length(unique(sign(zz$logFC))) > 1){
                  # cat(zz$ensembl_gene_id[1], " - discordant n.s. logFC direction\n") # n.s. = non significant
                  return()
                }
              }
              zz <- zz[which.min(zz$P.Value), ]
              zz
            })
            tmp.dupl.discarded <- tmp.dupl[!(tmp.dupl$ensembl_gene_id %in% tmp.dupl.filt$ensembl_gene_id), ]
            tmp <- rbind(tmp[!(tmp$ensembl_gene_id %in% duplGenes), ], tmp.dupl.filt) 
            stopifnot(sum(duplicated(tmp$ensembl_gene_id)) == 0)
          }
          
          
        }
        
        if(s == "logFC_ranking"){
          tmp             <- tmp[order(tmp$logFC, decreasing = T), ]
          geneList        <- tmp$logFC
          names(geneList) <- tmp$ensembl_gene_id
        } else{
          stop()
        }
        cols2select     <- c("gs_name", "ensembl_gene")
        
        
        out.list <- parallel::mclapply(mc.cores = length(names(msigdbr.geneSets)), X = names(msigdbr.geneSets), FUN = function(x){
          # x <- names(msigdbr.geneSets)[1]
          
          
          
          
          gsea <- clusterProfiler::GSEA(geneList      = geneList,
                                        pvalueCutoff  = 1,
                                        pAdjustMethod = "BH",
                                        TERM2GENE     = dplyr::select(msigdbr.geneSets[[x]], cols2select))
          gsea_df   <- gsea@result
          if(nrow(gsea_df) == 0){
            return(list(geneList = geneList,
                        gsea     = gsea,
                        gsea_df  = NULL))
          }
          colsSelect <- c("gs_name", "gs_exact_source", "gs_id", "gs_description")
          colsSelect <- colsSelect[colsSelect %in% colnames(msigdbr.geneSets[[x]])]
          termsMeta  <- unique(dplyr::select(msigdbr.geneSets[[x]], colsSelect))
          gsea_df    <- merge(gsea_df, termsMeta, by.x = "ID", by.y = "gs_name", all.x = T, all.y = F)
          gsea_df    <- gsea_df[order(gsea_df$qvalue), ]
          
          gsea_df$coreSize               <- unlist(lapply(gsea_df$core_enrichment, function(x) length(unlist(strsplit(x, "/")))))
          gsea_df$geneRatio              <- gsea_df$coreSize / gsea_df$setSize
          gsea_df$MSigDB_collection      <- x
          gsea_df$contrastID             <- k
          gsea_df$contrast_type          <- z
          gsea_df$dea_type               <- w
          gsea_df$ranking_type           <- s
          
          gsea_df
          
        })
        
        gsea_df.k <- as.data.frame(data.table::rbindlist(out.list, use.names = T, fill = T), stringsAsFactors = F)
        gsea_df.k <- type.convert(gsea_df.k, as.is = T)
        
        gsea_df.k$Description          <- gsub("^GOBP_|^GOCC_|^GOMF_|^HALLMARK_|^HP_|^KEGG_|^REACTOME_|^WP_|^BIOCARTA_|^PID_", "", gsea_df.k$Description)
        tmpLogi                        <- gsea_df.k$MSigDB_collection %in% c("CELLTYPE", "CGP", "IMMUNESIGDB")
        gsea_df.k$Description[tmpLogi] <- gsub("^.+?_(.*)", "\\1", gsea_df.k$Description[tmpLogi])
        gsea_df.k$Description          <- gsub("_", " ", gsea_df.k$Description)
        
        gsea_df.k <- gsea_df.k[order(gsea_df.k$qvalue, gsea_df.k$pvalue, gsea_df.k$Description), ]
        
        # NA
        gsea_df.k
      })
      names(out.list) <- names(top_tables.list[[w]][[z]])
      saveRDS(out.list, out_file_rds)
    })
  })
})



