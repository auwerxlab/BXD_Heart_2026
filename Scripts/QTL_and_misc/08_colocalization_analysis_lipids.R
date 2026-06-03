
library(data.table)


## (0) define analysis parameters
cis_thresh_Mb          <- 2
ci_overlap_thresh_perc <- 50
max_peaks_dist_Mb      <- 2


# gene table
geneTable                           <- readRDS("./Data/RNAseq_processing/geneConversionTables/geneConversionTable_GRCm38_release-102_complete_gtf.RDS")
geneTable                           <- type.convert(as.data.frame(geneTable), as.is = T)
geneTable                           <- geneTable[geneTable$seqnames %in% c(1:20, "X", "Y", "MT"), ]
geneTable                           <- dplyr::select(geneTable, dplyr::all_of(c("gene_id", "gene_name", "seqnames", "start", "end", "strand", "type")))
colnames(geneTable)                 <- c("gene_id_ensembl", "gene_symbol", "chromosome", "chromosome_start_position", "chromosome_end_position", "strand", "type")
geneTable                           <- geneTable[geneTable$type == "gene", ]
geneTable                           <- unique(geneTable)
geneTable                           <- dplyr::select(geneTable, -dplyr::all_of(c("strand", "type")))
geneTable$chromosome_start_position <- as.numeric(geneTable$chromosome_start_position)
geneTable$chromosome_end_position   <- as.numeric(geneTable$chromosome_end_position)
geneTable$posStart_Mb               <- geneTable$chromosome_start_position / 1e6
geneTable$posEnd_Mb                 <- geneTable$chromosome_end_position / 1e6
stopifnot(sum(duplicated(geneTable$gene_id_ensembl)) == 0)


## (2) load lQTL data 
lqtl_peaks.df <- lapply(c("avg", "ind"), function(x){
  # x <- "avg"
  if(x == "avg"){
    files.x        <- list.files("./Data/QTL_mapping/qtl_outputs_formatted", recursive = T, pattern = "average_lipidomic_derived_avgAggr", full.names = T)
    files.x        <- files.x[grepl("qtl_peaks_table", files.x) & grepl("RDS$", files.x)]
    names(files.x) <- gsub(".*qtl_peaks_table_average_lipidomic_derived_avgAggr_|\\.RDS$", "", files.x)
    names(files.x)[names(files.x) == "data"] <- "avg_data"
  } else if(x == "ind"){
    files.x        <- list.files("./Data/QTL_mapping/qtl_outputs_formatted", recursive = T, pattern = "individuals_lipidomic_Filt_MassNorm_BatchNorm_Data_NArm", full.names = T)
    files.x        <- files.x[grepl("qtl_peaks_table", files.x) & grepl("RDS$", files.x)]
    names(files.x) <- gsub(".*qtl_peaks_table_individuals_lipidomic_derived_|.*qtl_peaks_table_individuals_lipidomic_|\\.RDS$", "", files.x)
    names(files.x)[names(files.x) == "Filt_MassNorm_BatchNorm_Data_NArm"] <- "ind_data"
  }
  lqtl.df.x                  <- do.call(rbind, lapply(files.x, readRDS))
  lqtl.df.x$aggregation_type <- ifelse(x == "avg", "average", "individuals")
  lqtl.df.x                  <- lqtl.df.x[!is.na(lqtl.df.x$pos_Mb), ]
  lqtl.df.x
})
lqtl_peaks.df             <- do.call(rbind, lqtl_peaks.df)
lqtl_peaks.df             <- lqtl_peaks.df[!is.na(lqtl_peaks.df$ci_lo_Mb) & !is.na(lqtl_peaks.df$ci_hi_CM) & !is.na(lqtl_peaks.df$chr), ]
lqtl_peaks.df$analysis_id <- paste0(lqtl_peaks.df$condition, "__", lqtl_peaks.df$interaction, "__", lqtl_peaks.df$aggregation_type)
lqtl_peaks.df$qtl_type    <- "lQTL"
lqtl_peaks.df$pheno_id    <- lqtl_peaks.df$id


lqtl_peaks.species.df <- lqtl_peaks.df[lqtl_peaks.df$inputType %in% c("average_lipidomic_derived_avgAggr_data", "individuals_lipidomic_Filt_MassNorm_BatchNorm_Data_NArm"), ]
lqtl_peaks.derived.df <- lqtl_peaks.df[grepl("lipid_class_avg|lipid_class_ratio|lipid_indexes|lipid_PCA|lipid_class_fraction|lipid_class_sum|lipid_chain_length_class_fraction|lipid_nb_chains_class_fraction|lipid_nb_double_bonds_class_fraction", lqtl_peaks.df$inputType), ]
lqtl_peaks.other.df   <- lqtl_peaks.df[!(lqtl_peaks.df$inputType %in% c("average_lipidomic_derived_avgAggr_data", "individuals_lipidomic_Filt_MassNorm_BatchNorm_Data_NArm")), ]


## (3) load eQTL data 
eqtl_peaks.df              <- readRDS("./Data/QTL_mapping/qtl_outputs_formatted/RNAseq_cf_pca/qtl_peaks_table_RNAseq_cf_pca.RDS")
eqtl_peaks.df              <- eqtl_peaks.df[!is.na(eqtl_peaks.df$gene_posStart_Mb), ]
eqtl_peaks.df$isCisQTL     <- (((abs(eqtl_peaks.df$gene_posStart_Mb - eqtl_peaks.df$pos_Mb) < cis_thresh_Mb) | (abs(eqtl_peaks.df$gene_posEnd_Mb - eqtl_peaks.df$pos_Mb) < cis_thresh_Mb) | (eqtl_peaks.df$gene_posStart_Mb <= eqtl_peaks.df$pos_Mb & eqtl_peaks.df$gene_posEnd_Mb >= eqtl_peaks.df$pos_Mb)) & (eqtl_peaks.df$chr == eqtl_peaks.df$chr_gene))
eqtl_peaks.df              <- eqtl_peaks.df[!is.na(eqtl_peaks.df$ci_lo_Mb) & !is.na(eqtl_peaks.df$ci_hi_CM) & !is.na(eqtl_peaks.df$chr), ]
eqtl_peaks.df$analysis_id  <- paste0(eqtl_peaks.df$condition, "__", eqtl_peaks.df$interaction)
eqtl_peaks.df$qtl_type     <- ifelse(eqtl_peaks.df$isCisQTL, "cis_eQTL", "trans_eQTL")
eqtl_peaks.df$pheno_id     <- eqtl_peaks.df$id
stopifnot(sum(grepl(";", eqtl_peaks.df$id)) == 0 & sum(grepl(";", eqtl_peaks.df$gene_symbol)) == 0 & sum(grepl("ENS", eqtl_peaks.df$gene_symbol)) == 0)
sum(is.na(eqtl_peaks.df$gene_symbol) | eqtl_peaks.df$gene_symbol == "")

## (4) load pQTL data 
pqtl_peaks.df <- lapply(c("avg", "ind"), function(x){
  # x <- "avg"
  
    pqtl_peaks_file.x <- ifelse(x == "avg",
                              "./Data/QTL_mapping/qtl_outputs_formatted/average_coonProteomic_cf_pca/qtl_peaks_table_average_coonProteomic_cf_pca.RDS",
                              "./Data/QTL_mapping/qtl_outputs_formatted/individuals_coonProteomic_cf_pca/qtl_peaks_table_individuals_coonProteomic_cf_pca.RDS")
    
    pqtl.df.x                  <- readRDS(pqtl_peaks_file.x)
    pqtl.df.x                  <- pqtl.df.x[!is.na(pqtl.df.x$gene_posStart_Mb), ]
    pqtl.df.x$isCisQTL         <- (((abs(pqtl.df.x$gene_posStart_Mb - pqtl.df.x$pos_Mb) < cis_thresh_Mb) | (abs(pqtl.df.x$gene_posEnd_Mb - pqtl.df.x$pos_Mb) < cis_thresh_Mb) | (pqtl.df.x$gene_posStart_Mb <= pqtl.df.x$pos_Mb & pqtl.df.x$gene_posEnd_Mb >= pqtl.df.x$pos_Mb)) & (pqtl.df.x$chr == pqtl.df.x$chr_gene))
    pqtl.df.x$aggregation_type <- ifelse(x == "avg", "average", "individuals")
    pqtl.df.x
})
pqtl_peaks.df              <- do.call(rbind, pqtl_peaks.df)
pqtl_peaks.df              <- pqtl_peaks.df[!is.na(pqtl_peaks.df$ci_lo_Mb) & !is.na(pqtl_peaks.df$ci_hi_CM) & !is.na(pqtl_peaks.df$chr), ]
pqtl_peaks.df$analysis_id  <- paste0(pqtl_peaks.df$condition, "__", pqtl_peaks.df$interaction, "__", pqtl_peaks.df$aggregation_type)
pqtl_peaks.df$qtl_type     <- ifelse(pqtl_peaks.df$isCisQTL, "cis_pQTL", "trans_pQTL")
pqtl_peaks.df$pheno_id     <- pqtl_peaks.df$gene_id_ensembl
tmp.df <- unique(dplyr::select(eqtl_peaks.df, dplyr::all_of(c("id", "gene_symbol"))))
pqtl_peaks.df$gene_symbol  <- ifelse(pqtl_peaks.df$gene_id_ensembl %in% tmp.df$id,
                                     plyr::mapvalues(pqtl_peaks.df$gene_id_ensembl, from = tmp.df$id, to = tmp.df$gene_symbol, warn_missing = F),
                                     pqtl_peaks.df$gene_symbol)
pqtl_peaks.df$gene_symbol  <- ifelse(!(pqtl_peaks.df$gene_id_ensembl %in% tmp.df$id) & pqtl_peaks.df$gene_id_ensembl %in% geneTable$gene_id_ensembl,
                                     plyr::mapvalues(pqtl_peaks.df$gene_id_ensembl, from = geneTable$gene_id_ensembl, to = geneTable$gene_symbol, warn_missing = F),
                                     pqtl_peaks.df$gene_symbol)


stopifnot(sum(grepl(";", pqtl_peaks.df$gene_id_ensembl)) == 0 & sum(grepl(";", pqtl_peaks.df$gene_symbol)) == 0 & sum(grepl("ENS", pqtl_peaks.df$gene_symbol)) == 0)
sum(is.na(pqtl_peaks.df$gene_symbol) | pqtl_peaks.df$gene_symbol == "")

## (5) load pheQTL data 
pheno_peaks.df <- lapply(c("avg", "ind"), function(x){
  # x <- "avg"
  print(x)
  
  if(x == "avg"){
    peaks.pheno.df.fed             <- rbind(readRDS("./Data/QTL_mapping/qtl_outputs_formatted/average_pheno_mean/qtl_peaks_table_average_pheno_mean.RDS"),
                                            readRDS("./Data/QTL_mapping/qtl_outputs_formatted/average_pheno_derived_traits_1_mean/qtl_peaks_table_average_pheno_derived_traits_1_mean.RDS")
                                            # readRDS("./Data/QTL_mapping/qtl_outputs_formatted/average_pheno_meanRobust/qtl_peaks_table_average_pheno_meanRobust.RDS"),
                                            # readRDS("./Data/QTL_mapping/qtl_outputs_formatted/average_pheno_derived_traits_1_meanRobust/qtl_peaks_table_average_pheno_derived_traits_1_meanRobust.RDS")
                                            )
    peaks.pheno.df.fed$id          <- paste0(peaks.pheno.df.fed$id, " (", ifelse(peaks.pheno.df.fed$inputType == "average_pheno_mean", "avg", "avgR"), ")")
    
    
    peaks.pheno.df.fasted          <- rbind(readRDS("./Data/QTL_mapping/qtl_outputs_formatted/average_fasted_pheno_mean_commonStrains/qtl_peaks_table_average_fasted_pheno_mean_commonStrains.RDS"),
                                            readRDS("./Data/QTL_mapping/qtl_outputs_formatted/average_fasted_pheno_derived_traits_1_mean_commonStrains/qtl_peaks_table_average_fasted_pheno_derived_traits_1_mean_commonStrains.RDS")
                                            # readRDS("./Data/QTL_mapping/qtl_outputs_formatted/average_fasted_pheno_meanRobust_commonStrains/qtl_peaks_table_average_fasted_pheno_meanRobust_commonStrains.RDS"),
                                            # readRDS("./Data/QTL_mapping/qtl_outputs_formatted/average_fasted_pheno_derived_traits_1_meanRobust_commonStrains/qtl_peaks_table_average_fasted_pheno_derived_traits_1_meanRobust_commonStrains.RDS")
                                            )
    peaks.pheno.df.fasted$id       <- paste0(peaks.pheno.df.fasted$id, " (fasted_", ifelse(peaks.pheno.df.fasted$inputType == "average_pheno_mean", "avg", "avgR"), ")")
    peaks.pheno.df                 <- rbind(peaks.pheno.df.fed, peaks.pheno.df.fasted)
    
  } else if(x == "ind"){
    peaks.pheno.df.fed             <- rbind(readRDS("./Data/QTL_mapping/qtl_outputs_formatted/individuals_pheno/qtl_peaks_table_individuals_pheno.RDS"),
                                            readRDS("./Data/QTL_mapping/qtl_outputs_formatted/individuals_pheno_derived_traits_1/qtl_peaks_table_individuals_pheno_derived_traits_1.RDS"))
    peaks.pheno.df.fasted          <- rbind(readRDS("./Data/QTL_mapping/qtl_outputs_formatted/individuals_fasted_pheno_commonStrains/qtl_peaks_table_individuals_fasted_pheno_commonStrains.RDS"),
                                            readRDS("./Data/QTL_mapping/qtl_outputs_formatted/individuals_fasted_pheno_derived_traits_1_commonStrains/qtl_peaks_table_individuals_fasted_pheno_derived_traits_1_commonStrains.RDS"))
    peaks.pheno.df.fasted$id       <- paste0(peaks.pheno.df.fasted$id, "(fasted)")
    peaks.pheno.df                 <- rbind(peaks.pheno.df.fed, peaks.pheno.df.fasted)
  }
  
  peaks.pheno.df$aggregation_type <- ifelse(x == "avg", "average", "individuals")
  if(any(grepl("MV_E/A", peaks.pheno.df$id, fixed = T)) & any(grepl("MV_E/A_correct", peaks.pheno.df$id, fixed = T))){
    peaks.pheno.df <- peaks.pheno.df[!(grepl("MV_E/A", peaks.pheno.df$id, fixed = T) & !grepl("MV_E/A_correct", peaks.pheno.df$id, fixed = T)), ]
  }
  # peaks.pheno.df                 <- peaks.pheno.df[!grepl("caecum|week|^tibia|fc_|eWAT|Atria|Blood|BW|Gastro|Heart_mass|Intestine|Kidney|Liver|LV_mass|RV_mass|scWAT|Soleus|Spleen", peaks.pheno.df$id, ignore.case = T), ]
  print(sort(unique(peaks.pheno.df$id)))
  peaks.pheno.df
})
pheno_peaks.df              <- do.call(rbind, pheno_peaks.df)
pheno_peaks.df              <- pheno_peaks.df[!is.na(pheno_peaks.df$ci_lo_Mb) & !is.na(pheno_peaks.df$ci_hi_CM) & !is.na(pheno_peaks.df$chr), ]
pheno_peaks.df$analysis_id  <- paste0(pheno_peaks.df$condition, "__", pheno_peaks.df$interaction, "__", pheno_peaks.df$aggregation_type)
pheno_peaks.df$qtl_type     <- "pheQTL"
pheno_peaks.df$pheno_id     <- pheno_peaks.df$id




## mRNA DEA results - only significant genes
mRNA_DEA_results.sig                <- readRDS("./Data/RNAseq_processing/limma_DEA/limma_outputs/DEA_transcriptome_limma_table_all_contrasts.RDS")
mRNA_DEA_results.sig                <- mRNA_DEA_results.sig[mRNA_DEA_results.sig$contrastID == "HFD_vs_CD_diet", ]
mRNA_DEA_results.sig                <- mRNA_DEA_results.sig[mRNA_DEA_results.sig$adj.P.Val < 0.05 & !is.na(mRNA_DEA_results.sig$adj.P.Val), ]
mRNA_DEA_results.sig                <- dplyr::select(mRNA_DEA_results.sig,  c("gene_id", "logFC", "adj.P.Val"))
mRNA_DEA_results.sig                <- mRNA_DEA_results.sig[!is.na(mRNA_DEA_results.sig$gene_id) & !(mRNA_DEA_results.sig$gene_id == ""), ]
colnames(mRNA_DEA_results.sig)      <- c("gene_ensembl_id", "logFC", "adj.P.Val")
mRNA_DEA_results.sig$mRNA_gene_isDE <- ifelse(mRNA_DEA_results.sig$logFC > 0, "Up", "Down")
mRNA_DEA_results.sig$mRNA_gene_isDE <- paste0(mRNA_DEA_results.sig$mRNA_gene_isDE, " (", round(mRNA_DEA_results.sig$logFC, 2), "__", formatC(mRNA_DEA_results.sig$adj.P.Val, format = "e", digits = 2), ")")
mRNA_DEA_results.sig                <- dplyr::select(mRNA_DEA_results.sig, -dplyr::all_of(c("logFC", "adj.P.Val")))
stopifnot(all(!duplicated(mRNA_DEA_results.sig$gene_ensembl_id)))
# head(mRNA_DEA_results)

## mRNA DEA results - all genes
mRNA_DEA_results           <- readRDS("./Data/RNAseq_processing/limma_DEA/limma_outputs/DEA_transcriptome_limma_table_all_contrasts.RDS")
mRNA_DEA_results           <- mRNA_DEA_results[mRNA_DEA_results$contrastID == "HFD_vs_CD_diet", ]
mRNA_DEA_results           <- dplyr::select(mRNA_DEA_results,  c("gene_id", "logFC", "P.Value", "adj.P.Val"))
mRNA_DEA_results           <- mRNA_DEA_results[!is.na(mRNA_DEA_results$gene_id) & !(mRNA_DEA_results$gene_id == ""), ]
colnames(mRNA_DEA_results) <- c("gene_id", "DEA_logFC_mRNA", "DEA_pval_mRNA", "DEA_adjpval_mRNA")


## proteins DEA results - Coon data - only significant proteins
prot_coon_DEA_results.sig                <- readRDS("./Data/proteomics/coon_DEA/DEA_proteome_limma_table_all_contrasts__protein_level.RDS")
prot_coon_DEA_results.sig                <- prot_coon_DEA_results.sig[prot_coon_DEA_results.sig$contrastID == "HFD_vs_CD_diet", ]
prot_coon_DEA_results.sig                <- prot_coon_DEA_results.sig[prot_coon_DEA_results.sig$adj.P.Val < 0.05 & !is.na(prot_coon_DEA_results.sig$adj.P.Val), ]
prot_coon_DEA_results.sig                <- dplyr::select(prot_coon_DEA_results.sig,  c("ensembl_gene_id", "logFC", "adj.P.Val"))
prot_coon_DEA_results.sig                <- prot_coon_DEA_results.sig[!is.na(prot_coon_DEA_results.sig$ensembl_gene_id) & prot_coon_DEA_results.sig$ensembl_gene_id != "", ]
colnames(prot_coon_DEA_results.sig)      <- c("gene_ensembl_id", "logFC", "adj.P.Val")
prot_coon_DEA_results.sig                <- tidyr::separate_rows(prot_coon_DEA_results.sig, c("gene_ensembl_id"), sep = ";")
prot_coon_DEA_results.sig$prot_gene_isDE <- ifelse(prot_coon_DEA_results.sig$logFC > 0, "Up", "Down")
prot_coon_DEA_results.sig$prot_gene_isDE <- paste0(prot_coon_DEA_results.sig$prot_gene_isDE, " (", round(prot_coon_DEA_results.sig$logFC, 2), "__", formatC(prot_coon_DEA_results.sig$adj.P.Val, format = "e", digits = 2), ")")
prot_coon_DEA_results.aggrByGene.sig     <- data.table(prot_coon_DEA_results.sig)[, list(prot_isDE = paste(prot_gene_isDE, collapse = ";")), by = c("gene_ensembl_id")]
stopifnot(all(!duplicated(prot_coon_DEA_results.aggrByGene.sig$gene_ensembl_id)))
# head(prot_coon_DEA_results.aggrByGene.sig)


## proteins DEA results - Coon data - all proteins
prot_coon_DEA_results           <- readRDS("./Data/proteomics/coon_DEA/DEA_proteome_limma_table_all_contrasts__protein_level.RDS")
prot_coon_DEA_results           <- prot_coon_DEA_results[prot_coon_DEA_results$contrastID == "HFD_vs_CD_diet", ]
prot_coon_DEA_results           <- dplyr::select(prot_coon_DEA_results,  c("uniprot_id", "ensembl_gene_id", "logFC", "P.Value", "adj.P.Val"))
prot_coon_DEA_results           <- prot_coon_DEA_results[!is.na(prot_coon_DEA_results$ensembl_gene_id) & prot_coon_DEA_results$ensembl_gene_id != "", ]
prot_coon_DEA_results           <- tidyr::separate_rows(prot_coon_DEA_results, c("ensembl_gene_id"), sep = ";")
colnames(prot_coon_DEA_results) <- c("uniprot_id", "gene_id", "DEA_logFC_prot", "DEA_pval_prot", "DEA_adjpval_prot")

# for protein isoforms take the most significant one
prot_coon_DEA_results           <- plyr::ddply(prot_coon_DEA_results, c("gene_id"), function(x){
  x[which.min(x$DEA_adjpval_prot), ]
})

stopifnot(sum(duplicated(mRNA_DEA_results$gene_id)) == 0)
stopifnot(sum(duplicated(prot_coon_DEA_results$gene_id)) == 0)


## merge gene annotation tables
geneAnnot.df <- list(mRNA_DEA_results.sig,
                     prot_coon_DEA_results.aggrByGene.sig)
stopifnot(all(unlist(lapply(geneAnnot.df, function(x) "gene_ensembl_id" %in% colnames(x)))))
geneAnnot.df <- type.convert(Reduce(function(...) merge(..., all = T, by = c("gene_ensembl_id")), geneAnnot.df), as.is = T)
geneAnnot.df <- geneAnnot.df[!is.na(geneAnnot.df$gene_ensembl_id), ]
remove(mRNA_DEA_results.sig, prot_coon_DEA_results.aggrByGene.sig)


############
## colocalization analysis
############

testRun           <- F
overwrite         <- T
nbCores           <- 14
run_individuals   <- F
run_CD_HFD_noInt  <- F
conditions.run    <- unique(lqtl_peaks.species.df$analysis_id)
if(!run_individuals){
  conditions.run <- conditions.run[!grepl("individuals", conditions.run)]
}
if(!run_CD_HFD_noInt){
  conditions.run <- conditions.run[!grepl("CD_HFD__noInteraction", conditions.run)]
}
saveDir <- "./Data/gene_prioritization/prioritization_output_lipids"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
overlap.df      <- lapply(conditions.run[1], function(xx){
  # xx <- conditions.run[1]
  # xx <- "CD_HFD__dietInteraction__average"
  
  cat("Colocalization running for:", xx, "...\n")
  
  out_files <- c(paste0(saveDir, "/lipids_prioritized___", xx, "___", format(Sys.time(), "%d_%m_%y"), ".RDS"),
                 paste0(saveDir, "/lipids_all_genes_criteria___", xx, "___", format(Sys.time(), "%d_%m_%y"), ".RDS"))
  if(all(file.exists(out_files)) & !overwrite){
    cat("\tResults already saved... skipping!\n")
    return()
  }
  
  df.lqtl.species <- lqtl_peaks.species.df[lqtl_peaks.species.df$analysis_id == xx, ]
  df.lqtl.others  <- lqtl_peaks.other.df[lqtl_peaks.other.df$analysis_id == xx, ]
  df.lqtl.others  <- df.lqtl.others[!(df.lqtl.others$id %in% df.lqtl.species$id), ]
  df.eqtl         <- eqtl_peaks.df[eqtl_peaks.df$analysis_id == gsub("__([^_]*?)$", "", xx), ]
  df.pqtl         <- pqtl_peaks.df[pqtl_peaks.df$analysis_id == xx, ]
  df.pheqtl       <- pheno_peaks.df[pheno_peaks.df$analysis_id == xx, ]
  df.search       <- as.data.frame(rbindlist(list(df.lqtl.others, df.eqtl, df.pqtl, df.pheqtl), use.names = T, fill = T))
  idx.query       <- which(!grepl("unknown", df.lqtl.species$id, ignore.case = T))
  
  if(testRun){
    idx.query <- idx.query[1:min(c(20, length(idx.query)))]
  }                       
  cat("\t--> testing", length(idx.query), "lipid species QTL peaks...\n")
  
  
  overlap.df.xx <- pbmcapply::pbmclapply(mc.cores = nbCores, X = idx.query, FUN = function(kk){
    # kk <- which(df.lqtl.species$id == "PC 42:10_RT_8.033"  & grepl("63.133", df.lqtl.species$pos_Mb))
    # kk <- 171
    # kk <- 6

    return_error_warning_if_generated <- function(expr) {
      warning_message <- NULL
      
      result <- withCallingHandlers(
        tryCatch(expr, error = function(e) paste0(paste0("Error at idx ", kk, " -- ", paste(e, collapse = "___")))),
        warning = function(w) {
          warning_message <<- paste0(paste0("Warning at idx ", kk, " -- ", paste(w, collapse = "___")))
          invokeRestart("muffleWarning")  # suppress further handling
        }
      )
      
      if (!is.null(warning_message)) {
        return(warning_message)
      }
      
      return(result)
    }
    

    return_error_warning_if_generated({
      
      dt.search.kk <- as.data.table(rbindlist(list(df.search, df.lqtl.species[-kk, ]), use.names = T, fill = T))
      dt.ref.kk    <- as.data.table(df.lqtl.species[kk, ])
      
      setkey(dt.search.kk, chr, ci_lo_Mb, ci_hi_Mb)
      setkey(dt.ref.kk,    chr, ci_lo_Mb, ci_hi_Mb)
      
      overlap.dt.kk <- as.data.frame(data.table::foverlaps(dt.search.kk, dt.ref.kk,
                                                           by.x = c("chr", "ci_lo_Mb", "ci_hi_Mb"),
                                                           by.y = c("chr", "ci_lo_Mb", "ci_hi_Mb"),
                                                           type = "any",
                                                           mult = "all"))
      overlap.dt.kk <- overlap.dt.kk[!is.na(overlap.dt.kk$id) & overlap.dt.kk$i.id != dt.ref.kk$id[1], ]
      

      # get all genes under current reference QTL peak
      genes.underpk                  <- geneTable[(geneTable$posStart_Mb <= dt.ref.kk$ci_hi_Mb & geneTable$posEnd_Mb >= dt.ref.kk$ci_lo_Mb & geneTable$chromosome == dt.ref.kk$chr) | geneTable$gene_id_ensembl %in% overlap.dt.kk$i.pheno_id, ]
      rownames(genes.underpk)        <- NULL
      genes.underpk$is_under_ref_qtl <- genes.underpk$posStart_Mb <= dt.ref.kk$ci_hi_Mb & genes.underpk$posEnd_Mb >= dt.ref.kk$ci_lo_Mb & genes.underpk$chromosome == dt.ref.kk$chr
      genes.underpk                 <- dplyr::select(genes.underpk, -dplyr::all_of(c("chromosome_start_position", "chromosome_end_position")))
      colnames(genes.underpk)       <- c("gene_id", "gene_symbol", "chr", "gene_posStart_Mb", "gene_posEnd_Mb", "is_under_ref_qtl")
      for(i in unique(c("id", "chr", "pos_Mb", "ci_lo_Mb", "ci_hi_Mb", "LOD_threshold", "LOD", "pVal", "peak_type", "inputType", "aggregation_type"))){
        genes.underpk[[paste0("ref_qtl_", i)]] <- dt.ref.kk[[i]][1]
      }
      for(i in c("condition", "interaction", "inputType", "analysis_id")){
        genes.underpk[[i]] <- dt.ref.kk[[i]][1]
      }
      genes.underpk$gene_start_dist_from_pk_Mb <- ifelse(genes.underpk$chr == genes.underpk$ref_qtl_chr, abs(genes.underpk$gene_posStart_Mb - genes.underpk$ref_qtl_pos_Mb), NA) 
      
      
      # filter by QTL peaks distance and QTL peak CI overlap size
      # Note: pmax(overlap.dt.kk$ci_lo_Mb, overlap.dt.kk$i.ci_lo_Mb) computes the start position of the overlap
      #       pmin(overlap.dt.kk$ci_hi_Mb, overlap.dt.kk$i.hi_lo_Mb) computes the end position of the overlap
      overlap.dt.kk$peak_dist_Mb          <- overlap.dt.kk$pos_Mb - overlap.dt.kk$i.pos_Mb
      overlap.dt.kk$overlap_size          <- pmin(overlap.dt.kk$ci_hi_Mb, overlap.dt.kk$i.ci_hi_Mb) - pmax(overlap.dt.kk$ci_lo_Mb, overlap.dt.kk$i.ci_lo_Mb)
      overlap.dt.kk$perc_ci_overlap_ref   <- (overlap.dt.kk$overlap_size / (overlap.dt.kk$ci_hi_Mb - overlap.dt.kk$ci_lo_Mb)) * 100
      overlap.dt.kk$perc_ci_overlap_match <- (overlap.dt.kk$overlap_size / (overlap.dt.kk$i.ci_hi_Mb - overlap.dt.kk$i.ci_lo_Mb)) * 100
      # TODO: for some reason, we had a case where the ci_lo, ci_hi was the same (in Mb, not the same in CM - returning NA as division by 0 (0 interval length))
      overlap.dt.kk$perc_ci_overlap_match[is.na(overlap.dt.kk$perc_ci_overlap_match)] <- 100
      overlap.dt.kk$perc_ci_overlap_ref[is.na(overlap.dt.kk$perc_ci_overlap_ref)]     <- 100

      # !!!! Only keep QTL peaks that (1) have overlapping QTL peaks CI and (2) the peaks are within 2Mb (max_peaks_dist_Mb)
      overlap.dt.kk <- overlap.dt.kk[abs(overlap.dt.kk$peak_dist_Mb) <= max_peaks_dist_Mb, ]
      overlap.dt.kk <- overlap.dt.kk[overlap.dt.kk$perc_ci_overlap_ref > ci_overlap_thresh_perc | overlap.dt.kk$perc_ci_overlap_match > ci_overlap_thresh_perc, ]
      
      # check if there is at least an overlap with a (e/p)gene
      has_egene_overlap <- any(grepl("ENSMUSG", overlap.dt.kk$i.pheno_id) & !is.na(overlap.dt.kk$isCisQTL) & overlap.dt.kk$isCisQTL)
      # has_egene_overlap <- any(grepl("ENSMUSG", overlap.dt.kk$i.id) & !is.na(overlap.dt.kk$isCisQTL) & overlap.dt.kk$isCisQTL)
      if(!has_egene_overlap){
        return()
      }
      
      # update genes.underpk table (add (e/p)QTL overlap)
      # remove genes not under pheQL peak with no overlapping QTL peak
      genes.underpk                            <- genes.underpk[genes.underpk$is_under_ref_qtl | genes.underpk$gene_id %in% overlap.dt.kk$i.pheno_id, ]
      genes.underpk$cis_eQTL_overlap           <- unlist(lapply(genes.underpk$gene_id, function(zz) any(overlap.dt.kk$i.pheno_id == zz & overlap.dt.kk$i.qtl_type == "cis_eQTL")))
      genes.underpk$cis_pQTL_overlap           <- unlist(lapply(genes.underpk$gene_id, function(zz) any(overlap.dt.kk$i.pheno_id == zz & overlap.dt.kk$i.qtl_type == "cis_pQTL")))
      genes.underpk$cis_pQTL_overlap_feature   <- NA
      genes.underpk$cis_pQTL_overlap_feature[genes.underpk$cis_pQTL_overlap] <- unlist(lapply(genes.underpk$gene_id[genes.underpk$cis_pQTL_overlap], function(zz) paste(sort(unique(overlap.dt.kk$i.id[overlap.dt.kk$i.pheno_id == zz & overlap.dt.kk$i.qtl_type == "cis_pQTL"])), collapse = ";")))
      genes.underpk$trans_eQTL_overlap         <- unlist(lapply(genes.underpk$gene_id, function(zz) any(overlap.dt.kk$i.pheno_id == zz & overlap.dt.kk$i.qtl_type == "trans_eQTL")))
      genes.underpk$trans_pQTL_overlap         <- unlist(lapply(genes.underpk$gene_id, function(zz) any(overlap.dt.kk$i.pheno_id == zz & overlap.dt.kk$i.qtl_type == "trans_pQTL")))
      genes.underpk$trans_pQTL_overlap_feature <- NA
      genes.underpk$trans_pQTL_overlap_feature[genes.underpk$trans_pQTL_overlap] <- unlist(lapply(genes.underpk$gene_id[genes.underpk$trans_pQTL_overlap], function(zz) paste(sort(unique(overlap.dt.kk$i.id[overlap.dt.kk$i.pheno_id == zz & overlap.dt.kk$i.qtl_type == "trans_pQTL"])), collapse = ";")))
      genes.underpk$eGenes_peaks_table         <- lapply(genes.underpk$gene_id, function(zz){
        # zz <- genes.underpk$gene_id[284]
        df.zz <- as.data.frame(overlap.dt.kk[overlap.dt.kk$i.pheno_id == zz, ])
        
        if(nrow(df.zz) == 0){ return() }
        # print(df.zz$pos_Mb - df.zz$i.pos_Mb)
        df.zz <- dplyr::select(df.zz, dplyr::all_of(c("i.id", "i.pheno_id", "i.qtl_type", "chr", "chr_gene", "i.LOD_threshold", "i.LOD", "i.pVal", 
                                                      "i.peak_type", "i.pos_Mb", "i.ci_lo_Mb", "i.ci_hi_Mb", "i.inputType",
                                                      "i.condition", "i.interaction", "gene_posStart_Mb", "gene_posEnd_Mb")))
        colnames(df.zz) <- gsub("i\\.", "", colnames(df.zz))
        rownames(df.zz) <- NULL
        df.zz
      })
      genes.underpk$ref_qtl_other_lQTL_overlap <- any(overlap.dt.kk$i.qtl_type == "lQTL" & overlap.dt.kk$i.pheno_id != dt.ref.kk$id)
      if(any(genes.underpk$ref_qtl_other_lQTL_overlap)){
        genes.underpk$ref_qtl_other_lQTL_overlap_ids <- paste(sort(unique(overlap.dt.kk$i.pheno_id[overlap.dt.kk$i.qtl_type == "lQTL" & overlap.dt.kk$i.pheno_id != dt.ref.kk$id])), collapse = ";")
        ref_qtl_other_lQTL_table                     <- as.data.frame(dplyr::select(overlap.dt.kk[overlap.dt.kk$i.qtl_type == "lQTL" & overlap.dt.kk$i.pheno_id != dt.ref.kk$id, ], 
                                                                                    dplyr::all_of(c("i.id", "i.pheno_id", "i.qtl_type", "chr", "i.LOD_threshold", "i.LOD", 
                                                                                                    "i.pVal", "i.peak_type", "i.pos_Mb", "i.ci_lo_Mb", "i.ci_hi_Mb",
                                                                                                    "i.inputType", "i.condition", "i.interaction"))))
        colnames(ref_qtl_other_lQTL_table)           <- gsub("i\\.", "", colnames(ref_qtl_other_lQTL_table))
        rownames(ref_qtl_other_lQTL_table)           <- NULL
        genes.underpk$ref_qtl_other_lQTL_table       <- lapply(1:nrow(genes.underpk), function(zz) ref_qtl_other_lQTL_table)
      } else{
        genes.underpk$ref_qtl_other_lQTL_overlap_ids <- NA
        genes.underpk$ref_qtl_other_lQTL_table       <- lapply(1:nrow(genes.underpk), function(zz) NULL)
      }
      genes.underpk$ref_qtl_other_pheQTL_overlap       <- any(overlap.dt.kk$i.qtl_type == "pheQTL" & overlap.dt.kk$i.pheno_id != dt.ref.kk$id)
      if(any(genes.underpk$ref_qtl_other_pheQTL_overlap)){
        genes.underpk$ref_qtl_other_pheQTL_overlap_ids <- paste(sort(unique(overlap.dt.kk$i.pheno_id[overlap.dt.kk$i.qtl_type == "pheQTL" & overlap.dt.kk$i.pheno_id != dt.ref.kk$id])), collapse = ";")
        ref_qtl_other_pheQTL_table                     <- as.data.frame(dplyr::select(overlap.dt.kk[overlap.dt.kk$i.qtl_type == "pheQTL" & overlap.dt.kk$i.pheno_id != dt.ref.kk$id, ], 
                                                                                      dplyr::all_of(c("i.id", "i.pheno_id", "i.qtl_type", "chr", "i.LOD_threshold", "i.LOD", 
                                                                                                      "i.pVal", "i.peak_type", "i.pos_Mb", "i.ci_lo_Mb", "i.ci_hi_Mb",
                                                                                                      "i.inputType", "i.condition", "i.interaction"))))
        colnames(ref_qtl_other_pheQTL_table)           <- gsub("i\\.", "", colnames(ref_qtl_other_pheQTL_table))
        rownames(ref_qtl_other_pheQTL_table)           <- NULL
        genes.underpk$ref_qtl_other_pheQTL_table       <- lapply(1:nrow(genes.underpk), function(zz) ref_qtl_other_pheQTL_table)
      } else{
        genes.underpk$ref_qtl_other_pheQTL_overlap_ids <- NA
        genes.underpk$ref_qtl_other_pheQTL_table       <- lapply(1:nrow(genes.underpk), function(zz) NULL)
      }
      
     
      
      
      # if a gene has both a cis and trans-(e/p)QTL overlapping in the lipidomic QTL peak, keep only the cis-(e/p)QTL
      tmp.df        <- data.table(overlap.dt.kk)[, list(cis_trans_eQTL = all(c("cis_eQTL", "trans_eQTL") %in% i.qtl_type),
                                                        cis_trans_pQTL = all(c("cis_pQTL", "trans_pQTL") %in% i.qtl_type)), by = c("i.pheno_id")]
      entries.rm    <- c(paste0(tmp.df$i.pheno_id, "_trans_eQTL")[tmp.df$cis_trans_eQTL],
                         paste0(tmp.df$i.pheno_id, "trans_pQTL")[tmp.df$cis_trans_pQTL])
      overlap.dt.kk <- overlap.dt.kk[!(paste0(overlap.dt.kk$i.pheno_id, "_", overlap.dt.kk$i.qtl_type) %in% entries.rm), ]
      
      # paste together ensembl gene id and gene symbol
      tmpLogi                           <- grepl("^ENS", overlap.dt.kk$i.pheno_id)
      overlap.dt.kk$i.pheno_id[tmpLogi] <- paste0(overlap.dt.kk$i.pheno_id[tmpLogi], "__", overlap.dt.kk$gene_symbol[tmpLogi])
      
      dcast_formula                  <- as.formula(paste0(paste(colnames(dt.ref.kk), collapse = "+"), "~i.qtl_type"))
      overlap.dt.kk.dcast.eg         <- reshape2::dcast(overlap.dt.kk, dcast_formula, value.var = "i.pheno_id", fun.aggregate = function(x) paste(sort(unique(x[!is.na(x)])), collapse = ";"))
      
      if(!any(overlap.dt.kk$i.qtl_type == "cis_pQTL")){
        overlap.dt.kk.dcast.eg$cis_pQTL_feature <- NA
      } else{
        cis_pQTL_features <- paste0(overlap.dt.kk$i.pheno_id[overlap.dt.kk$i.qtl_type == "cis_pQTL"], "__", overlap.dt.kk$i.id[overlap.dt.kk$i.qtl_type == "cis_pQTL"])
        cis_pQTL_features[!grepl("ENS", cis_pQTL_features)] <- gsub(".*__", "", cis_pQTL_features)
        overlap.dt.kk.dcast.eg$cis_pQTL_feature             <- paste(unique(cis_pQTL_features), collapse = ";")
      }
      if(!any(overlap.dt.kk$i.qtl_type == "trans_pQTL")){
        overlap.dt.kk.dcast.eg$trans_pQTL_feature <- NA
      } else{
        trans_pQTL_features <- paste0(overlap.dt.kk$i.pheno_id[overlap.dt.kk$i.qtl_type == "trans_pQTL"], "__", overlap.dt.kk$i.id[overlap.dt.kk$i.qtl_type == "trans_pQTL"])
        trans_pQTL_features[!grepl("ENS", trans_pQTL_features)] <- gsub(".*__", "", trans_pQTL_features)
        overlap.dt.kk.dcast.eg$trans_pQTL_feature               <- paste(unique(trans_pQTL_features), collapse = ";")
      }
      stopifnot(nrow(overlap.dt.kk.dcast.eg) == 1)
      
      # update genes.underpk: paste together ensembl gene id and gene symbol
      genes.underpk$gene_id_full <- paste0(genes.underpk$gene_id, "__", genes.underpk$gene_symbol)
      
      
      # update genes.underpk table: add mRNA and prot mRNA
      genes.underpk <- merge(genes.underpk, mRNA_DEA_results, by = c("gene_id"), all.x = T, all.y = F)
      genes.underpk <- merge(genes.underpk, dplyr::select(prot_coon_DEA_results, -dplyr::all_of("uniprot_id")), by = c("gene_id"), all.x = T, all.y = F)
      genes.underpk <- merge(genes.underpk, dplyr::select(geneAnnot.df, -dplyr::all_of(c("mRNA_gene_isDE", "prot_isDE"))), by.x = "gene_id", by.y = "gene_ensembl_id", all.x = T, all.y = F)
      
      
      # get genes to annotate
      tmp.df                    <- dplyr::select(genes.underpk[!(genes.underpk$gene_id %in% gsub("__.*", "", overlap.dt.kk$i.pheno_id)), ], dplyr::all_of(c("gene_id_full", "gene_symbol")))
      colnames(tmp.df)          <- c("i.pheno_id", "gene_symbol")
      gene_map                  <- unique(rbind(dplyr::select(overlap.dt.kk, dplyr::all_of(c("i.pheno_id", "gene_symbol"))),
                                                tmp.df))
      stopifnot(sum(duplicated(gene_map$i.pheno_id)) == 0)
      gene_map                  <- gene_map[grepl("^ENSMUS", gene_map$i.pheno_id), ]
      stopifnot(nrow(gene_map) == nrow(genes.underpk) & all(gsub("__.*", "", gene_map$i.pheno_id) %in% genes.underpk$gene_id))
      gene_map$i.pheno_id       <- gsub("__.*", "", gene_map$i.pheno_id)
      cols_annotate             <- colnames(overlap.dt.kk.dcast.eg)[grepl("cis_|trans_", colnames(overlap.dt.kk.dcast.eg)) & !grepl("feature", colnames(overlap.dt.kk.dcast.eg))]
      genes.annnot              <- unique(unlist(lapply(cols_annotate, function(x) unlist(strsplit(overlap.dt.kk.dcast.eg[[x]], ";")))))
      genes.annnot              <- gsub("__.*", "", genes.annnot[!is.na(genes.annnot)])
      names(genes.annnot)       <- paste0(genes.annnot, "__", plyr::mapvalues(genes.annnot, from = gene_map$i.pheno_id, to = gene_map$gene_symbol, warn_missing = F))
      
      
      # get gene annotations for cis and trans (e/p)genes
      annot.df                <- geneAnnot.df[geneAnnot.df$gene_ensembl_id %in% genes.annnot, ]
      if(nrow(annot.df) > 0){
        annot.df$gene_symbol    <- plyr::mapvalues(annot.df$gene_ensembl_id, from = gene_map$i.pheno_id, to = gene_map$gene_symbol, warn_missing = F)
        annot.df$gene_id        <- paste0(annot.df$gene_ensembl_id, "__", annot.df$gene_symbol)
        annot.df                <- dplyr::select(annot.df, -dplyr::all_of(c("gene_ensembl_id", "gene_symbol")))
        annot.df                <- na.omit(reshape2::melt(annot.df, id.vars = c("gene_id")))
        annot.df.aggr           <- as.data.frame(data.table(annot.df)[, list(annot = paste(paste0(gene_id, "[", value, "]"), collapse = ";")), by = c("variable")])
        rownames(annot.df.aggr) <- annot.df.aggr$variable
        annot.df.aggr           <- annot.df.aggr[, -1, drop = F]
        annot.df.aggr           <- t(annot.df.aggr)
        
        annot.df.aggr           <- as.data.frame(data.table(annot.df)[, list(annot = paste(paste0(gene_id, "[", value, "]"), collapse = ";")), by = c("variable")])
        rownames(annot.df.aggr) <- annot.df.aggr$variable
        annot.df.aggr           <- annot.df.aggr[, -1, drop = F]
        annot.df.aggr           <- t(annot.df.aggr)
      } else{
        annot.df.aggr <- NULL
      }
      
      if(!is.null(annot.df.aggr)){
        overlap.dt.kk.dcast.eg <- cbind(overlap.dt.kk.dcast.eg, annot.df.aggr)
      }
      
      
      list(overlapping_features_table = as.data.frame(overlap.dt.kk.dcast.eg),
           all_genes_table            = as.data.frame(genes.underpk))
      
    })
    
  })
  tmpLogi.error_warnings <- unlist(lapply(overlap.df.xx, function(zz) !is.null(zz) & all(class(zz) == "character")))
  warning_errors         <- unlist(overlap.df.xx[tmpLogi.error_warnings])
  errors_idx             <- gsub("Error at idx | --.*", "", warning_errors[grepl("Error", warning_errors)])
  warnings_idx           <- gsub("Warning at idx | --.*", "", warning_errors[grepl("Warning", warning_errors)])
  if(length(errors_idx) > 0){
    cat("\nThe following idx generated errors - interrupting code execution:\n", paste0("\t", paste(errors_idx, collapse = ";")), "\n")
    stop()
  }
  if(length(warnings_idx) > 0){
    cat("\nThe following idx generated warnings - interrupting code execution:\n", paste0("\t", paste(warnings_idx, collapse = ";")), "\n")
    stop()
  }
  
  genes.underpk.xx <- rbindlist(lapply(overlap.df.xx, function(zz) zz$all_genes_table), use.names = T, fill = T)
  overlap.df.xx    <- rbindlist(lapply(overlap.df.xx, function(zz) zz$overlapping_features_table), use.names = T, fill = T)
  
  genes.underpk.xx <- as.data.frame(genes.underpk.xx)
  overlap.df.xx    <- as.data.frame(overlap.df.xx)
  
  # add potentially missing columns to avoid downstream code to fail
  cols_check    <- c("cis_eQTL", "cis_pQTL", 
                     "trans_eQTL", "trans_pQTL")
  for(i in cols_check){
    if(!(i %in% colnames(overlap.df.xx))){
      overlap.df.xx[[i]] <- NA
    }
  }
  
  # get genes with both cis-(e/p)QTLs
  e_p_cisQTL_list <- lapply(1:nrow(overlap.df.xx), function(x){
    # x <- 17
    if(is.na(overlap.df.xx$cis_eQTL[x]) | is.na(overlap.df.xx$cis_pQTL[x])){ return(list(out_eg = NA, out_gn = NA)) }
    common_genes.eg <- sort(unique(intersect(unlist(strsplit(overlap.df.xx$cis_eQTL[x], ";")), unlist(strsplit(overlap.df.xx$cis_pQTL[x], ";")))))
    list(out_eg = ifelse(length(common_genes.eg) == 0, NA, paste(common_genes.eg, collapse = ";")))
  })
  overlap.df.xx$e_p_cisQTL <- unlist(lapply(e_p_cisQTL_list, function(x) x$out_eg))
  # overlap.df.xx$e_p_cisQTL_gene_symbol <- unlist(lapply(e_p_cisQTL_list, function(x) x$out_gn))
  
  # get genes with both trans-(e/p)QTLs
  e_p_transQTL_list <- lapply(1:nrow(overlap.df.xx), function(x){
    # x <- 1
    if(is.na(overlap.df.xx$trans_eQTL[x]) | is.na(overlap.df.xx$trans_pQTL[x])){ return(list(out_eg = NA, out_gn = NA)) }
    common_genes.eg <- sort(unique(intersect(unlist(strsplit(overlap.df.xx$trans_eQTL[x], ";")), unlist(strsplit(overlap.df.xx$trans_pQTL[x], ";")))))
    list(out_eg = ifelse(length(common_genes.eg) == 0, NA, paste(common_genes.eg, collapse = ";")))
  })
  overlap.df.xx$e_p_transQTL <- unlist(lapply(e_p_transQTL_list, function(x) x$out_eg))

  # get genes with both (e/p)QTLs (cis or trans)
  # paste gene_name(cis_eQTL;trans_pQTL) for example
  e_p_QTL_list <- lapply(1:nrow(overlap.df.xx), function(x){
   
    if((is.na(overlap.df.xx$cis_eQTL[x]) & is.na(overlap.df.xx$trans_eQTL[x])) | (is.na(overlap.df.xx$cis_pQTL[x]) & is.na(overlap.df.xx$trans_pQTL[x])) | (!is.na(overlap.df.xx$e_p_cisQTL[x]) | !is.na(overlap.df.xx$e_p_transQTL[x]))){ return(list(out_eg = NA, out_gn = NA)) }
    
    labels.x <- lapply(c("cis_eQTL", "trans_eQTL", "cis_pQTL", "trans_pQTL"), function(y){
      if(!is.na(overlap.df.xx[[y]][x])){
        data.frame(gene_id = unlist(strsplit(overlap.df.xx[[y]][x], ";")), gene_label = y, stringsAsFactors = F)
      } else{
        NULL
      }
    })
    labels.x <- do.call(rbind, labels.x)
    
    labels.x      <- labels.x[!is.na(labels.x$gene_id), ]
    labels.x      <- labels.x[labels.x$gene_id %in% labels.x$gene_id[duplicated(labels.x$gene_id)], ]
    if(nrow(labels.x) == 0){ return(list(out_eg = NA, out_gn = NA)) }
    labels.x.aggr <- data.table(labels.x)[, list(labels_all = paste(sort(unique(gene_label)), collapse = ";")), by = c("gene_id")]
    out.eg        <- paste(sort(paste0(labels.x.aggr$gene_id, "(", labels.x.aggr$labels_all, ")")), collapse = ";")
    
    
    list(out_eg = out.eg)
  })
  overlap.df.xx$e_p_QTL <- unlist(lapply(e_p_QTL_list, function(x) x$out_eg))

  
  # rename and reorder columns of overlapping features table
  colnames(overlap.df.xx) <- gsub("_gene_symbol$", "_gs", colnames(overlap.df.xx))
  colnames(overlap.df.xx) <- gsub("_gene_synonyms$", "_gs_synonyms", colnames(overlap.df.xx))
  cols_order              <- c(colnames(df.lqtl.species), "data_type",
                               "pheQTL", "lQTL", "eigQTL", "cis_eQTL", "cis_pQTL", "cis_pQTL_feature", "trans_eQTL", "trans_pQTL","trans_pQTL_feature",
                               "e_p_cisQTL", "e_p_transQTL", "e_p_QTL",
                               colnames(geneAnnot.df)[colnames(geneAnnot.df) %in% colnames(overlap.df.xx)])
  cols_order <- cols_order[cols_order %in% colnames(overlap.df.xx)]
  stopifnot(all(colnames(overlap.df.xx) %in% cols_order))

  overlap.df.xx <- dplyr::select(overlap.df.xx, dplyr::all_of(cols_order))
  
  # rename and reorder columns of all genes table
  genes.underpk.xx <- dplyr::select(genes.underpk.xx, -dplyr::all_of(c("inputType")))
  colnames(genes.underpk.xx)[colnames(genes.underpk.xx) == "chr"]                    <- "gene_chr"
  colnames(genes.underpk.xx)[colnames(genes.underpk.xx) == "is_under_ref_qtl"]       <- "gene_is_under_ref_qtl"
  
                                 

  cols_order <- c("ref_qtl_id", "ref_qtl_chr", "ref_qtl_pos_Mb", "ref_qtl_ci_lo_Mb", "ref_qtl_ci_hi_Mb",
                  "ref_qtl_LOD_threshold", "ref_qtl_LOD", "ref_qtl_pVal", "ref_qtl_peak_type",
                  
                  "ref_qtl_inputType", "ref_qtl_data_type", "ref_qtl_aggregation_type", 
                  
                  "condition", "interaction", "analysis_id",
                  
                  "ref_qtl_other_lQTL_overlap", "ref_qtl_other_lQTL_overlap_ids", "ref_qtl_other_lQTL_table",
                  "ref_qtl_other_pheQTL_overlap", "ref_qtl_other_pheQTL_overlap_ids", "ref_qtl_other_pheQTL_table",
                  
                  "gene_id", "gene_symbol", "gene_id_full", "gene_chr", "gene_posStart_Mb", "gene_posEnd_Mb",
                  "gene_start_dist_from_pk_Mb", "gene_is_under_ref_qtl",

                  "cis_eQTL_overlap", "cis_pQTL_overlap", "cis_pQTL_overlap_feature", "trans_eQTL_overlap", "trans_pQTL_overlap", "trans_pQTL_overlap_feature",
    
                  "eGenes_peaks_table",
                  
                  "DEA_logFC_mRNA", "DEA_pval_mRNA", "DEA_adjpval_mRNA",
                  "DEA_logFC_prot", "DEA_pval_prot", "DEA_adjpval_prot")
  
  cols_order <- cols_order[cols_order %in% colnames(genes.underpk.xx)]
  stopifnot(all(colnames(genes.underpk.xx) %in% cols_order))
  # colnames(genes.underpk.xx)[!(colnames(genes.underpk.xx) %in% cols_order)]
  
  genes.underpk.xx <- dplyr::select(genes.underpk.xx, dplyr::all_of(cols_order))
  # table(genes.underpk.xx$ref_qtl_id)
  
  cols_char <- colnames(overlap.df.xx)[unlist(lapply(overlap.df.xx, function(zz) any(class(zz) == "character")))]
  for(i in cols_char){
    overlap.df.xx[[i]] <- gsub('"', "", overlap.df.xx[[i]])
    overlap.df.xx[[i]] <- gsub('\n', "", overlap.df.xx[[i]], fixed = T)
    overlap.df.xx[[i]] <- gsub('\t', "", overlap.df.xx[[i]], fixed = T)
  }
  
  cols_char <- colnames(genes.underpk.xx)[unlist(lapply(genes.underpk.xx, function(zz) any(class(zz) == "character")))]
  for(i in cols_char){
    genes.underpk.xx[[i]] <- gsub('"', "", genes.underpk.xx[[i]])
    genes.underpk.xx[[i]] <- gsub('\n', "", genes.underpk.xx[[i]], fixed = T)
    genes.underpk.xx[[i]] <- gsub('\t', "", genes.underpk.xx[[i]], fixed = T)
  }
  
  saveRDS(overlap.df.xx, paste0(saveDir, "/lipids_prioritized___", xx, ".RDS"))
  saveRDS(genes.underpk.xx, paste0(saveDir, "/lipids_all_genes_criteria___", xx, ".RDS"))
  remove(overlap.df.xx, genes.underpk.xx)
  cat("\n")
  NA
  # list(overlap.df.xx, genes.underpk.xx)
})
names(overlap.df) <- conditions.run



