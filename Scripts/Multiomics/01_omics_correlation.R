
library(data.table)
library(parallel)


################################################################################
## omic-omic correlation
################################################################################


saveDir        <- "./Data/lipidomics/omics_correlation"
saveDir_data   <- "./Data/lipidomics/omics_correlation_data"
for(i in c(saveDir, saveDir_data)){
  if(!dir.exists(i)){
    dir.create(i, recursive = T)
  }
}
overwrite      <- F
overwrite_data <- F
tt             <- lapply(c("pearson"), function(cc){
  # cc <- "pearson"
  
  lipid.derived.list  <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")
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
  
  
  
  lipid.avg.df.list <- lipid.derived.list[grepl("derived_avgAggr", names(lipid.derived.list))]
  lipid.avg.df.list <- lapply(names(lipid.avg.df.list), function(x){
    
    tryCatch({
      df.x <- lipid.avg.df.list[[x]]
      df.x <- dplyr::select(df.x, dplyr::all_of(colnames(df.x)[grepl("BXD|C57|DBA", colnames(df.x), ignore.case = T)]))
      colnames(df.x) <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", colnames(df.x)))
      df.x <- as.data.frame(t(df.x))
      df.x$strainDiet <- rownames(df.x)
      df.x[df.x == Inf] <- NA
      df.x[df.x == -Inf] <- NA
      
      cat("\t--> Transformation for data type", x, ":\n\t\t log2\n")
      df.x[, cols_data.logi] <- log2(df.x[, cols_data.logi] + 1)
      
    }, warning = function(w){
      print(x)
    })
    
    
    df.x
  })
  lipid.avg.df <- Reduce(function(...) merge(..., all = T, by = c("strainDiet")), lipid.avg.df.list)
  rownames(lipid.avg.df) <- lipid.avg.df$strainDiet
  lipid.avg.df <- dplyr::select(lipid.avg.df, -dplyr::all_of("strainDiet"))
  lipid.avg.df <- lipid.avg.df[rownames(pheno.avg.df), ]
  
  
  metadata                 <- readRDS("./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")
  rnaseq.df                <- readRDS("./Data/RNAseq_processing/countMatrix/geneLevel/STAR_RNAseq_count_matrix.RDS")
  metadata$strainDiet      <- paste0(metadata$strain, "_", metadata$diet)
  metadata                 <- metadata[metadata$sample_name %in% colnames(rnaseq.df), ]
  metadata                 <- metadata[match(colnames(rnaseq.df), metadata$sample_name), ]
  rownames(metadata)       <- metadata$sample_name
  colnames(metadata)       <- gsub("\\.", "_", gsub("\\[", "", gsub("\\]\\.|\\/", "_", colnames(metadata))))
  colnames(metadata)       <- make.names(colnames(metadata))
  stopifnot(identical(colnames(rnaseq.df), colnames(rnaseq.df)))
  rnaseq.df                <- as.data.frame(edgeR::cpm(rnaseq.df))
  rnaseq.df                <- rnaseq.df[rowSums(rnaseq.df > 0.01) > 0.5 * ncol(rnaseq.df), ]
  rnaseq.df                <- as.data.frame(t(rnaseq.df))
  rnaseq.df                <- cbind(strainDiet = gsub("CD.*", "CD", gsub("HFD.*", "HFD", rownames(rnaseq.df))), rnaseq.df)
  rnaseq.df.aggr           <- as.data.frame(data.table(rnaseq.df)[, lapply(.SD, function(zz) mean(zz, na.rm = T)), by = c("strainDiet")])
  rownames(rnaseq.df.aggr) <- rnaseq.df.aggr$strainDiet
  rnaseq.df.aggr           <- dplyr::select(rnaseq.df.aggr, -dplyr::all_of("strainDiet"))
  rnaseq.df.aggr           <- rnaseq.df.aggr[rownames(pheno.avg.df), ]

  
  
  
  prot.df           <- readRDS("./Data/input_data/proteomics/coon_formatted/raw_intensity_table_samples_excluded.RDS")
  samples.meta.df   <- readRDS("./Data/input_data/proteomics/coon_formatted/sample_metadata.RDS")
  prot.meta.df      <- readRDS("./Data/input_data/proteomics/coon_formatted/proteins_metadata.RDS")
  prot.meta.df.aggr <- data.table(prot.meta.df)[, list(gene_ids = ifelse(all(is.na(ensembl_gene_id)),
                                                                         as.character(NA),
                                                                         paste(unique(ensembl_gene_id[!is.na(ensembl_gene_id)]), collapse = "_"))), by = c("uniprot_id")]
  prot.df           <- type.convert(reshape2::melt(as.matrix(prot.df)), as.is = T)
  prot.df           <- merge(prot.df, samples.meta.df, by.x = "Var2", by.y = "sample", all.x = T, all.y = F)
  prot.df.aggr      <- data.table(prot.df)[, list(value = mean(value, na.rm = T)),
                                           by = c("strain", "diet", "Var1")]
  prot.df.aggr$strainDiet      <- paste0(prot.df.aggr$strain, "_", prot.df.aggr$diet)
  prot.df.aggr.meta            <- unique(dplyr::select(prot.df.aggr, c("strain", "diet", "strainDiet")))
  prot.df.aggr.meta$sample     <- prot.df.aggr.meta$strainDiet
  prot.df.aggr.meta            <- type.convert(as.data.frame(prot.df.aggr.meta), as.is = T)
  prot.df.aggr                 <- dcast(prot.df.aggr, Var1~strainDiet, value.var = "value")
  prot.df.aggr                 <- as.data.frame(prot.df.aggr)
  prot.df.aggr[is.na(prot.df.aggr)] <- NA
  prot.df.aggr                 <- cbind(ensembl_gene_id = plyr::mapvalues(prot.df.aggr$Var1, from = prot.meta.df.aggr$uniprot_id, to = prot.meta.df.aggr$gene_ids, warn_missing = F), prot.df.aggr)
  prot.df.aggr                 <- cbind(feature_id = ifelse((is.na(prot.df.aggr$ensembl_gene_id) | !grepl("ENS", prot.df.aggr$ensembl_gene_id)), 
                                                            prot.df.aggr$Var1, 
                                                            paste0(prot.df.aggr$ensembl_gene_id, "__", prot.df.aggr$Var1)), prot.df.aggr)
  rownames(prot.df.aggr)      <- prot.df.aggr$feature_id
  prot.df.aggr                <- as.data.frame(t(dplyr::select(prot.df.aggr, -dplyr::all_of(c("feature_id", "ensembl_gene_id", "Var1")))))
  prot.df.aggr                <- prot.df.aggr[, colSums(is.na(prot.df.aggr)) < 0.5 * nrow(prot.df.aggr)]
  prot.df.aggr                <- prot.df.aggr[match(rownames(pheno.avg.df), rownames(prot.df.aggr)), ]
  
  
  
  stopifnot(identical(rownames(pheno.avg.df), rownames(lipid.avg.df)))
  stopifnot(identical(rownames(pheno.avg.df), rownames(rnaseq.df.aggr)))
  stopifnot(identical(rownames(pheno.avg.df), rownames(prot.df.aggr)))
  
  
  rnaseq.df.aggr    <- log2(rnaseq.df.aggr + 1)
  prot.df.aggr      <- log2(prot.df.aggr + 1)
  
  df_list <- list(lipid           = lipid.avg.df,
                  phenotype       = pheno.avg.df,
                  gene            = rnaseq.df.aggr,
                  protein         = prot.df.aggr)
  
  if(!(file.exists(paste0(saveDir_data, "/omics_pheno_avg_data__for_pearson__omics_log2.RDS")) & !overwrite_data)){
    saveRDS(df_list, paste0(saveDir_data, "/omics_pheno_avg_data__for_pearson__omics_log2.RDS"))
  } else{
    cat("Data not saved\n")
  }
  corr_pairs <- c("lipid_phenotype", "lipid_gene", "lipid_protein", 
                  "phenotype_gene", "phenotype_protein", 
                  "lipid_lipid", "gene_gene", "protein_protein", "phenotype_phenotype")
  
  corr.list  <- lapply(corr_pairs, function(x){
    # x <- "lipid_phenotype"
    
    elements.x         <- unlist(strsplit(x, "_"))
    conditions.present <- c(rownames(df_list[[elements.x[1]]]), rownames(df_list[[elements.x[2]]]))
    conditions.do      <- c()
    if(any(grepl("CD", conditions.present))){ conditions.do <- c(conditions.do, "CD") }
    if(any(grepl("HFD", conditions.present))){ conditions.do <- c(conditions.do, "HFD") }
    if(any(grepl("CD", conditions.present)) & any(grepl("HFD", conditions.present))){ conditions.do <- c(conditions.do, "CD_HFD") }
    if(length(conditions.do) == 0){ conditions.do <- "all" }
    
    cat("Computing correlation pairs for (", which(x == corr_pairs), "/", length(corr_pairs), "):\t", paste(sort(conditions.do), collapse = ", "), "\t", gsub("_", "_vs_", x), "...\n")
    
    parallel::mclapply(mc.cores = 3, X = conditions.do, FUN = function(y){
      # y <- "CD_HFD"
      # y <- "CD"
      # y <- "HFD"
      # y <- "all"
      
      out_file <- paste0(saveDir, "/corr_table__", cc, "__", gsub("_", "_vs_", x), ifelse(y == "all", "", paste0("__", y)), ".RDS")
      if(file.exists(out_file) & !overwrite){
        return()
      }
      
      elements.x <- unlist(strsplit(x, "_"))
      df.1       <- df_list[[elements.x[1]]]
      df.2       <- df_list[[elements.x[2]]]
      if(y != "CD_HFD" & y != "all"){
        df.1 <- df.1[grepl(y, rownames(df.1)), ]
        df.2 <- df.2[grepl(y, rownames(df.2)), ]
      }
      cor_pval   <- WGCNA::corAndPvalue(df.1, df.2, method = cc)
      
      df.cor  <- cor_pval$cor
      df.pval <- cor_pval$p
      if(elements.x[1] == elements.x[2]){
        df.cor[upper.tri(df.cor, diag = T)]   <- NA
        df.pval[upper.tri(df.pval, diag = T)] <- NA
        # rename for columns naming
        elements.x <- paste0(elements.x, "_", 1:length(elements.x))
      }
      
      df.cor            <- reshape2::melt(df.cor)
      colnames(df.cor)  <- c(elements.x, "cor")
      df.pval           <- reshape2::melt(df.pval)
      colnames(df.pval) <- c(elements.x, "pval")
      df.cor$pval       <- df.pval$pval
      df.cor            <- df.cor[!is.na(df.cor$cor) & !is.na(df.cor$pval), ]
      df.cor$condition  <- ifelse(y == "all", NA, y)
      df.cor            <- df.cor[order(df.cor$cor, df.cor$pval), ]
      
      df.cor[, 1] <- as.character(df.cor[, 1])
      df.cor[, 2] <- as.character(df.cor[, 2])
      
      # df.cor            <- type.convert(df.cor, as.is = T)
      
      saveRDS(df.cor, out_file)
      
    })
  })
  
  
  ## also save a version only with all correlation pairs but no unknown lipids
  ## also save a version only with significant correlation pairs and no unknown lipids
  
  # tt <- lapply(c("pearson", "spearman", "kendall"), function(cc){
  # cc <- "spearman"
  corr_files <- list.files("./Data/lipidomics/omics_correlation", full.names = T, pattern = "RDS$")
  corr_files <- corr_files[!grepl("_filt_", corr_files)]
  corr_files <- corr_files[grepl(paste0("__", cc, "__"), corr_files)]
  corr_files <- corr_files[grepl(paste(gsub("_", "_vs_", corr_pairs), collapse = "|"), corr_files)]
  tt         <- parallel::mclapply(mc.cores = length(corr_files), X = corr_files, FUN = function(x){
    # x <- corr_files[15]
    
    out_files <- c(gsub("\\.RDS", "__filt__no_unknowns.RDS", x),
                   gsub("\\.RDS", "__pval_filt__no_unknowns.RDS", x))
    if(all(file.exists(out_files)) & !overwrite){
      return()
    }
    
    df.x <- readRDS(x)
    
    df.x.sub <- df.x[!grepl("unknown", df.x[, 1], ignore.case = T) & !grepl("unknown", df.x[, 2], ignore.case = T), ]
    saveRDS(df.x.sub, gsub("\\.RDS", "__filt__no_unknowns.RDS", x))
    
    df.x.sub <- df.x[df.x$pval < 0.05 & !grepl("unknown", df.x[, 1], ignore.case = T) & !grepl("unknown", df.x[, 2], ignore.case = T), ]
    saveRDS(df.x.sub, gsub("\\.RDS", "__pval_filt__no_unknowns.RDS", x))
    NA
  })
  # })
})


## also save a version only with all correlation pairs but no unknown lipids
## also save a version only with significant correlation pairs and no unknown lipids
tt <- lapply(c("pearson"), function(cc){
  # cc <- "spearman"
  cat(cc, "...\n")
  corr_files <- list.files("./Data/lipidomics/omics_correlation", full.names = T, pattern = "RDS$")
  corr_files <- corr_files[!grepl("_filt_", corr_files)]
  corr_files <- corr_files[grepl(paste0("__", cc, "__"), corr_files)]
  corr_files <- corr_files[grepl(paste(gsub("_", "_vs_", corr_pairs), collapse = "|"), corr_files)]
  tt         <- parallel::mclapply(mc.cores = length(corr_files), X = corr_files, FUN = function(x){
    # x <- corr_files[15]
    
    out_files <- c(gsub("\\.RDS", "__filt__no_unknowns.RDS", x),
                   gsub("\\.RDS", "__pval_filt__no_unknowns.RDS", x))
    if(all(file.exists(out_files)) & !overwrite){
      return()
    }
    
    df.x <- readRDS(x)
    
    df.x.sub <- df.x[!grepl("unknown", df.x[, 1], ignore.case = T) & !grepl("unknown", df.x[, 2], ignore.case = T), ]
    saveRDS(df.x.sub, gsub("\\.RDS", "__filt__no_unknowns.RDS", x))
    
    df.x.sub <- df.x[df.x$pval < 0.05 & !grepl("unknown", df.x[, 1], ignore.case = T) & !grepl("unknown", df.x[, 2], ignore.case = T), ]
    saveRDS(df.x.sub, gsub("\\.RDS", "__pval_filt__no_unknowns.RDS", x))
    NA
  })
})
