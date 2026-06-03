
library(data.table)
library(parallel)



genes_list <- c("Usp46" = "ENSMUSG00000054814",
                "Uchl1" = "ENSMUSG00000029223")

cat("\n Gene-gene correlation-ranked GSEA run for the following genes:\n",
    paste(paste0("\t--> ", paste0(names(genes_list), " (", unname(genes_list), ")")), collapse = "\n"), "\n",
    "!!! Modify gene selection is mediation scripts to run for other genes\n\n")

################################################################################
## GSEA using as ranking candidate gene - genes correlation
################################################################################


corr_mat.list     <- list("gene"    = c("CD"        = "./Data/lipidomics/omics_correlation/corr_table__pearson__gene_vs_gene__CD__filt__no_unknowns.RDS",
                                        "HFD"       = "./Data/lipidomics/omics_correlation/corr_table__pearson__gene_vs_gene__HFD__filt__no_unknowns.RDS",
                                        "CD_HFD"    = "./Data/lipidomics/omics_correlation/corr_table__pearson__gene_vs_gene__CD_HFD__filt__no_unknowns.RDS"),
                          
                          "protein" = c("CD"        = "./Data/lipidomics/omics_correlation/corr_table__pearson__protein_vs_protein__CD__filt__no_unknowns.RDS",
                                        "HFD"       = "./Data/lipidomics/omics_correlation/corr_table__pearson__protein_vs_protein__HFD__filt__no_unknowns.RDS",
                                        "CD_HFD"    = "./Data/lipidomics/omics_correlation/corr_table__pearson__protein_vs_protein__CD_HFD__filt__no_unknowns.RDS"))
stopifnot(all(file.exists(unlist(corr_mat.list))))

msigdbr.geneSets <- list(BP           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP"),
                         CC           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC"),
                         MF           = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "MF"),
                         REACTOME     = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"),
                         HPO          = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "HPO"),
                         HALLMARK     = msigdbr::msigdbr(species = "Mus musculus", category = "H"),
                         KEGG         = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG"),
                         WIKIPATHWAYS = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS"))



# load and subset correlation matrices just onces to improve speed
tmp_dir <- "./tmp_corr_mat_filt"
if(!dir.exists(tmp_dir)){
  dir.create(tmp_dir, recursive = T)
}
overwrite <- F
tt        <- parallel::mclapply(mc.cores = length(unlist(corr_mat.list)), X = unlist(corr_mat.list), FUN = function(x){
  # x <- unlist(corr_mat.list)[4]
  out_file <- paste0(tmp_dir, "/", gsub("corr_table", "tmp_filt_corr_table", gsub(".*\\/", "", x)))
  if(file.exists(out_file) & !overwrite){
    return()
  }
  df.x <- readRDS(x)
  if(grepl("gene_vs_gene", x)){
    df.x <- df.x[df.x$gene_1 %in% genes_list, ]
  } else if(grepl("protein_vs_protein", x)){
    gene_vec <- gsub("__.*", "", df.x$protein_1)
    df.x     <- df.x[gene_vec %in% unique(gsub("__.*", "", genes_list)), ]
  } else{
    stop("ERR")
  }
  saveRDS(df.x, out_file)
})

corr_mat.list     <- list("gene"    = c("CD"        = "./tmp_corr_mat_filt/tmp_filt_corr_table__pearson__gene_vs_gene__CD__filt__no_unknowns.RDS",
                                        "HFD"       = "./tmp_corr_mat_filt/tmp_filt_corr_table__pearson__gene_vs_gene__HFD__filt__no_unknowns.RDS",
                                        "CD_HFD"    = "./tmp_corr_mat_filt/tmp_filt_corr_table__pearson__gene_vs_gene__CD_HFD__filt__no_unknowns.RDS"),
                          
                          "protein" = c("CD"        = "./tmp_corr_mat_filt/tmp_filt_corr_table__pearson__protein_vs_protein__CD__filt__no_unknowns.RDS",
                                        "HFD"       = "./tmp_corr_mat_filt/tmp_filt_corr_table__pearson__protein_vs_protein__HFD__filt__no_unknowns.RDS",
                                        "CD_HFD"    = "./tmp_corr_mat_filt/tmp_filt_corr_table__pearson__protein_vs_protein__CD_HFD__filt__no_unknowns.RDS"))
stopifnot(all(file.exists(unlist(corr_mat.list))))


nbCores   <- 10
# nbCores   <- 20
overwrite <- F
set.seed(123)
gsea      <- lapply(names(corr_mat.list), function(y){
  # y <- names(corr_mat.list)[1]
  
  tt <- lapply(names(corr_mat.list[[y]]), function(w){
    # w <- names(corr_mat.list[[y]])[1]
    
    df.corr           <- readRDS(corr_mat.list[[y]][[w]])
    genes.run         <- unique(df.corr[, 1])
    genes.run         <- genes.run[genes.run %in% genes_list]
    colnames(df.corr) <- c("feat_1", "feat_2", "correlation", "pval", "condition")
    df.corr$feat_1    <- gsub("_.*", "", df.corr$feat_1)
    df.corr$feat_2    <- gsub("_.*", "", df.corr$feat_2)
    
    stopifnot(all(unique(df.corr$condition) == w))
    
    # before aggregation
    df.corr <- plyr::ddply(df.corr, c("feat_1"), function(kk){
      kk$p.adjust <- p.adjust(kk$pval, method = "BH") 
      kk
    })
    
    if(grepl("protein_vs_protein", corr_mat.list[[y]][[w]])){
      # keep feature with highest correlation
      df.corr.aggr <- data.table(df.corr)[, list(max_cor          = correlation[which.max(abs(correlation))],
                                                 max_cor_adj_pval = p.adjust[which.max(abs(correlation))]), by = c("feat_1", "feat_2", "condition")]
      df.corr.aggr <- as.data.frame(df.corr.aggr)
    } else{
      df.corr.aggr              <- df.corr
      colnames(df.corr.aggr)[3] <- "max_cor"
      colnames(df.corr.aggr)[6] <- "max_cor_adj_pval"
    }
    
    cat("Computing GSEA for correlation type:\t", y, "-", w, "-", length(genes.run), "features (genes)...\n")
    
    saveDir <- paste0("./Data/lipidomics/gene_gene_correlation_GSEA/full_tables/", y, "__", w, "__pearson")
    if(!dir.exists(saveDir)){
      dir.create(saveDir, recursive = T)
    }
    
    out.list <- parallel::mclapply(mc.cores = nbCores, X = genes.run, FUN = function(k){
      # k <- genes.run[1]
      
      
      out_file <- paste0(saveDir, "/", y, "__", w, "__pearson_corr_GSEA_gene_vs_gene_results__", k, ".RDS")
      if(file.exists(out_file) & !overwrite){
        cat("\t--> Results files already present. Skipping...\n")
        return()
      }
      
      tmp <- df.corr.aggr[df.corr.aggr$feat_1 == k, ]
      
      sig_thresh <- c(0.05, 0.1)
      genes.sig  <- lapply(sig_thresh, function(w){
        tmp$feat_2[tmp$max_cor_adj_pval < w]
      })
      names(genes.sig) <- paste0("sig_", gsub("\\.", "_", sig_thresh))
      
      
      
      tmp             <- tmp[order(tmp$max_cor, decreasing = T), ]
      geneList        <- tmp$max_cor
      names(geneList) <- tmp$feat_2
      
      remove(tmp)
      
      out.list <- lapply(names(msigdbr.geneSets), function(x){
        # x <- names(msigdbr.geneSets)[1]
        # x <- names(msigdbr.geneSets)[9]
        # print(x)
        
        set.seed(123)
        gsea <- clusterProfiler::GSEA(geneList      = geneList,
                                      # nPermSimple   = 100000,
                                      minGSSize     = 4,
                                      pvalueCutoff  = 1,
                                      pAdjustMethod = "BH",
                                      TERM2GENE     = dplyr::select(msigdbr.geneSets[[x]], dplyr::all_of(c("gs_name", "ensembl_gene"))))
        gsea_df <- gsea@result
        if(nrow(gsea_df) == 0){
          return()
        }
        remove(gsea)
        termsMeta  <- unique(dplyr::select(msigdbr.geneSets[[x]], dplyr::all_of(c("gs_name"))))
        gsea_df    <- merge(gsea_df, termsMeta, by.x = "ID", by.y = "gs_name", all.x = T, all.y = F)
        gsea_df    <- gsea_df[order(gsea_df$pvalue), ]
        
        rownames(gsea_df)         <- gsea_df$ID
        gsea_df$geneSetCollection <- x
        gsea_df$gene_id           <- k
        gsea_df$condition         <- w
        gsea_df$gene_feature_type <- y
        
        for(i in names(genes.sig)){
          if(length(genes.sig[[i]]) > 0){
            
            tmpLogi.genes                               <- lapply(genes.sig[[i]], function(g) grepl(g, gsea_df$core_enrichment))
            names(tmpLogi.genes)                        <- genes.sig[[i]]
            gsea_df[[paste0("core_has_corr_gene_", i)]] <- Reduce("|", tmpLogi.genes)
            gsea_df[[paste0("core_corr_gene_", i)]]     <- apply(as.matrix(Reduce("cbind", tmpLogi.genes)), 1, function(g) ifelse(all(!g), NA, paste(genes.sig[[i]][g], collapse = "/")))
            
          } else{
            
            gsea_df[[paste0("core_has_corr_gene_", i)]] <- F
            gsea_df[[paste0("core_corr_gene_", i)]]     <- NA
            
          }
        }
        
        gsea_df
        
      })
      if(all(unlist(lapply(out.list, is.null)))){
        return()
      }
      out.df           <- unique(rbindlist(out.list, use.names = T, fill = T))
      out.df$coreSize  <- unlist(lapply(out.df$core_enrichment, function(x) length(unlist(strsplit(x, "/")))))
      out.df$geneRatio <- out.df$coreSize / out.df$setSize
      
      
      saveRDS(out.df, out_file)
      NA
    })
  })
})

