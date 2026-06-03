

library(WGCNA)
library(ggplot2)
library(data.table)

source("./Scripts/Lipidomics/lipidomics_WGCNA_fun.R")
source("./Scripts/Lipidomics/lipidomics_fun.R")


################################################################################
## load required data
################################################################################

lipid.derived.list      <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")
lipid.avg.df            <- lipid.derived.list$derived_avgAggr_data[!lipid.derived.list$derived_avgAggr_data$isUnknown, ]
lipid.avg.df            <- dplyr::select(lipid.avg.df, dplyr::all_of(colnames(lipid.avg.df)[grepl("BXD|C57|DBA", colnames(lipid.avg.df), ignore.case = T)]))
colnames(lipid.avg.df)  <- gsub("C57BL6", "C57BL/6J", gsub("DBA2J", "DBA/2J", colnames(lipid.avg.df)))
lipid.avg.df            <- as.data.frame(t(lipid.avg.df))
lipid.avg.df[lipid.avg.df == Inf]  <- NA
lipid.avg.df[lipid.avg.df == -Inf] <- NA
lipid.avg.df            <- log2(lipid.avg.df + 1)
remove(lipid.derived.list)

#############################################
## complete WGCNA analysis
#############################################

saveDir <- "./Data/lipidomics/wgcna_coexpression_analysis"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}

# Note: perform analysis on avg-aggregated data to be able to do downstream association analyses with RNAseq data (pooled)
lipid.avg.list <- split(lipid.avg.df, gsub(".*_", "", rownames(lipid.avg.df)))
min_mod_size   <- 5
wgcna_res_list <- run_wgcna_analysis(omics_df_list         = lipid.avg.list,
                                     r2_thresh             = 0.85,
                                     min_mod_size          = min_mod_size,
                                     pam_respects_dendro   = T,
                                     nb_preservation_perm  = 500,
                                     tom_dir               = paste0(saveDir, "/tom_matrices"),
                                     do_stability_analysis = T,
                                     nbCores               = 40,
                                     verbose               = T)
saveRDS(wgcna_res_list, paste0(saveDir, "/wgcna_results_list__avg_data__min_mod_size_", min_mod_size, ".RDS"))


# re-run modules preservation analysis with NetRep
# Adaped from ModulePreservationNetRep implementation at:
# https://github.com/smorabit/hdWGCNA/blob/HEAD/R/ModulePreservation.R

wgcna_res_list <- readRDS(paste0("./Data/lipidomics/wgcna_coexpression_analysis/wgcna_results_list__avg_data__min_mod_size_", min_mod_size, ".RDS"))
obj_list       <- lapply(names(wgcna_res_list$wgcna_results_list), function(x){
  # x <- "CD"
  
  stopifnot(all(names(wgcna_res_list$wgcna_results_list[[x]]$wgcna_net$colors) == colnames(wgcna_res_list$wgcna_results_list[[x]]$wgcna_net$data_obj)))
  
  TOM_file    <- list.files("./Data/lipidomics/wgcna_coexpression_analysis/tom_matrices/", full.names = T)
  TOM_file    <- TOM_file[grepl(x, TOM_file)]
  TOM_mat     <- as.matrix(get(load(TOM_file)))
  rownames(TOM_mat) <- colnames(TOM_mat) <- names(wgcna_res_list$wgcna_results_list[[x]]$wgcna_net$colors)
  data_mat    <- as.matrix(wgcna_res_list$wgcna_results_list[[x]]$wgcna_net$data_obj)
  adj_mat     <- WGCNA::adjacency(datExpr = data_mat, 
                                  power   = wgcna_res_list$wgcna_results_list[[x]]$wgcna_net$power,
                                  type    = "signed",
                                  corFnc  = "bicor")
  cor_mat     <- WGCNA::bicor(data_mat, use = "pairwise.complete.obs")
  modules_vec <- wgcna_res_list$wgcna_results_list[[x]]$wgcna_net$colors
  
  # NetRep::modulePreservation doesn't allow missing or infinite data - impute as done internally by WGCNA
  data_mat[is.infinite(data_mat)] <- NA
  if(any(is.na(data_mat))){
    stopifnot(all(names(modules_vec) == colnames(data_mat)))
    data_mat <- WGCNA::imputeByModule(data_mat, labels = modules_vec, scale = F) 
  }
  
  
  list(data_mat    = data_mat,
       TOM_mat     = TOM_mat,
       adj_mat     = adj_mat,
       modules_vec = modules_vec)
})
names(obj_list) <- names(wgcna_res_list$wgcna_results_list)

nbCores   <- 70
nb_perm   <- 100000
overwrite <- T
tt        <- lapply(names(wgcna_res_list$wgcna_results_list), function(x){
  # x <- "CD"
  
  cat("Running analysis with reference diet:", x, "...\n")
  
  out_file <- paste0("./Data/lipidomics/wgcna_coexpression_analysis/NetRep_modules_preservation__", x, "__min_mod_size_", min_mod_size, ".RDS")
  if(file.exists(out_file) & !overwrite){
    return()
  }
  
  net_list   <- lapply(obj_list, function(y) y$TOM_mat)
  data_list  <- lapply(obj_list, function(y) y$data_mat)
  adj_list   <- lapply(obj_list, function(y) y$adj_mat)
  mod_labels <- obj_list[[x]]$modules_vec # reference modules assignment
  mod_labels <- paste0("M", mod_labels)
  names(mod_labels) <- names(obj_list[[x]]$modules_vec)
  
  stopifnot(all(unlist(lapply(net_list, function(x) sum(is.na(x) | is.infinite(x)))) == 0))
  stopifnot(all(unlist(lapply(data_list, function(x) sum(is.na(x) | is.infinite(x)))) == 0))
  stopifnot(all(unlist(lapply(adj_list, function(x) sum(is.na(x) | is.infinite(x)))) == 0))
  
  cat("\t--> Running Module preservation across diets...\n")
  set.seed(1234)
  pres_res <- NetRep::modulePreservation(network           = net_list,
                                         data              = data_list,
                                         correlation       = adj_list,
                                         moduleAssignments = mod_labels,
                                         backgroundLabel   = "M0",
                                         discovery         = x,                              # reference
                                         test              = ifelse(x == "CD", "HFD", "CD"), # query
                                         selfPreservation  = F,
                                         nPerm             = nb_perm,
                                         nThreads          = nbCores,
                                         verbose           = F)
  
  cat("\t--> Saving results...\n")
  pres_res <- list(pres_res      = pres_res,
                   ref_diet      = x,
                   test_diet     = ifelse(x == "CD", "HFD", "CD"),
                   nb_perm       = nb_perm,
                   net_list      = net_list,
                   data_list     = data_list,
                   adj_list      = adj_list,
                   mod_labels    = mod_labels)
  saveRDS(pres_res, out_file)
  
  cat("\t--> Done!\n\n")
})

remove(lipid.avg.df, tt)




#############################################
## Functional enrichment (module lipids)
#############################################

min_mod_size   <- 5
wgcna_res_list <- readRDS(paste0("./Data/lipidomics/wgcna_coexpression_analysis/wgcna_results_list__avg_data__min_mod_size_", min_mod_size, ".RDS"))
lipids_list    <- lapply(wgcna_res_list$wgcna_results_list, function(x){
  out <- data.frame(lipid  = names(x$wgcna_net$colors),
                    module = paste0("M", unname(x$wgcna_net$colors)), stringsAsFactors = F)
  split(out$lipid, out$module)
})
lipid_sets     <- readRDS("./Data/lipidomics/GSEA/lipid_gene_sets.RDS")
remove(wgcna_res_list)


saveDir <- paste0("./Data/lipidomics/wgcna_coexpression_functional_enrichment__lipid_annotations/min_mod_size_", min_mod_size)
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

overwrite <- F
nbCores   <- 26
tt        <- lapply(names(lipids_list), function(x){
  # x <- names(lipids_list)[1]
  cat("Running ORA for lipid modules in:", x, "...\n")
  
  out_file <- paste0(saveDir, "/ORA_results__lipid_annotations__", x, "__min_mod_size_", min_mod_size, ".RDS")
  if(file.exists(out_file) & !overwrite){
    return()
  }
  
  universe_lipids <- unname(unique(unlist(lipids_list[[x]])))
  ora_list        <- pbmcapply::pbmclapply(mc.cores = nbCores, X = names(lipids_list[[x]]), FUN = function(y){
    # y <- names(lipids_list[[x]])[1]
    
    ora_list <- lapply(names(lipid_sets), function(k){
      # k <- names(lipid_sets)[1]
      # print(k)
      
      ora_list <- lapply(names(lipid_sets[[k]]), function(z){
        # z <- names(lipid_sets[[k]])[1]
        # print(z)
        
        if(length(lipids_list[[x]][[y]]) == 0){
          return()
        }
        
        ora   <- clusterProfiler::enricher(gene          = unique(lipids_list[[x]][[y]]),
                                           pvalueCutoff  = 1,
                                           pAdjustMethod = "BH",
                                           universe      = universe_lipids,
                                           qvalueCutoff  = 1,
                                           minGSSize     = 4,
                                           TERM2GENE     = dplyr::select(lipid_sets[[k]][[z]], c("gs_name", "gene_symbol")))
        
        if(is.null(ora)){
          return()
        }
        ora_df   <- ora@result
        if(nrow(ora_df) == 0){
          return()
        }
        termsMeta  <- unique(dplyr::select(lipid_sets[[k]][[z]], dplyr::all_of(c("gs_name"))))
        ora_df     <- merge(ora_df, termsMeta, by.x = "ID", by.y = "gs_name", all.x = T, all.y = F)
        ora_df     <- ora_df[order(ora_df$pvalue), ]
        
        rownames(ora_df)         <- ora_df$ID
        ora_df$geneSetCollection <- z
        ora_df$genes_set_type    <- k
        ora_df$condition         <- x
        ora_df$wgcna_module      <- y
        ora_df$nb_tested_lipids  <- length(unique(lipids_list[[x]][[y]]))
        return(ora_df)
      })
      if(all(unlist(lapply(ora_list, is.null)))){
        return()
      }
      out.df <- unique(rbindlist(ora_list, use.names = T, fill = T))
      out.df
    })
    if(all(unlist(lapply(ora_list, is.null)))){
      return()
    }
    out.df <- unique(rbindlist(ora_list, use.names = T, fill = T))
    out.df
  })
  if(all(unlist(lapply(ora_list, is.null)))){
    return()
  }
  ora_df_all <- unique(rbindlist(ora_list, use.names = T, fill = T))
  
  ora_df_all$GeneRatio_num <- unlist(lapply(ora_df_all$GeneRatio, function(zz) eval(parse(text = zz))))
  ora_df_all$BgRatio_num   <- unlist(lapply(ora_df_all$BgRatio, function(zz) eval(parse(text = zz))))
  
  saveRDS(ora_df_all, out_file)
  
})



