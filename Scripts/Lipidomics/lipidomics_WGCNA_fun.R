


run_wgcna_analysis <- function(omics_df_list,
                               r2_thresh             = 0.85,
                               min_mod_size          = 3,
                               pam_respects_dendro   = T,
                               nb_preservation_perm  = 500,
                               tom_dir               = "./tmp_tom",
                               do_stability_analysis = T,
                               nbCores               = 40,
                               verbose               = F){
  

  
  if(!dir.exists(tom_dir)){
    dir.create(tom_dir, recursive = T)
  }
  
  
  library(ggplot2)
  library(WGCNA)
  enableWGCNAThreads(nThreads = nbCores)
  
  
  # explanations about parameters and WGCNA procedure:
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC12457846/
  # https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/wgcna.html#gsc.tab=0
  # https://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html
  # Differential expression across preserved modules (not yet implemented here):
  # https://bioc.r-universe.dev/articles/multiWGCNA/autism_full_workflow.html
  
  wgcna_results_list <- lapply(names(omics_df_list), function(omic_layer){
    # omic_layer <- names(omics_df_list)[1]

    
    cat("*******************************************************************************\n\t>Analyizing omic layer", 
        which(omic_layer == names(omics_df_list)), "/", length(omics_df_list) , "-", omic_layer, "...\n",
        "*******************************************************************************\n")
    
    cat("\t--> Identifying soft threshold...\n")
    
    # 1) Choose powers
    # identify the most suitable β power that makes the network satisfy the scale-free topology
    # ChatGPT:
    # > corFnc = "bicor" — biweight midcorrelation (robust) rather than Pearson. Why bicor? Lipidomics 
    # often has outliers and non-normal distributions. bicor down-weights outliers and produces more stable 
    # co-expression estimates. If your data is extremely clean and normally distributed, Pearson (cor) can 
    # be used, but bicor is generally safer.
    sft <- WGCNA::pickSoftThreshold(omics_df_list[[omic_layer]], 
                                    powerVector = 1:20, 
                                    verbose     = ifelse(verbose, 5, 0), 
                                    networkType = "signed", 
                                    blockSize   = 25000, 
                                    corFnc      = "bicor",
                                    corOptions  = list(use = 'pairwise.complete.obs'),
                                    RsquaredCut = r2_thresh,
                                    moreNetworkConcepts = T)
    
    
    cat("\n\t--> Plotting soft threshold selection...\n")
    dfPlot <- data.frame(power = sft$fitIndices$Power,
                         fit   = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
                         meank = sft$fitIndices$mean.k.,
                         stringsAsFactors = F)
    dfPlot <- dfPlot[dfPlot$fit > 0, ]
    
    pl.scale.ind <- ggplot(dfPlot, aes(x = power, y = fit)) +
      geom_point(color = "gray10") +
      geom_hline(yintercept = r2_thresh, color = "brown3") +
      ggrepel::geom_text_repel(aes(label = power), size = 2.3) +
      xlab("Soft threshold (power)") + ylab(expression(paste("Scale-free topology fit - ", R^{2}))) +
      ggtitle("Scale independence") +
      theme_classic() +
      theme(legend.position  = "none",
            panel.grid.major = element_line(),
            panel.grid.minor = element_line())
    pl.scale.ind
    
    
    pl.mean.conn <- ggplot(dfPlot, aes(x = power, y = meank)) +
      geom_point(color = "gray10") +
      ggrepel::geom_text_repel(aes(label = power), size = 2.3) +
      xlab("Soft threshold (power)") + ylab("Mean connectivity (k)") +
      ggtitle("Scale independence") +
      theme_classic() +
      theme(legend.position  = "none",
            panel.grid.major = element_line(),
            panel.grid.minor = element_line())
    pl.mean.conn
    
    wgcna_soft_power <- sft$powerEstimate
    
    cat("\n\t--> Building WGCNA network...\n")
    # > Meaning: when refining module membership with PAM (Partitioning Around Medoids), 
    # keep PAM decisions consistent with the hierarchical dendrogram — i.e., don't reassign 
    # genes across distant branches.
    # Why TRUE? It preserves the dendrogram structure and increases interpretability. 
    # FALSE allows more aggressive reassignment but may produce less interpretable modules.
    net <- WGCNA::blockwiseModules(omics_df_list[[omic_layer]], 
                                   power             = wgcna_soft_power,
                                   networkType       = "signed",
                                   TOMType           = "signed", 
                                   minModuleSize     = min_mod_size,
                                   numericLabels     = TRUE, 
                                   pamRespectsDendro = pam_respects_dendro, 
                                   saveTOMs          = T,
                                   saveTOMFileBase   = paste0(tom_dir, "/TOM_mat__", omic_layer),
                                   deepSplit         = 2, # the default (previously tested with different values)
                                   maxBlockSize      = 50000, 
                                   corType           = "bicor", 
                                   verbose           = ifelse(verbose, 3, 0), 
                                   randomSeed        = 123,
                                   nThreads          = 0)
  
    
    net$soft_threshold     <- sft
    net$power              <- wgcna_soft_power
    net$minModuleSize      <- min_mod_size
    net$corType            <- "bicor"
    net$data_obj           <- omics_df_list[[omic_layer]]
    

    cat("\n\t--> Plotting number of lipids per cluster...\n")
    
    dfPlot <- data.frame(lipid    = names(net$colors),
                         cluster  = as.character(unname(net$colors)),
                         stringsAsFactors = F)
    
    dfPlot           <- data.table(dfPlot)[, list(nb_lipids = length(unique(lipid))), by = c("cluster")]
    dfPlot           <- type.convert(as.data.frame(dfPlot), as.is = T)
    dfPlot           <- dfPlot[order(dfPlot$nb_lipids, -1*as.integer(dfPlot$cluster), decreasing = T), ]
    dfPlot$cluster   <- factor(dfPlot$cluster, levels = dfPlot$cluster)
    dfPlot$label_pos <- dfPlot$nb_lipids + max(dfPlot$nb_lipids) * 0.03
    
    pl.nbLipid.cluster <- ggplot(dfPlot, aes(x = cluster, y = nb_lipids)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = nb_lipids, y = label_pos), size = 2.2) +
      xlab("Nb. genes") + ylab("Cluster") +
      scale_y_discrete(expand = expand_scale(mult = c(0, 0.05))) +
      coord_flip() +
      theme_classic()
    # pl.nbLipid.cluster
    
    cat("\n\t--> Extracting lipids module membership...\n")
    
    ## compute kME, or signed eigengene-based connectivity, also known as module 
    ## membership
    ## highest membership <-> higest "level" of being a hub gene
    kMEs <- WGCNA::signedKME(datExpr = omics_df_list[[omic_layer]],
                             datME   = net$MEs,
                             corFnc  = "bicor")
    
    modules                    <- as.character(unique(unname(net$colors)))
    # modules                    <- modules[modules != "0"]
    mod_membership_list        <- lapply(modules, function(x){
      # x <- modules[1]
      genes.mod.x       <- names(net$colors[net$colors == x])
      kMEs.x            <- kMEs[rownames(kMEs) %in% genes.mod.x, ]
      kMEs.x.vec        <- kMEs.x[, paste0('kME', x)]
      names(kMEs.x.vec) <- rownames(kMEs.x)
      kMEs.x.vec        <- kMEs.x.vec[order(kMEs.x.vec, decreasing = T)]
      kMEs.x.vec
    })
    names(mod_membership_list) <- modules
    
    net$module_membership <- mod_membership_list
    
    # perform module stability analysis
    # note: sampledBlockwiseModules can do it but is very slow. I reimplemented a 
    # simple loop doing it faster in parallel.
    # Inspired from:
    # https://github.com/cran/WGCNA/blob/c8c3d8bc1204d6866c1157ff12898d73098a031f/R/sampledModules.R#L9
    if(do_stability_analysis){
      
      strains_vec   <- unique(gsub("_.*", "", rownames(omics_df_list[[omic_layer]])))
      
      cat("\n\t--> Performing clusters stability analysis...\n")
      cat("\n\t\t--> Computing leave-one-strain-out permutations ( N = ", length(strains_vec), ")...\n")
      
      nets.list     <- parallel::mclapply(mc.cores = min(nbCores, length(strains_vec)), X = strains_vec, FUN = function(x){
        # nets.list     <- lapply(strains_vec, function(x){
        
        set.seed(which(x == strains_vec))
        omics_df.x <- omics_df_list[[omic_layer]][x != strains_vec, ]
        stopifnot(nrow(omics_df.x) == (nrow(omics_df_list[[omic_layer]]) - 1))
        
        net.x <- WGCNA::blockwiseModules(omics_df_list[[omic_layer]], 
                                         power             = wgcna_soft_power,
                                         networkType       = "signed",
                                         TOMType           = "signed", 
                                         minModuleSize     = min_mod_size,
                                         numericLabels     = TRUE, 
                                         pamRespectsDendro = pam_respects_dendro, 
                                         deepSplit         = 2, # the default (previously tested with different values)
                                         maxBlockSize      = 50000, 
                                         corType           = "bicor", 
                                         verbose           = 0, 
                                         randomSeed        = which(x == strains_vec),
                                         nThreads          = 0)
        
        
        list(net_obj         = net.x,
             strain_removed  = x)
      })
      
      cat("\n\t\t--> Building lipid clusters map...\n")
      
      # Create a matrix of labels for the original and all resampling runs
      labs.df.ref <- data.frame(lipid                    = names(net$colors),
                                matched_cluster_original = unname(net$colors),
                                stringsAsFactors = F)
      
      ## Relabel modules in each of the resampling runs
      labs.list <- parallel::mclapply(mc.cores = nbCores, X = 1:length(nets.list), FUN = function(x){
        # x <- 1
        df.x <- data.frame(lipid           = names(nets.list[[x]]$net_obj$colors),
                           matched_cluster = WGCNA::matchLabels(nets.list[[x]]$net_obj$colors, net$colors, pThreshold = 1e-2),
                           stringsAsFactors = F)
        colnames(df.x)[2] <- ifelse(length(nets.list[[x]]$strain_removed) == 1, 
                                    paste0(colnames(df.x)[2], "_", x, "___", nets.list[[x]]$strain_removed),
                                    paste0(colnames(df.x)[2], "_", x))
        
        df.x
      })
      labs.df.all <- Reduce(function(...) merge(..., by = "lipid", all = T), labs.list)
      labs.df.all <- merge(labs.df.ref, labs.df.all, by = "lipid")
      
      cat("\n\t\t--> Building resampling overlap table...\n")
      
      overlap.df <- lapply(colnames(labs.df.all)[grepl("matched_cluster_", names(labs.df.all))], function(x){
        # x <- colnames(labs.df.all)[grepl("matched_cluster_", names(labs.df.all))][1]
        
        out <- lapply(colnames(labs.df.all)[grepl("matched_cluster_", names(labs.df.all))], function(y){
          # y <- colnames(labs.df.all)[grepl("matched_cluster_", names(labs.df.all))][2]
          
          overlap.perc <- dplyr::select(labs.df.all, c(x, y))
          overlap.perc <- (sum(overlap.perc[[x]] == overlap.perc[[y]]) / nrow(overlap.perc)) * 100
          overlap.perc <- data.frame(cluster1     = x,
                                     cluster2     = y,
                                     perc_overlap = overlap.perc,
                                     stringsAsFactors = F)
          overlap.perc
        })
        do.call(rbind, out)
        
      })
      overlap.df <- do.call(rbind, overlap.df)
      
      
      cat("\n\t\t--> Plotting resampling clusters stability heatmap...\n")
      
      dfPlot <- labs.df.all
      dfPlot <- reshape2::melt(dfPlot, id.vars = c("lipid"))
      dfPlot <- type.convert(dfPlot, as.is = T)
      dfPlot$resampling_id <- gsub("matched_cluster_original", "Original", dfPlot$variable)
      dfPlot$resampling_id <- gsub("matched_cluster_", "Resampling ", dfPlot$resampling_id)
      dfPlot$resampling_id[grepl("___", dfPlot$resampling_id)] <- paste0(dfPlot$resampling_id[grepl("___", dfPlot$resampling_id)], ")")
      dfPlot$resampling_id <- gsub("___", " (-", dfPlot$resampling_id)
      resampling_id_order  <- unique(dfPlot$resampling_id)
      resampling_id_order  <- c("Original", stringr::str_sort(resampling_id_order[resampling_id_order != "Original"], numeric = T))
      dfPlot$resampling_id <- factor(dfPlot$resampling_id, levels = resampling_id_order)
      tmp.df               <- dfPlot[dfPlot$resampling_id == "Original", ]
      tmp.df               <- tmp.df[order(tmp.df$value), ]
      dfPlot$lipid         <- factor(dfPlot$lipid, levels = tmp.df$lipid)
      dfPlot$cluster       <- paste0("cluster_", dfPlot$value)
      dfPlot$cluster       <- factor(dfPlot$cluster, levels = stringr::str_sort(unique(dfPlot$cluster), numeric = T))
      
      pl.stability.all <- ggplot(dfPlot, aes(x = lipid, y = resampling_id, fill = cluster)) +
        geom_tile() +
        scale_fill_manual(values = pals::brewer.dark2(length(unique(dfPlot$cluster)))) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_classic() +
        theme(legend.position = "none",
              axis.text.x     = element_blank(),
              axis.ticks.x    = element_blank(),
              axis.title.y    = element_blank())
      
      # tmp.df <- data.table(dfPlot)[, list(nb_gene = length(unique(gene))), by = "cluster"]
      # tmp.df <- tmp.df[order(tmp.df$nb_gene, decreasing = T), ]
      
      cat("\n\t\t--> Plotting resampling clusters overlap heatmap...\n")
      
      dfPlot <- overlap.df
      dfPlot$cluster1 <- gsub("matched_cluster_original", "Original", dfPlot$cluster1)
      dfPlot$cluster1 <- gsub("matched_cluster_", "Resampling ", dfPlot$cluster1)
      dfPlot$cluster1[grepl("___", dfPlot$cluster1)] <- paste0(dfPlot$cluster1[grepl("___", dfPlot$cluster1)], ")")
      dfPlot$cluster1 <- gsub("___", " (-", dfPlot$cluster1)
      
      dfPlot$cluster2 <- gsub("matched_cluster_original", "Original", dfPlot$cluster2)
      dfPlot$cluster2 <- gsub("matched_cluster_", "Resampling ", dfPlot$cluster2)
      dfPlot$cluster2[grepl("___", dfPlot$cluster2)] <- paste0(dfPlot$cluster2[grepl("___", dfPlot$cluster2)], ")")
      dfPlot$cluster2 <- gsub("___", " (-", dfPlot$cluster2)
      
      tmp.df <- dfPlot[dfPlot$cluster1 != "Original" & dfPlot$cluster2 != "Original", ]
      tmp.df <- reshape2::dcast(tmp.df, cluster1 ~ cluster2, value.var = "perc_overlap")
      rownames(tmp.df) <- tmp.df[, 1]
      tmp.df <- tmp.df[, -1]
      hc     <- hclust(as.dist(tmp.df), method = "ward.D2")
      
      resampling_id_order <- c("Original", hc$labels[hc$order])
      dfPlot$cluster1     <- factor(dfPlot$cluster1, levels = resampling_id_order)
      dfPlot$cluster2     <- factor(dfPlot$cluster2, levels = resampling_id_order)
      dfPlot              <- dfPlot[dfPlot$cluster1 != dfPlot$cluster2, ]
      
      pl.overlap <- ggplot(dfPlot, aes(x = cluster1, y = cluster2, fill = perc_overlap)) +
        geom_tile() +
        scale_fill_gradientn(colors = pals::coolwarm(200)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title  = element_blank()) +
        coord_equal() +
        labs(fill = "% overlap")
      
    } else{
      
      cat("\n\t--> Clusters stability analysis skipped...\n")
      
      nets.list        <- NULL
      labs.list        <- NULL
      overlap.df       <- NULL
      pl.stability.all <- NULL
      pl.overlap       <- NULL
      
    }
    
    
    
    
    # cat("\n\tDone!\n")
    
    list(wgcna_net                               = net,
         stability_analysis_nets_list            = nets.list,
         stability_analysis_matched_clust        = labs.list,
         stability_analysis_overlap_table        = overlap.df,
         
         sft_scale_topology_plot                 = pl.scale.ind,
         sft_mean_connectivity_plot              = pl.mean.conn,
         # nb_genes_cluster_plot                   = pl.nbGene.cluster,
         stability_analysis_matched_heatmap_plot = pl.stability.all,
         stability_analysis_overlap_plot         = pl.overlap)
    
  })
  names(wgcna_results_list) <- names(omics_df_list)

  
  cat("\n\t--> Building modules correspondance tables...\n")
  
  if(length(wgcna_results_list) > 1){
    
    ref_modules_vec <- wgcna_results_list[[1]]$wgcna_net$colors
    match_table     <- lapply(names(wgcna_results_list)[-1], function(zz){
      # zz <- names(wgcna_results_list)[-1][1]
      
      out.zz <- WGCNA::matchLabels(wgcna_results_list[[zz]]$wgcna_net$colors, ref_modules_vec, pThreshold = 1e-2)
      out.zz <- as.data.frame(out.zz)
      colnames(out.zz) <- paste0("module_source__", zz)
      out.zz
    })
    tmp.df <- data.frame(lipid  = names(ref_modules_vec),
                         module = unname(ref_modules_vec), stringsAsFactors = F)
    colnames(tmp.df)[2] <- paste0("module_ref__", names(wgcna_results_list)[1])
    
    match_table <- cbind(tmp.df, do.call(cbind, match_table))
    
  } else{
    
    match_table <- NULL
    
  }
  
  
  
  
  
  
  cat("\n\t--> Modules preservation analysis...\n")
  
  # > Z-summary score:  measures how well modules from the reference dataset are preserved 
  #                     in the test dataset. Higher Z-summary scores indicate stronger preservation.
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC12457846/
  
  # compute using as reference each of the inpput omic layers
  preservation_results.list <- lapply(names(wgcna_results_list), function(omic_layer){
    # zz <- names(wgcna_results_list)[1]
    
    cat("\n\t\t--> Running for omic layer as reference:", omic_layer, "...\n")
    
    ref_omic_name.use    <- omic_layer
    preservation_results <- WGCNA::modulePreservation(multiData           = lapply(wgcna_results_list, function(x) list(data = x$wgcna_net$data_obj)),
                                                      # multiColor          = lapply(wgcna_results_list, function(x) labels2colors(x$wgcna_net$colors)),
                                                      multiColor          = lapply(wgcna_results_list, function(x) x$wgcna_net$colors),
                                                      referenceNetworks   = which(names(wgcna_results_list) == ref_omic_name.use),
                                                      networkType         = "signed",
                                                      corFnc              = "bicor",
                                                      corOptions          = "use = 'pairwise.complete.obs'",
                                                      goldName            = "random_sampling_control",
                                                      greyName            = "unassigned_features_module",
                                                      nPermutations       = nb_preservation_perm,
                                                      randomSeed          = 1234,
                                                      parallelCalculation = F,
                                                      verbose             = ifelse(verbose, 5, 0))
    
    preservation_stats <- preservation_results$preservation$Z
    # (1) Remove the "gold" module, renamed here as "random_sampling_control".
    # The "gold" module is an artificial reference module automatically generated by 
    # WGCNA's modulePreservation() function as a technical control. It is not a biological 
    # module identified from your data. Instead, it contains randomly selected genes used to 
    # establish baseline preservation statistics and represents the null expectation for module 
    # preservation. This module serves as a statistical control, showing what preservation 
    # statistics would look like for a randomly assembled set of genes with no genuine biological 
    # relationships. By design, the gold module should show poor preservation because it lacks 
    # biological coherence. Since it is artificially constructed rather than biologically derived, 
    # it provides no meaningful biological insights and must be removed before proceeding with 
    # further analysis.
    # Z-summary < 2: No preservation.
    # 2 ≤ Z-summary < 10: Moderate preservation.
    # Z-summary ≥ 10: High preservation.
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC12457846/
    # (2) then plot
    preservation_stats_df <- lapply(names(preservation_stats), function(x){
      # x <- names(preservation_stats)[1]
      idx_test    <- which(names(preservation_stats[[x]]) != paste0("inColumnsAlsoPresentIn.", ref_omic_name.use))
      out         <- lapply(idx_test, function(y){
        # y <- idx_test[1]
        
        test_id.y   <- gsub("inColumnsAlsoPresentIn.", "", names(preservation_stats[[x]])[y], fixed = T)
        ref_id.y    <- gsub("ref\\.", "", x)
        stats.y     <- preservation_stats[[x]][[y]]
        stats.y     <- stats.y[!(rownames(stats.y) %in% c("random_sampling_control")), ]
        tmp.df      <- data.frame(module = unname(wgcna_results_list[[ref_omic_name.use]]$wgcna_net$colors),
                                  color  = labels2colors(wgcna_results_list[[ref_omic_name.use]]$wgcna_net$colors))
        tmp.df.aggr <- data.table(tmp.df)[, list(mod_size_ref = length(color)), by = "module"]
         
        stats.y     <- cbind("ref_module" = rownames(stats.y),
                             "ref_id"     = ref_id.y,
                             "test_id"    = test_id.y,
                             stats.y)
        
        pval.df.y  <- preservation_results$preservation$log.p[[x]][[y]]
        stats.y    <- merge(stats.y, 
                            dplyr::select(pval.df.y, -dplyr::all_of(colnames(pval.df.y)[colnames(pval.df.y) %in% colnames(stats.y)])),
                            by.x = "ref_module", by.y = "row.names")
        # rows:    test module names
        # columns: reference module names
        match.df.y    <- preservation_results$accuracy$observedCounts[[x]][[y]]
        nb_mod_ref    <- length(unique(wgcna_results_list[[ref_id.y]]$wgcna_net$colors))
        nb_mod_test   <- length(unique(wgcna_results_list[[test_id.y]]$wgcna_net$colors))
        stopifnot(nrow(match.df.y) == nb_mod_test & ncol(match.df.y) == nb_mod_ref)
        idx.best_ov   <- apply(match.df.y, 2, function(zz) unname(which(zz == max(zz))))
        size.best_of  <- apply(match.df.y, 2, function(zz) unname(zz[zz == max(zz)]))
        match.df.y    <- data.table(ref_module   = colnames(match.df.y),
                                    test_module  = unlist(lapply(idx.best_ov, function(zz) paste(rownames(match.df.y)[zz], collapse = ";"))),
                                    size_overlap = unlist(lapply(size.best_of, function(zz) paste(zz, collapse = ";"))), stringsAsFactors = F)
        # match.df.y$test_module <- plyr::mapvalues(match.df.y$test_id, from = test_col$color, to = test_col$module, warn_missing = F)
        stats.y       <- merge(stats.y, match.df.y, by = "ref_module")
        
        stopifnot(all(stats.y$moduleSize == plyr::mapvalues(stats.y$ref_module, from = tmp.df.aggr$module, to = tmp.df.aggr$mod_size_ref, warn_missing = F)))
        colnames(stats.y)[colnames(stats.y) == "mod_color"]  <- "ref_color"
        colnames(stats.y)[colnames(stats.y) == "moduleSize"] <- "ref_module_size"
        
        # cols_order <- c("ref_id", "test_id", "ref_color", "ref_module", "ref_module_size", "test_color", "test_module", "size_overlap")
        cols_order <- c("ref_id", "test_id", "ref_module", "ref_module_size", "test_module", "size_overlap")
        cols_order <- c(cols_order, colnames(stats.y)[!(colnames(stats.y) %in% cols_order)])
        stats.y    <- dplyr::select(stats.y, dplyr::all_of(cols_order))
        stats.y    <- type.convert(stats.y, as.is = T)
        
        stats.y
      })
      do.call(rbind, out)
    })
    preservation_stats_df <- do.call(rbind, preservation_stats_df)
    
   
    
    list(preservation_results = preservation_results,
         preservation_stats   = preservation_stats_df,
         preservation_plot    = pl.preserv)
  })
  names(preservation_results.list) <- paste0(names(wgcna_results_list), "_reference")
  
  preservation_plots <- lapply(preservation_results.list, function(x) x$preservation_plot)
  preservation_res   <- lapply(preservation_results.list, function(x) x[!(names(x) %in% c("preservation_plot"))])
  
  cat("\n\t--> Consensus modules analysis...\n")
  
  net_consensus <- WGCNA::blockwiseConsensusModules(multiExpr         = lapply(wgcna_results_list, function(x) list(data = x$wgcna_net$data_obj)),
                                                    power             = lapply(wgcna_results_list, function(x) x$wgcna_net$power),
                                                    networkType       = "signed",
                                                    TOMType           = "signed", 
                                                    minModuleSize     = min_mod_size,
                                                    numericLabels     = TRUE, 
                                                    pamRespectsDendro = pam_respects_dendro, 
                                                    deepSplit         = 2, # the default (previously tested with different values)
                                                    maxBlockSize      = 50000, 
                                                    corType           = "bicor", 
                                                    verbose           = ifelse(verbose, 3, 0), 
                                                    randomSeed        = 1234,
                                                    nThreads          = 0)
  
  consTree             <- net_consensus$dendrograms[[1]]
  mod_cd_match         <- matchLabels(wgcna_results_list$CD$wgcna_net$colors, net_consensus$colors)
  mod_hfd_match        <- matchLabels(wgcna_results_list$HFD$wgcna_net$colors, net_consensus$colors)
  mod_cd_match_colors  <- labels2colors(mod_cd_match)
  mod_hfd_match_colors <- labels2colors(mod_hfd_match)
  mod_cons_colors      <- labels2colors(net_consensus$colors)
  colors_df            <- data.frame(mod_cons_colors, mod_cd_match_colors, mod_hfd_match_colors)[net_consensus$blockGenes[[1]], ]
  
  
  cat("\n\tDone!\n\n")
  
  return(list(wgcna_results_list = wgcna_results_list,
              modules_match_df   = match_table,
              preservation_res   = preservation_res,
              preservation_plots = preservation_plots,
              consensus_net      = net_consensus))
  
}