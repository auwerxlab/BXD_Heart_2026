
library(data.table)
library(limma)


source("./Scripts/Lipidomics/lipidomics_fun.R")


################################################################################
## DEA with limma
################################################################################


###########################
## load and prepare data
###########################

lipid.derived.list  <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_with_derived_features.RDS")


features.list  <- lipid.derived.list[c("Filt_MassNorm_BatchNorm_Data_NArm")]
rownames(features.list$Filt_MassNorm_BatchNorm_Data_NArm) <- features.list$Filt_MassNorm_BatchNorm_Data_NArm$Identifier
features.list  <- lapply(features.list, function(x){
  x[!grepl("unknown", rownames(x), ignore.case = T), !grepl("PWT|Blank", colnames(x), ignore.case = T) & grepl("BXD|C57|DBA", colnames(x))]
})

meta.df            <- lipid.derived.list$MetaData
meta.df$strainDiet <- paste0(meta.df$Strain, "_", meta.df$Diet)


meta.lipids.df <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_metadata.RDS")
meta.lipids.df <- meta.lipids.df[!grepl("unknown", meta.lipids.df$Identifier, ignore.case = T), ]

###########################
## run DEA
###########################

saveDir <- "./Data/lipidomics/limma_DEA"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}

overwrite         <- F
allTopTables.list <- lapply(names(features.list), function(type){
  # type <- names(features.list)[1]

  cat("\nPerforming lipidomics DEA for input data table:\t", type, "...\n")
  
  out_file <- paste0(saveDir, "/limma_lipidomics_DEA_results___input_table_", type, ".RDS")
  if(file.exists(out_file) & !overwrite){
    cat("\t--> Results file already present. Skipping...\n")
    return()
  }
  
  expr_table            <- features.list[[type]]
  meta_table            <- meta.df[meta.df$SampleName %in% colnames(expr_table), ]
  meta_table            <- meta_table[match(colnames(expr_table), meta_table$SampleName), ]
  meta_table$Strain     <- gsub("\\/", "_", meta_table$Strain)
  meta_table$strainDiet <- gsub("\\/", "_", meta_table$strainDiet)
  colnames(expr_table)  <- gsub("\\/", "_", colnames(expr_table))
  rownames(meta_table)  <- meta_table$SampleName
  stopifnot(identical(rownames(meta_table), colnames(expr_table)))
  
  contrast_types        <- c("population_diet_effect", "strain_diet_effect")
  
  
  # !!!! 
  # !!!! log2 normalization
  # !!!! 
  expr_table <- log2(expr_table + 1)
  
  
  names(contrast_types) <- contrast_types
  top_tables.list       <- lapply(contrast_types, function(x){
    # x <- contrast_types[1]
    
    cat("\t--> Running contrast type:\t", x, "...\n")
    
    if(x == "population_diet_effect"){
      
      # Which lipids respond to HFD on average?
      designFormula  <- as.formula(paste0("~0+Diet"))
      mm             <- model.matrix(designFormula, data = meta_table)
      colnames(mm)   <- gsub("Diet", "", colnames(mm))
      contrasts_def  <- list("HFD_vs_CD_diet" = "(HFD)-(CD)")
      
    } else if(x == "strain_diet_effect"){
      
      # How does strain X respond to diet?
      designFormula        <- as.formula(paste0("~0+strainDiet"))
      mm                   <- model.matrix(designFormula, data = meta_table)
      colnames(mm)         <- gsub("strainDiet", "", colnames(mm))
      contrasts_def        <- lapply(unique(meta_table$Strain), function(y) paste0("(", y, "_HFD)-(", y, "_CD)"))
      names(contrasts_def) <- paste0("HFD_vs_CD_strain_", unique(meta_table$Strain))
      
    } 
    
    contrasts <- lapply(contrasts_def, function(k){
      # k <- contrasts_def[1]
      eval(parse(text = paste0("out <- makeContrasts(", unname(k), ", levels = mm)")))
      colnames(out) <- names(k)
      out
    })
    
    corrfit    <- duplicateCorrelation(expr_table, block = meta_table$Strain, design = mm)
    fitModel   <- limma::lmFit(expr_table, design = mm, correlation = corrfit$consensus.correlation, block = meta_table$Strain)
    
    top_tables <- lapply(1:length(contrasts), function(k){
      # k <- 1
      table.k             <- limma::topTable(limma::eBayes(limma::contrasts.fit(fitModel, contrasts = contrasts[[k]])), sort.by = "none", n = Inf)
      table.k$feature_id  <- rownames(table.k)
      table.k$contrast_id <- contrasts_def[[k]]
      table.k$contrast    <- names(contrasts_def)[k]
      table.k             <- table.k[!is.na(table.k$logFC), ]
      if(any(table.k$feature_id %in% meta.lipids.df$Identifier)){
        table.k <- merge(table.k, meta.lipids.df, by.x = "feature_id", by.y = "Identifier", all.x = T, all.y = F)
      }
      table.k             <- table.k[order(table.k$adj.P.Val, -1*abs(table.k$logFC)), ]
      print(table(table.k$adj.P.Val < 0.05))
      table.k
    })
    names(top_tables) <- names(contrasts)
    top_tables
  })
  
  saveRDS(top_tables.list, out_file)
  
})



################################################################################
## GSEA using lipid DEA logFC as ranking
################################################################################

top_table_files        <- list.files("./Data/lipidomics/limma_DEA", full.names = T)
names(top_table_files) <- gsub(".*input_table_|\\.RDS", "", top_table_files, ignore.case = T)
top_tables.list        <- lapply(top_table_files, readRDS)


saveDir <- "./Data/lipidomics/GSEA"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}

if(!file.exists("./Data/lipidomics/GSEA/lipid_gene_sets.RDS")){
  
  # define lipid-class gene sets
  meta.lipids.df           <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_metadata.RDS")
  meta.lipids.df           <- meta.lipids.df[!grepl("unknown", meta.lipids.df$Identifier, ignore.case = T), ]
  meta.lipids.df$nb_chains <- paste0(meta.lipids.df$nb_chains, "_side_chains")
  meta.lipids.df$sat_ratio <- paste0("C", meta.lipids.df$cLength, ":", meta.lipids.df$nbDoubleBonds)
  
  cols_gs <- c("lipid_grandparent_class" = "Lipid.Grandparent.Class",
               "lipid_parent_class"      = "Lipid.Parent.Class", 
               "lipid_class"             = "Lipid.Class", 
               "lipid_specie_nb_chain"   = "nb_chains")
  
  
  gene_sets <- lapply(names(cols_gs), function(x){
    # x <- names(cols_gs)[1]
    
    lipids.gs.x           <- unique(dplyr::select(meta.lipids.df, dplyr::all_of(c(cols_gs[[x]], "Identifier"))))
    lipids.gs.x           <- unique(lipids.gs.x)
    colnames(lipids.gs.x) <- c("gs_name", "gene_symbol")
    lipids.gs.x$gs_class  <- x
    rownames(lipids.gs.x) <- NULL
    lipids.gs.x
  })
  names(gene_sets) <- names(cols_gs)
  
  # Remove duplicated items
  # As we specify first smaller gene set collections in cols_gs, the duplicated items will be kept in those smaller classes and removed from the largest classes.
  # If identical in all the subclasses, those features are indeed part of the superclass as well and should be kept there
  gene_sets <- do.call(rbind, gene_sets)
  tmp.df    <- data.table(gene_sets)[, list(features = paste(sort(unique(gene_symbol)), collapse = ";")), by = c("gs_name", "gs_class")]
  tmp.df    <- tmp.df[!duplicated(tmp.df$features), ]
  gene_sets <- gene_sets[paste0(gene_sets$gs_name, "_", gene_sets$gs_class) %in% paste0(tmp.df$gs_name, "_", tmp.df$gs_class), ]
  gene_sets <- split(gene_sets, gene_sets$gs_class)
  
  
  breaks <- seq(2, max(meta.lipids.df$cLength) + 10, by = 10)
  bins   <- cut(meta.lipids.df$cLengt, breaks = breaks, right = TRUE, include.lowest = TRUE)
  meta.lipids.df$cLength_binned <- sapply(bins, function(b) {
    rng <- as.numeric(gsub("\\(|\\[|\\]|\\)", "", unlist(strsplit(as.character(b), ","))))
    paste(rng[1] + 1, rng[2], sep = "_")
  })
  
  breaks <- c(seq(1, max(meta.lipids.df$nbDoubleBonds) + 2, by = 2))
  bins   <- cut(meta.lipids.df$nbDoubleBonds, breaks = breaks, right = TRUE, include.lowest = TRUE)
  meta.lipids.df$nbDoubleBonds_binned <- sapply(bins, function(b) {
    rng <- as.numeric(gsub("\\(|\\[|\\]|\\)", "", unlist(strsplit(as.character(b), ","))))
    paste(rng[1] + 1, rng[2], sep = "_")
  })
  meta.lipids.df$nbDoubleBonds_binned[meta.lipids.df$nbDoubleBonds == 0] <- "0"
  table(meta.lipids.df$nbDoubleBonds_binned)
  
  gene_sets.2 <- lapply(c("cLength", "cLength_binned", "nbDoubleBonds", "nbDoubleBonds_binned", "sat_ratio"), function(x){
    # x <- "cLength"
    # x <- "nbDoubleBonds"
    # x <- "sat_ratio"
    
    df.x.overall           <- unique(dplyr::select(meta.lipids.df, dplyr::all_of(c(x, "Identifier"))))
    colnames(df.x.overall) <- c("gs_name", "gene_symbol")
    df.x.overall$gs_class  <- paste0(x, "__all_lipids")
    rownames(df.x.overall) <- NULL
    df.x.overall           <- df.x.overall[order(df.x.overall$gs_name, df.x.overall$gene_symbol), ]
    
    df.x.classes <- lapply(c("Lipid.Grandparent.Class", "Lipid.Parent.Class", "Lipid.Class"), function(y){
      # y <- "Lipid.Grandparent.Class"
      
      df.y           <- meta.lipids.df
      df.y$gs_name   <- paste0(df.y[[x]], "__", df.y[[y]], "__", gsub("\\.", "_", tolower(y)))          
      df.y           <- unique(dplyr::select(df.y, dplyr::all_of(c("gs_name", "Identifier"))))
      colnames(df.y) <- c("gs_name", "gene_symbol")
      df.y$gs_class  <- paste0(x, "__", gsub("\\.", "_", tolower(y)))
      rownames(df.y) <- NULL
      df.y
    })
    df.x.classes <- do.call(rbind, df.x.classes)
    
    df.x <- unique(rbind(df.x.overall, df.x.classes))
    
    if(x == "cLength"){
      df.x$gs_name <- paste0("c_length_C", df.x$gs_name)
    } else if(x == "cLength_binned"){
      df.x$gs_name <- paste0("c_length_binned__C", df.x$gs_name)
    } else if(x == "nbDoubleBonds"){
      df.x$gs_name <- paste0("nb_double_bonds_", df.x$gs_name)
    } else if(x == "nbDoubleBonds_binned"){
      df.x$gs_name <- paste0("nb_double_bonds_binned_", df.x$gs_name)
    }
    df.x
  })
  # Remove duplicated items
  # As we specify first smaller gene set collections in cols_gs, the duplicated items will be kept in those smaller classes and removed from the largest classes.
  # If identical in all the subclasses, those features are indeed part of the superclass as well and should be kept there
  gene_sets.2 <- do.call(rbind, gene_sets.2)
  tmp.df      <- data.table(gene_sets.2)[, list(features = paste(sort(unique(gene_symbol)), collapse = ";")), by = c("gs_name", "gs_class")]
  tmp.df      <- tmp.df[!duplicated(tmp.df$features), ]
  gene_sets.2 <- gene_sets.2[paste0(gene_sets.2$gs_name, "_", gene_sets.2$gs_class) %in% paste0(tmp.df$gs_name, "_", tmp.df$gs_class), ]
  gene_sets.2 <- split(gene_sets.2, gsub("__.*", "", gene_sets.2$gs_class))
  
  
  lipid_fun_meta <- readRDS("./Data/lipidomics/data_processing/BXD_heart_lipidomics_metadata__with_LipidSigR_lipid_meta.RDS")
  lipid_fun_meta <- dplyr::select(lipid_fun_meta, dplyr::all_of(c("Identifier", "Cellular.Component", "Function", "Headgroup.Charge", "Lateral.Diffusion", 
                                                                  "Bilayer.Thickness", "Intrinsic.Curvature", "Transition.Temperature")))
  lipid_fun_meta           <- type.convert(reshape2::melt(lipid_fun_meta, id.vars = c("Identifier")), as.is = T)
  lipid_fun_meta           <- tidyr::separate_rows(lipid_fun_meta, c("value"), sep = "\\|")
  lipid_fun_meta           <- lipid_fun_meta[!is.na(lipid_fun_meta$value), ]
  lipid_fun_meta$variable  <- tolower(gsub("\\.", "_", lipid_fun_meta$variable))
  colnames(lipid_fun_meta) <- c("gene_symbol", "gs_class", "gs_name")
  lipid_fun_meta           <- unique(dplyr::select(lipid_fun_meta, dplyr::all_of(c("gs_name", "gene_symbol", "gs_class"))))
  # split(lipid_fun_meta, lipid_fun_meta$gs_class)
  
  gene_sets <- list(lipid_composition   = gene_sets,
                    lipid_composition_2 = gene_sets.2,
                    lipid_properties    = list(lipid_properties = lipid_fun_meta))
  
  saveRDS(gene_sets, paste0(saveDir, "/lipid_gene_sets.RDS"))
  
} else{
  
  gene_sets <- readRDS("./Data/lipidomics/GSEA/lipid_gene_sets.RDS")
  
}




nbCores   <- 60
overwrite <- F
gsea      <- lapply(names(gene_sets), function(nn){
  # nn <- names(gene_sets)[2]
  lapply(c("logFC_ranking"), function(s){
    # s <- "logFC_ranking"
    
    lapply(names(top_tables.list), function(y){
      # y <- names(top_tables.list)[1]
      # y <- "Filt_MassNorm_BatchNorm_Data_NArm"
      
      cat("Computing GSEA for:\n\t> data type DEA results:\t", y, "\n\t> ranking type:\t\t\t", s, "\n\t> gene sets type:\t\t", nn, "...\n")
      
      # lipid properties added later. To avoid re-running all the GSEA, save them in a different format
      saveDir <- paste0("./Data/lipidomics/GSEA/", s, ifelse(nn == "lipid_composition", "", paste0("/", nn)))
      if(!dir.exists(saveDir)){
        dir.create(saveDir)
      }
      
      out_file <- paste0(saveDir, "/GSEA_lipidomics_results__input_table_", y, ".RDS")
      if(file.exists(out_file) & !overwrite){
        cat("\t--> Results files already present. Skipping...\n")
        return()
      }
      
      out.list <- lapply(names(top_tables.list[[y]]), function(z){
        # z <- names(top_tables.list[[y]])[1]
        # print(z)
        cat("\t--> Analysing DEA results for contrast type:\t", z, "...\n")
        
        out.list <- parallel::mclapply(mc.cores = nbCores, X = top_tables.list[[y]][[z]], FUN = function(k){
          # k <- top_tables.list[[y]][[z]][[1]]
          
          out.list <- lapply(names(gene_sets[[nn]]), function(x){
            # x <- names(gene_sets[[nn]])[1]
            
            tmp             <- k
            if(sum(tmp$feature_id %in% meta.lipids.df$Identifier) == 0){
              return()
            }
            tmp             <- tmp[!duplicated(tmp$feature_id) & !is.na(tmp$feature_id), ]
            tmp             <- tmp[!is.na(tmp$logFC), ]
            
            tmp             <- tmp[order(tmp$logFC, decreasing = T), ]
            geneList        <- tmp$logFC
            names(geneList) <- tmp$feature_id
            
            gsea <- clusterProfiler::GSEA(geneList      = geneList,
                                          nPermSimple   = 100000,
                                          pvalueCutoff  = 1,
                                          maxGSSize     = 800,
                                          pAdjustMethod = "BH",
                                          TERM2GENE     = dplyr::select(gene_sets[[nn]][[x]], dplyr::all_of(c("gs_name", "gene_symbol"))))
            gsea_df <- gsea@result
            if(nrow(gsea_df) == 0){
              return()
            }
            termsMeta  <- unique(dplyr::select(gene_sets[[nn]][[x]], dplyr::all_of(c("gs_name"))))
            gsea_df    <- merge(gsea_df, termsMeta, by.x = "ID", by.y = "gs_name", all.x = T, all.y = F)
            gsea_df    <- gsea_df[order(gsea_df$pvalue), ]
            
            rownames(gsea_df)         <- gsea_df$ID
            gsea_df$geneSetCollection <- x
            gsea_df$contrast_id       <- tmp$contrast_id[1]
            gsea_df$contrast          <- tmp$contrast[1]
            gsea_df$ranking_type      <- s
            
            gsea_df
            
          })
          if(all(unlist(lapply(out.list, is.null)))){
            return()
          }
          out.df           <- unique(rbindlist(out.list, use.names = T, fill = T))
          out.df$coreSize  <- unlist(lapply(out.df$core_enrichment, function(x) length(unlist(str_split(x, "/")))))
          out.df$geneRatio <- out.df$coreSize / out.df$setSize
          out.df           <- as.data.frame(out.df)
          out.df
        })
      })
      # tt <- out.list[[1]]$HFD_vs_CD_diet
      names(out.list) <- names(top_tables.list[[y]])
      saveRDS(out.list, out_file)
      
    })
  })
})







