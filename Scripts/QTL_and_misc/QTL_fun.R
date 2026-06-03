
library(parallel)


transform2normal_byCol_INT <- function(df){
  # df <- df.intensity
  # df <- df.counts

  df.tmp <- df
  transf <- lapply(1:ncol(df.tmp), function(x){
    # x <- 1
    vec   <- as.numeric(df.tmp[, x])
    vec.t <- stats::qnorm(base::rank(vec, na.last = "keep", ties.method = "average")/(length(vec[!is.na(vec)]) + 1))
    out   <- list(transf     = vec.t,
                  colName    = colnames(df.tmp)[x],
                  transfType = "INT")
  })

  transf.df           <- as.data.frame(do.call(cbind, lapply(transf, function(x) x$transf)))
  colnames(transf.df) <-  unlist(lapply(transf, function(x) x$colName))
  if(!all(colnames(transf.df) == colnames(df))){
    stop("ERROR: colnames mismatch")
  }
  transfType.df <- data.frame(colname = colnames(transf.df), transfType = unlist(lapply(transf, function(x) x$transfType)), stringsAsFactors = F)
  list(transf.df     = transf.df,
       transfType.df = transfType.df)
}


get_QTL_run_file_traits <- function(qtl_results_file){
  # qtl_results_file <- y
  diet            <- gsub("__.*", "", gsub(".*__HFD", "HFD", gsub(".*__CD", "CD", qtl_results_file)))
  pheno_file_name <- paste0(ifelse(grepl("batch", qtl_results_file), paste0("BXD_pheno__batch_", gsub(".*__batch|\\.RDS", "", qtl_results_file)), "BXD_pheno"), ".csv")
  pheno_file_path <- gsub("\\/qtl2_results__.*", paste0("/", diet, "/", pheno_file_name), gsub("qtl_outputs", "qtl_processed_inputs", qtl_results_file))
  as.character(fread(pheno_file_path, skip = 3, nrows = 1, header = F)[1, ])[-1]
}


get_list_QTL_file_traits <- function(QTL_results_file_list){
  # QTL_results_file_list <- p_QTL_result_files
  out <- lapply(QTL_results_file_list, function(x){ 
    out.x        <- lapply(x, get_QTL_run_file_traits)
    names(out.x) <- x
    out.x
  })
  out
}



estimate_QTLpeak_pValue <- function(peaksDf, permDf){
  # peaksDf <- peaks.df

  operm1  <- permDf
  pVal.df <- lapply(1:nrow(peaksDf), function(x){
    # x <- 1
    # print(x)
    id       <- peaksDf$id[x]
    LOD_peak <- peaksDf$LOD[x]
    chrSet   <- peaksDf$chrSet[x]
    
    if(chrSet == "X"){
      permutations <- operm1$X[, id]
    } else if(chrSet == "notX"){
      permutations <- operm1$A[, id]
    }
    counts   <- sum(permutations >= LOD_peak)
    pVal     <- counts / length(permutations)
    data.frame(id = id, LOD = LOD_peak, chrSet = chrSet, pVal = pVal, stringsAsFactors = F)
  })
  type.convert(rbindlist(pVal.df), as.is = T)
}

getQTL_results_all <- function(qtlResultPaths, 
                               mapTable                      = NULL, 
                               mapTable_joinBy               = "gene_id_ensembl",
                               nCores                        = 10,
                               verbose                       = T){
  
  

  
  results.list.batches <- lapply(qtlResultPaths, function(batch_path){
    # batch_path <- qtlResultPaths[1]
    
    qtlResultsList <- readRDS(batch_path)
    
    if(verbose){
      cat("\n***************************************************************************************************************************************************************************************\n** Collecting and formatting QTL data for condition: ", 
          qtlResultsList$data.type, ",", qtlResultsList$condition, ",", qtlResultsList$interaction, "\t\t- batch ", which(batch_path == qtlResultPaths), "/", length(qtlResultPaths),
          "\n***************************************************************************************************************************************************************************************\n")
    }
    
    
    isX_present <- "X" %in% qtl2::chr_names(qtlResultsList$cross)
    isA_present <- any(unlist(lapply(c(1:19, "Y"), function(x) x %in% qtl2::chr_names(qtlResultsList$cross))))
    if(!isX_present & isA_present){
      qtlResultsList$threshold <- list("A" = qtlResultsList$threshold, "X" = NULL)
      qtlResultsList$operm     <- list("A" = qtlResultsList$operm, "X" = NULL)
    } else if(!isA_present & isX_present){
      # TODO
    }
    
    gmap <- rbindlist(lapply(names(qtlResultsList$gmap), function(x) data.frame(marker = names(qtlResultsList$gmap[[x]]), 
                                                                                pos_CM = unname(qtlResultsList$gmap[[x]]), 
                                                                                chr    = x, stringsAsFactors = F)))
    pmap <- rbindlist(lapply(names(qtlResultsList$pmap), function(x) data.frame(marker = names(qtlResultsList$pmap[[x]]), 
                                                                                pos_Mb = unname(qtlResultsList$pmap[[x]]), 
                                                                                chr    = x, stringsAsFactors = F)))
    
    
    ######
    # extract peaks of significant QTLs (exceeding LOD threshold from permutation)
    ######
    if("qtl_peaks_summary" %in% results_list){
      peaks_A_unit <- c("peaks.g"     = "CM", "peaks.p"     = "Mb")
      peaks_A      <- lapply(names(peaks_A_unit), function(x){
        
        
        # main peaks
        peaks.df          <- qtlResultsList[[x]]
        # secondary peaks
        peaks.df.peakdrop <- qtlResultsList[[paste0(x, ".sec")]]
        if(is.null(peaks.df)){
          return(NULL)
        }
        if(nrow(peaks.df) == 0){
          return(NULL)
        }
        if(verbose){
          cat("> Extracting main and secondary chr QTL peaks for chromosome set: autosomal - ", peaks_A_unit[[x]], "...\n")
        }
        
        # concatenate the peaks data.frames
        if(nrow(peaks.df.peakdrop) > 0){
          peaks.df$peakType          <- "topPeak_chr"
          peaks.df$uniqueID          <- paste0(peaks.df$lodcolumn, "_chr", peaks.df$chr, "@", peaks.df$pos)
          
          peaks.df.peakdrop$peakType <- "minorPeak_chr"
          peaks.df.peakdrop$uniqueID <- paste0(peaks.df.peakdrop$lodcolumn, "_chr", peaks.df.peakdrop$chr, "@", peaks.df.peakdrop$pos)
          peaks.df.peakdrop          <- peaks.df.peakdrop[!(peaks.df.peakdrop$uniqueID %in% peaks.df$uniqueID), ]
          peaks.df                   <- rbind(peaks.df, peaks.df.peakdrop)
        }
        
        
        
        peaks.df           <- dplyr::select(peaks.df, -c("lodindex", "uniqueID"))
        colnames(peaks.df) <- c("id", "chr", paste0("pos_", peaks_A_unit[x]), "LOD", paste0(c("ci_lo_", "ci_hi_"), peaks_A_unit[x]), "peak_type")
        peaks.df$chrSet    <- "notX"
        peaks.df
      })
      peaks_A <- merge(peaks_A[[1]], peaks_A[[2]], by = intersect(colnames(peaks_A[[1]]), colnames(peaks_A[[2]])))
      
      if(nrow(peaks_A) == 0){
        cat("> 0 QTL peaks found for autosomal chr...\n")
        peaks_A <- NULL
      }
      
      peaks_X_unit <- c("peaks.X.g" = "CM", "peaks.X.p" = "Mb")
      peaks_X      <- lapply(names(peaks_X_unit), function(x){
        
        
        # x <- names(peaks_X_unit)[1]
        peaks.df <- qtlResultsList[[x]]
        if(is.null(peaks.df)){
          return(NULL)
        }
        if(nrow(peaks.df) == 0){
          return(NULL)
        }
        if(verbose){
          cat("> Extracting main and secondary chr QTL peaks for chromosome set: X - ", peaks_X_unit[[x]], "...\n")
        }
        
        # in addition to the highest peak per chromosome already calculated, 
        # compute all the other small peaks. For this, use a very small peakdrop 
        # (without the drop parameter to estimate Bayes CI instead of 
        # drop-intervals)
        # --> not yet implemented for chrX! See above for autosomal chromosomes
        
        peaks.df <- peaks.df[, -1]
        colnames(peaks.df) <- c("id", "chr", paste0("pos_", peaks_A_unit[x]), "LOD", paste0(c("ci_lo_", "ci_hi_"), peaks_A_unit[x]))
        peaks.df$chrSet    <- "X"
        peaks.df
      })
      peaks_X <- merge(peaks_X[[1]], peaks_X[[2]], by = intersect(colnames(peaks_X[[1]]), colnames(peaks_X[[2]])))
      if(nrow(peaks_X) == 0){
        cat("> 0 QTL peaks found for chr X...\n")
        peaks_X <- NULL
      }
      
      LOD_thresh.df <- type.convert(rbindlist(lapply(names(qtlResultsList$threshold), function(x){
        if(is.null(qtlResultsList$threshold[[x]])){
          return(NULL)
        }
        data.frame(id = colnames(qtlResultsList$threshold[[x]]), LOD_threshold = as.numeric(qtlResultsList$threshold[[x]]), chrSet = ifelse(x == "A", "notX", "X"), stringsAsFactors = F)
      })), as.is = T)
      
      peaks.df <- rbindlist(list(peaks_A, peaks_X))
      if(nrow(peaks.df) > 0){
        est.pVal <- estimate_QTLpeak_pValue(peaks.df, qtlResultsList$operm)
        peaks.df <- type.convert(peaks.df, as.is = T)
        peaks.df <- merge(peaks.df, est.pVal, by = c("id", "LOD", "chrSet"))
        peaks.df <- merge(peaks.df, LOD_thresh.df, by = c("id", "chrSet"), all = F)
        peaks.df <- dplyr::select(peaks.df, c("id", "chr", "chrSet", "LOD_threshold", "LOD", "pVal", "peak_type", "pos_CM", "pos_Mb", "ci_lo_CM", "ci_hi_CM", "pos_Mb", "ci_lo_Mb", "ci_hi_Mb"))
        peaks.df <- type.convert(as.data.frame(peaks.df), as.is = T)
        
        allIDs <- unique(peaks.df$id)
        allIDs <- allIDs[!is.na(allIDs) & !(allIDs %in% c("", " ", "NA"))]
        if(!is.null(mapTable) &  any(allIDs %in% mapTable[[mapTable_joinBy]])){
          # mapTable_tmp                      <- dplyr::select(mapTable, -c("chromosome_start_position", "chromosome_end_position"))
          # colnames(mapTable_tmp)            <- c("gene_id_ensembl", "gene_symbol", "chr_gene", "gene_posStart_Mb", "gene_posEnd_Mb")
          colsSelect                        <- unique(c(mapTable_joinBy, "gene_id_ensembl", "gene_symbol", "chromosome", "posStart_Mb", "posEnd_Mb"))
          colsSelect                        <- colsSelect[colsSelect %in% colnames(mapTable)]
          mapTable_tmp                      <- unique(dplyr::select(mapTable, colsSelect))
          colnames(mapTable_tmp)            <- plyr::mapvalues(colnames(mapTable_tmp), from = c("chromosome", "posStart_Mb", "posEnd_Mb"), to = c("chr_gene", "gene_posStart_Mb", "gene_posEnd_Mb"), warn_missing = T)
          mapTable_tmp$id                   <- mapTable_tmp[[mapTable_joinBy]]
          mapTable_tmp                      <- mapTable_tmp[!is.na(mapTable_tmp$id) & mapTable_tmp$id != "", ]
          
          peaks.df                           <- plyr::join(peaks.df, mapTable_tmp[mapTable_tmp$id %in% qtl_markers$id, ], by = c("id"), type = "full")
          peaks.df$dist_from_geneStart_Mb    <- peaks.df$pos_Mb - peaks.df$gene_posStart_Mb
          peaks.df$dist_from_geneStart_Mb[peaks.df$chr != peaks.df$chr_gene] <- NA
        }
        
        peaks.df$closestMarkers <- unlist(lapply(1:nrow(peaks.df), function(x){
          out      <- pmap[pmap$chr == peaks.df$chr[x], ]
          out$dist <- abs(out$pos_Mb - peaks.df$pos_Mb[x])
          paste(out$marker[out$dist == min(out$dist)], collapse = ";")
        }))
        pmap.tmp <- pmap[!grepl("^c", pmap$marker) & !grepl("loc", pmap$marker), ]
        peaks.df$closestMarkers_noPseudo <- unlist(lapply(1:nrow(peaks.df), function(x){
          out      <- pmap.tmp[pmap.tmp$chr == peaks.df$chr[x], ]
          out$dist <- abs(out$pos_Mb - peaks.df$pos_Mb[x])
          paste(out$marker[out$dist == min(out$dist)], collapse = ";")
        }))
      } else{
        peaks.df <- NULL
      }
      
      
    } else{
      peaks.df <- NULL
    }
    
    if(verbose){
      cat("> Done!\n")
    }
    
    list(peaks.df       = peaks.df)
    
  })
  
  
  # merge results for all batches
  results_names    <- unique(unlist(lapply(results.list.batches, names)))
  results.list.all <- lapply(results_names, function(x){
    # x <- results_names[1]
    # print(x)
    all.x    <- lapply(results.list.batches, function(kk) kk[[x]])
    tmpLogi  <- unlist(lapply(all.x, function(kk) all(class(kk) == "list" | is.null(kk))))
    if(all(tmpLogi)){
      return(NULL)
    }
    all.x    <- do.call(rbind, all.x[!tmpLogi])
    chrs_all <- stringr::str_sort(unique(all.x$chr), numeric = T)
    all.x    <- all.x[order(factor(all.x$chr, levels = chrs_all), all.x$pos_Mb, all.x$id), ]
    all.x
  })
  names(results.list.all) <- results_names
  
  results.list.all
  
}



