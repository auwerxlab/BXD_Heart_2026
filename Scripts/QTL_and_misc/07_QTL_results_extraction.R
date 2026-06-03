
library(qtl2)
library(data.table)


################################################################################
## load and aggregate qtl results
################################################################################

dirQTL_outputs <- "./Data/QTL_mapping/qtl_outputs"
resultFiles    <- sort(list.files(dirQTL_outputs, pattern = "qtl2_results", recursive = T, full.names = T))

geneTable                           <- readRDS("./Data/RNAseq_processing/geneConversionTables/geneConversionTable_GRCm38_release-102_complete_gtf.RDS")
geneTable                           <- type.convert(as.data.frame(geneTable), as.is = T)
geneTable                           <- geneTable[geneTable$seqnames %in% c(1:20, "X", "Y", "MT"), ]
geneTable                           <- dplyr::select(geneTable, c("gene_id", "gene_name", "seqnames", "start", "end", "strand", "type"))
colnames(geneTable)                 <- c("gene_id_ensembl", "gene_symbol", "chromosome", "chromosome_start_position", "chromosome_end_position", "strand", "type")
geneTable                           <- geneTable[geneTable$type == "gene", ]
geneTable                           <- unique(geneTable)
geneTable                           <- dplyr::select(geneTable, -c("strand", "type"))
geneTable$chromosome_start_position <- as.numeric(geneTable$chromosome_start_position)
geneTable$chromosome_end_position   <- as.numeric(geneTable$chromosome_end_position)
geneTable$posStart_Mb               <- geneTable$chromosome_start_position / 1e6
geneTable$posEnd_Mb                 <- geneTable$chromosome_end_position / 1e6


protTable_Coon                          <- readRDS("./Data/input_data/proteomics/coon_formatted/proteins_metadata.RDS")
colnames(protTable_Coon)                <- c("protein_id_uniprot", "peptide_id_ensembl", "gene_id_entrez", "gene_id_ensembl", "transcript_id_ensembl", "gene_symbol", "gene_symbol_uniprot")
protTable_Coon$protein_id_iso_uniprot   <- protTable_Coon$protein_id_uniprot
protTable_Coon$protein_id_noIso_uniprot <- gsub("-.", "", protTable_Coon$protein_id_uniprot)
protTable_Coon                          <- unique(dplyr::select(protTable_Coon, c("protein_id_uniprot", "protein_id_iso_uniprot", "protein_id_noIso_uniprot", "gene_id_ensembl")))
protTable_Coon                          <- unique(merge(protTable_Coon, geneTable, by = "gene_id_ensembl", all.x = T, all.y = F))
protTable_Coon                          <- protTable_Coon[!is.na(protTable_Coon$gene_id_ensembl) & !(protTable_Coon$gene_id_ensembl %in% c("", " ", "NA")) & !is.na(protTable_Coon$posStart_Mb), ]
protTable_Coon                          <- dplyr::select(protTable_Coon, c("gene_id_ensembl", "protein_id_uniprot", "protein_id_iso_uniprot", "protein_id_noIso_uniprot", "gene_symbol", "chromosome", "chromosome_start_position", "chromosome_end_position", "posStart_Mb", "posEnd_Mb"))


# Split result paths by input type and run type (nested list)
types              <- dir(dirQTL_outputs)
resultFiles_byType <- unlist(lapply(resultFiles, function(x) types[unlist(lapply(types, function(y) grepl(paste0("\\/", y, "\\/"), x)))]))
resultFiles_byType <- split(resultFiles, resultFiles_byType)
resultFiles_byType <- rev(resultFiles_byType)
resultFiles_byType <- lapply(resultFiles_byType, function(x){
  # x <- resultFiles_byType$individuals_proteomic_Coon
  names(x) <- gsub(".*qtl2_results__|__\\.RDS|\\.RDS", "", x)
  split(unname(x), names(x))
})

nCores     <- 40
resultsAll <- lapply(names(resultFiles_byType), function(kk){
  # kk <- names(resultFiles_byType)[1]
  
  cat("\n\nExtracting formatted results for data type: ", kk, "\n")
  saveDir <- paste0(dataPath, "/QTL_mapping/qtl_outputs_formatted/", kk)
  if(!dir.exists(saveDir)){
    dir.create(saveDir, recursive = T)
  }
  results <- lapply(names(resultFiles_byType[[kk]]), function(x){
    # x <- names(resultFiles_byType[[kk]])[2]
    
    if(grepl("RNAseq", x)){
      mapperTable.x        <- geneTable
      mapperJoinBy.x       <- "gene_id_ensembl"
    } else if(grepl("proteo", x, ignore.case = T)){
      mapperTable.x        <- protTable_Coon
      mapperJoinBy.x       <- "protein_id_uniprot"
    } else if(grepl("pheno", x)){
      mapperTable.x        <- NULL
      mapperJoinBy.x       <- NULL
    } else if(grepl("lipidomic", x)){
      mapperTable.x        <- NULL
      mapperJoinBy.x       <- NULL
    } 
    
    # get results by looping over all result files (batches, if multiple were used to speed up QTL mapping) 
    # for a given input type and analysis condition 
    res.list <- getQTL_results_all(resultFiles_byType[[kk]][[x]], 
                                   mapperTable.x, 
                                   results_list                  = resultsList.x,
                                   get_LOD_scores_element_idList = lod_profile_traits.x,
                                   mapTable_joinBy               = mapperJoinBy.x,
                                   nCores                        = nCores)
    
    res.list <- lapply(res.list, function(kk){
      if(!is.null(kk)){
        if("inputType" %in% colnames(kk)){
          out           <- kk
          print(unique(out$inputType))
          out$inputType <- gsub("__", "", gsub("_derived_traits_\\d+|_fasted|_allStrains|_commonStrains|_$", "", out$inputType))
          print(unique(out$inputType))
          return(out)
        }
      }
      kk
    })
    
    res.list
  })
  
  peaks.df       <- type.convert(as.data.frame(rbindlist(lapply(results, function(x) x$peaks.df), use.names = T, fill = T)), as.is = T)

  saveRDS(peaks.df, paste0(saveDir, "/qtl_peaks_table_", kk, ".RDS"))
  write.table(peaks.df, paste0(saveDir, "/qtl_peaks_table_", kk, ".csv"), row.names = F, sep = ",")
  
  NULL
})
names(resultsAll) <- names(resultFiles_byType)


