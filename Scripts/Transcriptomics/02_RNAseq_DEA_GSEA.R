

library(qtl)
library(limma)
library(edgeR)
library(ggplot2)
library(data.table)


###############
## perform DEA
###############

metadata            <- readRDS("./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")
countMatrix         <- readRDS("./Data/RNAseq_processing/countMatrix/geneLevel/STAR_RNAseq_count_matrix.RDS")
metadata$strainDiet <- paste0(metadata$strain, "_", metadata$diet)
metadata            <- metadata[metadata$sample_name %in% colnames(countMatrix), ]
metadata            <- metadata[match(colnames(countMatrix), metadata$sample_name), ]
rownames(metadata)  <- metadata$sample_name
colnames(metadata)  <- gsub("\\.", "_", gsub("\\[", "", gsub("\\]\\.|\\/", "_", colnames(metadata))))
colnames(metadata)  <- make.names(colnames(metadata))
countMatrix         <- as.matrix(countMatrix)
stopifnot(identical(rownames(metadata), colnames(countMatrix)))

y               <- edgeR::DGEList(counts = countMatrix, samples = metadata)
keep.exprs      <- edgeR::filterByExpr(y, group = "strain_diet")
y               <- y[keep.exprs, , keep.lib.sizes = FALSE]
hasCounts       <- rowSums(is.na(y$counts)) == 0
y               <- y[hasCounts, , keep.lib.sizes = FALSE]
y               <- edgeR::calcNormFactors(y)

y$samples$strain        <- gsub("\\/", "_", y$samples$strain)
y$samples$strainDiet    <- gsub("\\/", "_", y$samples$strainDiet)
geneTable               <- readRDS("./Data/RNAseq_processing/geneConversionTables/geneConversionTable_GRCm38_release-102.RDS")

# prepare contrasts
globalDietContrast <- c("HFD_vs_CD_diet" = "(HFD)-(CD)")

# define covariates
covariates <-  c("concentration_[ng/uL].BGI", "RIN.BGI", "28S/18S.BGI")
covariates <- gsub("\\.", "_", gsub("\\[", "", gsub("\\]\\.|\\/", "_", covariates)))
covariates <- make.names(covariates)


designFormula  <- reformulate(response = NULL, termlabels = c(0, "diet", covariates))
mm             <- model.matrix(designFormula, data = y$samples)
yv             <- voom(y, mm, plot = T)
dc             <- duplicateCorrelation(yv, design = mm, block = yv$targets$correctStrain)
yv             <- voom(y, mm, block = yv$targets$correctStrain, correlation = dc$consensus.correlation,  plot = T)
colnames(mm)   <- gsub("diet", "", colnames(mm))
allContrasts   <- globalDietContrast

contrasts    <- lapply(allContrasts, function(kk){
  # kk <- allContrasts[1]
  commandstr <- paste0("out <- makeContrasts(",unname(kk),", levels = mm)")
  eval(parse(text = commandstr))
  colnames(out) <- names(kk)
  out
})
fitModel <- limma::lmFit(yv, design = mm, block = yv$targets$correctStrain, correlation = dc$consensus.correlation)
allFits  <- lapply(contrasts, function(x){
  limma::eBayes(limma::contrasts.fit(fitModel, contrasts = x))
})
allTopTables <- lapply(names(allFits), function(x){
  # x <- names(allFits)[1]
  out                    <- limma::topTable(allFits[[x]], sort.by = "none", n = Inf)
  out$gene_id            <- rownames(out)
  cNames                 <- colnames(out)
  out                    <- merge(out, unique(dplyr::select(geneTable, c("gene_id", "gene_name"))), by = "gene_id", all.x = T, all.y = F)
  out                    <- dplyr::select(out, c(cNames, "gene_name"))
  tmpLogi                <- is.na(out$gene_name)
  out$gene_name[tmpLogi] <- out$gene_id[tmpLogi]
  out$contrastID         <- x
  out$contrast           <- unname(allContrasts[x])
  out
})
allTopTables <- do.call(rbind, allTopTables)

saveDir      <- "./Data/RNAseq_processing/limma_DEA/limma_outputs"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}
write.table(allTopTables, paste0(saveDir, "/DEA_transcriptome_limma_table_all_contrasts.csv"), row.names = F, sep = ",")
saveRDS(allTopTables, paste0(saveDir, "/DEA_transcriptome_limma_table_all_contrasts.RDS"))



############
# Download mitocarta gene sets
############

saveDir <- "./Data/input_data/mitochondrial_ressources_geneSets"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
system("wget -P ./Data/input_data/mitochondrial_ressources_geneSets https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls")


############
# Functional enrichment - GSEA
############

geneTable         <- readRDS("./Data/RNAseq_processing/geneConversionTables/geneConversionTable_GRCm38_release-102.RDS")
tmp.df            <- readRDS("./Data/RNAseq_processing/limma_DEA/limma_outputs/DEA_transcriptome_limma_table_all_contrasts.RDS")
tmp.df            <- split(tmp.df, tmp.df$contrastID)
dea.df.list       <- list("population_diet_effect" = tmp.df)

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


overwrite <- F
gsea      <- lapply(c("logFC_ranking"), function(s){
  # s <- "logFC_ranking"
  
  
  lapply(names(dea.df.list), function(x){
    # x <- names(dea.df.list)[2]
    
    saveDir <- paste0("./Data/RNAseq_processing/limma_DEA/GSEA/", s)
    if(!dir.exists(saveDir)){
      dir.create(saveDir)
    }
    
    out_file <- paste0(saveDir, "/GSEA_trancriptomic__contrast_type_", x, ".RDS")
    if(file.exists(out_file) & !overwrite){
      cat("\t--> Results files already present. Skipping...\n")
      return()
    }
    
    out.list <- lapply(names(dea.df.list[[x]]), function(y){
      # y <- names(dea.df.list[[x]])[1]
      
      cat(paste0("Running GSEA on:\n\t> contrast type:\t", x, "\n\t> contrast:\t\t", y, "\n\t> ranking type:\t\t", s, "\n"))
      
      tmp <- dea.df.list[[x]][[y]]
      stopifnot(length(unique(tmp$contrast)) == 1)
      if(nrow(tmp) == 0){
        return(NULL)
      }
      
      gene_id_col     <- ifelse("gene_id" %in% colnames(tmp), "gene_id", "ensembl_gene_id")
      tmp             <- tmp[!duplicated(tmp[[gene_id_col]]), ]
      tmp             <- tmp[!is.na(tmp$logFC), ]
      
     
      if(s == "logFC_ranking"){
        tmp             <- tmp[order(tmp$logFC, decreasing = T), ]
        geneList        <- tmp$logFC
        names(geneList) <- tmp[[gene_id_col]]
      }
      
      out <- mclapply(mc.cores = length(names(msigdbr.geneSets)), X = names(msigdbr.geneSets), FUN = function(k){
        # k <- names(msigdbr.geneSets)[1]
        
        
        gsea <- clusterProfiler::GSEA(geneList      = geneList,
                                      nPermSimple   = 10000,
                                      pvalueCutoff  = 1,
                                      pAdjustMethod = "BH",
                                      TERM2GENE     = dplyr::select(msigdbr.geneSets[[k]], c("gs_name", "ensembl_gene")))
        gsea_df   <- gsea@result
        if(nrow(gsea_df) == 0){
          return(list(geneList = geneList,
                      gsea     = gsea,
                      gsea_df  = NULL))
        }
        colsSelect <- c("ensembl_gene", "gs_exact_source", "gs_id", "gs_description")
        colsSelect <- colsSelect[colsSelect %in% colnames(msigdbr.geneSets[[k]])]
        termsMeta  <- unique(dplyr::select(msigdbr.geneSets[[k]], dplyr::all_of(colsSelect)))
        gsea_df    <- merge(gsea_df, termsMeta, by.x = "ID", by.y = "ensembl_gene", all.x = T, all.y = F)
        gsea_df    <- gsea_df[order(gsea_df$pvalue), ]
        
        rownames(gsea_df)         <- gsea_df$ID
        gsea_df$MSigDB_collection <- k
        # gsea_df$LFCdirection      <- k
        gsea_df$contrastID        <- y
        gsea_df$contrast          <- tmp$contrast[1]
        
        list(geneList = geneList,
             gsea     = gsea,
             gsea_df  = gsea_df)
        
      })
      
      
      gsea_df           <- as.data.frame(data.table::rbindlist(lapply(out, function(z) z$gsea_df), use.names = T, fill = T), stringsAsFactors = F)
      gsea_df           <- type.convert(gsea_df, as.is = T)
      gsea_df$coreSize  <- unlist(lapply(gsea_df$core_enrichment, function(x) length(unlist(str_split(x, "/")))))
      gsea_df$geneRatio <- gsea_df$coreSize / gsea_df$setSize
      
      
      gsea_df$Description          <- gsub("^GOBP_|^GOCC_|^GOMF_|^HALLMARK_|^HP_|^KEGG_|^REACTOME_|^WP_|^BIOCARTA_|^PID_", "", gsea_df$Description)
      tmpLogi                      <- gsea_df$MSigDB_collection %in% c("CELLTYPE", "CGP", "IMMUNESIGDB")
      gsea_df$Description[tmpLogi] <- gsub("^.+?_(.*)", "\\1", gsea_df$Description[tmpLogi])
      gsea_df$Description          <- gsub("_", " ", gsea_df$Description)
      gsea_df
    })
    names(out.list) <- names(dea.df.list[[x]])
    
    saveRDS(out.list, out_file)
  })
})






