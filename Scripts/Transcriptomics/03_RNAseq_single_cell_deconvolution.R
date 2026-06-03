
# MuSiC forked repo as mentioned in manuscript's methods
remotes::install_github("vonalven/MuSiC")

library(MuSiC)
library(Seurat)
library(Biobase)
library(data.table)
library(SummarizedExperiment)


##################
## scRNA-seq data download
## GSE109774 (Tabula Muris)
##################


## First download manually files from
## https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687/7
## https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_microfluidic_emulsion_v2_/5968960/2
## in folder ./Data/input_data/heart_singleCell_datasets/GSE109774/
# and deflate them. Recenty Figshare disabled programmatic download.

required_files <- c("./Data/input_data/heart_singleCell_datasets/GSE109774/annotations_FACS.csv",
                    "./Data/input_data/heart_singleCell_datasets/GSE109774/annotations_droplet.csv",
                    "./Data/input_data/heart_singleCell_datasets/GSE109774/droplet_Heart_seurat_tiss.Robj",
                    "./Data/input_data/heart_singleCell_datasets/GSE109774/facs_Heart_seurat_tiss.Robj")
tmpLogi <- file.exists(required_files)
if(any(!tmpLogi)){
  cat("\nManually download the following scRNA-seq files from Figshare, as detailed in the single-cell deconvolution scrip data download's section:\n",
      paste(paste0("\t--> ", required_files[!tmpLogi]), collapse = "\n"), "\n\n")
}

##################
## scRNA-seq data preparation
## GSE109774 (Tabula Muris)
##################

## check supplementary data
## https://www.nature.com/articles/s41586-018-0590-4
## https://figshare.com/articles/dataset/Robject_files_for_tissues_processed_by_Seurat/5821263/1
## https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687/7
## https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_microfluidic_emulsion_v2_/5968960/2


dirSave <- "./Data/input_data/heart_singleCell_datasets/GSE109774/formatted"
if(!dir.exists(dirSave)){
  dir.create(dirSave)
}
facs.annotation    <- fread("./Data/input_data/heart_singleCell_datasets/GSE109774/annotations_FACS.csv")
droplet.annotation <- fread("./Data/input_data/heart_singleCell_datasets/GSE109774/annotations_droplet.csv")
annotation.all     <- rbindlist(list(facs.annotation, droplet.annotation), use.names = T, fill = T)

# ## From the paper: 4 mice used in total for liver (4 samples)
# ## Four 10–15 week old male and four virgin female C57BL/6JN mice
# ## Isolating viable single cells from both the pancreas and the liver of the 
# ## same mouse was not possible; therefore, two males and two females were used 
# ## for each.

file_paths <- c("droplet" = "./Data/input_data/heart_singleCell_datasets/GSE109774/droplet_Heart_seurat_tiss.Robj",
                "facs"    = "./Data/input_data/heart_singleCell_datasets/GSE109774/facs_Heart_seurat_tiss.Robj")

overwrite <- F
if(any(grepl("GSE109774__gene_conversion_table", list.files(dirSave))) & !overwrite){
  
  cat("Gene conversion table already existing!\n")
  
} else{
  
  genes.all <- lapply(file_paths, function(x){
    seurat.obj <- get(load(x))
    rownames(seurat.obj@raw.data)
  })
  genes.all <- sort(unique(unlist(genes.all)))
  
  mart           <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  filters.df     <- biomaRt::listFilters(mart)
  attributes.df  <- biomaRt::listAttributes(mart)
  attr.vec       <- c("external_gene_name", "mgi_symbol", "hgnc_symbol", "wikigene_name", "uniprot_gn_id", "uniprot_gn_symbol", "entrezgene_accession", "external_synonym",
                      "external_transcript_name", "entrezgene_trans_name", "mgi_trans_name", "hgnc_trans_name", "protein_id")
  stopifnot(all(attr.vec %in% filters.df$name))
  stopifnot(all(attr.vec %in% attributes.df$name))
  gene_info.list <- lapply(attr.vec, function(x){
                               cat("\nQuerying:", x, "...\n\n")
                               out <- biomaRt::getBM(attributes = c("ensembl_gene_id", x), 
                                                     filters = x, 
                                                     values  = genes.all, 
                                                     mart    = mart,
                                                     useCache = F)
                               colnames(out) <- c("ensembl_gene_id", "gene_symbol")
                               out
                             })
  names(gene_info.list) <- attr.vec
  gene_info.df <- unique(rbindlist(gene_info.list, use.names = T))
  table(genes.all %in% gene_info.df$gene_symbol)
  gene.missing <- genes.all[!(genes.all %in% gene_info.df$gene_symbol)]
  gene.missing
  
  geneConversionTable <- readRDS("./Data/RNAseq_processing/geneConversionTables/geneConversionTable_GRCm38_release-102.RDS")
  genes.all[!(genes.all %in% geneConversionTable$gene_name)]
  table(genes.all %in% geneConversionTable$gene_name)
  table(genes.all %in% gene_info.df$gene_symbol)
  
  tmp.df <- geneConversionTable[geneConversionTable$gene_name %in% gene.missing, ]
  tmp.df <- dplyr::select(tmp.df, c("gene_id", "gene_name"))
  colnames(tmp.df) <- c("ensembl_gene_id", "gene_symbol")
  
  gene_info.df <- unique(rbind(gene_info.df, tmp.df))
  gene.missing <- genes.all[!(genes.all %in% gene_info.df$gene_symbol)]
  gene.missing
  
  gene_info.df <- type.convert(as.data.frame(gene_info.df), as.is = T)
  
  # for the duplciated genes with multiple ensembl gene IDs, keep the ensembl gene ID present in our data
  dupl.genes <- unique(gene_info.df$gene_symbol[duplicated(gene_info.df$gene_symbol)])
  tmp.df.1   <- gene_info.df[!(gene_info.df$gene_symbol %in% dupl.genes), ]
  tmp.df.2   <- gene_info.df[(gene_info.df$gene_symbol %in% dupl.genes), ]
  tmp.df.2   <- plyr::ddply(tmp.df.2, "gene_symbol", function(kk){
    tmpLogi <- kk$ensembl_gene_id %in% geneConversionTable$gene_id
    if(any(tmpLogi)){
      kk[tmpLogi, ]
    } else{
      kk
    }
  })
  gene_info.df <- rbind(tmp.df.1, tmp.df.2)
  dupl.genes   <- unique(gene_info.df$gene_symbol[duplicated(gene_info.df$gene_symbol)])
  length(unique(dupl.genes))
  
  gene_info.df <- gene_info.df[gene_info.df$gene_symbol %in% genes.all, ]
  
  saveRDS(gene_info.df, paste0(dirSave, "/GSE109774__gene_conversion_table.RDS"))
  
}


gene_info.df        <- readRDS("./Data/input_data/heart_singleCell_datasets/GSE109774/formatted/GSE109774__gene_conversion_table.RDS")
dupl.genes.symbol   <- unique(gene_info.df$gene_symbol[duplicated(gene_info.df$gene_symbol)])
gene_info.df        <- gene_info.df[!(gene_info.df$gene_symbol %in% dupl.genes.symbol), ]
dupl.genes.id       <- unique(gene_info.df$ensembl_gene_id[duplicated(gene_info.df$ensembl_gene_id)])
gene_info.df        <- gene_info.df[!(gene_info.df$ensembl_gene_id %in% dupl.genes.id), ]

raw.UMI.list <- lapply(c("droplet", "facs"), function(x){
  # x <- "facs"
  print(x)
  filePath   <- file_paths[[x]]
  # Deconvolution requires raw UMI
  seurat.obj <- get(load(filePath))
  raw.UMI    <- seurat.obj@raw.data # the UMI counts
  print(dim(raw.UMI))
  
  raw.UMI    <- raw.UMI[rownames(raw.UMI) %in% gene_info.df$gene_symbol, ]
  print(dim(raw.UMI))
  rownames(raw.UMI) <- plyr::mapvalues(rownames(raw.UMI), from = gene_info.df$gene_symbol, to = gene_info.df$ensembl_gene_id, warn_missing = F)
  data       <- seurat.obj@data     # the processed data
  # keep only filtered cells (more than 500 genes and 50,000 reads), according
  # to original paper methods:
  # https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0590-4/MediaObjects/41586_2018_590_MOESM4_ESM.pdf
  raw.UMI    <- dplyr::select(as.data.frame(raw.UMI), colnames(data))
  raw.UMI
})
names(raw.UMI.list) <- c("droplet", "facs")
stopifnot(length(intersect(colnames(raw.UMI.list$droplet), colnames(raw.UMI.list$facs))) == 0)
unlist(lapply(raw.UMI.list, nrow))

allCellIDs     <- c(colnames(raw.UMI.list$droplet), colnames(raw.UMI.list$facs))
stopifnot(!any(!(allCellIDs %in% c(facs.annotation$cell, droplet.annotation$cell))))
annotation.all <- annotation.all[annotation.all$cell %in% allCellIDs, ]
table(annotation.all$mouse.sex)


# prepare objects for deconvolution, use separated object for droplet and facs.
# from original paper:
# Two distinct technical approaches were used for most organs: one approach, 
# microfluidic droplet-based 3′-end counting, enabled the survey of thousands 
# of cells at relatively low coverage, whereas the other, full-length transcript 
# analysis based on fluorescence-activated cell sorting, enabled the 
# characterization of cell types with high sensitivity and coverage. The 
# cumulative data provide the foundation for an atlas of transcriptomic cell 
# biology.

# https://tabula-sapiens-portal.ds.czbiohub.org/organs
# cell_ontology_class: Cell type annotations using the Cell Ontology.
# free_annotation: Cell type annotations using free text.

expr.list.all        <- readRDS("./Data/external_RNAseq_datasets_formatted/formatted_gene_expression_datasets/expression_tables_mouse_human_filtered_lowly_expressed_genes.RDS")
expr.mouse.bulk.filt <- expr.list.all$`Mouse LV CD_HFD (CPM)`

tt <- lapply(c("all_cell_types", "selected_cell_types", "selected_cell_types_2", "selected_cell_types_3"), function(ct){
  # ct <- "all_cell_types"
  lapply(names(raw.UMI.list), function(x){
    cat("Saving object for:", x, ",", ct, "...\n")
    
    expr <- raw.UMI.list[[x]]
    meta <- as.data.frame(annotation.all[annotation.all$cell %in% colnames(expr), ])
    meta <- meta[meta$tissue == "Heart", ]
    
    if(ct == "selected_cell_types"){
      meta <- meta[meta$cell_ontology_class %in% c("cardiac muscle cell", "endothelial cell", "fibroblast"), ]
    } else if(ct == "selected_cell_types_2"){
      meta <- meta[meta$cell_ontology_class %in% c("cardiac muscle cell", "endothelial cell", "fibroblast", "myofibroblast cell"), ]
    } else if(ct == "selected_cell_types_3"){
      meta <- meta[meta$cell_ontology_class %in% c("cardiac muscle cell", "endothelial cell", "fibroblast", "myofibroblast cell", "endocardial cell"), ]
    }
    if(nrow(meta) == 0){
      return()
    }
    expr <- dplyr::select(expr, meta$cell)
    stopifnot(ncol(expr) == nrow(meta))
    stopifnot(all(meta$cell == colnames(expr)))
    print(dim(expr))
    
    meta$cellID    <- meta$cell
    meta$sampleID  <- meta$mouse.id
    meta$celltype1 <- ifelse(!is.na(meta$cell_ontology_class) & !(meta$cell_ontology_class %in% c("", " ")), meta$cell_ontology_class, NA)
    rownames(meta) <- meta$cell
    meta           <- meta[colnames(expr), ]
    scExp          <- SingleCellExperiment::SingleCellExperiment(list(counts = as.matrix(expr)), colData = meta)
    
    saveRDS(scExp, paste0(dirSave, "/GSE109774_mouse_scExp__", x, "__", ct, "__all_genes.RDS"))
    NA
  })
})


##############
## Prepare bulk sequencing data objects
## MuSiC works best with raw counts
##############

saveDir <- "./Data/singleCell_deconvolution/bulkInput"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

##############
## mouse data
##############

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


expr.mouse <- as.matrix(expr.mouse)

metadata                              <- dplyr::select(metadata, c("sample_name", "strainDiet", "strain", "diet"))
metadata$sampleID                     <- metadata$sample_name
rownames(metadata)                    <- metadata$sampleID
expr.mouse                            <- expr.mouse[, colnames(expr.mouse) %in% metadata$sampleID]
colnames(expr.mouse)[!(colnames(expr.mouse) %in% metadata$sampleID)]
metadata                              <- metadata[colnames(expr.mouse), ]


scExp.mouse <- SingleCellExperiment::SingleCellExperiment(list(counts = as.matrix(expr.mouse)), colData = metadata)
saveRDS(scExp.mouse, paste0(saveDir, "/mouse_bulk_scExp.RDS"))


# split mouse eSet by diet
dirSave.tmp <- paste0(saveDir, "/mouse_single_diet")
if(!dir.exists(dirSave.tmp)){
  dir.create(dirSave.tmp)
}
tt <- lapply(unique(metadata$diet), function(aa){
  # aa <- unique(metadata$strainDiet)[1]
  print(aa)
  
  meta.aa  <- metadata[metadata$diet == aa, ]
  expr.aa  <- expr.mouse[, colnames(expr.mouse) %in% rownames(meta.aa)]
  meta.aa  <- meta.aa[colnames(expr.aa), ]
  all(colnames(expr.aa) == rownames(meta.aa))
  
  scExp.mouse.aa <- SingleCellExperiment::SingleCellExperiment(list(counts = as.matrix(expr.aa)), colData = meta.aa)
  saveRDS(scExp.mouse.aa, paste0(dirSave.tmp, "/diet_", gsub("\\/", "___", aa), "_bulk_scExp.RDS"))
})


################################################################################
## Run Single Cell Deconvolution with MuSiC
################################################################################

################################################################################
## Run music deconvolution
## MuSiC works best with raw UMI counts
################################################################################

saveDir <- "./Data/singleCell_deconvolution/deconvolutionResults"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
dir.exists(saveDir)

scDatasets        <- list.files("./Data/input_data/heart_singleCell_datasets/GSE109774/formatted/", full.names = T, pattern = "GSE109774_mouse_scExp__")
scDatasets        <- scDatasets[grepl("RDS$", scDatasets)]
names(scDatasets) <- paste0("GSE109774__", gsub(".*GSE109774_mouse_scExp__|\\.RDS", "", scDatasets))
stopifnot(length(unique(names(scDatasets))) == length(scDatasets))
scDatasets        <- lapply(scDatasets, function(x){
  list(esetPath    = unname(x),
       clusterCols = c("celltype1"))
})
scDatasets        <- list("mouse" = scDatasets)

mouse.bulk.list        <- list.files("./Data/singleCell_deconvolution/bulkInput/mouse_single_diet", full.names = T, pattern = "_bulk_scExp")
listNames              <- gsub("_bulk_scExp\\.RDS", "", gsub(".*diet_", "", mouse.bulk.list))
mouse.bulk.list        <- lapply(mouse.bulk.list, function(x) list(esetPath = x))
names(mouse.bulk.list) <- listNames
bulkDatasets           <- c(mouse = list(mouse.bulk.list))

celltypes.grouped <- list("cardiomyocytes"      = c("cardiac muscle cell", "cardiac muscle cell_M", "cardiac muscle cell_F"),
                          "endothelial_cells"   = c("endothelial cell", "endothelial cell_F", "endocardial cell_M", 
                                                    "endothelial cells_1", "endothelial cells_2", "endothelial cells_3", 
                                                    "endothelial cells_1_M", "endothelial cells_2_M", "endothelial cells_3_M", 
                                                    "endothelial cells_1_F", "endothelial cells_2_F", "endothelial cells_3_F"),
                          "fibroblasts"         = c("fibroblast", "fibroblast_F", "fibroblast_M", "fibroblasts", "fibroblasts_M"),
                          "myofibroblasts"      = c("myofibroblast cell", "myofibroblast cell_F", "myofibroblast cell_M", 
                                                    "Myofibroblast", "Myofibroblast_F", "Myofibroblast_M"),
                          "endocardial_cells"   = c("endocardial cell", "endocardial cell_F", "endothelial cell_M"),
                          "smooth_muscle_cells" = c("smooth muscle cell", "smooth muscle cell_M", "smooth muscle cell_F", 
                                                    "Smooth_Muscle_Cells", "Smooth_Muscle_Cells_M", "Smooth_Muscle_Cells_F"),
                          "erythrocytes"        = c("erythrocyte", "erythrocyte_M", "erythrocyte_F", 
                                                    "red blood cells", "red blood cells_M", "red blood cells_F"),
                          "leukocytes"          = c("leukocyte", "leukocyte_F", "leukocyte_M"),
                          "professional_APCs"   = c("professional antigen presenting cell", "professional antigen presenting cell_M", 
                                                    "professional antigen presenting cell_F", "antigen presenting cells",
                                                    "antigen presenting cells_M", "antigen presenting cells_F"),
                          "cardiac neurons"     = c("cardiac neuron", "cardiac neuron_M", "cardiac neuron_F", 
                                                  "conduction cells", "conduction cells_M", "conduction cells_F"))

#######################
## Run using MuSiC 2 ##
#######################

saveDir_music2_table <- "./Data/singleCell_deconvolution/deconvolutionResults_music2/prop_tables"
saveDir_music2_list  <- "./Data/singleCell_deconvolution/deconvolutionResults_music2/iteration_results"
for(i in c(saveDir_music2_table, saveDir_music2_list)){
  if(!dir.exists(i)){
    dir.create(i, recursive = T)
  }
}


mouse.bulk.list        <- list.files("./Data/singleCell_deconvolution/bulkInput", full.names = T, pattern = "_bulk_scExp")
listNames              <- "CD_HFD"
mouse.bulk.list        <- lapply(mouse.bulk.list, function(x) list(esetPath = x))
names(mouse.bulk.list) <- listNames
bulkDatasets           <- c(mouse = list(mouse.bulk.list))

overwrite        <- F
org              <- "mouse"
nb_retry         <- 20
nbCores.external <- 10
nbCores.internal <- 9
tt               <- lapply(names(bulkDatasets[[org]]), function(bkID){
  # bkID <- names(bulkDatasets[[org]])[1]
  
  if(is.na(bulkDatasets[[org]][[bkID]])){
    return(NULL)
  }
  
  lapply(names(scDatasets[[org]]), function(scID){
    # scID <- names(scDatasets[[org]])[1]
    
    lapply(scDatasets[[org]][[scID]]$clusterCols[1], function(clCol){
      # clCol <- scDatasets[[org]][[scID]]$clusterCols[1]
      
      cat("\nComputing scDec for conditions: \tsc dataset id:\t\t", scID, "\n\t\t\t\t\tbulk condition:\t\t", bkID, "\n\t\t\t\t\tcell annotations:\t", clCol, "\n")
      file_1 <- paste0(saveDir_music2_table, "/", org, "___bulkRNAseq_", bkID, "___scRNAseq_", scID, "___", clCol, "___inferred_cell_size.RDS")
      file_2 <- paste0(saveDir_music2_list, "/", org, "___bulkRNAseq_", bkID, "___scRNAseq_", scID, "___", clCol, "___inferred_cell_size___iterations_result_list.RDS")
      if(file.exists(file_1) & file.exists(file_2) & !overwrite){
        cat("\t\t\t\t\t--> output file already present, skipping...\n")
        return()
      }
      
      tryCatch({
        bulk.scexp <- readRDS(bulkDatasets[[org]][[bkID]]$esetPath)
        sc.scexp   <- readRDS(scDatasets[[org]][[scID]]$esetPath)
        
        # remove genes with some NA in some cells. This is not compatible with 
        # MuSiC::music_basis non.zero = T, that is called by default inside 
        ## the MuSiC::music_prop function
        tmpLogi    <- unname(apply(SingleCellExperiment::counts(bulk.scexp), 1, function(qq) any(is.na(qq))))
        bulk.scexp <- bulk.scexp[!tmpLogi, ]
        tmpLogi    <- unname(apply(SingleCellExperiment::counts(sc.scexp), 1, function(qq) any(is.na(qq))))
        sc.scexp   <- sc.scexp[!tmpLogi, ]
        
        
        genes.bulk   <- rownames(SingleCellExperiment::counts(bulk.scexp))
        genes.sc     <- rownames(SingleCellExperiment::counts(sc.scexp))
        commonGenes  <- intersect(genes.bulk, genes.sc)
        commonGenes  <- commonGenes[!(commonGenes %in% c("", " ")) & !is.na(commonGenes)]
        meta.bulk    <- as.data.frame(SingleCellExperiment::colData(bulk.scexp))
        meta.sc      <- as.data.frame(SingleCellExperiment::colData(sc.scexp))
        meta.sc      <- meta.sc[!is.na(meta.sc[[clCol]]) & !(meta.sc[[clCol]] %in% c("", " ")), ]
        bulk.scexp   <- bulk.scexp[commonGenes, , drop = F]
        sc.scexp     <- sc.scexp[commonGenes, colnames(sc.scexp) %in% rownames(meta.sc), drop = F]
        # identical(colnames(sc.scexp), rownames(meta.sc))
        
        
        if(length(unique(meta.sc$sampleID)) <= 1){
          cat("\t\t\t\t\t--> less than 1 sample available, skipping...\n")
          return(NULL)
        }
        
        ########################################################################
        ## Estimation of cell type proportions
        ########################################################################
        
        bulk.mat         <- SingleCellExperiment::counts(bulk.scexp)
        bulk.control.mat <- bulk.mat[, grepl("_CD", colnames(bulk.mat))]
        bulk.case.mat    <- bulk.mat[, grepl("_HFD", colnames(bulk.mat))]
        
        
        deconv.obj.list <- parallel::mclapply(mc.cores = nbCores.external, X = 1:nb_retry, FUN = function(x){
          deconv.obj.x <- music2_prop(bulk.control.mtx = bulk.control.mat,
                                      bulk.case.mtx    = bulk.case.mat,
                                      sc.sce           = sc.scexp,
                                      clusters         = clCol,
                                      samples          = "sampleID",
                                      select.ct        = NULL,
                                      cell_size        = NULL,
                                      method           = "t_stats",
                                      maxiter          = 200,
                                      sample_prop      = 0.5,
                                      n_resample       = 40,
                                      eps_c            = 0.01,
                                      eps_r            = 0.005,
                                      nb_cores         = nbCores.internal)
          deconv.obj.x
        })
        tmpLogi <- unlist(lapply(deconv.obj.list, function(x) x$convergence))
        
        prop_mat.list <- lapply(deconv.obj.list, function(x) x$Est.prop)
        rows_order    <- rownames(prop_mat.list[[1]])
        cols_order    <- colnames(prop_mat.list[[1]])
        prop_mat.list <- lapply(prop_mat.list, function(x) x[match(rows_order, rownames(x)), match(cols_order, colnames(x))])
        
        
        prop_array  <- array(unlist(prop_mat.list), dim = c(nrow(prop_mat.list[[1]]), ncol(prop_mat.list[[1]]), length(prop_mat.list)))
        prop_mean   <- as.data.frame(apply(prop_array, c(1, 2), mean))
        prop_sd     <- as.data.frame(apply(prop_array, c(1, 2), sd))
        
        colnames(prop_mean)     <- colnames(prop_sd) <- cols_order
        rownames(prop_mean)     <- rownames(prop_sd) <- rows_order
        prop_mean$bulk_sampleID <- rownames(prop_mean)
        prop_sd$bulk_sampleID   <- rownames(prop_sd)
        
        prop_mean <- type.convert(reshape2::melt(prop_mean, id.vars = "bulk_sampleID", variable.name = "celltype", value.name = "proportion"), as.is = T)
        prop_sd   <- type.convert(reshape2::melt(prop_sd, id.vars = "bulk_sampleID", variable.name = "celltype", value.name = "poportion_sd"), as.is = T)
        deconv.df <- merge(prop_mean, prop_sd, by = c("bulk_sampleID", "celltype"))
        deconv.df <- merge(deconv.df, meta.bulk, by.x = "bulk_sampleID", by.y = "row.names", all.x = T, all.y = F)
        
        deconv.df$nbCells                <- nrow(meta.sc)
        deconv.df$nb_iterations_converse <- sum(tmpLogi)
        
        saveRDS(deconv.df, paste0(saveDir_music2_table, "/", org, "___bulkRNAseq_", bkID, "___scRNAseq_", scID, "___", clCol, "___inferred_cell_size.RDS"))
        saveRDS(deconv.obj.list, paste0(saveDir_music2_list, "/", org, "___bulkRNAseq_", bkID, "___scRNAseq_", scID, "___", clCol, "___inferred_cell_size___iterations_result_list.RDS"))
        
        remove(deconv.obj.list, deconv.df, genes.bulk, genes.sc, meta.bulk, meta.sc, bulk.scexp, sc.scexp, commonGenes)
      }, error = function(e){
        print(paste0("ERROR for analysis:", org, "___bulkRNAseq_", bkID, "___scRNAseq_", scID, "___", clCol, "___inferred_cell_size"))
      })
      NA
    }) 
    NA
  })
  NA
})




################
## PCA
################

saveDir.pca <- "./Data/singleCell_deconvolution/cells_proportion_pca"
if(!dir.exists(saveDir.pca)){
  dir.create(saveDir.pca)
}

df_prop     <- readRDS("./Data/singleCell_deconvolution/deconvolutionResults_music2/prop_tables/mouse___bulkRNAseq_CD_HFD___scRNAseq_GSE109774__facs__heart__all_cell_types__all_genes___celltype1___inferred_cell_size.RDS")
save_suffix <- paste0("___nbIt_", df_prop$nb_iterations_converse[1])


pca.obj.list <- lapply(c("CD", "HFD", "CD_HFD"), function(dd){
  # dd <- "CD"
  
  # Compute pca for average strain_diet conditions (few strains have duplicated RNAseq samples)
  if(dd != "CD_HFD"){
    df_pca <- reshape2::dcast(df_prop[df_prop$diet == dd, ], sampleID~celltype, value.var = "proportion")
  } else{
    df_pca <- reshape2::dcast(df_prop, sampleID~celltype, value.var = "proportion")
  }
  df_pca$sampleID       <- gsub("_HFD.*", "_HFD", gsub("_CD.*", "_CD", df_pca$sampleID))
  # sum(duplicated(df_pca$sampleID))
  df_pca.aggr           <- data.table(df_pca)[, lapply(.SD, mean), by = "sampleID"]
  df_pca.aggr           <- type.convert(as.data.frame(df_pca.aggr), as.is = T)
  rownames(df_pca.aggr) <- df_pca.aggr$sampleID
  df_pca.aggr           <- df_pca.aggr[, -1]
  
  pca.obj          <- prcomp(df_pca.aggr, center = T, scale. = T)
  pca.obj.summary  <- summary(pca.obj)
  pca.explPerc     <- pca.obj.summary$importance[2, ] * 100
  cumsum(pca.explPerc)
  
  saveRDS(pca.obj$x, paste0(saveDir.pca, "/avg_pca_projections___music2___mouse___bulkRNAseq_", dd, "___scRNAseq_GSE109774__facs__heart__all_cell_types__all_genes___celltype1___inferred_cell_size.RDS"))
  saveRDS(pca.obj$rotation, paste0(saveDir.pca, "/avg_pca_eigenvectors___music2___mouse___bulkRNAseq_", dd, "___scRNAseq_GSE109774__facs__heart__all_cell_types__all_genes___celltype1___inferred_cell_size.RDS"))
  saveRDS(pca.explPerc, paste0(saveDir.pca, "/avg_pca_explained_var___music2___mouse___bulkRNAseq_", dd, "___scRNAseq_GSE109774__facs__heart__all_cell_types__all_genes___celltype1___inferred_cell_size.RDS"))
  saveRDS(pca.obj, paste0(saveDir.pca, "/avg_pca_object___music2___mouse___bulkRNAseq_", dd, "___scRNAseq_GSE109774__facs__heart__all_cell_types__all_genes___celltype1___inferred_cell_size.RDS"))
})
names(pca.obj.list) <- c("CD", "HFD", "CD_HFD")






