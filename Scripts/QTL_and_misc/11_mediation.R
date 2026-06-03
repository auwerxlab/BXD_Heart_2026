
library(qtl2)
library(data.table)


source("./Scripts/QTL_and_misc/QTL_fun.R")

################################################################################
## select lipids to run for full mediation scan testing at given lQTL locus
## every possible gene (forward) or lipid (reverse) mediator for the lQTL lipid
## or colocalizing (e/p)QTL gene target respectively
## As the full scan is computationally intensive, we also perform a reduced scan
## testing only specific colocalizing lipid/gene QTL pairs only at the lQTL locus
## Therefore this script does:
## (1) Reduced mediation scan for all colocalizing lipid/gene QTL pairs
## (2) For specific lipids, a full mediation scan testing all possible mediators 
##     expressed in LV - seleted lipids specified here
################################################################################

lipids_run <- c("HexCer[NS] d18:0_22:1_RT_12.391", "HexCer[NS] d41:1_RT_12.925", "HexCer[NS] d42:2_RT_12.286")

cat("\n Full mediation scan analysis run for all the possible gene (forward) and lipid (reverse) mediators lQTL peaks for the following lipids and their colocalizing (e/p)QTL gene targets respectively:\n",
    paste(paste0("\t--> ", lipids_run), collapse = "\n"), "\n",
    "!!! Modify lipid selection is mediation scripts to run for other lipids.\n\n")




################################################################################
## qtl2 LOD drop "mediation" (conditioning) analysis - reduced scan limited at 
## lipid/gene features colocalizing through QTL peak at lQTL locus
## genotype > gene > lipid
## or
## genotype > lipid > gene
## mediation model for each eQTL/lQTL colocalization
## always at lQTL locus
## For now exclude proteins (anyway much less pQTL peaks than eQTL peaks)
################################################################################

lipid_QTL_result_files       <- list.files("./Data/QTL_mapping/qtl_outputs", recursive = T, full.names = T, pattern = "average_lipidomic_derived|individuals_lipidomic_derived|individuals_lipidomic_Filt")
lipid_QTL_result_files       <- lipid_QTL_result_files[!grepl("previous_runs_no_diet_stratification", lipid_QTL_result_files)]
lipid_QTL_result_files       <- lipid_QTL_result_files[grepl("qtl2_results__", lipid_QTL_result_files)]
diet_vec                     <- gsub(".*Interaction__|__.*", "", lipid_QTL_result_files)
interaction_vec              <- gsub("__.*" , "", gsub(".*noInt", "noInt", gsub(".*dietInt", "dietInt", lipid_QTL_result_files)))
data_type_vec                <- ifelse(grepl("average", lipid_QTL_result_files), "average", "individuals")
lipid_QTL_result_files       <- split(lipid_QTL_result_files, paste0(diet_vec, "__", interaction_vec, "__", data_type_vec))
lipid_QTL_result_files       <- lipid_QTL_result_files[sort(names(lipid_QTL_result_files))]
lipid_QTL_result_files       <- lipid_QTL_result_files[!grepl("indiv", names(lipid_QTL_result_files))]
lipid_QTL_result_file_traits <- get_list_QTL_file_traits(lipid_QTL_result_files)
remove(lipid_QTL_result_files)

e_QTL_pheno_files        <- list.files("./Data/QTL_mapping/qtl_processed_inputs/RNAseq_cf_pca", recursive = T, full.names = T, pattern = "pheno_original_distribution\\.")
names(e_QTL_pheno_files) <- gsub(".*\\/", "", gsub("\\/BXD.*|\\/pheno_origi.*", "", e_QTL_pheno_files))

e_QTL_cross_files        <- list.files("./Data/QTL_mapping/qtl_processed_inputs/RNAseq_cf_pca", recursive = T, full.names = T, pattern = "BXD\\.json")
names(e_QTL_cross_files) <- gsub(".*\\/", "", gsub("\\/BXD.*|\\/pheno_origi.*", "", e_QTL_cross_files))

p_QTL_cross_files        <- list.files("./Data/QTL_mapping/qtl_processed_inputs/average_coonProteomic_cf_pca", recursive = T, full.names = T, pattern = "BXD\\.json")
names(p_QTL_cross_files) <- gsub(".*\\/", "", gsub("\\/BXD.*|\\/pheno_origi.*", "", p_QTL_cross_files))

l_QTL_cross_files        <- list.files("./Data/QTL_mapping/qtl_processed_inputs/average_lipidomic_derived_avgAggr_data", recursive = T, full.names = T, pattern = "BXD\\.json")
names(l_QTL_cross_files) <- gsub(".*\\/", "", gsub("\\/BXD.*|\\/pheno_origi.*", "", l_QTL_cross_files))

# load original gene count data used for QTL mapping (prior RINT normalization, then mean-aggregate and RINT-transform)
gene_pheno_list <- lapply(e_QTL_pheno_files, function(x){
  # x <- e_QTL_pheno_files[[2]]
  
  gene_pheno.x           <- as.data.frame(readRDS(x))
  gene_pheno.x           <- cbind(strain_diet = gsub("[0-9]$", "", rownames(gene_pheno.x)), gene_pheno.x)
  gene_pheno.aggr.x      <- data.table(gene_pheno.x)[, lapply(.SD, function(zz) mean(zz, na.rm = T)), by = c("strain_diet")]
  gene_pheno.x           <- as.data.frame(gene_pheno.aggr.x)
  rownames(gene_pheno.x) <- gene_pheno.x$strain_diet
  gene_pheno.x           <- gene_pheno.x[, -1]
  strain_vec             <- gsub("_.*", "", rownames(gene_pheno.x))
  # remove genes measured in less than 50% of the strains. 
  # Note: this is a slightly different version compared to the filtering applied prior eQTL mapping, but it's a better filter
  tmpLogi                <- apply(gene_pheno.x, 2, function(zz) length(unique(strain_vec[!is.na(zz) & zz > 0]))) > length(unique(strain_vec)) * 0.5
  gene_pheno.x           <- gene_pheno.x[, tmpLogi]
  gene_pheno.x
})

protein_pheno_list <- lapply(p_QTL_cross_files, function(x){
  # x <- p_QTL_cross_files[[1]]
  
  protein_pheno.x <- qtl2::read_cross2(x)
  protein_pheno.x <- as.data.frame(protein_pheno.x$pheno)
  protein_pheno.x
})

protein_addcovar_list <- lapply(p_QTL_cross_files, function(x){
  # x <- p_QTL_cross_files[[1]]
  
  protein_covar.x <- qtl2::read_cross2(x)
  protein_covar.x <- as.data.frame(protein_covar.x$covar)
  protein_covar.x <- type.convert(protein_covar.x[, grepl("cf__", colnames(protein_covar.x))], as.is = T)
  protein_covar.x
})

lipid_pheno_list <- lapply(l_QTL_cross_files, function(x){
  # x <- l_QTL_cross_files[[1]]
  
  lipid_pheno.x <- qtl2::read_cross2(x)
  lipid_pheno.x <- as.data.frame(lipid_pheno.x$pheno)
  lipid_pheno.x <- lipid_pheno.x[, !grepl("Unknown", colnames(lipid_pheno.x), ignore.case = T)]
  lipid_pheno.x
})


# load lipid/genes colocalization tables
all_genes_criteria_files        <- list.files("./Data/gene_prioritization/prioritization_output_lipids", full.names = T, pattern = "all_genes_criteria")
all_genes_criteria_files        <- all_genes_criteria_files[!grepl("derived", all_genes_criteria_files)]
names(all_genes_criteria_files) <- gsub(".*criteria___|__average.*", "", all_genes_criteria_files)
stopifnot(length(all_genes_criteria_files) == 3)
all_genes_criteria.list         <- lapply(all_genes_criteria_files, readRDS)

# get all overlapping lQTL/(e/p)QTL couples
nbCores       <- 10
overlap_pairs <- lapply(all_genes_criteria.list, function(x){
  # x <- all_genes_criteria.list[[1]]
  out <- x[!unlist(lapply(x$eGenes_peaks_table, is.null)), ]
  # filter only genes with cis or trans eQTLs as for now we don't care about pQTLs
  out <- out[out$cis_eQTL_overlap | out$trans_eQTL_overlap, ]
  out <- parallel::mclapply(mc.cores = nbCores, X = 1:nrow(out), FUN = function(y){
    # y <- 113
    df.y.1           <- dplyr::select(out[y, ], dplyr::all_of(c("ref_qtl_id", "ref_qtl_chr", "ref_qtl_inputType", "ref_qtl_data_type", "ref_qtl_pos_Mb", "ref_qtl_pVal",  "condition", "interaction", "analysis_id", 
                                                                "gene_id_full", "gene_id", "gene_chr", "gene_posStart_Mb", "gene_posEnd_Mb", "gene_start_dist_from_pk_Mb", "gene_is_under_ref_qtl", "gene_symbol")))
    colnames(df.y.1) <- gsub("ref_qtl_", "", colnames(df.y.1))
    df.y.2           <- unique(dplyr::select(out$eGenes_peaks_table[[y]], dplyr::all_of(c("id", "qtl_type"))))
    df.y.2$feature_type  <- ifelse(grepl("eQTL", df.y.2$qtl_type), "mRNA", "protein")
    colnames(df.y.2) <- c("feature_id", "qtl_type", "feature_type")
    df.y.2.aggr      <- data.table(df.y.2)[, list(eGene_qtl_types = paste(sort(unique(qtl_type)), collapse = ";")), by = c("feature_id", "feature_type")]
    cbind(df.y.1, df.y.2.aggr)
  })
  out <- rbindlist(out)
  # remove proteins as for now we don't care about proteins
  out <- out[out$feature_type == "mRNA", ]
  out
})
overlap_pairs <- unique(as.data.frame(rbindlist(overlap_pairs)))
if(nrow(overlap_pairs) == 0){
  stop("No lQTL/(e/p)QTL overlapping pairs found for the selected lipids!!!")
}




nbCores    <- 20
overwrite  <- F
tt         <- parallel::mclapply(mc.cores = nbCores, X = 1:nrow(overlap_pairs), FUN = function(x){
  # x <- 1
  
  saveDir.qtldrop.all <- paste0("./Data/lipidomics/mediation_all_eGene_lQTL_colocalized/qtl_drop_method/", overlap_pairs$analysis_id[x])
  if(!dir.exists(saveDir.qtldrop.all)){
    dir.create(saveDir.qtldrop.all, recursive = T)
  }
  
  
  cat("Running mediation analysis for entry", x, "/", nrow(overlap_pairs), "...\n")
  
  unique_id.x <- paste0(overlap_pairs$id[x], "__", overlap_pairs$analysis_id[x], "__chr", 
                        overlap_pairs$chr[x], "@", round(overlap_pairs$pos_Mb[x], 2), "__", 
                        overlap_pairs$peak_type[x], "__", overlap_pairs$gene_id[x])
  out_file    <-  paste0(saveDir.qtldrop.all, "/qtl_lod_drop_meditation___lQTL_locus___", gsub("\\.|\\:|\\-| ", "_",  unique_id.x), ".RDS")
  if(file.exists(out_file) & !overwrite){
    cat("\t--> results already found! Skipping iteration...\n")
    return()
  }
  
  lipid_pheno <- lipid_pheno_list[[overlap_pairs$condition[x]]]
  lipid_pheno <- lipid_pheno[, overlap_pairs$id[x], drop = F]
  # keep only strain/diet conditions without missing lipid data
  lipid_pheno <- lipid_pheno[!is.na(lipid_pheno[, 1]), , drop = F]
  
  # get original gene count data, mean-aggregate by strain/diet condition, filter only strain/diet conditions without missing lipid data
  gene_pheno           <- gene_pheno_list[[overlap_pairs$condition[x]]]
  if(!(overlap_pairs$gene_id[x] %in% colnames(gene_pheno))){
    # Some genes removed upstream when initializing gene_pheno_list (only genes with non-zero expression in at least half of strains retained). Better not to run mediation as "low-quality" genes
    return()
  }
  gene_pheno           <- gene_pheno[rownames(lipid_pheno), overlap_pairs$gene_id[x], drop = F]
  gene_pheno           <- gene_pheno[!is.na(gene_pheno[, 1]), , drop = F]
  gene_norm            <- transform2normal_byCol_INT(gene_pheno)
  gene_norm            <- gene_norm$transf.df
  rownames(gene_norm)  <- rownames(gene_pheno)
  gene_pheno           <- gene_norm
  lipid_pheno          <- lipid_pheno[rownames(gene_pheno), , drop = F]
  
  
  lQTL.obj    <- lipid_QTL_result_file_traits[[overlap_pairs$analysis_id[x]]]
  tmpLogi     <- unlist(lapply(lQTL.obj, function(y) overlap_pairs$id[x] %in% y))
  lQTL.obj    <- readRDS(names(lQTL.obj)[tmpLogi])
  stopifnot(sum(tmpLogi) == 1)
  
  int_covar <- lQTL.obj$intcovar
  if(!is.null(int_covar)){
    int_covar <- int_covar[names(int_covar) %in% rownames(lipid_pheno)]
    int_covar <- int_covar[match(rownames(lipid_pheno), names(int_covar))]
  }
  
  pr              <- lQTL.obj$pr
  pr              <- lapply(pr, function(zz) zz[rownames(lipid_pheno), , ])
  attributes(pr)  <- attributes(lQTL.obj$pr)
  
  kin             <- lQTL.obj$kinship
  kin             <- lapply(kin, function(zz) zz[rownames(lipid_pheno), rownames(lipid_pheno)])
  
  pk_idx           <- lQTL.obj$pmap[[overlap_pairs$chr[x]]]
  dist             <- abs(pk_idx - overlap_pairs$pos_Mb[x])
  pk_pos           <- unname(pk_idx)[which.min(dist)]
  pk_idx           <- names(pk_idx)[which.min(dist)]
  stopifnot(sum(dist < 1e-7) == 1) 
  
  # scan only at peak location
  # adapted from qtl2::pull_genoprobint source code
  pr.use      <- pr[, overlap_pairs$chr[x]]
  pr.use[[1]] <- pr.use[[1]][, , pk_idx, drop = F]
  kin.use     <- kin[names(kin) == overlap_pairs$chr[x]]
  
  
  add_covar <- protein_addcovar_list[[overlap_pairs$condition[x]]]
  add_covar <- add_covar[rownames(lipid_pheno), ]
  
  out_list.x <- lapply(c("direct", "reverse"), function(dd){
    # dd <- "direct"
    # dd <- "reverse"
    
    if(dd == "direct"){
      target_pheno   <- lipid_pheno
      mediator_pheno <- gene_pheno
      target_id      <- "lipid_id"
      mediator_id    <- "gene_id"
    } else if(dd == "reverse"){
      target_pheno   <- gene_pheno
      mediator_pheno <- lipid_pheno
      target_id      <- "gene_id"
      mediator_id    <- "lipid_id"
    }
    
    # conditional additive covariate matrix (with mediator added as column)
    if(!is.null(add_covar)){
      stopifnot(identical(rownames(add_covar), rownames(mediator_pheno)))
      add_covar_con  <- cbind(add_covar, mediator_pheno)
    } else{
      add_covar_con  <- mediator_pheno
    }
    
    stopifnot(identical(rownames(target_pheno), rownames(mediator_pheno)))
    stopifnot(identical(rownames(target_pheno), rownames(add_covar_con)))
    if(!is.null(add_covar)){
      stopifnot(identical(rownames(target_pheno), rownames(add_covar)))
    }
    if(!is.null(int_covar)){
      stopifnot(all(rownames(target_pheno) == names(int_covar)))
    }
    stopifnot(all(unlist(lapply(kin, function(y) identical(rownames(target_pheno), colnames(y))))))
    stopifnot(all(unlist(lapply(pr, function(y) identical(rownames(target_pheno), rownames(y))))))
    
    # Unconditioned (we already have oscan1_unc from the saved QTL mapping results, but recompute on samples with lipidomic data (NA-excluded) for fairness)
    oscan1_unc <- qtl2::scan1(genoprobs = pr.use, pheno = target_pheno, kinship = kin.use, addcovar = add_covar, intcovar = int_covar, cores = 1)
    lod_pk_unc <- oscan1_unc[pk_idx, ]
    
    # Conditional on mediator
    oscan1_con <- qtl2::scan1(genoprobs = pr.use, pheno = target_pheno, kinship = kin.use, addcovar = add_covar_con, intcovar = int_covar, cores = 1)
    lod_pk_con <- oscan1_con[pk_idx, ]
    
    lod_drop.df <- data.frame(mediator_id          = colnames(mediator_pheno),
                              target_id            = colnames(target_pheno),
                              mediation_direction  = dd,
                              lod_unc              = lod_pk_unc,
                              lod_cond             = lod_pk_con, 
                              lod_drop             = lod_pk_con - lod_pk_unc,
                              peak_chr             = overlap_pairs$chr[x],
                              peak_pos             = pk_pos,
                              peak_marker          = pk_idx,
                              lipid_input_type     = overlap_pairs$inputType[x],
                              lipid_data_type      = overlap_pairs$data_type[x],
                              eGene_type           = overlap_pairs$eGene_qtl_types[x],
                              eGene_chr            = overlap_pairs$gene_chr[x],
                              eGene_start          = overlap_pairs$gene_posStart_Mb[x],
                              eGene_end            = overlap_pairs$gene_posEnd_Mb[x],
                              eGene_lQTL_peak_dist = overlap_pairs$gene_start_dist_from_pk_Mb[x],
                              eGene_under_lQTL     = overlap_pairs$gene_is_under_ref_qtl[x],
                              eGene_id             = overlap_pairs$gene_id[x],
                              eGene_symbol         = overlap_pairs$gene_symbol[x],
                              condition            = overlap_pairs$condition[x],
                              interaction          = overlap_pairs$interaction[x],
                              analysis_id          = overlap_pairs$analysis_id[x],
                              has_cf_covar         = T,
                              stringsAsFactors = F)
    
    lod_drop.df
  })
  lod_drop.df.all <- do.call(rbind, out_list.x)
  saveRDS(lod_drop.df.all, out_file)
})




################################################################################
## qtl2 "mediation"/conditioning analysis - full scan as indicated in the comment 
## above
## genotype > gene (all expressed ones) > lipid
## or
## genotype > lipid (all measured ones) > gene
## always at lQTL locus
## For every mediation scan, test all the models (CD, HFD, CD_HFD, CD_HFD GxD) 
## even if the initial QTL overlap was found just in one - to have the results in
## case of downstream need
################################################################################


geneTable                           <- readRDS("./Data/RNAseq_processing/geneConversionTables/geneConversionTable_GRCm38_release-102_complete_gtf.RDS")
geneTable                           <- type.convert(as.data.frame(geneTable), as.is = T)
geneTable                           <- geneTable[geneTable$seqnames %in% c(1:20, "X", "Y", "MT"), ]
geneTable                           <- geneTable[geneTable$type == "gene", ]
geneTable                           <- dplyr::select(geneTable, c("gene_id", "gene_name", "seqnames", "start", "end"))
colnames(geneTable)                 <- c("gene_id_ensembl", "gene_symbol", "gene_chr", "chromosome_start_position", "chromosome_end_position")
geneTable                           <- unique(geneTable)
geneTable$chromosome_start_position <- as.numeric(geneTable$chromosome_start_position)
geneTable$chromosome_end_position   <- as.numeric(geneTable$chromosome_end_position)
geneTable$posStart_Mb               <- geneTable$chromosome_start_position / 1e6
geneTable$posEnd_Mb                 <- geneTable$chromosome_end_position / 1e6
geneTable                           <- dplyr::select(geneTable, -dplyr::all_of(c("chromosome_start_position", "chromosome_end_position")))
stopifnot(sum(duplicated(geneTable$gene_id_ensembl)) == 0)

protTable_Coon                          <- readRDS("./Data/input_data/proteomics/coon_formatted/proteins_metadata.RDS")
colnames(protTable_Coon)                <- c("protein_id_uniprot", "peptide_id_ensembl", "gene_id_entrez", "gene_id_ensembl", "transcript_id_ensembl", "gene_symbol", "gene_symbol_uniprot")
protTable_Coon$protein_id_iso_uniprot   <- protTable_Coon$protein_id_uniprot
protTable_Coon$protein_id_noIso_uniprot <- gsub("-.", "", protTable_Coon$protein_id_uniprot)
protTable_Coon                          <- unique(dplyr::select(protTable_Coon, c("protein_id_uniprot", "protein_id_iso_uniprot", "protein_id_noIso_uniprot", "gene_id_ensembl")))
protTable_Coon                          <- unique(merge(protTable_Coon, geneTable, by = "gene_id_ensembl", all.x = T, all.y = F))
protTable_Coon                          <- protTable_Coon[!is.na(protTable_Coon$gene_id_ensembl) & !(protTable_Coon$gene_id_ensembl %in% c("", " ", "NA")) & !is.na(protTable_Coon$posStart_Mb), ]
protTable_Coon                          <- unique(dplyr::select(protTable_Coon, c("gene_id_ensembl", "protein_id_uniprot")))
protTable_Coon$gene_symbol              <- plyr::mapvalues(protTable_Coon$gene_id_ensembl, from = geneTable$gene_id_ensembl, to = geneTable$gene_symbol, warn_missing = F)
protTable_Coon$gene_symbol[grepl("ENS", protTable_Coon$gene_symbol)] <- NA
protTable_Coon                          <- merge(protTable_Coon, dplyr::select(geneTable, -dplyr::all_of("gene_symbol")), by = "gene_id_ensembl", all.x = T, all.y = F)
protTable_Coon.aggr                     <- data.table(protTable_Coon)[, list(gene_id_ensembl_all = paste(gene_id_ensembl, collapse = ";"),
                                                                             gene_symbol_all     = paste(gene_symbol, collapse = ";"),
                                                                             gene_chr_all        = paste(gene_chr, collapse = ";"),
                                                                             posStart_Mb_all     = paste(posStart_Mb, collapse = ";"),
                                                                             posEnd_Mb_all       = paste(posEnd_Mb, collapse = ";")), by = c("protein_id_uniprot")]
protTable_Coon                          <- as.data.frame(protTable_Coon.aggr)
colnames(protTable_Coon)                <- gsub("_all", "", colnames(protTable_Coon))
stopifnot(sum(duplicated(protTable_Coon$protein_id_uniprot)) == 0)


lipid_QTL_result_files       <- list.files("./Data/QTL_mapping/qtl_outputs", recursive = T, full.names = T, pattern = "average_lipidomic_derived|individuals_lipidomic_derived|individuals_lipidomic_Filt")
lipid_QTL_result_files       <- lipid_QTL_result_files[grepl("qtl2_results__", lipid_QTL_result_files)]
diet_vec                     <- gsub(".*Interaction__|__.*", "", lipid_QTL_result_files)
interaction_vec              <- gsub("__.*" , "", gsub(".*noInt", "noInt", gsub(".*dietInt", "dietInt", lipid_QTL_result_files)))
data_type_vec                <- ifelse(grepl("average", lipid_QTL_result_files), "average", "individuals")
lipid_QTL_result_files       <- split(lipid_QTL_result_files, paste0(diet_vec, "__", interaction_vec, "__", data_type_vec))
lipid_QTL_result_files       <- lipid_QTL_result_files[sort(names(lipid_QTL_result_files))]
lipid_QTL_result_files       <- lipid_QTL_result_files[!grepl("indiv", names(lipid_QTL_result_files))]
lipid_QTL_result_file_traits <- get_list_QTL_file_traits(lipid_QTL_result_files)

e_QTL_pheno_files        <- list.files("./Data/QTL_mapping/qtl_processed_inputs/RNAseq_cf_pca", recursive = T, full.names = T, pattern = "pheno_original_distribution\\.")
names(e_QTL_pheno_files) <- gsub(".*\\/", "", gsub("\\/BXD.*|\\/pheno_origi.*", "", e_QTL_pheno_files))

e_QTL_cross_files        <- list.files("./Data/QTL_mapping/qtl_processed_inputs/RNAseq_cf_pca", recursive = T, full.names = T, pattern = "BXD\\.json")
names(e_QTL_cross_files) <- gsub(".*\\/", "", gsub("\\/BXD.*|\\/pheno_origi.*", "", e_QTL_cross_files))

p_QTL_cross_files        <- list.files("./Data/QTL_mapping/qtl_processed_inputs/average_coonProteomic_cf_pca", recursive = T, full.names = T, pattern = "BXD\\.json")
names(p_QTL_cross_files) <- gsub(".*\\/", "", gsub("\\/BXD.*|\\/pheno_origi.*", "", p_QTL_cross_files))

l_QTL_cross_files        <- list.files("./Data/QTL_mapping/qtl_processed_inputs/average_lipidomic_derived_avgAggr_data", recursive = T, full.names = T, pattern = "BXD\\.json")
names(l_QTL_cross_files) <- gsub(".*\\/", "", gsub("\\/BXD.*|\\/pheno_origi.*", "", l_QTL_cross_files))


# load lipid/genes colocalization tables
all_genes_criteria_files        <- list.files("./Data/gene_prioritization/prioritization_output_lipids", full.names = T, pattern = "all_genes_criteria")
all_genes_criteria_files        <- all_genes_criteria_files[!grepl("derived", all_genes_criteria_files)]
names(all_genes_criteria_files) <- gsub(".*criteria___|__average.*", "", all_genes_criteria_files)
stopifnot(length(all_genes_criteria_files) == 3)
all_genes_criteria.list         <- lapply(all_genes_criteria_files, readRDS)

# get all overlapping lQTL/(e/p)QTL couples
nbCores       <- 10
overlap_pairs <- lapply(all_genes_criteria.list, function(x){
  # x <- all_genes_criteria.list[[1]]
  out <- x[!unlist(lapply(x$eGenes_peaks_table, is.null)), ]
  out <- parallel::mclapply(mc.cores = nbCores, X = 1:nrow(out), FUN = function(y){
    # y <- 113
    df.y.1           <- dplyr::select(out[y, ], dplyr::all_of(c("ref_qtl_id", "ref_qtl_chr", "ref_qtl_pos_Mb", "ref_qtl_pVal", "analysis_id", "gene_id_full")))
    colnames(df.y.1) <- gsub("ref_qtl_", "", colnames(df.y.1))
    df.y.2           <- unique(dplyr::select(out$eGenes_peaks_table[[y]], dplyr::all_of(c("id", "qtl_type"))))
    df.y.2$qtl_type  <- ifelse(grepl("eQTL", df.y.2$qtl_type), "mRNA", "protein")
    colnames(df.y.2) <- c("feature_id", "feature_type")
    cbind(df.y.1, df.y.2)
  })
  rbindlist(out)
})
overlap_pairs <- unique(as.data.frame(rbindlist(overlap_pairs)))
overlap_pairs <- overlap_pairs[overlap_pairs$id %in% lipids_run, ]
if(nrow(overlap_pairs) == 0){
  stop("No lQTL/(e/p)QTL overlapping pairs found for the selected lipids!!!")
}


# load original gene count data used for QTL mapping (prior RINT normalization, then mean-aggregate and RINT-transform)
gene_pheno_list <- lapply(e_QTL_pheno_files, function(x){
  # x <- e_QTL_pheno_files[[1]]
  
  gene_pheno.x           <- as.data.frame(readRDS(x))
  gene_pheno.x           <- cbind(strain_diet = gsub("[0-9]$", "", rownames(gene_pheno.x)), gene_pheno.x)
  gene_pheno.aggr.x      <- data.table(gene_pheno.x)[, lapply(.SD, function(zz) mean(zz, na.rm = T)), by = c("strain_diet")]
  gene_pheno.x           <- as.data.frame(gene_pheno.aggr.x)
  rownames(gene_pheno.x) <- gene_pheno.x$strain_diet
  gene_pheno.x           <- gene_pheno.x[, -1]
  strain_vec             <- gsub("_.*", "", rownames(gene_pheno.x))
  tmpLogi                <- apply(gene_pheno.x, 2, function(zz) length(unique(strain_vec[!is.na(zz) & zz > 0]))) > length(unique(strain_vec)) * 0.5
  gene_pheno.x           <- gene_pheno.x[, tmpLogi]
  gene_pheno.x
})

protein_pheno_list <- lapply(p_QTL_cross_files, function(x){
  # x <- p_QTL_cross_files[[1]]
  
  protein_pheno.x <- qtl2::read_cross2(x)
  protein_pheno.x <- as.data.frame(protein_pheno.x$pheno)
  protein_pheno.x
})

protein_addcovar_list <- lapply(p_QTL_cross_files, function(x){
  # x <- p_QTL_cross_files[[1]]
  
  protein_covar.x <- qtl2::read_cross2(x)
  protein_covar.x <- as.data.frame(protein_covar.x$covar)
  protein_covar.x <- type.convert(protein_covar.x[, grepl("cf__", colnames(protein_covar.x))], as.is = T)
  protein_covar.x
})

lipid_pheno_list <- lapply(l_QTL_cross_files, function(x){
  # x <- l_QTL_cross_files[[1]]
  
  lipid_pheno.x <- qtl2::read_cross2(x)
  lipid_pheno.x <- as.data.frame(lipid_pheno.x$pheno)
  lipid_pheno.x <- lipid_pheno.x[, !grepl("Unknown", colnames(lipid_pheno.x), ignore.case = T)]
  lipid_pheno.x
})


saveDir.qtldrop  <- "./Data/lipidomics/mediation/qtl_drop_method"
for(i in c(saveDir.qtldrop)){
  if(!dir.exists(i)){
    dir.create(i, recursive = T)
  }
}

nbCores_mediatiors_scan <- 40
overwrite               <- F
qtl_scan_only_at_pk     <- T # speed up computations
tt                      <- lapply(1:nrow(overlap_pairs), function(x){
  # x <- 1
  
  cat("\nRunning mediation analysis for entry", x, "/", nrow(overlap_pairs), "...\n")
  
  
  #######################################
  ## Direct (forward) mediation:
  ## Genotype > Gene/Protein > Lipid
  ## Reverse (backward) mediation:
  ## Genotype > Lipid > Gene/Protein
  #######################################
  
  
  tt <- lapply(c("direct", "reverse"), function(dd){
    # dd <- "direct"
    # dd <- "reverse"
    
    lQTL.obj <- lapply(lipid_QTL_result_file_traits, function(kk){
      # kk <- lipid_QTL_result_file_traits[[1]]
      tmpLogi  <- unlist(lapply(kk, function(y) overlap_pairs$id[x] %in% y))
      stopifnot(sum(tmpLogi) <= 1)
      if(sum(tmpLogi) == 0){ return() }
      readRDS(names(kk)[tmpLogi])
    })
    
    
    lipid_pheno <- lipid_pheno_list
    if(dd == "direct"){
      lipid_pheno <- lapply(lipid_pheno, function(kk){
        # kk <- lipid_pheno$CD
        out <- kk[, overlap_pairs$id[x], drop = F]
        # keep only strain/diet conditions without missing lipid data
        out <- out[!is.na(out[, 1]), , drop = F]
      })
    }
    
    if(overlap_pairs$feature_type[x] == "mRNA"){
      
      # get original gene count data, mean-aggregate by strain/diet condition, filter only strain/diet conditions without missing lipid data
      gene_pheno <- gene_pheno_list
      gene_pheno <- lapply(names(gene_pheno), function(kk){
        # kk <- names(gene_pheno)[1]
        out <- gene_pheno[[kk]][match(rownames(lipid_pheno[[kk]]), rownames(gene_pheno[[kk]])), ]
        if(dd == "reverse"){
          if(!(overlap_pairs$feature_id[x] %in% colnames(out))){ return() }
          out <- out[, overlap_pairs$feature_id[x], drop = F]
          out <- out[!is.na(out[, 1]), , drop = F]
        }
        gene_norm            <- transform2normal_byCol_INT(out)
        gene_norm            <- gene_norm$transf.df
        rownames(gene_norm)  <- rownames(out)
        out                  <- gene_norm
        out
      }) 
      names(gene_pheno) <- names(gene_pheno_list)
      
    } else{
      
      gene_pheno <- protein_pheno_list
      gene_pheno <- lapply(names(gene_pheno), function(kk){
        # kk <- names(gene_pheno)[1]
        out <- gene_pheno[[kk]][match(rownames(lipid_pheno[[kk]]), rownames(gene_pheno[[kk]])), ]
        if(dd == "reverse"){
          if(!(overlap_pairs$feature_id[x] %in% colnames(out))){ return() }
          out <- out[, overlap_pairs$feature_id[x], drop = F]
          out <- out[!is.na(out[, 1]), , drop = F]
        }
        out
      })
      names(gene_pheno) <-  names(protein_pheno_list)
      
    }
    
    lipid_pheno <- lapply(names(lipid_pheno), function(kk){
      lipid_pheno[[kk]][rownames(lipid_pheno[[kk]]) %in% rownames(gene_pheno[[kk]]), , drop = F]
    })
    names(lipid_pheno) <- names(lipid_pheno_list)
    
    
    if(dd == "direct"){
      target_pheno   <- lipid_pheno
      mediator_pheno <- gene_pheno
      target_id      <- "lipid_id"
      mediator_id    <- "gene_id"
    } else if(dd == "reverse"){
      target_pheno   <- gene_pheno
      mediator_pheno <- lipid_pheno
      target_id      <- "gene_id"
      mediator_id    <- "lipid_id"
    }
    
    
    int_covar <- lapply(names(lQTL.obj), function(kk){
      # kk <- names(lQTL.obj)[2]
      out <- lQTL.obj[[kk]]$intcovar
      if(!is.null(out)){
        out <- out[names(out) %in% rownames(lipid_pheno[[gsub("__.*", "", kk)]])]
        out <- out[match(rownames(lipid_pheno[[gsub("__.*", "", kk)]]), names(out))]
      }
    })
    names(int_covar) <- names(lQTL.obj)
    
    add_covar <- lapply(names(protein_addcovar_list), function(kk){
      protein_addcovar_list[[kk]][match(rownames(lipid_pheno[[gsub("__.*", "", kk)]]), rownames(protein_addcovar_list[[kk]])), ]
    })
    names(add_covar) <- names(protein_addcovar_list)
    
    if(!all(unlist(lapply(add_covar, is.null)))){
      # to be implemented
      add_covar_con        <- lapply(names(add_covar), function(kk){
        stopifnot(identical(rownames(add_covar[[kk]]), rownames(mediator_pheno[[kk]])))
        cbind(add_covar[[kk]], mediator_pheno[[kk]])
      })
      names(add_covar_con) <- names(add_covar)
      add_covar_cols       <- lapply(add_covar, colnames)
    } else{
      add_covar_con  <- mediator_pheno
      add_covar_cols <- lapply(add_covar, function(kk) NULL)
    }
    
    
    pr <- lapply(names(lQTL.obj), function(kk){
      # kk <- names(lQTL.obj)[1]
      out             <- lQTL.obj[[kk]]$pr
      out             <- lapply(out, function(zz) zz[rownames(lipid_pheno[[gsub("__.*", "", kk)]]), , ])
      attributes(out) <- attributes(lQTL.obj[[kk]]$pr)
      out
    })
    names(pr) <- names(lQTL.obj)
    
    
    kin <-lapply(names(lQTL.obj), function(kk){
      # kk <- names(lQTL.obj)[1]
      out <- lQTL.obj[[kk]]$kinship
      out <- lapply(out, function(zz) zz[rownames(lipid_pheno[[gsub("__.*", "", kk)]]), rownames(lipid_pheno[[gsub("__.*", "", kk)]])])
      out
    })
    names(kin) <- names(lQTL.obj)
    
    
    conditions.keep <- lapply(list(target_pheno, mediator_pheno, add_covar_con), function(kk){
      # kk <- target_pheno
      kk      <- kk[!unlist(lapply(kk, is.null))]
      tmpLogi <- unlist(lapply(kk, function(qq) nrow(qq) > 0 & ncol(qq) > 0))
      names(kk)[tmpLogi]
    })
    conditions.keep <- unlist(unique(conditions.keep))
    
    target_pheno   <- target_pheno[names(target_pheno) %in% conditions.keep]
    mediator_pheno <- mediator_pheno[names(mediator_pheno) %in% conditions.keep]
    add_covar_con  <- add_covar_con[names(add_covar_con) %in% conditions.keep]
    lQTL.obj       <- lQTL.obj[grepl(paste(paste0("^", conditions.keep, "__"), collapse = "|"), names(lQTL.obj))]
    kin            <- kin[grepl(paste(paste0("^", conditions.keep, "__"), collapse = "|"), names(kin))]
    pr             <- pr[grepl(paste(paste0("^", conditions.keep, "__"), collapse = "|"), names(pr))]
    int_covar      <- int_covar[grepl(paste(paste0("^", conditions.keep, "__"), collapse = "|"), names(int_covar))]
    
    tt <- lapply(names(target_pheno), function(kk){
      # kk <- names(target_pheno)[3]
      stopifnot(identical(rownames(target_pheno[[kk]]), rownames(mediator_pheno[[kk]])))
      stopifnot(identical(rownames(target_pheno[[kk]]), rownames(add_covar_con[[kk]])))
      stopifnot(all(rownames(target_pheno[[kk]]) == names(int_covar[[which(grepl(paste0("^", kk, "__"), names(int_covar)))[1]]])))
      stopifnot(all(unlist(lapply(kin[[which(grepl(paste0("^", kk, "__"), names(kin)))[1]]], function(y) identical(rownames(target_pheno[[kk]]), colnames(y))))))
      stopifnot(all(unlist(lapply(pr[[which(grepl(paste0("^", kk, "__"), names(pr)))[1]]], function(y) identical(rownames(target_pheno[[kk]]), rownames(y))))))
    })
    
    # perform mediation at lQTL peak
    pk_idx           <- lQTL.obj[[overlap_pairs$analysis_id[x]]]$pmap[[overlap_pairs$chr[x]]]
    dist             <- abs(pk_idx - overlap_pairs$pos_Mb[x])
    pk_pos           <- unname(pk_idx)[which.min(dist)]
    pk_idx           <- names(pk_idx)[which.min(dist)]
    stopifnot(sum(dist < 1e-7) == 1)

  
    ##########################
    ## qtl LOD drop method at lQTL peak
    ##########################
    
    tt <- lapply(names(lQTL.obj), function(qq){
      # qq <- names(lQTL.obj)[1]
      
      out_file <-  paste0(saveDir.qtldrop, "/", dd, "_qtl_lod_drop_meditation", 
                          "___", qq, 
                          "___cf_covar",
                          "___", overlap_pairs$feature_type[x],
                          "___mediation_lQTL_locus_", gsub("\\[|\\]| |\\:|\\.", "_", overlap_pairs$id[x]),
                          "_chr", overlap_pairs$chr[x], "@", gsub("\\.", "_", round(pk_pos, 2)), "_", pk_idx,
                          "___target_", gsub("\\[|\\]| |\\:|\\.", "_", ifelse(dd == "direct", overlap_pairs$id[x], overlap_pairs$value[x])), 
                          ".RDS")

      if(file.exists(out_file) & !overwrite){
        
        cat("\t-->", dd, "-", qq, "qtl lod drop mediation analysis already performed. Skipping...\n")
        
      } else{
        
        cat("\t--> Running", dd, "-", qq, "qtl lod drop mediation analysis...\n")
        
        if(qtl_scan_only_at_pk){
          
          # adapted from qtl2::pull_genoprobint source code
          pr.use      <- pr[[qq]][, overlap_pairs$chr[x]]
          pr.use[[1]] <- pr.use[[1]][, , pk_idx, drop = F]
          kin.use     <- kin[[qq]][names(kin[[qq]]) == overlap_pairs$chr[x]]
          
        } else{
          
          pr.use  <- pr[[qq]]
          kin.use <- kin[[qq]]
          
        }
        diet.qq <- gsub("__.*", "", qq)
        
        # Unconditioned (we already have oscan1_unc from the saved QTL mapping results, but recompute on samples with lipidomic data (NA-excluded) for fairness)
        oscan1_unc <- qtl2::scan1(genoprobs = pr.use, pheno = target_pheno[[diet.qq]], kinship = kin.use, addcovar = add_covar[[diet.qq]], intcovar = int_covar[[qq]], cores = 1)
        lod_pk_unc <- oscan1_unc[pk_idx, ]
        
        # Conditional on gene expression
        lod_pk_con <- pbmcapply::pbmclapply(mc.cores = nbCores_mediatiors_scan, X = colnames(mediator_pheno[[diet.qq]]), function(zz){
          # zz <- colnames(mediator_pheno[[diet.qq]])[2]
          oscan1_con <- qtl2::scan1(genoprobs = pr.use, pheno = target_pheno[[diet.qq]], kinship = kin.use, addcovar = add_covar_con[[diet.qq]][, c(zz, add_covar_cols[[diet.qq]]), drop = F], intcovar = int_covar[[qq]], cores = 1)
          lod_pk_con <- oscan1_con[pk_idx, ]
          data.frame(mediator_id = zz, lod_cond = lod_pk_con, stringsAsFactors = F)
        })
        
        
        lod_drop.df                     <- as.data.frame(rbindlist(lod_pk_con))
        # colnames(lod_drop.df)[1]    <- mediator_id
        lod_drop.df$target_id           <- ifelse(dd == "direct", colnames(target_pheno[[diet.qq]]), overlap_pairs$gene_id_full[x])
        lod_drop.df$lod_unc             <- lod_pk_unc
        lod_drop.df$peak_chr            <- overlap_pairs$chr[x]
        lod_drop.df$peak_pos_id         <- "lQTL_locus"
        lod_drop.df$peak_pos            <- pk_pos
        lod_drop.df$peak_marker         <- pk_idx
        lod_drop.df$lod_drop            <- lod_drop.df$lod_cond - lod_drop.df$lod_unc
        lod_drop.df$gene_feature_type   <- overlap_pairs$feature_type[x]
        lod_drop.df$analysis_id         <- qq
        lod_drop.df$has_cf_covar        <- T
        lod_drop.df$mediation_direction <- dd
        if(dd == "direct" & overlap_pairs$feature_type[x] == "mRNA"){
          lod_drop.df <- merge(lod_drop.df, geneTable, by.x = "mediator_id", by.y = "gene_id_ensembl", all.x = T, all.y = F)
        } else if(dd == "direct" & overlap_pairs$feature_type[x] == "protein"){
          lod_drop.df <- merge(lod_drop.df, protTable_Coon, by.x = "mediator_id", by.y = "protein_id_uniprot", all.x = T, all.y = F)
          lod_drop.df <- type.convert(tidyr::separate_rows(lod_drop.df, c("gene_id_ensembl", "gene_symbol", "gene_chr", "posStart_Mb", "posEnd_Mb"), sep = "\\;"), as.is = T)
        }
        if(dd == "reverse" & overlap_pairs$feature_type[x] == "mRNA"){
          lod_drop.df$gene_id_ensembl <- gsub("__.*", "", lod_drop.df$target_id)
          lod_drop.df                 <- merge(lod_drop.df, geneTable, by = "gene_id_ensembl", all.x = T, all.y = F)
        } else if(dd == "reverse" & overlap_pairs$feature_type[x] == "protein"){
          lod_drop.df$protein_id_uniprot <- gsub(".*__", "", lod_drop.df$target_id)
          lod_drop.df                    <- merge(lod_drop.df, protTable_Coon, by = "protein_id_uniprot", all.x = T, all.y = F)
          lod_drop.df                    <- type.convert(tidyr::separate_rows(lod_drop.df, c("gene_id_ensembl", "gene_symbol", "gene_chr", "posStart_Mb", "posEnd_Mb"), sep = "\\;"), as.is = T)
        }
        lod_drop.df <- as.data.frame(lod_drop.df)
        
        saveRDS(lod_drop.df, out_file)
        
      }
      
      
    })
  })
})










