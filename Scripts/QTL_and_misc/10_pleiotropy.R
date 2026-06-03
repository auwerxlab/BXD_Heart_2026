
library(qtl2pleio)
library(data.table)

source("./Scripts/QTL_and_misc/QTL_fun.R")



################################################################################
## select lipids to run
################################################################################

lipids_run <- c("HexCer[NS] d18:0_22:1_RT_12.391", "HexCer[NS] d41:1_RT_12.925", "HexCer[NS] d42:2_RT_12.286")

cat("\n Pleiotropy analysis run for all the lQTL/(e/p)QTL pairs identified for the following lipdis:\n",
    paste(paste0("\t--> ", lipids_run), collapse = "\n"), "\n",
      "!!! Modify lipid selection is pleiotropy scripts to run for other lipids.\n\n")

################################################################################
## qtl2 pleiotropy analysis
################################################################################


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


lipid_QTL_result_files       <- list.files("./Data/QTL_mapping/qtl_outputs", recursive = T, full.names = T, pattern = "average_lipidomic_derived|individuals_lipidomic_derived|individuals_lipidomic_Filt")
lipid_QTL_result_files       <- lipid_QTL_result_files[grepl("qtl2_results__", lipid_QTL_result_files)]
lipid_QTL_result_files       <- lipid_QTL_result_files[!grepl("00__previous_runs_no_diet_stratification|previous_runs_no_diet_stratification", lipid_QTL_result_files)]
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




saveDir  <- "./Data/lipidomics/qtl_colocalization_pleiotropy_test"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

# Otherwise qtl2pleio uses all available threads
# Note: this changes how the random state is internally handled. Results generated
#       when limiting the number of threads or with default behaviour are different
#       but consistent in terms of outcomes, even when specifying a seed (that makes
#       multiple runs reproducible - tested)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

overwrite               <- F
pk_region_half_size_Mb  <- 20
nbCores_bootstrap       <- 68
nb_bootstrap_iterations <- 1000
tt                      <- lapply(1:nrow(overlap_pairs), function(x){
  # x <- 1
  
  cat("\nRunning pleiotropy analysis for entry", x, "/", nrow(overlap_pairs), "...\n")
  
  tt <- lapply(c("whole_chr", "peak_region"), function(rr){
    # rr <- "whole_chr"
    # rr <- "peak_region"
    
    cat("\t--> Running qtl2pleio pleiotropy test for idx ", x, "/", nrow(overlap_pairs), " - ", rr, "...\n")
    
    lQTL.obj <- lipid_QTL_result_file_traits[[overlap_pairs$analysis_id[x]]]
    tmpLogi  <- unlist(lapply(lQTL.obj, function(y) overlap_pairs$id[x] %in% y))
    stopifnot(sum(tmpLogi) == 1)
    lQTL.obj <- readRDS(names(lQTL.obj)[tmpLogi])
    
    
    lipid_pheno <- lipid_pheno_list[[gsub("__.*", "", overlap_pairs$analysis_id[x])]]
    lipid_pheno <- lipid_pheno[, overlap_pairs$id[x], drop = F]
    # keep only strain/diet conditions without missing lipid data
    lipid_pheno <- lipid_pheno[!is.na(lipid_pheno[, 1]), , drop = F]
    
    
    if(overlap_pairs$feature_type[x] == "mRNA"){
      
      # get original gene count data, mean-aggregate by strain/diet condition, filter only strain/diet conditions without missing lipid data
      gene_pheno           <- gene_pheno_list[[gsub("__.*", "", overlap_pairs$analysis_id[x])]]
      gene_pheno           <- gene_pheno[match(rownames(lipid_pheno), rownames(gene_pheno)), ]
      gene_pheno           <- gene_pheno[, overlap_pairs$feature_id[x], drop = F]
      gene_pheno           <- gene_pheno[!is.na(gene_pheno[, 1]), , drop = F]
      gene_norm            <- transform2normal_byCol_INT(gene_pheno)
      gene_norm            <- gene_norm$transf.df
      rownames(gene_norm)  <- rownames(gene_pheno)
      gene_pheno           <- gene_norm
      
    } else{
      
      gene_pheno <- protein_pheno_list[[gsub("__.*", "", overlap_pairs$analysis_id[x])]]
      gene_pheno <- gene_pheno[match(rownames(lipid_pheno), rownames(gene_pheno)), ]
      gene_pheno <- gene_pheno[, overlap_pairs$feature_id[x], drop = F]
      gene_pheno <- gene_pheno[!is.na(gene_pheno[, 1]), , drop = F]
      
    }
    lipid_pheno <- lipid_pheno[rownames(lipid_pheno) %in% rownames(gene_pheno), , drop = F]
    
    
    add_covar <- protein_addcovar_list[[gsub("__.*", "", overlap_pairs$analysis_id[x])]]
    if(grepl("CD_HFD", overlap_pairs$analysis_id[x])){
      diet_vec  <- ifelse(gsub(".*_", "", rownames(add_covar)) == "CD", 0, 1)
      stopifnot(length(unique(diet_vec)) == 2)
      add_covar <- cbind(diet = diet_vec, add_covar)
    }
    add_covar   <- add_covar[match(rownames(lipid_pheno), rownames(add_covar)), , drop = F]
    int_covar   <- lQTL.obj$intcovar
    
    pr          <- lQTL.obj$pr[[which(names(lQTL.obj$pr) == overlap_pairs$chr[x])]]
    kin         <- lQTL.obj$kinship[[which(names(lQTL.obj$kinship) == overlap_pairs$chr[x])]]
    all_mat     <- cbind(lipid_pheno, gene_pheno)
    
    kin         <- kin[rownames(all_mat), rownames(all_mat)]
    pr          <- pr[rownames(all_mat), , ]
    if(!is.null(add_covar)){
      add_covar   <- add_covar[rownames(all_mat), , drop = F]
    }
    if(!is.null(int_covar)){
      int_covar   <- int_covar[match(rownames(all_mat), names(int_covar))]
    }
    
    # get lipid and gene QTL peaks
    lQTl.pk.pos    <- overlap_pairs$pos_Mb[x]
    lQTl.pk.idx    <- lQTL.obj$pmap[[which(names(lQTL.obj$pmap) == overlap_pairs$chr[x])]]
    lQTl.pk.idx    <- names(lQTl.pk.idx)[which.min(abs(lQTl.pk.idx - lQTl.pk.pos))]
    geneQTL.pk.pos <- all_genes_criteria.list[[gsub("__average.*|__ind.*", "", overlap_pairs$analysis_id[x])]]
    geneQTL.pk.pos <- geneQTL.pk.pos[geneQTL.pk.pos$ref_qtl_id == overlap_pairs$id[x] & geneQTL.pk.pos$gene_id == gsub("__.*", "", overlap_pairs$gene_id_full[x]), ]$eGenes_peaks_table[[1]]
    geneQTL.pk.pos <- geneQTL.pk.pos[grepl(ifelse(overlap_pairs$feature_type[x] == "mRNA", "eQTL", "pQTL"), geneQTL.pk.pos$qtl_type), ]
    geneQTL.pk.pos <- geneQTL.pk.pos[geneQTL.pk.pos$id == overlap_pairs$feature_id[x] & geneQTL.pk.pos$chr == overlap_pairs$chr[x], ]
    if(nrow(geneQTL.pk.pos) == 0){
      cat("\t\t--> qtl for molecular feature not found...\n")
      return()
    }
    geneQTL.pk.pos$dist <- abs(geneQTL.pk.pos$pos_Mb - overlap_pairs$pos_Mb[x])
    geneQTL.pk.pos      <- geneQTL.pk.pos$pos_Mb[which.min(geneQTL.pk.pos$dist)]
    geneQTL.pk.idx      <- lQTL.obj$pmap[[which(names(lQTL.obj$pmap) == overlap_pairs$chr[x])]]
    geneQTL.pk.idx      <- names(geneQTL.pk.idx)[which.min(abs(geneQTL.pk.idx - geneQTL.pk.pos))]
    
    
    if(rr == "peak_region"){
      
      # only scan in the (+/-)pk_region_half_size_Mb region around the peaks, to speed up computations
      scan_range   <- sort(c(lQTl.pk.pos, geneQTL.pk.pos))
      
      scan_range   <- c(lQTl.pk.pos[1] - pk_region_half_size_Mb, lQTl.pk.pos[1] + pk_region_half_size_Mb)
      scan_markers <- lQTL.obj$pmap[[which(names(lQTL.obj$pmap) == overlap_pairs$chr[x])]]
      scan_markers <- scan_markers[scan_markers >= scan_range[1] & scan_markers <= scan_range[2]]
      sort(c(lQTl.pk.pos, geneQTL.pk.pos))
      
      pr <- pr[, , names(scan_markers), drop = F]
      
    } else if(rr == "whole_chr"){
      
      scan_range <- base::range(unname(lQTL.obj$pmap[[which(names(lQTL.obj$pmap) == overlap_pairs$chr[x])]]))
      
    } else{
      stop()
    }
    
    
    stopifnot(identical(rownames(lipid_pheno), rownames(gene_pheno)))
    stopifnot(identical(rownames(all_mat), rownames(kin)))
    stopifnot(identical(rownames(all_mat), colnames(kin)))
    stopifnot(identical(rownames(all_mat), rownames(pr)))
    if(!is.null(add_covar)){
      stopifnot(identical(rownames(all_mat), rownames(add_covar)))
    }
    if(!is.null(int_covar)){
      stopifnot(identical(rownames(all_mat), names(int_covar)))
    }
    
    
    
    
    # Note: to model the GxD interaction, using a matrix of genotyping prob. (BB, DD, BB * diet, DD * diet) will not work as qtl2pleio assumes a row-wise sum of genetic
    #       probabilities that sum up to 1. Adding (BB *diet, BB * diet) to the additive covariate matrix doesn't work as this matrix is multiplied to the genetic 
    #       probabilities matrix, meaning that G * C would result in terms where we have twice the same terms multiplied: BB * (BB * diet), DD * (DD * diet).
    #       For now we replace the GxD just by a D additive effect, and we save the names accordingly.
    
    
    out_file <- paste0(saveDir, "/qtl2pleio_results", 
                       "___", rr, ifelse(rr == "peak_region", paste0("_half_size_", pk_region_half_size_Mb, "_Mb"), ""),
                       "___cf_covar",
                       "___", gsub("__average.*|__ind.*", "", overlap_pairs$analysis_id[x]), 
                       "___", gsub("\\[|\\]| |\\:|\\.", "_", overlap_pairs$id[x]), 
                       "___", overlap_pairs$gene_id_full[x],
                       "___", overlap_pairs$chr[x], "@", gsub("\\.", "_", round(lQTl.pk.pos, 2)), "_", lQTl.pk.idx,
                       "___", overlap_pairs$feature_type[x], ".RDS")
    # Note: as the interaction effect in CD_HFD__dietInteraction is replaced with a diet additive effect, change saving id accordingly
    out_file <- gsub("CD_HFD__dietInteraction", "CD_HFD__dietAdditive", out_file)
    if(file.exists(out_file) & !overwrite){
      cat("\t--> Results already saved. Skipping...\n")
      return()
    }
    
    
    
    # Perform model fitting for all ordered pairs of markers in a genomic region of interest
    set.seed(123)
    pleio_scan <- qtl2pleio::scan_pvl(probs     = pr,
                                      pheno     = all_mat,
                                      addcovar  = add_covar,
                                      kinship   = kin,
                                      start_snp = 1,
                                      n_snp     = dim(pr)[3],
                                      cores     = 1)
    
    
    
    # calculate the likelihood ratio test statistic value for the specified traits and specified genomic region
    lrt <- qtl2pleio::calc_lrt_tib(pleio_scan)

    # find the marker index corresponding to the peak of the pleiotropy trace
    # identify the index (on the chromosome under study) of the marker that maximizes the likelihood under the pleiotropy constraint
    pleio_index  <- qtl2pleio::find_pleio_peak_tib(pleio_scan, start_snp = 1)
    pleio_pos    <- lQTL.obj$pmap[[which(names(lQTL.obj$pmap) == overlap_pairs$chr[x])]]
    pleio_pos    <- pleio_pos[names(pleio_pos) %in% dimnames(pr)[[3]]]
    pleio_marker <- names(pleio_pos)[pleio_index]
    pleio_pos    <- unname(pleio_pos)[pleio_index]
    
    pleio_lod     <- qtl2pleio::calc_profile_lods(pleio_scan)
    pleio_lod$trait[pleio_lod$trait != "pleiotropy"] <- plyr::mapvalues(pleio_lod$trait[pleio_lod$trait != "pleiotropy"], from = c("tr1", "tr2"), to = colnames(all_mat))
    pleio_lod     <- qtl2pleio::add_pmap(tib = pleio_lod, pmap = lQTL.obj$pmap[[which(names(lQTL.obj$pmap) == overlap_pairs$chr[x])]])
    pleio_lod     <- pleio_lod[order(pleio_lod$trait, pleio_lod$marker_position), ]
    pleio_lod$chr <- overlap_pairs$chr[x]
    
  
    # bootstrapping for pvalue
    set.seed(123)
    pleio_boot <- qtl2pleio::boot_pvl(probs            = pr,
                                      pheno            = all_mat,
                                      addcovar         = add_covar,
                                      kinship          = kin,
                                      pleio_peak_index = pleio_index,
                                      nboot            = nb_bootstrap_iterations,
                                      start_snp        = 1,
                                      n_snp            = dim(pr)[3],
                                      cores            = nbCores_bootstrap)
    
    
    pvalue <- mean(pleio_boot >= lrt)
    
    results_list <- list(probs               = pr,
                         pheno               = all_mat,
                         kinship             = kin,
                         addcovar            = add_covar,
                         nb_perm             = nb_bootstrap_iterations,
                         pleio_scan          = pleio_scan,
                         lrt                 = lrt,
                         pleio_index         = pleio_index,
                         pleio_lod           = pleio_lod,
                         pleio_marker        = pleio_marker,
                         pleio_pos           = pleio_pos,
                         pleio_boot          = pleio_boot,
                         pvalue              = pvalue,
                         scan_region         = rr,
                         scan_width_half_Mb  = ifelse(rr == "peak_region", pk_region_half_size_Mb, NA),
                         scan_range          = scan_range,
                         peak_chr            = overlap_pairs$chr[x],
                         lipid_peak_pos      = overlap_pairs$pos_Mb[x],
                         feature_peak_pos    = geneQTL.pk.pos,
                         analysis_id         = overlap_pairs$analysis_id[x],
                         lipid_id            = overlap_pairs$id[x],
                         gene_id             = overlap_pairs$gene_id_full[x],
                         feature_id          = overlap_pairs$feature_id[x],
                         feature_type        = overlap_pairs$feature_type[x],
                         has_cf_covar        = T)
    
    saveRDS(results_list, out_file)
  })
})