library(qtl2)
library(dplyr)
library(parallel)


################################################################################
## set analysis parameters
################################################################################

# Directory containing the R/qtl2 input files prepared as described here:
# https://kbroman.org/qtl2/assets/vignettes/input_files.html
input_dir       <- "./Data/QTL_mapping/qtl_processed_inputs/average_lipidomic/CD_HFD"
# data type
input_data_type <- "average_lipidomic"
# condition
input_condition <- "CD_HFD"
# nb. permutations
nbPerm          <- 5000
# significance threshold
signif_thresh   <- 0.05
# CI probability
ci_prob         <- 0.95
# output saving directory
dirQTL_outputs  <- "./Data/qtl_outputs"
# Nb. cores
nbCores         <- parallel::detectCores() - 2
# overwrite existing results
overwrite       <- F


################################################################################
## run analysis
################################################################################

set.seed(22)
qtlResults <- lapply(input_dir, function(dirQTL){
  # dirQTL <- input_dir[1]

  input_data_type      <- tail(unlist(strsplit(dirQTL, "/")), 2)[1]
  input_condition <- tail(unlist(strsplit(dirQTL, "/")), 1)
  
  dirSave <- paste0(dirQTL_outputs, "/", input_data_type)
  if(!dir.exists(dirSave)){
    dir.create(dirSave, recursive = T)
  }
  
  res <- lapply(c("noInteraction", "dietInteraction"), function(y){
    # y <- "noInteraction"
    
    cat(paste0("Running QTL analysis for...",
                "\n\t\t\t\t--> input type:\t", input_data_type,
                "\n\t\t\t\t--> condition:\t", input_condition,
                "\n\t\t\t\t--> interaction:\t", y, "\n"))
    
    input_file  <- paste0(dirQTL, "/BXD.json")
    output_file <- paste0(dirSave, "/qtl2_results__", input_data_type, "__", y, "__", input_condition, "__.RDS")
    if(file.exists(output_file) & !overwrite){
      cat("\t\t\t\t--> Analysis already done!\n")
      return(NA)
    }
    
    cat("\t\t\t\t--> Loading QTL objects...\n")
    a              <- Sys.time()
    qtl.cross      <- qtl2::read_cross2(input_file)

    if(y == "noInteraction"){
      cat("\t\t\t\t--> No interactive covariate used...\n")
      intercovar        <- NULL
      if(input_condition == "CD_HFD"){
        cat("\t\t\t\t--> Performing diet-stratified permutations...\n")
        strat_vec         <- as.numeric(factor(qtl.cross$covar$diet, levels = c("CD", "HFD"))) - 1
        names(strat_vec)  <- rownames(qtl.cross$covar)
      } else{
        cat("\t\t\t\t--> Performing non-stratified permutations...\n")
        strat_vec         <- NULL
      }
    } else if(y == "dietInteraction" & input_condition == "CD_HFD"){
      cat("\t\t\t\t--> Using diet interactive covariate...\n")
      intercovar        <- as.numeric(factor(qtl.cross$covar$diet, levels = c("CD", "HFD"))) - 1
      names(intercovar) <- rownames(qtl.cross$covar)
      cat("\t\t\t\t--> Performing diet-stratified permutations...\n")
      strat_vec         <- as.numeric(factor(qtl.cross$covar$diet, levels = c("CD", "HFD"))) - 1
      names(strat_vec)  <- rownames(qtl.cross$covar)
    } else{
      cat("\t\t\t\t--> Interaction not possible (only one diet); iteration skipped!\n")
      return(NA)
    }
    tmpLogi <- grepl("add_covar__", colnames(qtl.cross$covar))
    if(any(tmpLogi)){
      cat("\t\t\t\t--> Using additive covariates...\n")
      additivecovar <- as.matrix(qtl.cross$covar[, tmpLogi])
      cat("\t\t\t\t--> The following additive covariates will be used:\n", 
          paste(paste0("\t\t\t\t    (", 1:ncol(additivecovar), ") ", colnames(additivecovar)), collapse = "\n"), "\n")
    } else{
      cat("\t\t\t\t--> No additive covariate used...\n")
      additivecovar <- NULL
    }
    cat("\t\t\t\t--> Using", nbPerm, "permutations...\n")
    cat("\t\t\t\t--> Analyzing", ncol(qtl.cross$pheno), "phenotypes...\n")
    
    cat("\t\t\t\t--> Removing X and M chromosomes...\n")
    chrNames       <- qtl2::chr_names(qtl.cross)
    qtl.cross      <- subset(qtl.cross, chr = chrNames[!(chrNames %in% c("X", "M"))])
    perm_Xsp       <- F
    
    cat("\t\t\t\t--> Inserting pseudomarkers and interpolating new physical map...\n")
    gmap           <- qtl2::insert_pseudomarkers(map = qtl.cross$gmap, step = 0.2, stepwidth = "max")
    pmap           <- qtl2::interp_map(gmap, qtl.cross$gmap, qtl.cross$pmap)  # linear interpolation to get estimated physical positions for the inserted pseudomarkers.
    markersX       <- names(gmap$X)

    cat("\t\t\t\t--> Computing conditional genotype probilities and kinship matrix...\n")
    pr             <- qtl2::calc_genoprob(cross = qtl.cross, map = gmap, map_function = "c-f", cores = nbCores) 
    Xcovar         <- qtl2::get_x_covar(qtl.cross)
    kinship        <- qtl2::calc_kinship(pr, type = "loco", cores = nbCores) 
    kinship.o      <- qtl2::calc_kinship(pr, type = "overall", cores = nbCores)
    
    cat("\t\t\t\t--> Running scan1...\n")
    oscan1         <- qtl2::scan1(genoprobs = pr, pheno = qtl.cross$pheno, kinship = kinship, Xcovar = Xcovar, addcovar = additivecovar, intcovar = intercovar, cores = nbCores)

    cat("\t\t\t\t--> Running scan1perm...\n")
    if(perm_Xsp){
      chrLength <- qtl2::chr_lengths(gmap)
    } else{
      chrLength <- NULL
    }
    operm          <- qtl2::scan1perm(genoprobs = pr, pheno = qtl.cross$pheno, kinship = kinship, Xcovar = Xcovar, intcovar = intercovar, addcovar = additivecovar, n_perm = nbPerm, perm_Xsp = perm_Xsp, chr_lengths = chrLength, cores = nbCores, perm_strata = strat_vec)
    thr            <- summary(operm, alpha = signif_thresh)

    cat("\t\t\t\t--> Finding QTL Peaks - main peak per chr...\n")
    peaks.g        <- qtl2::find_peaks(scan1_output = oscan1[!(rownames(oscan1) %in% markersX), , drop = F], map = gmap, threshold = thr, prob = ci_prob, expand2markers = FALSE, cores = nbCores)
    peaks.p        <- qtl2::find_peaks(scan1_output = oscan1[!(rownames(oscan1) %in% markersX), , drop = F], map = pmap, threshold = thr, prob = ci_prob, expand2markers = FALSE, cores = nbCores)
    peaks.X.g      <- NULL
    peaks.X.p      <- NULL
    cat("\t\t\t\t--> Finding QTL Peaks - secondary peak per chr...\n")
    peaks.g.sec    <- qtl2::find_peaks(scan1_output = oscan1[!(rownames(oscan1) %in% markersX), , drop = F], map = gmap, threshold = thr, prob = ci_prob, expand2markers = FALSE, peakdrop = 0.05, cores = nbCores)
    peaks.p.sec    <- qtl2::find_peaks(scan1_output = oscan1[!(rownames(oscan1) %in% markersX), , drop = F], map = pmap, threshold = thr, prob = ci_prob, expand2markers = FALSE, peakdrop = 0.05, cores = nbCores)
    peaks.X.g.sec  <- NULL
    peaks.X.p.sec  <- NULL
    b              <- Sys.time()

    cat("\t\t\t\t--> Saving results...\n")
    results <- list(data.type      = input_data_type,
                    condition      = input_condition,
                    interaction    = y,
                    nbPerm         = nbPerm,
                    intcovar       = intercovar,
                    addcovar       = additivecovar,
                    perm_strata    = strat_vec,
                    droppedMarkers = markers2remove,
                    cross          = qtl.cross,
                    gmap           = gmap,
                    pmap           = pmap,
                    pr             = pr,
                    kinship        = kinship,
                    kinship.o      = kinship.o,
                    oscan1         = oscan1,
                    operm          = operm,
                    threshold      = thr,
                    peaks.g        = peaks.g,
                    peaks.p        = peaks.p,
                    peaks.X.g      = peaks.X.g,
                    peaks.X.p      = peaks.X.p,
                    peaks.g.sec    = peaks.g.sec,
                    peaks.p.sec    = peaks.p.sec,
                    peaks.X.g.sec  = peaks.X.g.sec,
                    peaks.X.p.sec  = peaks.X.p.sec,
                    elapsedTime    = as.numeric(difftime(b, a, units = "secs")))
    
    saveRDS(kinship, paste0(dirSave, "/kinship_loco__", input_data_type, "__", y, "__", input_condition, ".RDS"))
    saveRDS(kinship.o, paste0(dirSave, "/kinship_overall__", input_data_type, "__", y, "__", input_condition, ".RDS"))
    saveRDS(results, paste0(dirSave, "/qtl2_results__", input_data_type, "__", y, "__", input_condition, ".RDS"))

    remove(a, b, qtl.cross, intercovar, nbPerm, chrNames, perm_Xsp,
           markers2remove, allMarkers, markers2keep, gmap, pmap, markersX,
           pr, Xcovar, kinship, kinship.o, oscan1, chrLength, operm, thr,
           peaks.g, peaks.p, peaks.g.sec, peaks.p.sec, results)
    
    NA
  })
  NA
})

