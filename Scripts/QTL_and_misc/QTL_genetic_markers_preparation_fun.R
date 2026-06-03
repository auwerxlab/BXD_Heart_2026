
getMarkers2remove <- function(geno.chr, gmap.chr, min_distance = 0.01, H_code = 0, cores = 1){

  # remove marker closer than 0.01 cM if they have the same ghenotype across all strains
  # Notes: (1) qtl2::reduce_markers remove markers closer than X cM but doesn't check if the genotype at those markers is identical for all strains
  #        (2) if 2 markers are adjacent and one has a missing genotype H = 0 (not D or B) for some strains and if for all other strains the 
  #            genotype is identical (and if at for the strain with H all other markers have same non-H allele, otherwise keep all possible alleles 0/1/2 that are present), 
  #            the marker with no H or with less H is retained
  
  markers2remove.list <- mclapply(mc.cores = cores, X = names(geno.chr), FUN = function(ee){
    # ee <- names(geno.chr)[1]
    print(ee)
    tmp.geno      <- geno.chr[[ee]]
    tmp.geno      <- as.data.frame(t(tmp.geno))
    tmp.gmap      <- data.frame(marker = names(gmap.chr[[ee]]), pos = as.numeric(unname(gmap.chr[[ee]])), stringsAsFactors = F)
    tmp           <- merge(tmp.gmap, tmp.geno, by.x = "marker", by.y = "row.names")
    tmp           <- tmp[order(tmp$pos), ]
    tmp           <- cbind(dist = c(1, diff(tmp$pos)), tmp) 
    tmp           <- cbind(logiDist = tmp$dist <= min_distance, tmp)
    rownames(tmp) <- NULL
    tmpIdx        <- which(tmp$logiDist)
    tt <- tmp[tmpIdx, ]
    tmpIdx        <- split(tmpIdx, cumsum(c(1, diff(tmpIdx) != 1)))
    if(length(unlist(tmpIdx)) == 0){
      return(NULL)
    }
    markersRm     <- unlist(lapply(tmpIdx, function(uu){
      # uu <- tmpIdx[[55]]
      geno.uu             <- dplyr::select(tmp[c(min(uu) - 1, uu), ], -c("logiDist", "dist"))
      geno.uu$hasZeroGeno <- apply(geno.uu[, -c(1, 2)], 1, function(qq) any(qq == H_code))
      geno.uu$nbZeroGeno  <- apply(geno.uu[, -c(1, 2, ncol(geno.uu))], 1, function(qq) sum(qq == H_code))
      
      geno.uu <- geno.uu[order(!geno.uu$hasZeroGeno, geno.uu$nbZeroGeno), ]
      geno.uu <- dplyr::select(geno.uu, -c("hasZeroGeno", "pos", "nbZeroGeno"))
      tmpLogi <- which(c(F, (apply(geno.uu[, -1], 2, function(qq) any(qq == H_code) & !all(qq == H_code)))))
      for(i in tmpLogi){
        # i <- tmpLogi[1]
        nonZero <- unique(geno.uu[[i]][geno.uu[[i]] != H_code])
        if(length(nonZero) == 1){
          geno.uu[[i]] <- nonZero
        }
      }
      mark.uu <- geno.uu$marker
      geno.uu <- geno.uu[!duplicated(geno.uu[, -1]), ]
      mark.uu[!(mark.uu %in% geno.uu$marker)]
    }))
    unname(markersRm)
  })
  names(markers2remove.list) <- names(geno.chr)
  markers2remove.list
}




