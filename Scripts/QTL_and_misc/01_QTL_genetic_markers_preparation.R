
library(reshape2)
library(data.table)




source("./Scripts/QTL_and_misc/QTL_genetic_markers_preparation_fun.R")


################################################################################
## QTL preparation
################################################################################

################################################################################
## (1a) prepare BXD markers, genetic and physical maps from files downloaded from 
##     geneNetwork (GN).
##     Adapted from: https://github.com/rqtl/qtl2data/tree/master/BXD
##     https://gn1.genenetwork.org/webqtl/main.py?FormID=sharinginfo&GN_AccessionId=600
##     BXD genotype files downloaded from:
##     https://gn1.genenetwork.org/webqtl/main.py?FormID=sharinginfo&GN_AccessionId=600
##
##     Also fix switched markers (by position)
################################################################################


# Download genetic markers
saveDir <- "./Data/input_data/BXD_markers"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
if(!file.exists(paste0(saveDir, "/BXD_Geno_2017.xlsx"))){
  system(paste0("wget -P ", saveDir, " http://files.genenetwork.org/current/GN600/BXD_Geno_2017.xlsx"))
}


dirQTL_input <- "./Data/QTL_mapping/qtl_processed_inputs"
if(!dir.exists(dirQTL_input)){
  dir.create(dirQTL_input, recursive = T)
}

geno                   <- openxlsx::read.xlsx("./Data/input_data/BXD_markers/BXD_geno_2017.xlsx", colNames = F, rowNames = F)
rowStart               <- which(unname(apply(geno, 1, function(x) any(grepl("^chr", x, ignore.case = T)) & any(grepl("^locus", x, ignore.case = T)))))
geno                   <- geno[rowStart:nrow(geno), ]
colnames(geno)         <- geno[1, ]
geno                   <- geno[-1, ]
rownames(geno)         <- NULL
geno                   <- type.convert(geno, as.is = T)
col_cM                 <- "cM_BXD"
col_Mb                 <- "Mb_mm10"
drop_switched_markers  <- F # switch them

cols.header    <- c("Chr", "Locus", "Mb_mm9", "Mb_mm10", "Mb", "cM_BXD", "cM")
unique.alleles <- unique(do.call(c, dplyr::select(geno, -cols.header[cols.header %in% colnames(geno)])))
cat("The following alleles have been detected:", paste(sort(unique.alleles), collapse = ", "), "\n")

geno           <- geno[order(factor(geno$Chr, c(1:19,"X", "Y", "M")), as.numeric(geno[[col_Mb]])),]
rownames(geno) <- geno$Locus
geno.raw       <- geno

# check that physical and genetic maps match
allChromosomes <- unique(geno$Chr)
isOrderCorrect <- unlist(lapply(allChromosomes, function(x){
  tmp.g   <- geno[[col_cM]][geno$Chr == x]
  tmp.p   <- geno[[col_Mb]][geno$Chr == x]
  all(order(tmp.g) == order(tmp.p))
}))
names(isOrderCorrect) <- allChromosomes
# print(knitr::kable(isOrderCorrect))

if(!drop_switched_markers & any(!isOrderCorrect)){
  cat("Switching inverted markers...\n")
  
  # fix differences in order; take physical map as correct and swap genetic positions
  locusSwaps     <- lapply(allChromosomes, function(x){
    tmp.chr           <- geno[geno$Chr == x, c("Locus", col_cM, col_Mb)]
    tmp.chr           <- tmp.chr[order(tmp.chr[[col_Mb]]), ]
    rownames(tmp.chr) <- NULL
    tmpLogi           <- tmp.chr[[col_cM]][-1] < tmp.chr[[col_cM]][-nrow(tmp.chr)]
    if(any(tmpLogi)){
      orderCM  <- order(tmp.chr[[col_cM]])
      # tt <- tmp.chr[orderCM[orderCM != 1:length(orderCM)], ]
      switches <- match(orderCM, 1:nrow(tmp.chr))
      from     <- (1:nrow(tmp.chr))[which(switches != 1:nrow(tmp.chr))]
      to       <- orderCM[which(orderCM != 1:nrow(tmp.chr))]
      out      <- data.frame(chr         = x,
                             fromMarker  = tmp.chr$Locus[from],
                             toMarker    = tmp.chr$Locus[to],
                             fromPos_CM  = tmp.chr[[col_cM]][from],
                             toPos_CM    = tmp.chr[[col_cM]][to],
                             stringsAsFactors = F)
      out[!duplicated(unlist(lapply(1:nrow(out), function(y) paste(sort(as.character(out[y, 2:3])), collapse = ";")))), ]
    } else{
      NULL
    }
  })
  locusSwaps <- rbindlist(locusSwaps)
  print(knitr::kable(locusSwaps))
  cat("\t--> ", nrow(locusSwaps), "markers switched.\n")
  
  for(i in 1:nrow(locusSwaps)){
    geno[locusSwaps$fromMarker[i], ][[col_cM]] <- locusSwaps$toPos_CM[i]
    geno[locusSwaps$toMarker[i], ][[col_cM]]   <- locusSwaps$fromPos_CM[i]
  }
  
  
} else if(drop_switched_markers & any(!isOrderCorrect)){
  cat("Dropping inverted markers...\n")
  
  markers.rm  <- lapply(allChromosomes, function(x){
    tmp.chr           <- geno[geno$Chr == x, c("Locus", col_cM, col_Mb)]
    tmp.chr           <- tmp.chr[order(tmp.chr[[col_Mb]]), ]
    rownames(tmp.chr) <- NULL
    tmpLogi           <- tmp.chr[[col_cM]][-1] < tmp.chr[[col_cM]][-nrow(tmp.chr)]
    if(any(tmpLogi)){
      tmp.chr$Locus[which(tmpLogi)]
    } else{
      NULL
    }
  })
  markers.rm <- unlist(markers.rm)
  cat("\t--> ", length(markers.rm), "markers removed,", nrow(geno) - length(markers.rm), "left\n")
  cat("\tRemoved markers:\t", paste(markers.rm, collapse = "\n\t\t\t\t"), "\n")
  geno       <- geno[!(geno$Locus %in% markers.rm), ]
  
}

# check that physical and genetic maps match
isOrderCorrect <- unlist(lapply(allChromosomes, function(x){
  tmp.g   <- geno[[col_cM]][geno$Chr == x]
  tmp.p   <- geno[[col_Mb]][geno$Chr == x]
  all(order(tmp.g) == order(tmp.p))
}))
names(isOrderCorrect) <- allChromosomes
stopifnot(all(isOrderCorrect))

# genetic map
gmap           <- data.frame(geno[, c("Locus", "Chr", col_cM)])
colnames(gmap) <- c("marker", "chr", "pos")
rownames(gmap) <- gmap[, "marker"]

# physical map
pmap           <- data.frame(geno[, c("Locus", "Chr", col_Mb)])
colnames(pmap) <- c("marker", "chr", "pos")
rownames(pmap) <- pmap[, "marker"]

# check that physical and genetic maps match
isOrderCorrect <- unlist(lapply(allChromosomes, function(x){
  tmp.g   <- gmap[gmap$chr == x, ]
  tmp.p   <- pmap[pmap$chr == x, ]
  all(order(tmp.g$pos) == order(tmp.p$pos))
}))
names(isOrderCorrect) <- allChromosomes
stopifnot(all(isOrderCorrect))

# reorder genetic and physical map by physical order (check if they are already ordered)
isOrdered.pmap <- order(factor(pmap$chr, c(1:19,"X", "Y", "M")), as.numeric(pmap$pos)) == 1:nrow(pmap)
isOrdered.gmap <- order(factor(gmap$chr, c(1:19,"X", "Y", "M")), as.numeric(gmap$pos)) == 1:nrow(gmap)
stopifnot(all(isOrdered.pmap))
stopifnot(all(isOrdered.gmap))


# get genotypes
strains2keep   <- fread("./Data/RNAseq_metadata.tsv", data.table = F, stringsAsFactors = F)
strains2keep   <- stringr::str_sort(unique(strains2keep$strain), numeric = T)
tmpLogi        <- strains2keep %in% colnames(geno)
strains2remove <- unique(colnames(geno)[-(1:5)][!(colnames(geno)[-(1:5)] %in% strains2keep)])
stopifnot(all(tmpLogi))

saveRDS(geno, paste0(dirQTL_input, "/BXD_geno_all_BXD_strains.RDS"))

geno                                      <- geno[, !(colnames(geno) %in% c("Chr", "Mb_mm9", col_Mb, col_cM, strains2remove))]
colnames(geno)[colnames(geno) == "Locus"] <- "marker"
rownames(geno)                            <- geno[, "marker"]

# reorder markers as in the maps
geno <- geno[rownames(pmap), ]


# remove redundant markers
# First transform objects to qtl2cross object-like format we originally implemented function getMarkers2remove to work on data extracted from featurees extracted from qtl2 cross objects
geno.chr        <- geno[, -1]
# transform to numerical as done by qtl2
geno.chr[geno.chr == "B"] <- 1
geno.chr[geno.chr == "D"] <- 2
geno.chr[geno.chr == "H"] <- 0
geno.chr        <- split(type.convert(geno.chr, as.is = T), gmap$chr)
geno.chr        <- lapply(geno.chr, function(x) t(as.matrix(x)))
gmap.chr        <- gmap$pos
names(gmap.chr) <- gmap$marker
gmap.chr        <- split(gmap.chr, gmap$chr)
markers2remove  <- unlist(getMarkers2remove(geno.chr, gmap.chr, min_distance = 0.01, H_code = 0, cores = max(c(1, round(detectCores()) / 2))))

allMarkers     <- unname(unlist(lapply(gmap.chr, names)))
markers2keep   <- allMarkers[!(allMarkers %in% markers2remove)]
cat("\t> ", length(markers2remove), " redundant markers removed (", round((length(markers2remove)/length(allMarkers)) * 100, 1), "%); ", length(markers2keep), " markers retained\n")

gmap <- gmap[!(gmap$marker %in% markers2remove), ]
pmap <- pmap[!(pmap$marker %in% markers2remove), ]
geno <- geno[geno$marker %in% gmap$marker, ]
stopifnot(nrow(gmap) == nrow(geno) & nrow(pmap) == nrow(geno))

# write genetic and physical maps
qtl2convert::write2csv(gmap, paste0(dirQTL_input, "/BXD_gmap.csv"), comment = "Genetic map for BXD data", overwrite = TRUE)
qtl2convert::write2csv(pmap, paste0(dirQTL_input, "/BXD_pmap.csv"), comment = "Physical map (mm10 Mbp) for BXD data", overwrite=TRUE)

saveRDS(gmap, paste0(dirQTL_input, "/BXD_gmap.RDS"))
saveRDS(pmap, paste0(dirQTL_input, "/BXD_pmap.RDS"))

# write the genotypes
qtl2convert::write2csv(geno, paste0(dirQTL_input, "/BXD_geno.csv"), comment = "Genotypes for BXD data", overwrite = TRUE)
saveRDS(geno, paste0(dirQTL_input, "/BXD_geno.RDS"))

# write cross info (BxD for all strains)
crossinfo <- data.frame(id = names(geno)[-1], cross_direction = rep("BxD", ncol(geno) - 1))
qtl2convert::write2csv(crossinfo, paste0(dirQTL_input, "/BXD_crossinfo.csv"), overwrite = TRUE, comment = paste0("Cross info for BXD data\n#", "(all lines formed from cross between female B and male D)"))
saveRDS(crossinfo, paste0(dirQTL_input, "/BXD_crossinfo.RDS"))


na_strings <- c("-", "NA")      # the default (see ?qtl2::write_control_file)

qtl2::write_control_file(paste0(dirQTL_input,  "/BXD.json"),
                         description     = "BXD mouse data from GeneNetwork",
                         crosstype       = "risib",
                         geno_file       = "BXD_geno.csv",
                         geno_transposed = TRUE,
                         geno_codes      = list(B = 1, D = 2),
                         xchr            = "X",
                         gmap_file       = "BXD_gmap.csv",
                         pmap_file       = "BXD_pmap.csv",
                         crossinfo_file  = "BXD_crossinfo.csv",
                         crossinfo_codes = c("BxD" = 0),
                         alleles         = c("B", "D"),
                         na.strings      = na_strings,
                         overwrite       = TRUE)





