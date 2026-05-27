# clear working space
rm(list=ls())
# # clear all plots
# dev.off()
# clear console
cat("\014")

library(qtl)
library(qtl2)
library(ggplot2)
library(cowplot)
library(Rsamtools)
library(data.table)
library(rtracklayer)
library(GenomicFeatures)
library(VariantAnnotation)



dataPath                     <- "./Data"
plotsDirPath                 <- "./Plots"
fastQC_toolPath              <- "./Tools/FastQC/fastqc"
starAlignerPath              <- "./Tools/STAR-2.7.9a/bin/Linux_x86_64_static/STAR"
metadata_RNAseq_folder       <- "./Data/input_data/metadata_RNAseq"
reference_genome_indexPath   <- "./Data/input_data/reference_genome_indexed/Mus_musculus"


fq_dir_path <- "./Data/RNA_seq_fq_files"


###############
## download tools
###############

if(!dir.exists("./Data/Tools")){
  dir.create("./Data/Tools", recursive = T)
}
if(!file.exists("./Data/Tools/STAR_2.7.11b.zip")){
  system("wget -P ./Data/Tools https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.zip")
  system("unzip ./Data/Tools/STAR-2.7.9a -d ./Data/Tools/")
}
if(!file.exists("./Data/Tools/fastqc_v0.12.1.zip")){
  system("wget -P ./Data/Tools https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip")
  system("unzip ./Data/Tools/fastqc_v0.12.1.zip -d ./Data/Tools/")
}

###############
## download Mus musculus reference genome (last version - GRCm39.104)
###############

saveDir <- "./Data/input_data/reference_genome/Mus_musculus/GRCm38_release-102"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
gtf_ftp_url              <- "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
primary_assembly_ftp_url <- "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
if(!file.exists(paste0(saveDir, "/", gsub(".*\\/", "", gtf_ftp_url)))){
  system(paste0("wget -P ", saveDir, " ", gtf_ftp_url))
}
if(!file.exists(paste0(saveDir, "/", gsub(".*\\/", "", primary_assembly_ftp_url)))){
  system(paste0("wget -P ", saveDir, " ", primary_assembly_ftp_url))
}

###############
## prepare gene conversion table
###############

saveDir <- "./Data/RNAseq_processing/geneConversionTables"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}
gtfref <- list.files("./Data/input_data/reference_genome/Mus_musculus/GRCm38_release-102", recursive = T, full.names = T, pattern = "gtf.gz")
if(!file.exists(gsub(".gz$", "", gtfref))){
  R.utils::gunzip(gtfref, remove = F) # deflate
}
gtf    <- rtracklayer::import(str_replace(gtfref, ".gz", ""))
gtf_df <- as.data.frame(gtf, stringsAsFactors = F)

geneConversionTable <- unique(dplyr::select(gtf_df, c("gene_id", "gene_name", "gene_biotype")))
write.table(geneConversionTable, paste0(saveDir, "/geneConversionTable_GRCm38_release-102.csv"), row.names = F, sep = ",")
saveRDS(geneConversionTable, paste0(saveDir, "/geneConversionTable_GRCm38_release-102.RDS"))
write.table(gtf_df, paste0(saveDir, "/geneConversionTable_GRCm38_release-102_complete_gtf.csv"), row.names = F, sep = ",")
saveRDS(gtf_df, paste0(saveDir, "/geneConversionTable_GRCm38_release-102_complete_gtf.RDS"))


###############
## do QC of RNA-seq samples
###############

RNAseq_QC_Path <- "./Data/RNAseq_processing/RNAseq_QC/samples"
if(!dir.exists(RNAseq_QC_Path)){
  dir.create(RNAseq_QC_Path)
}
fqfiles         <- list.files(fq_dir_path, full.names = T, recursive = T, pattern = "\\.fq\\.gz$")
names(fqfiles)  <- gsub(".fq.gz", "", gsub("(.)+\\/", "", fqfiles))
done.dirs       <- list.dirs(RNAseq_QC_Path, recursive = F)
html            <- list.files(done.dirs, pattern = "html")
done.file       <- gsub("_fastqc.html", "", html)
remaining       <- setdiff(names(fqfiles), done.file)
remaining.files <- fqfiles[remaining]

tt <- parallel::mclapply(mc.cores = 20, X = remaining.files, FUN = function(f){
  s      <- gsub("(.)+\\/", "", f)
  fqcdir <- paste0(RNAseq_QC_Path, "/", s)
  dir.create(fqcdir)
  fqcommand <- paste0("./Tools/FastQC/fastqc ", unname(f)," --outdir ", RNAseq_QC_Path, "/", s)
  system(fqcommand)
  return(s)
})

#######
## aggregate multiqc results (multiqc was installed with miniconda)
#######

RNAseq_QC_Path_aggregated <- "./Data/RNAseq_processing/RNAseq_QC/summary"
if(!dir.exists(RNAseq_QC_Path_aggregated)){
  dir.create(RNAseq_QC_Path_aggregated)
}
system("docker pull multiqc/multiqc")
command <- paste0("docker run --rm",
                  " -v ", normalizePath(RNAseq_QC_Path), ":", normalizePath(RNAseq_QC_Path),
                  " -v ", normalizePath(RNAseq_QC_Path_aggregated), ":", normalizePath(RNAseq_QC_Path_aggregated),
                  " multiqc/multiqc multiqc ",  normalizePath(RNAseq_QC_Path), " --outdir ", normalizePath(RNAseq_QC_Path_aggregated))
system(command)



###############
## concatenate multiple samples (multiple lanes fasta files) for STAR mapping
###############

concatenatedDir <- "./Data/RNAseq_processing/concatenated_fastq"
if(!dir.exists(concatenatedDir)){
  dir.create(concatenatedDir)
}
metadata        <- readRDS("./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")
allFqFolders    <- dir(fq_dir_path, full.names = T, recursive = T)
tt              <- mclapply(mc.cores = 40, X = metadata$simplifiedSampleName, FUN = function(x){
  # x <- metadata$simplifiedSampleName[76]
  # x <- "BXD1_CD1"
  
  sample.x   <- gsub("C57BL\\/6J", "C57BL-6J", gsub("DBA\\/2J", "DBA-2J", gsub("\\_", "-", x)))
  if(x == "DBA/2J_HFD1"){
    sample.x <- "DBA2J-HFD1"
  }
  folder.x   <- allFqFolders[grepl(sample.x, allFqFolders)]
  fqfiles    <- list.files(folder.x, full.names = T)
  forward.fq <- c(fqfiles[grepl("_L1_", fqfiles) & grepl("_1\\.fq", fqfiles)], fqfiles[grepl("_L2_", fqfiles) & grepl("_1\\.fq", fqfiles)])
  reverse.fq <- c(fqfiles[grepl("_L1_", fqfiles) & grepl("_2\\.fq", fqfiles)], fqfiles[grepl("_L2_", fqfiles) & grepl("_2\\.fq", fqfiles)])
 
  if(file.exists(paste0(concatenatedDir, "/", gsub("\\/", "_", x), "_concatenated_1.fq.gz")) & file.exists(paste0(concatenatedDir, "/", gsub("\\/", "_", x), "_concatenated_2.fq.gz"))){
    cat("Sample already concatenated!\n")
    return(folder.x)
  } else{
    cat(paste0("Concatenating multiple lanes sequencing data for sample: ", x, "...\n"))
  }
  
  command.forward <- paste0("cat ", paste(forward.fq, collapse = " "), " > ", paste0(concatenatedDir, "/", gsub("\\/", "_", x), "_concatenated_1.fq.gz"))
  system(command.forward)
  command.reverse <- paste0("cat ", paste(reverse.fq, collapse = " "), " > ", paste0(concatenatedDir, "/", gsub("\\/", "_", x), "_concatenated_2.fq.gz"))
  system(command.reverse)
})


###############
## generate genome (build STAR indexed genome from reference genome) for 1st 
## pass STAR alignment.
## The reads have length 100, therefore sjdbOverhang is set to 100 - 1
###############

saveDir <- paste0(dataPath, "/RNAseq_processing/reference_genome_STAR_indexed_pass1")
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}
gtfref         <- "./Data/input_data/reference_genome/Mus_musculus/GRCm38_release-102/Mus_musculus.GRCm38.102.gtf.gz"
fastaref       <- "./Data/input_data/reference_genome/Mus_musculus/GRCm38_release-102/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
if(!file.exists(gsub(".gz$", "", gtfref))){
  R.utils::gunzip(gtfref, remove = F) # deflate
}
if(!file.exists(gsub(".gz$", "", fastaref))){
  R.utils::gunzip(fastaref, remove = F) # deflate
}
nThreads       <- max(c(1, round(parallel::detectCores() / 2)))
readLength     <- 100
command        <- paste0(starAlignerPath,
                         " --runThreadN ", nThreads, 
                         " --runMode genomeGenerate",
                         " --genomeDir ", saveDir, 
                         " --genomeFastaFiles ", gsub(".gz$", "", fastaref), 
                         " --sjdbGTFfile ", gsub(".gz$", "", gtfref), 
                         " --sjdbOverhang ", readLength - 1)
system(command)


###############
## STAR mapping / alignment - 1st pass
###############

metadata              <- readRDS("./Data/metadata_RNAseq.RDS")
tempLocalFolder_fasta <- "./tmp_fastaFiles"
tempLocalFolder_align <- "./tmp_starAlign"
starAlignedPath       <- "./Data/RNAseq_processing/star_alignment/GRCm38_release-102/pass1"
for(i in c(tempLocalFolder_fasta, tempLocalFolder_align, starAlignedPath)){
  if(!dir.exists(i)){
    dir.create(i, recursive = T)
  }
}
fqfiles        <- list.files("./Data/RNAseq_processing/concatenated_fastq", pattern = "\\.fq\\.gz$", full.names = T, recursive = T)
fqfiles        <- fqfiles[grepl("_concatenated_", fqfiles)]

nThreads <- 10
nn       <- mclapply(mc.cores = 5, X = metadata$simplifiedSampleName, FUN = function(s){
  # s <- metadata$simplifiedSampleName[2]
  
  sample.s         <- gsub("\\/", "_", s)
  fastaFiles       <- c(fqfiles[grepl(sample.s, fqfiles) & grepl("_1\\.fq", fqfiles)], fqfiles[grepl(sample.s, fqfiles) & grepl("_2\\.fq", fqfiles)])
  fastaFiles.local <- gsub(".*\\/", paste0(tempLocalFolder_fasta, "/"), fastaFiles)
  if(length(fastaFiles) != 2){
    return(NULL)
  }
  # skip if already mapped
  if(file.exists(paste0(starAlignedPath, "/refalign.", sample.s, ".ReadsPerGene.out.tab"))){
    cat(paste0("Sample ", sample.s, " already mapped succesfully\n\n"))
    return(NULL)
  }
  cat("Making local copy of fastq files...\n")
  system(paste0("cp ", paste(fastaFiles, collapse = " "), " ", tempLocalFolder_fasta))
  cat("Running STAR alignment...\n")
  alignCommand <- paste0(starAlignerPath, 
                         " --runThreadN ", nThreads, 
                         " --genomeDir ", paste0(dataPath, "/RNAseq_processing/reference_genome_STAR_indexed_pass1"),
                         " --quantMode GeneCounts TranscriptomeSAM",
                         " --readFilesCommand zcat",
                         " --outSAMtype BAM SortedByCoordinate",
                         " --readFilesIn ", paste(fastaFiles.local, collapse = " "),
                         " --outFileNamePrefix ", tempLocalFolder_align, "/refAlign.", sample.s, ".")
  system(alignCommand)
  cat("Alignment complete, transferring output to common server...\n")
  starOutputs  <- list.files(tempLocalFolder_align, pattern = sample.s, full.names = T)
  moveCommands <- paste0("mv ", starOutputs, " ", starAlignedPath)
  sapply(moveCommands, system)
  print(paste0("Cleaning temporary files"))
  system(paste0("rm ", paste(fastaFiles.local, collapse = " ")))
  cat(paste0("Sample ", sample.s, " mapped successfully\n\n"))
  sample.s
})
unlink(tempLocalFolder_fasta, recursive = T)
unlink(tempLocalFolder_align, recursive = T)


###############
## STAR alignment QC - 1st pass
###############

starAlignedPath    <- "./Data/RNAseq_processing/star_alignment/GRCm38_release-102/pass1"
star_qc_folderPath <- "./Data/RNAseq_processing/starAlignment_pass1_QC"
tmpFolderPath      <- paste0(folderPath, "/tmp_logs")
for(i in c(star_qc_folderPath, tmpFolderPath)){
  if(!dir.exists(i)){
    dir.create(i, recursive = T)
  }
}
logFiles <- list.files(starAlignedPath, pattern = "Log.final.out$", recursive = T, full.names = T)
system(paste0("cp ", paste(logFiles, collapse = " "), " ", tmpFolderPath))

system("docker pull multiqc/multiqc")
command <- paste0("docker run --rm",
                  " -v ", normalizePath(tmpFolderPath), ":", normalizePath(tmpFolderPath),
                  " -v ", normalizePath(star_qc_folderPath), ":", normalizePath(star_qc_folderPath),
                  " multiqc/multiqc multiqc ",  normalizePath(tmpFolderPath), " --outdir ", normalizePath(star_qc_folderPath))
system(command)


###############
## get STAR RNAseq gene-level counts
###############

dirSave <- "./Data/RNAseq_processing/countMatrix/geneLevel"
if(!dir.exists(dirSave)){
  dir.create(dirSave, recursive = T)
}

metadata   <- readRDS("./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")
rpg        <- list.files(paste0(dataPath, "/RNAseq_processing/star_alignment/GRCm38_release-102/pass1"), pattern = "ReadsPerGene", full.names = T, recursive = T)
names(rpg) <- gsub("^(.)+refAlign\\.", "", gsub("\\.ReadsPerGene(.)+$", "", rpg))

# Here we have have an unstranded design, therefore we extract the first two 
# columns of the "ReadsPerGene" files, which contain the gene ID and the number 
# of reads that map on both strands
allReads   <- lapply(names(rpg), FUN = function(x){
  # x <- names(rpg)[1]
  out           <- as.data.frame(fread(rpg[x], select = 1:2), stringsAsFactors = FALSE)
  colnames(out) <- c("gene", metadata$simplifiedSampleName[gsub("\\/", "_", metadata$simplifiedSampleName) == x])
  return(out)
})

allGenes    <- unique(unlist(lapply(allReads, function(x){as.character(x$gene)})))
countMatrix <- do.call(cbind, lapply(allReads, function(x){
  # x <- allReads[[1]]
  out <- x[match(allGenes, x$gene), 2, drop = F]
  if(!identical(x[match(allGenes, x$gene), 1], allGenes)){
    stop("ERROR")
  }
  out
}))
rownames(countMatrix) <- allGenes
countMatrix           <- countMatrix[-(1:4), ]

saveRDS(countMatrix, paste0(dirSave, "/STAR_RNAseq_count_matrix.RDS"))


###############
## perform DEA
###############

metadata            <- readRDS("./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")
metadata$strainDiet <- paste0(metadata$correctStrain, "_", metadata$diet)
countMatrix_qcpass  <- readRDS("./Data/RNAseq_processing/countMatrix/geneLevel/STAR_RNAseq_count_matrix.RDS")
metadata            <- metadata[metadata$simplifiedSampleName %in% colnames(countMatrix_qcpass), ]
metadata            <- metadata[match(colnames(countMatrix_qcpass), metadata$simplifiedSampleName), ]
rownames(metadata)  <- metadata$simplifiedSampleName
colnames(metadata)  <- gsub("\\.", "_", gsub("\\[", "", gsub("\\]\\.|\\/", "_", colnames(metadata))))
colnames(metadata)  <- make.names(colnames(metadata))
identical(rownames(metadata), colnames(countMatrix_qcpass))
countMatrix_qcpass  <- as.matrix(countMatrix_qcpass)


numericCov <- c("concentration_[ng/uL].BGI", "volume_[uL].BGI", "totalMass_[ug].BGI", "RIN.BGI", "28S/18S.BGI")
numericCov <- gsub("\\.", "_", gsub("\\[", "", gsub("\\]\\.|\\/", "_", numericCov)))
numericCov <- make.names(numericCov)
for(i in numericCov){
  metadata[[paste0(i, "_range")]] <- convert2factor(metadata[[i]], 4)
}


y               <- edgeR::DGEList(counts = countMatrix_qcpass, samples = metadata)
keep.exprs      <- edgeR::filterByExpr(y, group = "strain_diet")
y               <- y[keep.exprs, , keep.lib.sizes = FALSE]
hasCounts       <- rowSums(is.na(y$counts)) == 0
y               <- y[hasCounts, , keep.lib.sizes = FALSE]
y               <- edgeR::calcNormFactors(y)


y$samples$strain        <- gsub("\\/", "_", y$samples$strain)
y$samples$strainDiet    <- gsub("\\/", "_", y$samples$strainDiet)
geneTable               <- readRDS("./Data/RNAseq_processing/geneConversionTables/geneConversionTable_GRCm38_release-102.RDS")

# prepare contrasts
globalDietContrast                   <- c("HFD_vs_CD_diet" = "(HFD)-(CD)")

# define covariates
covariates <-  c("concentration_[ng/uL].BGI", "RIN.BGI", "28S/18S.BGI")
covariates <- gsub("\\.", "_", gsub("\\[", "", gsub("\\]\\.|\\/", "_", covariates)))
covariates <- make.names(covariates)
numericCov <- c("concentration_[ng/uL].BGI", "volume_[uL].BGI", "totalMass_[ug].BGI", "RIN.BGI", "28S/18S.BGI")
numericCov <- gsub("\\.", "_", gsub("\\[", "", gsub("\\]\\.|\\/", "_", numericCov)))
numericCov <- make.names(numericCov)

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

if(!file.exists("./Data/Tools/STAR_2.7.11b.zip")){
  
  system("unzip ./Data/Tools/STAR-2.7.9a -d ./Data/Tools/")
}

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






