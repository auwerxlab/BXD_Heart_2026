
library(qtl)
library(ggplot2)
library(data.table)


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
## save .txt metadata as RDS
###############

saveDir <- "./Data/input_data/metadata_RNAseq"
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

metadata <- fread("./Data/RNAseq_metadata.tsv", data.table = F, stringsAsFactors = F)
saveRDS(metadata, "./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")


###############
## check if all fq files exist
###############

metadata  <- readRDS("./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")
files_all <- paste0(fq_dir_path, "/", c(metadata$file_name_1, metadata$file_name_2, metadata$file_name_3, metadata$file_name_4))
tmpLogi   <- file.exists(files_all)
if(any(!tmpLogi)){
  cat("\n!!! The following .fq.gz files are missing. Please download first all the files from GEO under accession ID GSE330676 and save them in folder './Data/RNA_seq_fq_files' !!!\n",
      paste(paste0("\t--> ", gsub(".*\\/", "", files_all[!tmpLogi])), collapse = "\n"), "\n\n")
  stop()
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
tt              <- mclapply(mc.cores = 40, X = 1:nrow(metadata), FUN = function(x){
  # x <- 1

  forward.fq <- paste0(fq_dir_path, "/", c(metadata$file_name_1[x], metadata$file_name_2[x]))
  reverse.fq <- paste0(fq_dir_path, "/", c(metadata$file_name_3[x], metadata$file_name_4[x]))
  stopifnot(all(file.exists(c(forward.fq, reverse.fq))))
 
  if(file.exists(paste0(concatenatedDir, "/", gsub("\\/", "", metadata$sample_name[x]), "_concatenated_1.fq.gz")) & file.exists(paste0(concatenatedDir, "/",  gsub("\\/", "", metadata$sample_name[x]), "_concatenated_2.fq.gz"))){
    cat("Sample already concatenated!\n")
    return(folder.x)
  } else{
    cat(paste0("Concatenating multiple lanes sequencing data for sample: ", x, "...\n"))
  }
  
  command.forward <- paste0("cat ", paste(forward.fq, collapse = " "), " > ", paste0(concatenatedDir, "/", gsub("\\/", "", metadata$sample_name[x]), "_concatenated_1.fq.gz"))
  system(command.forward)
  command.reverse <- paste0("cat ", paste(reverse.fq, collapse = " "), " > ", paste0(concatenatedDir, "/", gsub("\\/", "", metadata$sample_name[x]), "_concatenated_2.fq.gz"))
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

metadata              <- readRDS("./Data/input_data/metadata_RNAseq/metadata_RNAseq.RDS")
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
nn       <- mclapply(mc.cores = 5, X = 1:nrow(metadata), FUN = function(s){
  # s <- 1
  
  sample.s         <- gsub("\\/", "_", metadata$sample_name[s])
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
  colnames(out) <- c("gene", metadata$sample_name[gsub("\\/", "_", metadata$sample_name) == x])
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

write.table(countMatrix, paste0(dirSave, "/STAR_RNAseq_count_matrix.txt"), sep = "\t", row.names = T, quote = T)
saveRDS(countMatrix, paste0(dirSave, "/STAR_RNAseq_count_matrix.RDS"))
