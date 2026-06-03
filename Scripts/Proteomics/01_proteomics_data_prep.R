


library(limma)
library(data.table)


###########
## save text files to RDS
###########

saveDir <- "./Data/input_data/proteomics/coon_formatted"
if(!dir.exists(saveDir)){
  dir.create(saveDir)
}

raw.data          <- fread("./Data/BXDMice_heart_proteomicsData.txt", data.table = F, stringsAsFactors = F)
meta.df           <- fread("./Data/BXDMice_heart_proteomics_samples_metadata.txt", data.table = F, stringsAsFactors = F)
meta.proteins.df  <- fread("./Data/BXDMice_heart_proteomics_proteins_metadata.txt", data.table = F, stringsAsFactors = F)

rownames(raw.data) <- raw.data[, 1]
raw.data           <- raw.data[, -1]

rownames(meta.df)  <- meta.df[, 1]
meta.df            <- meta.df[, -1]

rownames(meta.proteins.df) <- meta.proteins.df[, 1]
meta.proteins.df           <- meta.proteins.df[, -1]

saveRDS(raw.data, paste0(saveDir, "/raw_intensity_table_samples_excluded.RDS"))
saveRDS(meta.df, paste0(saveDir, "/sample_metadata.RDS"))
saveRDS(meta.proteins.df, paste0(saveDir, "/proteins_metadata.RDS"))


