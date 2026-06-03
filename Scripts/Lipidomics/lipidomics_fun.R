
color_pal <- c("ST"               = "#7644c3", # grandparent class; Sterol Lipids,
               "CE"               = "#be29ec", # parent class;      Cholesteryl Esters
               "CE_cl"            = "#be29ec", # lipid class;
               
               # "Cer"        = "#99badd",
               
               "SM"               = "#ADD8E6", # parent class;      Sphingomyelins
               "Cer"              = "#00b4d8", # parent class;      Ceramides
               "HexCer"           = "#4169E1", # parent class;      Hexosyl Ceramides
               "SP"               = "#0096FF", # parent class;      Sphingolipids
               "SP_gp"            = "#120a8f", # grandparent class; Sphingolipids
               
               "CerP_cl"          = "#00c6d8", # lipid class;
               "Cer[NS]_cl"       = "#00d8c6", # lipid class;
               "HexCer[NS]_cl"    = "#4169E1", # lipid class;
               
               "SM_cl"            = "#ADD8E6", # lipid class;
               "SP_cl"            = "#0096FF", # lipid class;
               
               
               "GP"               = "#F38701", # grandparent class; Glycerophospholipids (Phospholipids)
               
               "PC"               = "#E30B5C", # parent class;      Phosphatidylcholines
               "PE"               = "#a30000", # parent class;      Phosphatidylethanolamines
               "PI"               = "#f93339", # parent class;      Phosphatidylinositols
               "PS"               = "#FF5F1F", # parent class;      Phosphatidylserines
               "PG"               = "#F08000", # parent class;      Phosphatidylglycerols
               
               "PC_cl"            = "#e30b92", # lipid class;
               "PC[OH]_cl"        = "#e30b5c", # lipid class;
               "Plasmanyl-PC_cl"  = "#e30b26", # lipid class;
               "Plasmenyl-PC_cl"  = "#f75994", # lipid class;

               "PE_cl"            = "#2d0000", # lipid class;
               "PE-NMe_cl"        = "#680000", # lipid class;
               "PE-NMe2_cl"       = "#a30000", # lipid class;
               "Plasmanyl-PE_cl"  = "#904e4e", # lipid class;
               "Plasmenyl-PE_cl"  = "#6a3939", # lipid class;
               
               "PG_cl"            = "#F08000", # lipid class;
               
               "PI_cl"            = "#d15b5f", # lipid class;
               "PI[OH]_cl"        = "#f13b41", # lipid class;
               
               "PS_cl"            = "#FF5F1F", # lipid class;
               
               "LysoPC"           = "#FF69B4", # parent class;      LysoPhosphatidylcholines
               "LysoPE"           = "#F8C8DC", # parent class;      LysoPhosphatidylethanolamines
               "LysoPI"           = "#ff99cc", # parent class;      LysoPhosphatidylinositols
               "LysoPG"           = "#F3CFC6", # parent class;      LysoPhosphatidylglycerols
               
               "LysoPC_cl"        = "#FF69B4", # lipid class;
               "LysoPE_cl"        = "#F8C8DC", # lipid class;
               "LysoPI_cl"        = "#ff99cc", # lipid class;
               "LysoPG_cl"        = "#F3CFC6", # lipid class;
               
               "CL"               = "#B87333", # parent class;      Cardiolypins
               "MLCL"             = "#E97451", # parent class;      Monolysocardiolipins
          
               "CL_cl"            = "#B87333", # lipid class;      
               "MLCL_cl"          = "#E97451", # lipid class;
 
               "GL"               = "#097969", # grandparent class; Glycerolipids
               "TG_cl"            = "#32CD32",
               "TG"               = "#32CD32", # parent class;      Triacylglycerols
               "DG"               = "#AFE1AF", # parent class;      Diacylglycerols
               
               "TG_cl"            = "#1b6e1b", # lipid class;
               "Alkenyl-TG_cl"    = "#2aad2a", # lipid class;
               "DG_cl"            = "#daf1da", # lipid class;
               "Alkenyl-DG_cl"    = "#85d185", # lipid class;
               "Alkanyl-DG_cl"    = "#98be98", # lipid class;
               
               
               "FA"               = "#FFEA00", # grandparent class; Fatty Acyls
               "AC"               = "#FFFAA0", # parent class;      Acylcarnitines
               "AC_cl"            = "#FFFAA0", # lipid class;
               
               
               "rest"             = "grey")
# pals::pal.bands(color_pal)

# 
# |Lipid.Parent.Class |Lipid.Grandparent.Class |   avg_perc|
# |CE                 |ST                      |  0.0017780|
# |SP                 |SP                      |  0.0018229|
# |HexCer             |SP                      |  0.0050653|
# |Cer                |SP                      |  1.3258836|
# |SM                 |SP                      |  4.9780964|
# |LysoPG             |GP                      |  0.0035998|
# |LysoPI             |GP                      |  0.0085727|
# |LysoPE             |GP                      |  0.1218744|
# |MLCL               |GP                      |  0.2563993|
# |PG                 |GP                      |  0.9513828|
# |PS                 |GP                      |  1.1376475|
# |LysoPC             |GP                      |  1.3789241|
# |PI                 |GP                      |  2.2218376|
# |CL                 |GP                      |  5.3788114|
# |PE                 |GP                      | 14.9652156|
# |PC                 |GP                      | 42.0955800|
# |DG                 |GL                      |  1.6959124|
# |TG                 |GL                      | 23.4578758|
# |AC                 |FA                      |  0.0137204|


# |Lipid.Class  |Lipid.Parent.Class |Lipid.Grandparent.Class | avg_perc|
# |:------------|:------------------|:-----------------------|--------:|
# |AC           |AC                 |FA                      | 18.07848|
# |Alkanyl-DG   |DG                 |GL                      | 23.65417|
# |Alkenyl-DG   |DG                 |GL                      | 23.22150|
# |DG           |DG                 |GL                      | 23.09225|
# |Alkenyl-TG   |TG                 |GL                      | 21.07984|
# |TG           |TG                 |GL                      | 23.43835|
# |CL           |CL                 |GP                      | 21.95172|
# |LysoPC       |LysoPC             |GP                      | 22.09108|
# |LysoPE       |LysoPE             |GP                      | 21.13936|
# |LysoPG       |LysoPG             |GP                      | 17.54877|
# |LysoPI       |LysoPI             |GP                      | 18.11035|
# |MLCL         |MLCL               |GP                      | 20.62465|
# |PC           |PC                 |GP                      | 22.68066|
# |PC[OH]       |PC                 |GP                      | 20.44983|
# |Plasmanyl-PC |PC                 |GP                      | 22.36136|
# |Plasmenyl-PC |PC                 |GP                      | 21.08179|
# |PE           |PE                 |GP                      | 21.85448|
# |PE-NMe       |PE                 |GP                      | 21.59902|
# |PE-NMe2      |PE                 |GP                      | 21.00391|
# |Plasmanyl-PE |PE                 |GP                      | 22.04455|
# |Plasmenyl-PE |PE                 |GP                      | 23.24759|
# |PG           |PG                 |GP                      | 20.26303|
# |PI           |PI                 |GP                      | 22.44392|
# |PI[OH]       |PI                 |GP                      | 20.49351|
# |PS           |PS                 |GP                      | 23.77406|
# |CerP         |Cer                |SP                      | 19.06630|
# |Cer[NS]      |Cer                |SP                      | 22.24232|
# |HexCer[NS]   |HexCer             |SP                      | 18.79590|
# |SM           |SM                 |SP                      | 23.71957|
# |SP           |SP                 |SP                      | 19.59298|
# |CE           |CE                 |ST                      | 19.37499|


diagnosticPlots <- function(data.df, meta.df){
  
  # data.df <- lipids.df
  # meta.df <- meta
  
  # data.df <- lipid.list[[x]]
  # meta.df <- lipid.list$MetaData[lipid.list$MetaData$SampleName %in% colnames(lipid.list[[x]]), ]
  
  cat("\nCreating the following plots:\n")
  
  lipids.df            <- data.df
  lipids.df            <- lipids.df[!lipids.df$isUnknown, ]
  lipids.df            <- dplyr::select(lipids.df, c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min.", colnames(lipids.df)[grepl("BDX|BXD|C57|DBA", colnames(lipids.df))]))
  lipids.df            <- lipids.df[, !grepl("PWT|Blank", colnames(lipids.df))]
  tmpLogi              <- grepl("BXD|C57B|DBA2", colnames(lipids.df))
  lipids.df[, tmpLogi] <- log2(lipids.df[, tmpLogi])
  meta                 <- meta.df
  
  cat("\t--> Duplicated IDs retention time...\n")
  dfPlot.duplicatedID <- lipids.df[grepl("RT", lipids.df$Identification), ]
  dfPlot.duplicatedID <- dfPlot.duplicatedID[!dfPlot.duplicatedID$isUnknown, ]
  # dfPlot.duplicatedID <- dplyr::select(dfPlot.duplicatedID, c("Identification", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", colnames(dfPlot.duplicatedID)[grepl("BDX|BXD|C57|DBA", colnames(dfPlot.duplicatedID))]))
  dfPlot.duplicatedID <- dfPlot.duplicatedID[, !grepl("PWT|Blank", colnames(dfPlot.duplicatedID))]
  dfPlot.duplicatedID <- unique(dplyr::select(dfPlot.duplicatedID, "Identification", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min."))
  dfPlot.duplicatedID$simpleID <- gsub("_RT.*$", "", dfPlot.duplicatedID$Identification)
  tmp.avg <- data.table(dfPlot.duplicatedID)[, list(avg = mean(Retention.Time..min.)), by = c("simpleID")]
  tmp.avg <- tmp.avg[order(tmp.avg$avg, decreasing = T)]
  dfPlot.duplicatedID$simpleID <- factor(dfPlot.duplicatedID$simpleID, levels = unique(tmp.avg$simpleID))
  pl.scatter.duplicatedLipids.retTime <- ggplot(dfPlot.duplicatedID, aes(x = simpleID, y = Retention.Time..min., color = Lipid.Class)) +
    geom_point() +
    # scale_color_manual(values = colors) +
    scale_color_manual(values = unname(pals::alphabet2(length(unique(dfPlot.duplicatedID$Lipid.Class))))) +
    xlab(paste0("Duplicated Lipid Species [N=", length(unique(dfPlot.duplicatedID$simpleID)),"/", length(unique(gsub("_RT.*$", "", lipids.df$Identification))), "]")) +
    ylab("Retention Time [min]") +
    theme_classic() +
    theme(axis.text.x        = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(),
          panel.grid.major.x = element_line(),
          legend.key.size    = unit(0.2, "cm"))
  # pl.scatter.duplicatedLipids.retTime

  
  cat("\t--> Duplicated IDs amplitude range...\n")
  dfPlot.duplicatedID <- lipids.df[grepl("RT", lipids.df$Identification), ]
  dfPlot.duplicatedID <- dfPlot.duplicatedID[!dfPlot.duplicatedID$isUnknown, ]
  dfPlot.duplicatedID <- dfPlot.duplicatedID[, !grepl("PWT|Blank", colnames(dfPlot.duplicatedID))]
  dfPlot.duplicatedID <- dplyr::select(dfPlot.duplicatedID, -c("Identifier", "Retention.Time..min.", "Polarity", "isUnknown", "cLength", "nbDoubleBonds"))
  dfPlot.duplicatedID$simpleID <- gsub("_RT.*$", "", dfPlot.duplicatedID$Identification)
  dfPlot.duplicatedID <- reshape2::melt(dfPlot.duplicatedID, id.vars = c("Identification", "Lipid.Class", "simpleID"))

  dfRange <- data.table(dfPlot.duplicatedID)[, list(rangeDiff = diff(base::range(value))), by = c("simpleID", "variable", "Lipid.Class")]
  pl.scatter.duplLipids.range <- ggplot(dfRange, aes(x = simpleID, y = rangeDiff, color = Lipid.Class)) +
    geom_point() +
    # scale_color_manual(values = colors) +
    scale_color_manual(values = unname(pals::alphabet2(length(unique(dfRange$Lipid.Class))))) +
    xlab(paste0("Duplicated Lipid Species [N=", length(unique(dfPlot.duplicatedID$simpleID)),"/", length(unique(gsub("_RT.*$", "", lipids.df$Identification))), "]")) +
    ylab("log2(Intensity)\nrange withing samples") +
    theme_classic() +
    theme(axis.text.x        = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(),
          panel.grid.major.x = element_line(),
          legend.key.size    = unit(0.2, "cm"))
  # pl.scatter.duplLipids.range

  cat("\t--> Duplicated IDs amplitude range (sorted)...\n")
  tmp.avg <- data.table(dfRange)[, list(avg = mean(rangeDiff)), by = c("simpleID")]
  tmp.avg <- tmp.avg[order(tmp.avg$avg, decreasing = T), ]
  dfRange$simpleID <- factor(dfRange$simpleID, levels = unique(tmp.avg$simpleID))
  pl.scatter.duplLipids.range.ordered <- ggplot(dfRange, aes(x = simpleID, y = rangeDiff, color = Lipid.Class)) +
    geom_point() +
    scale_color_manual(values = unname(pals::alphabet2(length(unique(dfRange$Lipid.Class))))) +
    xlab(paste0("Duplicated Lipid Species [N=", length(unique(dfPlot.duplicatedID$simpleID)),"/", length(unique(gsub("_RT.*$", "", lipids.df$Identification))), "]")) +
    ylab("log2(Intensity)\nrange withing samples") +
    theme_classic() +
    theme(axis.text.x        = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(),
          panel.grid.major.x = element_line(),
          legend.key.size    = unit(0.2, "cm"))
  # pl.scatter.duplLipids.range

  # chain length and saturation distribution also with intensities
  # duplicated lipids: plot retention time, how many, ...
  
  cat("\t--> Nb. samples by batch...\n")
  dfPlot <- unique(dplyr::select(meta, c("SampleName", "Batch")))
  dfPlot <- data.table(dfPlot)[, list(nbSamples = length(unique(SampleName))), by = c("Batch")]
  dfPlot$Batch <- as.character(dfPlot$Batch)
  dfPlot <- dfPlot[order(as.numeric(dfPlot$Batch)), ]
  dfPlot$Batch <- factor(dfPlot$Batch, levels = unique(dfPlot$Batch))
  pl.nbSamples.byBatch <- ggplot(dfPlot, aes(x = Batch, y = nbSamples)) +
    geom_bar(stat = "Identity", fill = "#BEBEBE") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ylab("Nb. Samples") +
    theme_classic()
  # pl.nbSamples.byBatch

  
  cat("\t--> Samples Intensity boxplot...\n")
  dfPlot <- reshape2::melt(lipids.df, id.vars = c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min."))
  dfPlot <- type.convert(dfPlot, as.is = T)
  dfPlot <- merge(dfPlot, meta, by.x = "variable", by.y = "SampleName", all.x = T, all.y = F)
  dfPlot <- na.omit(dfPlot)

  dfPlot$samplePol <- paste0(dfPlot$variable, " (", dfPlot$Polarity, ")")

  pl.boxlot.byPolarization <- ggplot(dfPlot, aes(x = variable, y = value)) +
    geom_boxplot(size = 0.07, outlier.shape = NA) +
    xlab("Samples") + ylab("log2(intensity)") +
    facet_wrap(~Polarity, ncol = 1) +
    theme_classic() +
    theme(axis.text.x =  element_text(angle = 45, hjust = 1, size = 1.4),
          axis.ticks.x = element_line(size = 0.1))
  # pl.boxlot.byPolarization

  dfPlot$polarity_all <- "+/-"
  pl.boxlot.overall <- ggplot(dfPlot, aes(x = variable, y = value)) +
    geom_boxplot(size = 0.07, outlier.shape = NA) +
    xlab("Samples") + ylab("log2(intensity)") +
    facet_wrap(~polarity_all, ncol = 1) +
    theme_classic() +
    theme(axis.text.x =  element_text(angle = 45, hjust = 1, size = 1.4),
          axis.ticks.x = element_line(size = 0.1))
  # pl.boxlot.overall
  
  cat("\t--> C-length vs saturation (overview)...\n")
  dfPlot.carbChains     <- unique(dplyr::select(dfPlot, c("cLength", "nbDoubleBonds", "Polarity", "Identification")))
  dfPlot.carbChains$idx <- 1:nrow(dfPlot.carbChains)
  dfPlot.carbChains     <- data.table(dfPlot.carbChains)[, list(nb = length(unique(idx))), by = c("cLength", "nbDoubleBonds", "Polarity")]
  pl.dotplot.cLength.vs.saturation.byPolarization <- ggplot(dfPlot.carbChains, aes(x = cLength, y = nbDoubleBonds, size = nb, color = nb)) +
    geom_point() +
    scale_size(range = c(0.5, 4.5)) +
    scale_color_gradientn(colors = pals::plasma(100)) +
    xlab("Carbon Chain Length") + ylab("Nb. Double Bonds") +
    facet_wrap(~Polarity) +
    theme_classic() +
    theme(panel.grid.major = element_line(size = 0.5)) 
    # coord_fixed()
  # pl.dotplot.cLength.vs.saturation.byPolarization


  dfPlot.carbChains              <- unique(dplyr::select(dfPlot, c("cLength", "nbDoubleBonds", "Identification")))
  dfPlot.carbChains$idx          <- 1:nrow(dfPlot.carbChains)
  dfPlot.carbChains$polarity_all <- "+/-"
  dfPlot.carbChains              <- data.table(dfPlot.carbChains)[, list(nb = length(unique(idx))), by = c("cLength", "nbDoubleBonds", "polarity_all")]
  pl.dotplot.cLength.vs.saturation.overall <- ggplot(dfPlot.carbChains, aes(x = cLength, y = nbDoubleBonds, size = nb, color = nb)) +
    geom_point() +
    scale_size(range = c(0.5, 4.5)) +
    scale_color_gradientn(colors = pals::plasma(100)) +
    xlab("Carbon Chain Length") + ylab("Nb. Double Bonds") +
    facet_wrap(~polarity_all) +
    theme_classic() +
    theme(panel.grid.major = element_line(size = 0.5))
  # pl.dotplot.cLength.vs.saturation.overall
  
  cat("\t--> C-length vs saturation (by lipid class)...\n")
  dfPlot.carbChains     <- unique(dplyr::select(dfPlot, c("cLength", "nbDoubleBonds", "Polarity", "Identification", "Lipid.Class")))
  dfPlot.carbChains$idx <- 1:nrow(dfPlot.carbChains)
  dfPlot.carbChains     <- data.table(dfPlot.carbChains)[, list(nb = length(unique(idx))), by = c("cLength", "nbDoubleBonds", "Polarity", "Lipid.Class")]
  pl.dotplot.cLength.vs.saturation.byClass <- ggplot(dfPlot.carbChains, aes(x = cLength, y = nbDoubleBonds, size = nb, color = nb)) +
    geom_point() +
    scale_size(range = c(0.5, 4.5)) +
    scale_color_gradientn(colors = pals::plasma(100)) +
    xlab("Carbon Chain Length") + ylab("Nb. Double Bonds") +
    facet_grid(Lipid.Class~Polarity) +
    theme_classic() +
    theme(panel.grid.major = element_line(size = 0.5))
  # pl.dotplot.cLength.vs.saturation.byClass

  cat("\t--> Saturation histogram...\n")
  dfPlot.saturationRatio       <- unique(dplyr::select(dfPlot, c("cLength", "nbDoubleBonds", "Polarity", "Identification", "Lipid.Class")))
  dfPlot.saturationRatio$ratio <- dfPlot.saturationRatio$cLength / dfPlot.saturationRatio$nbDoubleBonds
  pl.histo.saturationRatio.byPolarization <- ggplot(dfPlot.saturationRatio, aes(x = ratio, fill = Polarity)) +
    geom_histogram(bins = 200, alpha = 0.9) +
    xlab("Ratio (Carbon Chain Length):(Nb. Double Bonds)") +
    ylab("Count") +
    scale_fill_manual(values = c("#6799e4", "#b0222b")) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    facet_wrap(~Polarity) +
    theme_classic()
  # pl.histo.saturationRatio.byPolarization

  dfPlot.saturationRatio$polarity_all <- "+/-"
  pl.histo.saturationRatio.overall <- ggplot(dfPlot.saturationRatio, aes(x = ratio)) +
    geom_histogram(bins = 200, alpha = 0.9) +
    xlab("Ratio (Carbon Chain Length):(Nb. Double Bonds)") +
    ylab("Count") +
    # scale_fill_manual(values = c("#6799e4", "#b0222b")) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    facet_wrap(~polarity_all) +
    theme_classic()
  # pl.histo.saturationRatio.overall
  
  cat("\t--> Nb. lipids by class...\n")
  dfPlot.nbByClass <- unique(dplyr::select(dfPlot, c("Polarity", "Identification", "Lipid.Class")))
  dfPlot.nbByClass <- data.table(dfPlot.nbByClass)[, list(nb = length(unique(Identification))), by = c("Polarity", "Lipid.Class")]
  dfPlot.nbByClass <- dfPlot.nbByClass[order(dfPlot.nbByClass$nb, decreasing = T), ]
  dfPlot.nbByClass$Lipid.Class <- factor(dfPlot.nbByClass$Lipid.Class, levels = unique(dfPlot.nbByClass$Lipid.Class))
  pl.histo.nbLipidByClass.byPolarization <- ggplot(dfPlot.nbByClass, aes(x = Lipid.Class, y = nb)) +
    geom_bar(stat = "Identity") +
    xlab("Lipid Class") + ylab("Nb. Identified Lipid Species") +
    facet_wrap(~Polarity, ncol = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # pl.histo.nbLipidByClass.byPolarization

  dfPlot.nbByClass                <- unique(dplyr::select(dfPlot, c("Identification", "Lipid.Class")))
  dfPlot.nbByClass$polarity_all   <- "+/-"
  dfPlot.nbByClass                <- data.table(dfPlot.nbByClass)[, list(nb = length(unique(Identification))), by = c("polarity_all", "Lipid.Class")]
  dfPlot.nbByClass                <- dfPlot.nbByClass[order(dfPlot.nbByClass$nb, decreasing = T), ]
  dfPlot.nbByClass$Lipid.Class    <- factor(dfPlot.nbByClass$Lipid.Class, levels = unique(dfPlot.nbByClass$Lipid.Class))
  pl.histo.nbLipidByClass.overall <- ggplot(dfPlot.nbByClass, aes(x = Lipid.Class, y = nb)) +
    geom_bar(stat = "Identity") +
    xlab("Lipid Class") + ylab("Nb. Identified Lipid Species") +
    facet_wrap(~polarity_all, ncol = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pl.histo.nbLipidByClass.overall
  
  cat("\t--> Density plots (non annotated and annotated)...\n")
  dfPlot.dens <- lapply(unique(dfPlot$samplePol), function(x){
    # x <- unique(dfPlot$samplePol)[1]
    # out <- density(na.omit(dfPlot$valueLog2[dfPlot$samplePol == x]))
    out <- density(na.omit(dfPlot$value[dfPlot$samplePol == x]))
    data.frame(x          = out$x, 
               y          = out$y, 
               samplePol  = x, 
               sampleName = dfPlot$variable[dfPlot$samplePol == x][1],
               polarity   = dfPlot$Polarity[dfPlot$samplePol == x][1],
               batch      = dfPlot$Batch[dfPlot$samplePol == x][1], 
               diet       = dfPlot$Diet[dfPlot$samplePol == x][1], 
               strain     = dfPlot$Strain[dfPlot$samplePol == x][1], 
               maxPeakX = out$x[which.max(out$y)], maxPeakY = max(out$y), 
               stringsAsFactors = F)
  })
  dfPlot.dens       <- as.data.frame(rbindlist(dfPlot.dens), stringsAsFactors = F)
  dfPlot.dens$batch <- as.character(dfPlot.dens$batch)
  
  set.seed(22)
  colorBy <- c("sampleName", "batch", "diet", "strain")
  set.seed(123)
  pl.density.list.byPolarization <- lapply(colorBy, function(cby){
    # cby <- "batch"
    if(cby == "diet"){
      colors <- c("CD" = "#6799e4", "HFD" = "#b0222b")
    } else{
      colors <- colors()
      set.seed(123)
      colors <- colors[sample(1:length(colors), length(unique(dfPlot.dens[[cby]])))]
    }
    pl.density <- ggplot(dfPlot.dens, aes_string(x = "x", y = "y", color = cby, group = "sampleName")) +
      geom_line(size = 0.2) +
      scale_color_manual(values = colors) +
      xlab("log2(intensity)") + ylab("Density") +
      theme_classic() +
      facet_wrap(~polarity) +
      guides(color = F) +
      theme(panel.grid.major = element_line(size = 0.5))
    # pl.density
  })
  names(pl.density.list.byPolarization) <- paste0("colorBy_", colorBy)

  
  peaksDf       <- unique(dplyr::select(dfPlot.dens, c("sampleName", "polarity", "maxPeakX", "maxPeakY", "batch", "diet", "strain")))
  peaksDf$label <- NA
  peaksDf       <- peaksDf[order(peaksDf$maxPeakY), ]
  peaksDf$label[base::range(which(peaksDf$polarity == "+"))] <- paste0(peaksDf$sampleName[base::range(which(peaksDf$polarity == "+"))], " - B", peaksDf$batch[base::range(which(peaksDf$polarity == "+"))])
  peaksDf$label[base::range(which(peaksDf$polarity == "-"))] <- paste0(peaksDf$sampleName[base::range(which(peaksDf$polarity == "-"))], " - B", peaksDf$batch[base::range(which(peaksDf$polarity == "-"))])
  peaksDf       <- peaksDf[order(peaksDf$maxPeakX), ]
  peaksDf$label[base::range(which(peaksDf$polarity == "+"))] <- paste0(peaksDf$sampleName[base::range(which(peaksDf$polarity == "+"))], " - B", peaksDf$batch[base::range(which(peaksDf$polarity == "+"))])
  peaksDf$label[base::range(which(peaksDf$polarity == "-"))] <- paste0(peaksDf$sampleName[base::range(which(peaksDf$polarity == "-"))], " - B", peaksDf$batch[base::range(which(peaksDf$polarity == "-"))])
  peaksDf$label <- gsub("\\.|_", " ", peaksDf$label)
  avgPeakX      <- mean(peaksDf$maxPeakX)
  peaksDf$posX  <- ifelse(peaksDf$maxPeakX < avgPeakX, peaksDf$maxPeakX - 10, peaksDf$maxPeakX + 10)
  
  # tt <- na.omit(peaksDf)
  
  colorBy <- c("sampleName", "batch", "diet", "strain")
  set.seed(123)
  pl.density.pkAnn.list.byPolarization <- lapply(colorBy, function(cby){
    # cby <- "diet"
    if(cby == "diet"){
      colors <- c("CD" = "#6799e4", "HFD" = "#b0222b")
    } else{
      colors <- colors()
      set.seed(123)
      colors <- colors[sample(1:length(colors), length(unique(dfPlot.dens[[cby]])))]
    }
    pl.density <- ggplot(dfPlot.dens, aes_string(x = "x", y = "y", color = cby, group = "sampleName")) +
      geom_line(size = 0.2) +
      geom_text_repel(data = peaksDf[!is.na(peaksDf$label) & peaksDf$maxPeakX > avgPeakX, ], aes_string(x = "maxPeakX", y = "maxPeakY", color = cby, label = "label"), nudge_x = 8, size = 2.5, segment.size = 0.2, min.segment.length = 0) +
      geom_text_repel(data = peaksDf[!is.na(peaksDf$label) & peaksDf$maxPeakX < avgPeakX, ], aes_string(x = "maxPeakX", y = "maxPeakY", color = cby, label = "label"), nudge_x = -8, size = 2.5, segment.size = 0.2, min.segment.length = 0) +
      scale_color_manual(values = colors) +
      xlab("log2(intensity)") + ylab("Density") +
      theme_classic() +
      facet_wrap(~polarity) +
      guides(color = F) +
      theme(panel.grid.major = element_line(size = 0.5))
    # pl.density
  })
  names(pl.density.pkAnn.list.byPolarization) <- paste0("colorBy_", colorBy)
  
  
  dfPlot.dens <- lapply(unique(dfPlot$sample), function(x){
    # x <- unique(dfPlot$sample)[1]
    out <- density(na.omit(dfPlot$value[dfPlot$sample == x]))
    data.frame(x          = out$x, 
               y          = out$y, 
               sample     = x, 
               sampleName = dfPlot$variable[dfPlot$sample == x][1],
               # polarity   = dfPlot$Polarity[dfPlot$sample == x][1],
               batch      = dfPlot$Batch[dfPlot$sample == x][1], 
               diet       = dfPlot$Diet[dfPlot$sample == x][1], 
               strain     = dfPlot$Strain[dfPlot$sample == x][1], 
               maxPeakX = out$x[which.max(out$y)], maxPeakY = max(out$y), 
               stringsAsFactors = F)
  })
  dfPlot.dens       <- as.data.frame(rbindlist(dfPlot.dens), stringsAsFactors = F)
  dfPlot.dens$batch <- as.character(dfPlot.dens$batch)
  dfPlot.dens$polarity_all <- "+/-"
  
  set.seed(22)
  colorBy <- c("sampleName", "batch", "diet", "strain")
  set.seed(123)
  pl.density.list.overall <- lapply(colorBy, function(cby){
    # cby <- "batch"
    if(cby == "diet"){
      colors <- c("CD" = "#6799e4", "HFD" = "#b0222b")
    } else{
      colors <- colors()
      set.seed(123)
      colors <- colors[sample(1:length(colors), length(unique(dfPlot.dens[[cby]])))]
    }
    pl.density <- ggplot(dfPlot.dens, aes_string(x = "x", y = "y", color = cby, group = "sample")) +
      geom_line(size = 0.2) +
      scale_color_manual(values = colors) +
      xlab("log2(intensity)") + ylab("Density") +
      theme_classic() +
      facet_wrap(~polarity_all) +
      guides(color = F) +
      theme(panel.grid.major = element_line(size = 0.5))
    # pl.density
  })
  names(pl.density.list.overall) <- paste0("colorBy_", colorBy)
  
  cat("\t--> Annotated density plots...\n")
  peaksDf       <- unique(dplyr::select(dfPlot.dens, c("sampleName", "polarity_all", "maxPeakX", "maxPeakY", "batch", "diet", "strain")))
  peaksDf$label <- NA
  peaksDf       <- peaksDf[order(peaksDf$maxPeakY), ]
  peaksDf$label <- paste0(peaksDf$sampleName, " - B", peaksDf$batch)
  peaksDf       <- peaksDf[order(peaksDf$maxPeakX), ]
  peaksDf$label <- paste0(peaksDf$sampleName, " - B", peaksDf$batch)
  peaksDf$label <- gsub("\\.|_", " ", peaksDf$label)
  avgPeakX      <- mean(peaksDf$maxPeakX)
  peaksDf$posX  <- ifelse(peaksDf$maxPeakX < avgPeakX, peaksDf$maxPeakX - 10, peaksDf$maxPeakX + 10)
  
  # tt <- na.omit(peaksDf)
  
  colorBy <- c("sampleName", "batch", "diet", "strain")
  set.seed(123)
  pl.density.pkAnn.list.overall <- lapply(colorBy, function(cby){
    # cby <- "diet"
    if(cby == "diet"){
      colors <- c("CD" = "#6799e4", "HFD" = "#b0222b")
    } else{
      colors <- colors()
      set.seed(123)
      colors <- colors[sample(1:length(colors), length(unique(dfPlot.dens[[cby]])))]
    }
    pl.density <- ggplot(dfPlot.dens, aes_string(x = "x", y = "y", color = cby, group = "sampleName")) +
      geom_line(size = 0.2) +
      geom_text_repel(data = peaksDf[!is.na(peaksDf$label) & peaksDf$maxPeakX > avgPeakX, ], aes_string(x = "maxPeakX", y = "maxPeakY", color = cby, label = "label"), nudge_x = 8, size = 2.5, segment.size = 0.2, min.segment.length = 0) +
      geom_text_repel(data = peaksDf[!is.na(peaksDf$label) & peaksDf$maxPeakX < avgPeakX, ], aes_string(x = "maxPeakX", y = "maxPeakY", color = cby, label = "label"), nudge_x = -8, size = 2.5, segment.size = 0.2, min.segment.length = 0) +
      scale_color_manual(values = colors) +
      xlab("log2(intensity)") + ylab("Density") +
      theme_classic() +
      facet_wrap(~polarity_all) +
      guides(color = F) +
      theme(panel.grid.major = element_line(size = 0.5))
    # pl.density
  })
  names(pl.density.pkAnn.list.overall) <- paste0("colorBy_", colorBy)
  
  
  
  
  cat("\t--> PCA (samples)...\n")
  tmp    <- lipids.df
  tmp    <- c(list("+/-" = tmp), split(tmp, tmp$Polarity))
  tmp    <- lapply(names(tmp), function(tmp.pol){
    # tmp.pol <- names(tmp)[1]
    
    df.pca            <- tmp[[tmp.pol]]
    rownames(df.pca)  <- df.pca$Identification
    df.pca            <- dplyr::select(df.pca, -c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min."))
    # df.pca            <- log2(df.pca)
    df.pca            <- df.pca[unname(apply(df.pca, 1, function(zz) sum(!is.na(zz)) > 0.5*ncol(df.pca))), ]
    pca               <- FactoMineR::PCA(t(df.pca), graph = F)
    df.pca            <- pca$ind$coord
    df.pca            <- merge(df.pca, meta, by.x = "row.names", by.y = "SampleName", all = F)
    df.pca$varDim1    <-  round(pca$eig[1, "percentage of variance"], 2)
    df.pca$varDim2    <-  round(pca$eig[2, "percentage of variance"], 2)
    df.pca$polarity   <- tmp.pol
    df.pca
  })
  dfPlot            <- type.convert(as.data.frame(rbindlist(tmp)), as.is = T)
  dfPlot$sampleName <- dfPlot$Sample.ID
  
  dfPlot$xlab <- paste0("(", dfPlot$polarity, ")", dfPlot$varDim1, "%")
  dfPlot$ylab <- paste0("(", dfPlot$polarity, ")", dfPlot$varDim2, "%")
  
  dfPlot$Batch <- as.character(dfPlot$Batch)
  
  xlab <- paste0("Dim.1: ", paste(unique(dfPlot$xlab), collapse = "; "))
  ylab <- paste0("Dim.2:\n", paste(unique(dfPlot$ylab), collapse = "; "))
  
  colorBy     <- c("Strain", "Batch", "Diet")
  set.seed(123)
  pl.pca.list.samples <- lapply(colorBy, function(cby){
    # cby <- "Diet"
    
    if(cby == "Diet"){
      colors <- c("CD" = "#6799e4", "HFD" = "#b0222b")
    } else{
      colors <- colors()
      colors <- colors[!grepl("grey|gray|white", colors)]
      set.seed(123)
      colors <- colors[sample(1:length(colors), length(unique(dfPlot[[cby]])))]
    }
    out.list <- lapply(c("by_polarity", "all_polarities"), function(zz){
      if(zz == "by_polarity"){
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+", "-"), ]
      } else{
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+/-"), ]
      }
      ggplot(dfPlot.zz, aes_string(x = "Dim.1", y = "Dim.2", color = cby)) +
        geom_point(size = 0.9) +
        scale_color_manual(values = colors) +
        xlab(xlab) +
        ylab(ylab) +
        theme_classic() +
        facet_wrap(~polarity) +
        theme(legend.spacing.y = unit(0.1, 'cm'),
              legend.text      = element_text(size = 7),
              legend.key.size = unit(0.1, "cm"),
              legend.margin = margin(0.02,0,0,0, unit="cm")) +
        labs(color = cby, shape = "Diet")
    })
    names(out.list) <- c("by_polarity", "all_polarities")
    out.list
  })
  names(pl.pca.list.samples) <- paste0("colorBy_", colorBy)
  
  # ggsave(paste0(saveFolder, "/pca_plot_colorBy_Strain.pdf"), plot = pl.pca.list$colorBy_Strain, width = 6, height = 2.4, useDingbats = F)
  
 
  
  cat("\t--> MDS (samples)...\n")
  tmp    <- lipids.df
  tmp    <- c(list("+/-" = tmp), split(tmp, tmp$Polarity))
  tmp    <- lapply(names(tmp), function(tmp.pol){
    # tmp.pol <- names(tmp)[1]
    
    df.mds            <- tmp[[tmp.pol]]
    rownames(df.mds)  <- df.mds$Identification
    df.mds            <- dplyr::select(df.mds, -c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min."))
    df.mds            <- df.mds[unname(apply(df.mds, 1, function(zz) sum(!is.na(zz)) > 0.5*ncol(df.mds))), ]
    # df.mds            <- log2(df.mds)
    df.mds            <- as.data.frame(cmdscale(dist(t(df.mds)), k = 2))
    df.mds            <- merge(df.mds, meta, by.x = "row.names", by.y = "SampleName", all = F)
    df.mds$polarity   <- tmp.pol
    df.mds
  })
  dfPlot <- type.convert(as.data.frame(rbindlist(tmp)), as.is = T)
  dfPlot$sampleName <- dfPlot$Sample.ID
  dfPlot$Batch      <- as.character(dfPlot$Batch)
  
  colorBy <- c("Strain", "Batch", "Diet")
  set.seed(123)
  pl.mds.list.samples <- lapply(colorBy, function(cby){
    # cby <- "Diet"
    
    if(cby == "Diet"){
      colors <- c("CD" = "#6799e4", "HFD" = "#b0222b")
    } else{
      colors <- colors()
      set.seed(123)
      colors <- colors[!grepl("grey|gray|white", colors)]
      colors <- colors[sample(1:length(colors), length(unique(dfPlot[[cby]])))]
    }
    out.list <- lapply(c("by_polarity", "all_polarities"), function(zz){
      if(zz == "by_polarity"){
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+", "-"), ]
      } else{
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+/-"), ]
      }
      ggplot(dfPlot.zz, aes_string(x = "V1", y = "V2", color = cby)) +
        geom_point(size = 0.9) +
        scale_color_manual(values = colors) +
        xlab("Dim.1") +
        ylab("Dim.2") +
        theme_classic() +
        facet_wrap(~polarity) +
        theme(legend.spacing.y = unit(0.1, 'cm'),
              legend.text      = element_text(size = 7),
              legend.key.size = unit(0.1, "cm"),
              legend.margin = margin(0.02,0,0,0, unit="cm")) +
        labs(color = cby, shape = "Diet")
    })
    names(out.list) <- c("by_polarity", "all_polarities")
    out.list
  })
  names(pl.mds.list.samples) <- paste0("colorBy_", colorBy)
  
  
  cat("\t--> UMAP (samples)...\n")
  tmp    <- lipids.df
  tmp    <- c(list("+/-" = tmp), split(tmp, tmp$Polarity))
  tmp    <- lapply(names(tmp), function(tmp.pol){
    # tmp.pol <- names(tmp)[1]

    df.umap            <- tmp[[tmp.pol]]
    rownames(df.umap)  <- df.umap$Identification
    df.umap            <- dplyr::select(df.umap, -c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min."))
    # df.umap            <- log2(df.umap)
    df.umap            <- df.umap[unname(apply(df.umap, 1, function(zz) all(!is.na(zz)))), ]
    df.umap            <- as.data.frame(umap::umap(t(df.umap))$layout)
    df.umap            <- merge(df.umap, meta, by.x = "row.names", by.y = "SampleName", all = F)
    df.umap$polarity   <- tmp.pol
    df.umap
  })
  dfPlot            <- type.convert(as.data.frame(rbindlist(tmp)), as.is = T)
  dfPlot$sampleName <- dfPlot$Sample.ID
  dfPlot$Batch      <- as.character(dfPlot$Batch)
  
  colorBy <- c("Strain", "Batch", "Diet")
  set.seed(123)
  pl.umap.list.samples <- lapply(colorBy, function(cby){
    # cby <- "Diet"
    
    if(cby == "Diet"){
      colors <- c("CD" = "#6799e4", "HFD" = "#b0222b")
    } else{
      colors <- colors()
      set.seed(123)
      colors <- colors[!grepl("grey|gray|white", colors)]
      colors <- colors[sample(1:length(colors), length(unique(dfPlot[[cby]])))]
    }
    out.list <- lapply(c("by_polarity", "all_polarities"), function(zz){
      if(zz == "by_polarity"){
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+", "-"), ]
      } else{
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+/-"), ]
      }
      ggplot(dfPlot.zz, aes_string(x = "V1", y = "V2", color = cby)) +
        geom_point(size = 0.9) +
        scale_color_manual(values = colors) +
        xlab("Dim.1") +
        ylab("Dim.2") +
        theme_classic() +
        facet_wrap(~polarity) +
        theme(legend.spacing.y = unit(0.1, 'cm'),
              legend.text      = element_text(size = 7),
              legend.key.size = unit(0.1, "cm"),
              legend.margin = margin(0.02,0,0,0, unit="cm")) +
        labs(color = cby, shape = "Diet")
    })
    names(out.list) <- c("by_polarity", "all_polarities")
    out.list
  })
  names(pl.umap.list.samples) <- paste0("colorBy_", colorBy)
  
  
  cat("\t--> PCA (lipids)...\n")
  tmp                  <- lipids.df
  tmp$saturation_ratio <- tmp$cLength / tmp$nbDoubleBonds
  tmp                  <- c(list("+/-" = tmp), split(tmp, tmp$Polarity))
  tmp                  <- lapply(names(tmp), function(tmp.pol){
    # tmp.pol <- names(tmp)[1]
    
    df.pca            <- tmp[[tmp.pol]]
    rownames(df.pca)  <- df.pca$Identification
    df.pca            <- dplyr::select(df.pca, -c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min.", "saturation_ratio"))
    # df.pca            <- log2(df.pca)
    df.pca            <- df.pca[unname(apply(df.pca, 1, function(zz) sum(!is.na(zz)) > 0.5*ncol(df.pca))), ]
    pca               <- FactoMineR::PCA(df.pca, graph = F)
    df.pca            <- pca$ind$coord
    df.pca            <- merge(df.pca, dplyr::select(tmp[[tmp.pol]], c("Identification", "Lipid.Class", "saturation_ratio", "cLength", "nbDoubleBonds")), by.x = "row.names", by.y = "Identification", all = F)
    df.pca$varDim1    <-  round(pca$eig[1, "percentage of variance"], 2)
    df.pca$varDim2    <-  round(pca$eig[2, "percentage of variance"], 2)
    df.pca$polarity   <- tmp.pol
    df.pca
  })
  dfPlot <- type.convert(as.data.frame(rbindlist(tmp)), as.is = T)

  dfPlot$xlab <- paste0("(", dfPlot$polarity, ")", dfPlot$varDim1, "%")
  dfPlot$ylab <- paste0("(", dfPlot$polarity, ")", dfPlot$varDim2, "%")
  
  xlab <- paste0("Dim.1: ", paste(unique(dfPlot$xlab), collapse = "; "))
  ylab <- paste0("Dim.2:\n", paste(unique(dfPlot$ylab), collapse = "; "))
  
  colorBy     <- c("Lipid.Class", "cLength", "nbDoubleBonds", "saturation_ratio")
  set.seed(123)
  pl.pca.list.lipids <- lapply(colorBy, function(cby){
    # cby <- "Lipid.Class"
    # cby <- "saturation_ratio"
    
    if(cby %in% c("cLength", "nbDoubleBonds", "saturation_ratio")){
      colors_layer <- scale_color_gradientn(colors = pals::viridis(200))
    } else{
      colors       <- colors()
      colors       <- colors[!grepl("grey|gray|white", colors)]
      set.seed(123)
      colors       <- colors[sample(1:length(colors), length(unique(dfPlot[[cby]])))]
      colors_layer <- scale_color_manual(values = colors)
    }
    out.list <- lapply(c("by_polarity", "all_polarities"), function(zz){
      if(zz == "by_polarity"){
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+", "-"), ]
      } else{
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+/-"), ]
      }
      ggplot(dfPlot.zz, aes_string(x = "Dim.1", y = "Dim.2", color = cby)) +
        geom_point(size = 0.9) +
        colors_layer +
        xlab(xlab) +
        ylab(ylab) +
        theme_classic() +
        facet_wrap(~polarity) +
        theme(legend.spacing.y = unit(0.1, 'cm'),
              legend.text      = element_text(size = 7),
              legend.key.size = unit(0.1, "cm"),
              legend.margin = margin(0.02,0,0,0, unit="cm")) +
        labs(color = cby, shape = "Diet")
    })
    names(out.list) <- c("by_polarity", "all_polarities")
    out.list
  })
  names(pl.pca.list.lipids) <- paste0("colorBy_", colorBy)
  

  
  
  cat("\t--> MDS (samples)...\n")
  tmp                  <- lipids.df
  tmp$saturation_ratio <- tmp$cLength / tmp$nbDoubleBonds
  tmp                  <- c(list("+/-" = tmp), split(tmp, tmp$Polarity))
  tmp                  <- lapply(names(tmp), function(tmp.pol){
    # tmp.pol <- names(tmp)[1]
    
    df.mds            <- tmp[[tmp.pol]]
    rownames(df.mds)  <- df.mds$Identification
    df.mds            <- dplyr::select(df.mds, -c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min.", "saturation_ratio"))
    df.mds            <- df.mds[unname(apply(df.mds, 1, function(zz) sum(!is.na(zz)) > 0.5*ncol(df.mds))), ]
    df.mds            <- as.data.frame(cmdscale(dist(df.mds), k = 2))
    df.mds            <- merge(df.mds, dplyr::select(tmp[[tmp.pol]], c("Identification", "Lipid.Class", "saturation_ratio", "cLength", "nbDoubleBonds")), by.x = "row.names", by.y = "Identification", all = F)
    df.mds$polarity   <- tmp.pol
    df.mds
  })
  dfPlot <- type.convert(as.data.frame(rbindlist(tmp)), as.is = T)
  
  dfPlot$xlab <- paste0("(", dfPlot$polarity, ")", dfPlot$varDim1, "%")
  dfPlot$ylab <- paste0("(", dfPlot$polarity, ")", dfPlot$varDim2, "%")
  
  xlab <- paste0("Dim.1: ", paste(unique(dfPlot$xlab), collapse = "; "))
  ylab <- paste0("Dim.2:\n", paste(unique(dfPlot$ylab), collapse = "; "))
  
  colorBy     <- c("Lipid.Class", "cLength", "nbDoubleBonds", "saturation_ratio")
  set.seed(123)
  pl.mds.list.lipids <- lapply(colorBy, function(cby){
    # cby <- "Lipid.Class"
    # cby <- "cLength"
    
    if(cby %in% c("cLength", "nbDoubleBonds", "saturation_ratio")){
      colors_layer  <- scale_color_gradientn(colors = pals::viridis(200))
    } else{
      colors       <- colors()
      colors       <- colors[!grepl("grey|gray|white", colors)]
      set.seed(123)
      colors       <- colors[sample(1:length(colors), length(unique(dfPlot[[cby]])))]
      colors_layer <- scale_color_manual(values = colors)
    }
    out.list <- lapply(c("by_polarity", "all_polarities"), function(zz){
      # zz <- "by_polarity"
      if(zz == "by_polarity"){
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+", "-"), ]
      } else{
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+/-"), ]
      }
      ggplot(dfPlot.zz, aes_string(x = "V1", y = "V2", color = cby)) +
        geom_point(size = 0.9) +
        colors_layer +
        xlab("Dim.1") +
        ylab("Dim.2") +
        theme_classic() +
        facet_wrap(~polarity) +
        theme(legend.spacing.y = unit(0.1, 'cm'),
              legend.text      = element_text(size = 7),
              legend.key.size = unit(0.1, "cm"),
              legend.margin = margin(0.02,0,0,0, unit="cm")) +
        labs(color = cby, shape = "Diet")
    })
    names(out.list) <- c("by_polarity", "all_polarities")
    out.list
  })
  names(pl.mds.list.lipids) <- paste0("colorBy_", colorBy)
  
  
  
  
  cat("\t--> UMAP (lipids)...\n")
  tmp                  <- lipids.df
  tmp$saturation_ratio <- tmp$cLength / tmp$nbDoubleBonds
  tmp                  <- c(list("+/-" = tmp), split(tmp, tmp$Polarity))
  tmp                  <- lapply(names(tmp), function(tmp.pol){
    # tmp.pol <- names(tmp)[1]
    
    df.umap            <- tmp[[tmp.pol]]
    rownames(df.umap)  <- df.umap$Identification
    df.umap            <- dplyr::select(df.umap, -c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min.", "saturation_ratio"))
    df.umap            <- df.umap[unname(apply(df.umap, 1, function(zz) all(!is.na(zz)))), ]
    df.umap            <- as.data.frame(umap::umap(df.umap)$layout)
    df.umap            <- merge(df.umap, dplyr::select(tmp[[tmp.pol]], c("Identification", "Lipid.Class", "saturation_ratio", "cLength", "nbDoubleBonds")), by.x = "row.names", by.y = "Identification", all = F)
    df.umap$polarity   <- tmp.pol
    df.umap
  })
  dfPlot <- type.convert(as.data.frame(rbindlist(tmp)), as.is = T)

  colorBy     <- c("Lipid.Class", "cLength", "nbDoubleBonds", "saturation_ratio")
  set.seed(123)
  pl.umap.list.lipids <- lapply(colorBy, function(cby){
    # cby <- "Lipid.Class"
    # cby <- "cLength"
    
    if(cby %in% c("cLength", "nbDoubleBonds", "saturation_ratio")){
      colors_layer  <- scale_color_gradientn(colors = pals::viridis(200))
    } else{
      colors       <- colors()
      colors       <- colors[!grepl("grey|gray|white", colors)]
      set.seed(123)
      colors       <- colors[sample(1:length(colors), length(unique(dfPlot[[cby]])))]
      colors_layer <- scale_color_manual(values = colors)
    }
    out.list <- lapply(c("by_polarity", "all_polarities"), function(zz){
      # zz <- "by_polarity"
      if(zz == "by_polarity"){
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+", "-"), ]
      } else{
        dfPlot.zz <- dfPlot[dfPlot$polarity %in% c("+/-"), ]
      }
      ggplot(dfPlot.zz, aes_string(x = "V1", y = "V2", color = cby)) +
        geom_point(size = 0.9) +
        colors_layer +
        xlab("Dim.1") +
        ylab("Dim.2") +
        theme_classic() +
        facet_wrap(~polarity) +
        theme(legend.spacing.y = unit(0.1, 'cm'),
              legend.text      = element_text(size = 7),
              legend.key.size = unit(0.1, "cm"),
              legend.margin = margin(0.02,0,0,0, unit="cm")) +
        labs(color = cby, shape = "Diet")
    })
    names(out.list) <- c("by_polarity", "all_polarities")
    out.list
  })
  names(pl.umap.list.lipids) <- paste0("colorBy_", colorBy)
  
  
  
  cat("\t--> Sample - blank intensity differences...\n")
  if(any(grepl("Blank", colnames(data.df)))){
    
    lipids.df            <- data.df
    lipids.df            <- lipids.df[!lipids.df$isUnknown, ]
    lipids.df            <- dplyr::select(lipids.df, c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min.", colnames(lipids.df)[grepl("BDX|BXD|C57|DBA", colnames(lipids.df))]))
    lipids.df            <- lipids.df[, !grepl("PWT", colnames(lipids.df))]
    tmpLogi              <- grepl("BXD|C57B|DBA2|Blank", colnames(lipids.df))
    lipids.df[, tmpLogi] <- log2(lipids.df[, tmpLogi])
    meta                 <- meta.df
    dfPlot               <- lapply(colnames(lipids.df)[grepl("BXD|C57B|DBA2", colnames(lipids.df)) & !grepl("Blank", colnames(lipids.df))], function(zz){
      # zz <- colnames(lipids.df)[grepl("BXD|C57B|DBA2", colnames(lipids.df))][1]
      # zz <- "BXD27_CD_1"
      # print(zz)
      batch      <- meta$Batch[meta$SampleName == zz]
      blank_col  <- meta$SampleName[meta$Batch == batch & meta$is_blank_control_sample]
      if(length(blank_col) != 1){
        return()
      }
      blank.vec  <- lipids.df[[blank_col]]
      sample.vec <- lipids.df[[zz]]
      diff.vec   <- sample.vec - blank.vec
      data.frame(sample_id         = zz, 
                 batch             = batch, 
                 sample_log2       = sample.vec,
                 blank_log2        = blank.vec,
                 sample_blank_diff = diff.vec, stringsAsFactors = F)
    })
    dfPlot <- as.data.frame(rbindlist(dfPlot))
    
    dfPlot$batch_label <- paste0("batch ", dfPlot$batch)
    dfPlot$batch_label <- factor(dfPlot$batch_label, levels = unique(dfPlot$batch_label))
    
    
    pl.sample.blank.diff <- ggplot(dfPlot, aes(x = sample_id, y = sample_blank_diff)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_jitter(aes(color = sample_log2), size = 0.2, width = 0.3) +
      geom_boxplot(outlier.shape = NA, fill = NA) +
      scale_color_gradientn(colors = pals::brewer.ylorrd(200)) +
      ylab("log2(sample intensity) - log2(batch blank intensity)") +
      theme_classic() +
      facet_wrap(~batch, scales = "free_x") +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      labs(color = "Sample\nlog2 intensity")
    # pl.sample.blank.diff
    
    tmpLogi        <- grepl("BXD|C57B|DBA2", colnames(lipids.df))
    samples_median <- apply(lipids.df[, tmpLogi], 1, function(zz) log2(median(2^zz, na.rm = T)))
    tmpLogi        <- grepl("Blank", colnames(lipids.df))
    blank_median   <- apply(lipids.df[, tmpLogi], 1, function(zz) log2(median(2^zz, na.rm = T)))
    
    dfPlot <- data.frame(Identifier          = lipids.df$Identifier,
                         log2_samples_median = samples_median,
                         log2_blank_median   = blank_median,
                         stringsAsFactors = F)
    dfPlot$med_blank_diff <- dfPlot$log2_samples_median - dfPlot$log2_blank_median
    
    # threshold of a two-fold difference between the sample and blank intensities used for peak filtering
    threshold <- log2(2)
    
    dfPlot$color <- ifelse(dfPlot$med_blank_diff > threshold, "high-quality peaks", "others")
    
    pl.median.sample.blank.diff <- ggplot(dfPlot, aes(x = log2_samples_median, y = med_blank_diff)) +
      geom_hline(yintercept = c(0, threshold), linetype = "dashed") +
      geom_point(aes(color = color)) +
      scale_color_manual(values = c("high-quality peaks" = "#b0222b", "others" = "#BEBEBE")) +
      xlab("log2 sample intensity") + ylab("log2 sample - log2 blank") +
      theme_classic() +
      labs(color = "Peak type")
    # pl.median.sample.blank.diff
    
    
    # dfPlot <- lipids.df[, !grepl("Blank", colnames(lipids.df))]
    # dfPlot <- reshape2::melt(dfPlot, id.vars = c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min."))
    # 
    # ggplot(dfPlot[1:300000, ], aes(x = Retention.Time..min., color = variable, group = variable)) +
    #   geom_density() +
    #   theme_classic() +
    #   theme(legend.position = "none")
    
  } else{
    
    pl.sample.blank.diff        <- NULL
    pl.median.sample.blank.diff <- NULL
    
  }
  
  
  cat("\t--> Sample - pooled controls intensity differences...\n")
  if(any(grepl("PWT", colnames(data.df)))){
    
    lipids.df            <- data.df
    lipids.df            <- lipids.df[!lipids.df$isUnknown, ]
    lipids.df            <- dplyr::select(lipids.df, c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min.", colnames(lipids.df)[grepl("BDX|BXD|C57|DBA", colnames(lipids.df))]))
    lipids.df            <- lipids.df[, !grepl("Blank", colnames(lipids.df))]
    tmpLogi              <- grepl("BXD|C57B|DBA2|PWT", colnames(lipids.df))
    lipids.df[, tmpLogi] <- log2(lipids.df[, tmpLogi])
    meta                 <- meta.df
    dfPlot               <- lapply(colnames(lipids.df)[grepl("BXD|C57B|DBA2", colnames(lipids.df)) & !grepl("Blank|PWT", colnames(lipids.df))], function(zz){
      # zz <- colnames(lipids.df)[grepl("BXD|C57B|DBA2", colnames(lipids.df))][1]
      # zz <- "BXD27_CD_1"
      # print(zz)
      batch       <- meta$Batch[meta$SampleName == zz]
      pooled_col  <- meta$SampleName[meta$Batch == batch & meta$is_pooled_control_sample]
      if(length(pooled_col) == 0){
        return()
      }
      pooled.vec <- unname(apply(dplyr::select(lipids.df, dplyr::all_of(pooled_col)), 1, function(zz) log2(mean(2^zz, na.rm = T))))
      sample.vec <- lipids.df[[zz]]
      diff.vec   <- sample.vec - pooled.vec
      data.frame(sample_id          = zz, 
                 batch              = batch, 
                 sample_log2        = sample.vec,
                 avg_pooled_log2    = pooled.vec,
                 sample_pooled_diff = diff.vec, stringsAsFactors = F)
    })
    dfPlot <- as.data.frame(rbindlist(dfPlot))
    
    dfPlot$batch_label <- paste0("batch ", dfPlot$batch)
    dfPlot$batch_label <- factor(dfPlot$batch_label, levels = unique(dfPlot$batch_label))
    
    
    pl.sample.pooled.diff <- ggplot(dfPlot, aes(x = sample_id, y = sample_pooled_diff)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_jitter(aes(color = sample_log2), size = 0.2, width = 0.3) +
      geom_boxplot(outlier.shape = NA, fill = NA) +
      scale_color_gradientn(colors = pals::brewer.ylorrd(200)) +
      ylab("log2(sample intensity) - log2(avg pooled control intensity)") +
      theme_classic() +
      facet_wrap(~batch, scales = "free_x") +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      labs(color = "Sample\nlog2 intensity")
    # pl.sample.pooled.diff
    
    
    tmpLogi        <- grepl("BXD|C57B|DBA2", colnames(lipids.df))
    samples_median <- apply(lipids.df[, tmpLogi], 1, function(zz) log2(median(2^zz, na.rm = T)))
    tmpLogi        <- grepl("PWT", colnames(lipids.df))
    pooled_median  <- apply(lipids.df[, tmpLogi], 1, function(zz) log2(median(2^zz, na.rm = T)))
    
    dfPlot <- data.frame(Identifier           = lipids.df$Identifier,
                         lipid_class          = lipids.df$Lipid.Class,
                         log2_samples_median  = samples_median,
                         log2_pooled_median   = pooled_median,
                         stringsAsFactors = F)
    dfPlot$med_pooled_diff <- dfPlot$log2_samples_median - dfPlot$log2_pooled_median
    
    
    pl.median.sample.pooled.diff <- ggplot(dfPlot, aes(x = log2_samples_median, y = med_pooled_diff)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point() +
      xlab("log2 sample intensity") + ylab("log2 sample - log2 pooled") +
      theme_classic() +
      labs(color = "Peak type")
    # pl.median.sample.pooled.diff
    
  } else{
    
    pl.sample.pooled.diff        <- NULL
    pl.median.sample.pooled.diff <- NULL
    
  }
  
  
  
  cat("\t--> RT vs m/z...\n")
  lipids.df            <- data.df
  # lipids.df            <- lipids.df[!lipids.df$isUnknown, ]
  lipids.df            <- dplyr::select(lipids.df, c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min.", "Quant.Ion", colnames(lipids.df)[grepl("BDX|BXD|C57|DBA", colnames(lipids.df))]))
  lipids.df            <- lipids.df[, !grepl("PWT|Blank", colnames(lipids.df))]
  tmpLogi              <- grepl("BXD|C57B|DBA2", colnames(lipids.df))
  lipids.df[, tmpLogi] <- log2(lipids.df[, tmpLogi])
  meta                 <- meta.df
  
  dfPlot          <- dplyr::select(lipids.df, c("Identifier", "Retention.Time..min.", "Quant.Ion", "Lipid.Class", "isUnknown"))

  pl.rt.vs.mz <- ggplot(dfPlot[!dfPlot$isUnknown, ], aes(x = Retention.Time..min., y = Quant.Ion, color = Lipid.Class)) +
    geom_point() +
    xlab("RT [min]") + ylab("m/z") +
    theme_classic() +
    labs(color = "Lipid class")
  # pl.rt.vs.mz
  
  pl.rt.vs.mz.withUnknown <- ggplot(dfPlot, aes(x = Retention.Time..min., y = Quant.Ion, color = Lipid.Class)) +
    geom_point() +
    xlab("RT [min]") + ylab("m/z") +
    theme_classic() +
    labs(color = "Lipid class")
  # pl.rt.vs.mz.withUnknown
  
  
  cat("\t--> Pooled controls distributions...\n")
  if(any(grepl("PWT", colnames(data.df)))){
    
    dfPlot     <- data.df
    colsSelect <- c("Identifier", "Retention.Time..min.", "Quant.Ion", "Polarity", "Area..max.", 
                    "Identification", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", 
                    "Features.Found", colnames(dfPlot)[grepl("BXDPWT", colnames(dfPlot))])
    dfPlot     <- dplyr::select(dfPlot, dplyr::all_of(colsSelect))
    dfPlot     <- reshape2::melt(dfPlot, id.vars = c("Identifier", "Retention.Time..min.", "Quant.Ion", "Polarity", "Area..max.", 
                                                     "Identification", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", 
                                                     "Features.Found"))
    dfPlot$value      <- log2(dfPlot$value)
    dfPlot$batch_nb   <- paste0("batch ", gsub("_.*", "", gsub("BXDPWT_", "", dfPlot$variable)))
    dfPlot$control_nb <- paste0("Pooled control ", gsub(".*_", "", dfPlot$variable))
    dfPlot$batch_nb   <- factor(dfPlot$batch_nb, levels = unique(dfPlot$batch_nb))
    
    pl.pooled_intensity.distrib.all <- ggplot(dfPlot, aes(x = value, color = control_nb)) +
      geom_density() +
      scale_color_manual(values = pals::jet(length(unique(dfPlot$control_nb)))) +
      xlab("log2(raw intensity) - all features") + ylab("Density") +
      theme_classic() +
      facet_wrap(~batch_nb, scales = "free") +
      labs(color = "Pooled control nb.")
    # pl.pooled_intensity.distrib.all
    
    pl.pooled_intensity.distrib.noUnknown <- ggplot(dfPlot[!dfPlot$isUnknown, ], aes(x = value, color = control_nb)) +
      geom_density() +
      scale_color_manual(values = pals::jet(length(unique(dfPlot$control_nb)))) +
      xlab("log2(raw intensity) - known features") + ylab("Density") +
      theme_classic() +
      facet_wrap(~batch_nb, scales = "free") +
      labs(color = "Pooled control nb.")
    # pl.pooled_intensity.distrib.noUnknown
    
    pl.pooled_intensity.distrib.unknown <- ggplot(dfPlot[dfPlot$isUnknown, ], aes(x = value, color = control_nb)) +
      geom_density() +
      scale_color_manual(values = pals::jet(length(unique(dfPlot$control_nb)))) +
      xlab("log2(raw intensity) - unknown features") + ylab("Density") +
      theme_classic() +
      facet_wrap(~batch_nb, scales = "free") +
      labs(color = "Pooled control nb.")
    # pl.pooled_intensity.distrib.unknown
    
  } else{
    
    pl.pooled_intensity.distrib.all       <- NULL
    pl.pooled_intensity.distrib.noUnknown <- NULL
    pl.pooled_intensity.distrib.unknown   <- NULL
    
  }
  
  
  cat("\t--> Blank controls distributions...\n")
  if(any(grepl("Blank", colnames(data.df)))){
    
    dfPlot     <- data.df
    colsSelect <- c("Identifier", "Retention.Time..min.", "Quant.Ion", "Polarity", "Area..max.", 
                    "Identification", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", 
                    "Features.Found", colnames(dfPlot)[grepl("BXDBlank", colnames(dfPlot))])
    dfPlot     <- dplyr::select(dfPlot, dplyr::all_of(colsSelect))
    dfPlot     <- reshape2::melt(dfPlot, id.vars = c("Identifier", "Retention.Time..min.", "Quant.Ion", "Polarity", "Area..max.", 
                                                     "Identification", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", 
                                                     "Features.Found"))
    dfPlot$value      <- log2(dfPlot$value)
    dfPlot$batch_nb   <- paste0("batch ", gsub("_.*", "", gsub("BXDBlank_", "", dfPlot$variable)))
    dfPlot$batch_nb   <- factor(dfPlot$batch_nb, levels = unique(dfPlot$batch_nb))
    
    pl.blank_intensity.distrib <- ggplot(dfPlot, aes(x = value, color = batch_nb)) +
      geom_density() +
      scale_color_manual(values = pals::jet(length(unique(dfPlot$batch_nb)))) +
      xlab("log2(raw intensity) - all features") + ylab("Density") +
      theme_classic() +
      labs(color = "Background control nb.")
    # pl.blank_intensity.distrib
    
  } else{
    
    pl.blank_intensity.distrib <- NULL
    
  }
  
  
  
  
  cat("\t--> Sample weight variation...\n")
  dfPlot <- meta.df
  dfPlot <- plyr::ddply(dfPlot[!grepl("PWT|Blank", dfPlot$SampleName), ], c("Batch"), function(x){
    x$weight_perc_first  <- (x$Weight / x$Weight[1]) * 100
    x$is_first_sample    <- F
    x$is_first_sample[1] <- T
    x$sample_nb <- 1:nrow(x)
    x
  })
  dfPlot$weight_perc_variation <- dfPlot$weight_perc_first - 100
  dfPlot                       <- dfPlot[!dfPlot$is_first_sample, ]
  dfPlot$batch_nb              <- paste0("batch ", dfPlot$Batch)
  dfPlot$perc_lab              <- paste0(ifelse(dfPlot$weight_perc_variation > 0, "+", ""), round(dfPlot$weight_perc_variation, 1), "%")
  dfPlot$lab_pos               <- dfPlot$weight_perc_variation + sign(dfPlot$weight_perc_variation) * 1
  dfPlot$lab_pos[dfPlot$weight_perc_variation == 0] <- 1
  
  pl.weight_variation.batch <- ggplot(dfPlot, aes(x = sample_nb, y = weight_perc_variation)) +
    geom_hline(yintercept = 0) +
    geom_bar(stat = "Identity") +
    geom_text(aes(label = perc_lab, y = lab_pos), size = 1.8) +
    xlab("Sample nb.") + ylab("Sample weight variation\ncompared to first batch sample") +
    facet_wrap(~batch_nb, scales = "free") +
    theme_classic()
  # pl.weight_variation.batch
  
  
  dfPlot <- meta.df
  dfPlot <- dfPlot[!grepl("PWT|Blank", dfPlot$SampleName), ]
  dfPlot$weight_perc_first  <- (dfPlot$Weight / dfPlot$Weight[1]) * 100
  dfPlot <- dfPlot[-1, ]
  dfPlot$sample_nb <- (1:nrow(dfPlot)) + 1
  dfPlot$weight_perc_variation <- 100 - dfPlot$weight_perc_first
  dfPlot$perc_lab              <- paste0(ifelse(dfPlot$weight_perc_variation > 0, "+", ""), round(dfPlot$weight_perc_variation, 1), "%")
  dfPlot$lab_pos               <- dfPlot$weight_perc_variation + sign(dfPlot$weight_perc_variation) * 1
  dfPlot$lab_pos[dfPlot$weight_perc_variation == 0] <- 1
  
  pl.weight_variation.all <- ggplot(dfPlot, aes(x = sample_nb, y = weight_perc_variation)) +
    geom_hline(yintercept = 0) +
    geom_bar(stat = "Identity") +
    geom_text(aes(label = perc_lab, y = lab_pos), size = 1.8) +
    xlab("Sample nb.") + ylab("Sample weight variation\ncompared to first sample") +
    theme_classic()
  pl.weight_variation.all
  
  
  
  
  
  
  
  
  cat("\t--> PCA - label experimental and control samples...\n")
  if(any(grepl("PWT", colnames(data.df)))){
    
    lipids.df            <- data.df
    lipids.df            <- lipids.df[!lipids.df$isUnknown, ]
    lipids.df            <- dplyr::select(lipids.df, c("Identification", "Identifier", "Polarity", "Lipid.Class", "isUnknown", "cLength", "nbDoubleBonds", "Retention.Time..min.", colnames(lipids.df)[grepl("BDX|BXD|C57|DBA", colnames(lipids.df))]))
    lipids.df            <- lipids.df[, !grepl("Blank", colnames(lipids.df))]
    tmpLogi              <- grepl("BXD|C57B|DBA2|PWT", colnames(lipids.df))
    lipids.df[, tmpLogi] <- log2(lipids.df[, tmpLogi])
    rownames(lipids.df)  <- lipids.df$Identifier
    
    df.pca            <- lipids.df[, tmpLogi]
    df.pca            <- df.pca[unname(apply(df.pca, 1, function(zz) sum(!is.na(zz)) > 0.5*ncol(df.pca))), ]
    pca               <- FactoMineR::PCA(t(df.pca), graph = F)
    df.pca            <- as.data.frame(pca$ind$coord)
    df.pca$label      <- ifelse(grepl("HFD", rownames(df.pca)), "HFD", "CD")
    df.pca$label[grepl("PWT", rownames(df.pca))] <- "Pooled control"
    df.pca$batch      <- plyr::mapvalues(rownames(df.pca), from = meta.df$SampleName, to = meta.df$Batch, warn_missing = F)
    df.pca$label_2    <- df.pca$label
    tmpLogi           <- df.pca$label_2 == "Pooled control"
    df.pca$label_2[tmpLogi] <- paste0("Pooled control (batch ", df.pca$batch[tmpLogi], ")")
    df.pca$label_2    <- factor(df.pca$label_2, levels = unique(df.pca$label_2))
    
    pl.pca.with_pooled_controls <- ggplot(df.pca, aes(x = Dim.1, y = Dim.2, color = label)) +
      geom_point(size = 0.9) +
      scale_color_manual(values = c("CD" = "#6799e4", "HFD" = "#b0222b", "Pooled control" = "#C9CC3F")) +
      xlab(paste0("Dim.1 (", round(pca$eig[1, "percentage of variance"], 2), "%)")) +
      ylab(paste0("Dim.2 (", round(pca$eig[2, "percentage of variance"], 2), "%)")) +
      theme_classic() +
      theme(legend.spacing.y = unit(0.5, 'cm'),
            legend.text      = element_text(size = 7),
            legend.key.size = unit(0.5, "cm"),
            legend.margin = margin(0.02,0,0,0, unit="cm")) +
      labs(color = "Group")
    # pl.pca.with_pooled_controls
    
    pooled_ctrls_labels <- as.character(unique(df.pca$label_2[grepl("Pooled", df.pca$label_2)]))
    colors_vec          <- c("CD" = "#6799e4", "HFD" = "#b0222b", pals::jet(length(pooled_ctrls_labels)))
    names(colors_vec)[3:length(colors_vec)] <- pooled_ctrls_labels
    
    pl.pca.with_pooled_controls.batches <- ggplot(df.pca, aes(x = Dim.1, y = Dim.2, color = label_2, shape = label)) +
      geom_point(size = 0.9) +
      scale_color_manual(values = colors_vec) +
      scale_shape_manual(values = c("CD" = 19, "HFD" = 19, "Pooled control" = 1)) +
      xlab(paste0("Dim.1 (", round(pca$eig[1, "percentage of variance"], 2), "%)")) +
      ylab(paste0("Dim.2 (", round(pca$eig[2, "percentage of variance"], 2), "%)")) +
      theme_classic() +
      theme(legend.spacing.y = unit(0.5, 'cm'),
            legend.text      = element_text(size = 7),
            legend.key.size = unit(0.5, "cm"),
            legend.margin = margin(0.02,0,0,0, unit="cm")) +
      labs(color = "Group", shape = "Type")
    # pl.pca.with_pooled_controls.batches
    
    
  } else{
    
    pl.pca.with_pooled_controls         <- NULL
    pl.pca.with_pooled_controls.batches <- NULL
    
  }
 
  
  
  
  cat("\nDone!\n\n")
  
  list(pl.scatter.duplicatedLipids.retTime               = pl.scatter.duplicatedLipids.retTime,
       pl.scatter.duplLipids.range                       = pl.scatter.duplLipids.range,
       pl.scatter.duplLipids.range.ordered               = pl.scatter.duplLipids.range.ordered,
       pl.nbSamples.byBatch                              = pl.nbSamples.byBatch,
       pl.boxlot.byPolarization                          = pl.boxlot.byPolarization,
       pl.boxlot.overall                                 = pl.boxlot.overall,
       pl.dotplot.cLength.vs.saturation.byPolarization   = pl.dotplot.cLength.vs.saturation.byPolarization,
       pl.dotplot.cLength.vs.saturation.overall          = pl.dotplot.cLength.vs.saturation.overall,
       pl.dotplot.cLength.vs.saturation.byClass          = pl.dotplot.cLength.vs.saturation.byClass,
       pl.histo.saturationRatio.byPolarization           = pl.histo.saturationRatio.byPolarization,
       pl.histo.saturationRatio.overall                  = pl.histo.saturationRatio.overall,
       pl.histo.nbLipidByClass.byPolarization            = pl.histo.nbLipidByClass.byPolarization,
       pl.histo.nbLipidByClass.overall                   = pl.histo.nbLipidByClass.overall,
       pl.density.list.byPolarization                    = pl.density.list.byPolarization,
       pl.density.pkAnn.list.byPolarization              = pl.density.pkAnn.list.byPolarization,
       pl.density.list.overall                           = pl.density.list.overall,
       pl.density.pkAnn.list.overall                     = pl.density.pkAnn.list.overall,
       pl.pca.list.samples                               = pl.pca.list.samples,
       pl.mds.list.samples                               = pl.mds.list.samples,
       pl.umap.list.samples                              = pl.umap.list.samples,
       pl.pca.list.lipids                                = pl.pca.list.lipids,
       pl.mds.list.lipids                                = pl.mds.list.lipids,
       pl.umap.list.lipids                               = pl.umap.list.lipids,
       pl.sample.blank.diff                              = pl.sample.blank.diff,
       pl.median.sample.blank.diff                       = pl.median.sample.blank.diff,
       pl.sample.pooled.diff                             = pl.sample.pooled.diff,
       pl.median.sample.pooled.diff                      = pl.median.sample.pooled.diff,
       pl.rt.vs.mz                                       = pl.rt.vs.mz,
       pl.rt.vs.mz.withUnknown                           = pl.rt.vs.mz.withUnknown,
       pl.pooled_intensity.distrib.all                   = pl.pooled_intensity.distrib.all,
       pl.pooled_intensity.distrib.noUnknown             = pl.pooled_intensity.distrib.noUnknown,
       pl.pooled_intensity.distrib.unknown               = pl.pooled_intensity.distrib.unknown,
       pl.blank_intensity.distrib                        = pl.blank_intensity.distrib,
       pl.weight_variation.batch                         = pl.weight_variation.batch,
       pl.weight_variation.all                           = pl.weight_variation.all,
       pl.pca.with_pooled_controls                       = pl.pca.with_pooled_controls,
       pl.pca.with_pooled_controls.batches               = pl.pca.with_pooled_controls.batches)
  
  
  
}


