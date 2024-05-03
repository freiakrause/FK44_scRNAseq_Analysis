#This is the script to perform QC on FK44.1 scRNAseq COUNT data provided by BSF

############################ Install Packages #################################
#install.packages("renv")
#renv::init()    #creates project library,  lockfile and Rprofile
#renv::snapshot() # updates lockfile abou current-used packages-> can be shared and reproduced by others when restore() is used
#renv::restore() # installs the exact same versions of packages determined in lockfile
#renv::update()) #
#renv::history() #
#renv::install("usethis")
#usethis::create_github_token()
#renv::install("gitcreds")
#gitcreds::gitcreds_set()
#renv::install("remotes")
#renv::install("ggplot2", "dplyr" ,"RColorBrewer") # package "SingleR" required by tutorial but not found when trying isntallation
#renv::install("igraph") # für das scheiß igraph (benötigt für seurat) nach tausend jahren troublehsooting gefunden: ich brauche: sudo apt install build-essential gfortran UND sudo apt install build-essential gfortran, dann gehts
#renv::install("Seurat")
#renv::install("SoupX")
#if (!require("BiocManager", quietly = TRUE))
#  renv::install("BiocManager")
#renv::install("gprofiler2")
#BiocManager::install(version = "3.18")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install('glmGamPoi')
#BiocManager::install('decontX')
#remotes::install_github("Moonerss/scrubletR")
#remotes::install_github('immunogenomics/presto')
#remotes::install_github("cysouw/qlcMatrix")
#remotes::install_github("Nawijn-Group-Bioinformatics/FastCAR")
############################# Load Libraries ###################################

library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(scrubletR)
library(gprofiler2)
#library(decontX)
library(SoupX)
library(celda)
library(scater)
library(stringr)
library(Matrix)
library(gridExtra)
library(qlcMatrix)
library(FastCAR)
library(pheatmap)
set.seed(42)

# # ################### Load Output from CellRanger and generate decontaminated Counts with DeconX ############
# # #deconX= ambient RNA is floating around and got incorporated into GEMs and now cells which should not expres certain RNA seem to express it. 
# # #e.g all cell types show Saa upregulation or whatever. but this is from dying heps that spread their RNA everywhere
# # #So deconX kinda filters 
# NPC_87 <- readRDS("./01_tidy_data/iAL87_FastCar.rds") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
# counts <- GetAssayData(object = NPC_87, layer = "counts")
# sce_87 <- SingleCellExperiment(list(counts = counts))
# sce_87 <- decontX(sce_87)
# NPC_87[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce_87))
# 
# NPC_88 <- readRDS("./01_tidy_data/iAL88_FastCar.rds") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
# counts <- GetAssayData(object = NPC_88, layer = "counts")
# sce_88 <- SingleCellExperiment(list(counts = counts))
# sce_88 <- decontX(sce_88)
# NPC_88[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce_88))
# 
# NPC_91 <- readRDS("./01_tidy_data/iAL91_FastCar.rds") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
# counts <- GetAssayData(object = NPC_91, layer = "counts")
# sce_91 <- SingleCellExperiment(list(counts = counts))
# sce_91 <- decontX(sce_91)
# NPC_91[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce_91))
# 
# 
# NPC_92 <- readRDS("./01_tidy_data/iAL92_FastCar.rds") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
# counts <- GetAssayData(object = NPC_92, layer = "counts")
# sce_92 <- SingleCellExperiment(list(counts = counts))
# sce_92 <- decontX(sce_92)
# NPC_92[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce_92))

# ###################################### Visualize the effect of DeconX #############################################################
# umap <- reducedDim(sce_87, "decontX_UMAP")
# plotDimReduceCluster(x = sce_87$decontX_clusters,dim1 = umap[, 1], dim2 = umap[, 2])
# plotDecontXContamination(sce_87)
# sce_87 <- logNormCounts(sce_87)
# plotDimReduceFeature(as.matrix(logcounts(sce_87)),
#                      dim1 = umap[, 1],
#                      dim2 = umap[, 2],
#                      features = c("Hpd","C1qc","Ptprb", "Cd79a","Skap1", "Ccl5", "Tyrobp", "Spp1","Dapp1","S100a9", "Dcn", "Adgb", "Rbms3", "Snca" , "Cadm1"),
#                      exactMatch = TRUE)
# markers <- list(Hep_Markers=c( "Hpd", "Alb"),
#                 T_Markers =c("Cd3e", "Cd3d"),
#                 NK_Markers =c("Cd56","Cd16"),
#                 B_Markers = c("Cd19, Cd79a"),
#                 M_Markers = c("C1qc"),
#                 SAA = c("Saa1"))
# #cellTypeMappings <- list(Tcells = 2, Bcells = 5, Hepatocytes = 1, NKcells = 6)
# #The percetage of cells within a cluster that have detectable expression of marker genes
# plotDecontXMarkerPercentage(sce_87,
#                             markers = markers,
#                             #groupClusters = cellTypeMappings,
#                             assayName = c("counts", "decontXcounts"))
# 
# 
# #Some helpful hints when using plotDecontXMarkerPercentage:
# #
# #  Cell clusters can be renamed and re-grouped using the groupCluster parameter, which also needs to be a named list.
# #If groupCluster is used, cell clusters not included in the list will be excluded in the barplot. For example, if we wanted to group T-cells and NK-cells together,
# #we could set cellTypeMappings <- list(NK_Tcells = c(2,6), Bcells = 5, Monocytes = 1)
# #The level a gene that needs to be expressed to be considered detected in a cell can be adjusted using the threshold parameter.
# #If you are not using a SingleCellExperiment, then you will need to supply the original counts matrix or the decontaminated counts matrix as the first argument to generate the barplots.
# 
# plotDecontXMarkerExpression(sce_87,
#                             markers = markers[["T_Markers"]],
#                             groupClusters = list(Tcells = 2, Bcells = 5, Hepatocytes = 1, NKcells = 6),
#                             ncol = 3)
# sce_87 <- logNormCounts(sce_87,
#                         exprs_values = "decontXcounts",
#                         name = "decontXlogcounts")
# 
# plotDecontXMarkerExpression(sce_87,
#                             markers = markers[["M_Markers"]],
#                             #groupClusters = cellTypeMappings,
#                             ncol = 3,
#                             assayName = c("logcounts", "decontXlogcounts"))
# 
# sce_87 <- logNormCounts(sce_87,
#                         exprs_values = "decontXcounts",
#                         name = "decontXlogcounts")
# 
# 
# 
# 
# 
# ###########################Load Data and Decontaminate with SoupX ##################
#  #### Soup Decont for iAL87----
# 
#  sc = load10X('./00_raw_data/Liver_NPC_iAL87_transcriptome')
#  sc = setClusters(sc, sc$metaData$clustersFine)
#  ##Vizuals----
# 
#  dd = sc$metaData[colnames(sc$toc), ]
#  mids = aggregate(cbind(tSNE1, tSNE2) ~ clustersFine, data = dd, FUN = mean)
#  gg = ggplot(dd, aes(tSNE1, tSNE2))+
#  geom_point(aes(colour = clustersFine), size = 0.2) +
#  geom_label(data = mids, aes(label = clustersFine)) + ggtitle("NPC 87") +
#  guides(colour = guide_legend(override.aes = list(size = 1)))
# 
#  png("./03_plots/QC_1_SoupX_87_ClustersFine.png")
#  plot(gg)
#  dev.off()
# 
#  dd$val = sc$toc["Saa1", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_87_ClustersSAA1.png")
#  plot(gg)
#  dev.off()
#  dd$val = sc$toc["Alb", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_87_ClustersAlb.png")
#  plot(gg)
#  dev.off()
#  png("./03_plots/QC_1_SoupX_87_ClustersAlb_SIG.png")
#   gg = plotMarkerMap(sc, "Alb")
#  dev.off()
#  plot(gg)
# 
#  dd$val = sc$toc["Cd3e", ]
# gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
# png("./03_plots/QC_1_SoupX_87_ClustersCD3e.png")
#  plot(gg)
# dev.off()
# gg = plotMarkerMap(sc, "Cd3e")
# png("./03_plots/QC_1_SoupX_87_ClustersCD3e_SIG.png")
# plot(gg)
# 
# 
#  dd$val = sc$toc["Cd19", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_87_ClustersCD19.png")
#   plot(gg)
# dev.off()
#  gg = plotMarkerMap(sc, "Cd19")
#  png("./03_plots/QC_1_SoupX_87_ClustersCD19_SIG.png")
#  plot(gg)
#  dev.off()
# 
#  #### Stuff----
# 
# sc = autoEstCont(sc)
# head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
# 
#  out = adjustCounts(sc)
#  png("./03_plots/QC_1_SoupX_87_ChangedSAA1.png")
#  plotChangeMap(sc, out, "Saa1")
#  dev.off()
# 
#  dd$Saa1 = sc$toc["Saa1", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = Saa1 > 0))
#  png("./03_plots/QC_1_SoupX_87_CLUSTERSSAA AFTER.png")
#  plot(gg)
#  dev.off()
# 
#  gg = plotMarkerMap(sc, "Saa1")
#  png("./03_plots/QC_1_SoupX_87_CLUSTERSSAA SIG AFTER.png")
#  plot(gg)
#  dev.off()
#  NPC_87 = CreateSeuratObject(out)
# 
# 
#  #### Soup Decont for iAL88----
# 
#  sc = load10X('./00_raw_data/Liver_NPC_iAL88_transcriptome')
#  sc = setClusters(sc, sc$metaData$clustersFine)
#  ##Vizuals----
# 
#  dd = sc$metaData[colnames(sc$toc), ]
#  mids = aggregate(cbind(tSNE1, tSNE2) ~ clustersFine, data = dd, FUN = mean)
#  gg = ggplot(dd, aes(tSNE1, tSNE2))+
#    geom_point(aes(colour = clustersFine), size = 0.2) +
#    geom_label(data = mids, aes(label = clustersFine)) + ggtitle("NPC 88") +
#    guides(colour = guide_legend(override.aes = list(size = 1)))
# 
#  png("./03_plots/QC_1_SoupX_88_ClustersFine.png")
#  plot(gg)
#  dev.off()
# 
#  dd$val = sc$toc["Saa1", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_88_ClustersSAA1.png")
#  plot(gg)
#  dev.off()
#  dd$val = sc$toc["Alb", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_88_ClustersAlb.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Alb")
#  png("./03_plots/QC_1_SoupX_88_ClustersAlb_SIG.png")
#  plot(gg)
#  dev.off()
# 
#  dd$val = sc$toc["Cd3e", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_88_ClustersCD3e.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Cd3e")
#  png("./03_plots/QC_1_SoupX_88_ClustersCD3e_SIG.png")
#  plot(gg)
#  dev.off()
# 
#  dd$val = sc$toc["Cd19", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_88_ClustersCD19.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Cd19")
#  png("./03_plots/QC_1_SoupX_88_ClustersCD19_SIG.png")
#  plot(gg)
#  dev.off()
# 
#  #### Stuff----
# 
#  sc = autoEstCont(sc)
#  head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
# 
#  out = adjustCounts(sc)
#  png("./03_plots/QC_1_SoupX_88_ChangedSAA1.png")
#  plotChangeMap(sc, out, "Saa1")
#  dev.off()
# 
#  dd$Saa1 = sc$toc["Saa1", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = Saa1 > 0))
#  png("./03_plots/QC_1_SoupX_88_CLUSTERSSAA AFTER.png")
#  plot(gg)
#  dev.off()
# 
#  gg = plotMarkerMap(sc, "Saa1")
#  png("./03_plots/QC_1_SoupX_88_CLUSTERSSAA SIG AFTER.png")
#  plot(gg)
#  dev.off()
#  NPC_88 = CreateSeuratObject(out)
# 
# 
#  #### Soup Decont for iAL91----
# 
#  sc = load10X('./00_raw_data/Liver_NPC_iAL91_transcriptome')
#  sc = setClusters(sc, sc$metaData$clustersFine)
#  ##Vizuals----
# 
#  dd = sc$metaData[colnames(sc$toc), ]
#  mids = aggregate(cbind(tSNE1, tSNE2) ~ clustersFine, data = dd, FUN = mean)
#  gg = ggplot(dd, aes(tSNE1, tSNE2))+
#    geom_point(aes(colour = clustersFine), size = 0.2) +
#    geom_label(data = mids, aes(label = clustersFine)) + ggtitle("NPC 91") +
#    guides(colour = guide_legend(override.aes = list(size = 1)))
# 
#  png("./03_plots/QC_1_SoupX_91_ClustersFine.png")
#  plot(gg)
#  dev.off()
# 
#  dd$val = sc$toc["Saa1", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_91_ClustersSAA1.png")
#  plot(gg)
#  dev.off()
#  dd$val = sc$toc["Alb", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_91_ClustersAlb.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Alb")
#  png("./03_plots/QC_1_SoupX_91_ClustersAlb_SIG.png")
# plot(gg)
#   dev.off()
#    dd$val = sc$toc["Cd3e", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_91_ClustersCD3e.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Cd3e")
#  png("./03_plots/QC_1_SoupX_91_ClustersCD3e_SIG.png")
#  plot(gg)
#  dev.off()
# 
#  dd$val = sc$toc["Cd19", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_91_ClustersCD19.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Cd19")
#  png("./03_plots/QC_1_SoupX_91_ClustersCD19_SIG.png")
#  plot(gg)
#  dev.off()
# 
#  #### Stuff----
# 
#  sc = autoEstCont(sc)
#  head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
# 
#  out = adjustCounts(sc)
#  png("./03_plots/QC_1_SoupX_91_ChangedSAA1.png")
#  plotChangeMap(sc, out, "Saa1")
#  dev.off()
# 
#  dd$Saa1 = sc$toc["Saa1", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = Saa1 > 0))
#  png("./03_plots/QC_1_SoupX_91_CLUSTERSSAA AFTER.png")
#  plot(gg)
#  dev.off()
# 
#  gg = plotMarkerMap(sc, "Saa1")
#  png("./03_plots/QC_1_SoupX_91_CLUSTERSSAA SIG AFTER.png")
#  plot(gg)
#  dev.off()
#  NPC_91 = CreateSeuratObject(out)
# 
# 
#  #### Soup Decont for iAL92----
# 
#  sc = load10X('./00_raw_data/Liver_NPC_iAL92_transcriptome')
#  sc = setClusters(sc, sc$metaData$clustersFine)
#  ##Vizuals----
# 
#  dd = sc$metaData[colnames(sc$toc), ]
#  mids = aggregate(cbind(tSNE1, tSNE2) ~ clustersFine, data = dd, FUN = mean)
#  gg = ggplot(dd, aes(tSNE1, tSNE2))+
#    geom_point(aes(colour = clustersFine), size = 0.2) +
#    geom_label(data = mids, aes(label = clustersFine)) + ggtitle("NPC 92") +
#    guides(colour = guide_legend(override.aes = list(size = 1)))
# 
#  png("./03_plots/QC_1_SoupX_92_ClustersFine.png")
#  plot(gg)
#  dev.off()
# 
#  dd$val = sc$toc["Saa1", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_92_ClustersSAA1.png")
#  plot(gg)
#  dev.off()
#  dd$val = sc$toc["Alb", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_92_ClustersAlb.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Alb")
#  png("./03_plots/QC_1_SoupX_92_ClustersAlb_SIG.png")
#  plot(gg)
#  dev.off()
#  dd$val = sc$toc["Cd3e", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_92_ClustersCD3e.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Cd3e")
#  png("./03_plots/QC_1_SoupX_92_ClustersCD3e_SIG.png")
#  plot(gg)
#  dev.off()
# 
#  dd$val = sc$toc["Cd19", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
#  png("./03_plots/QC_1_SoupX_92_ClustersCD19.png")
#  plot(gg)
#  dev.off()
#  gg = plotMarkerMap(sc, "Cd19")
#  png("./03_plots/QC_1_SoupX_92_ClustersCD19_SIG.png")
#  plot(gg)
#  dev.off()
# 
#  #### Stuff----
# 
#  sc = autoEstCont(sc)
#  head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
# 
#  out = adjustCounts(sc)
#  png("./03_plots/QC_1_SoupX_92_ChangedSAA1.png")
#  plotChangeMap(sc, out, "Saa1")
#  dev.off()
# 
#  dd$Saa1 = sc$toc["Saa1", ]
#  gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = Saa1 > 0))
#  png("./03_plots/QC_1_SoupX_92_CLUSTERSSAA AFTER.png")
#  plot(gg)
#  dev.off()
# 
#  gg = plotMarkerMap(sc, "Saa1")
#  png("./03_plots/QC_1_SoupX_92_CLUSTERSSAA SIG AFTER.png")
#  plot(gg)
#  dev.off()
#  NPC_92 = CreateSeuratObject(out)
# 

########################## Kategorien Stimlulation und Sex hinzufügen, evtl noch age? ####################
numbers = c("002", "150", "250", "350", "450", "550")

for ( i in numbers){
NPC_87 <- readRDS(paste0("./01_tidy_data/iAL87_FastCar_",i,".rds"))
NPC_88 <- readRDS(paste0("./01_tidy_data/iAL88_FastCar_",i,".rds")) 
NPC_91 <- readRDS(paste0("./01_tidy_data/iAL91_FastCar_",i,".rds")) 
NPC_82 <- readRDS(paste0("./01_tidy_data/iAL92_FastCar_",i,".rds")) 
NPC_88$stim <- "EtOH"
NPC_91$stim <- "TAM"
NPC_92$stim <- "EtOH"
NPC_87$sex <- "female"
NPC_88$sex <- "female"
NPC_91$sex <- "male"
NPC_92$sex <- "male"
NPC_87$sample <- "iAL87"
NPC_88$sample <- "iAL88"
NPC_91$sample <- "iAL91"
NPC_92$sample <- "iAL92"
NPC_88$sample <- "iAL88"
NPC_91$sample <- "iAL91"
NPC_92$sample <- "iAL92"
############################### After Applying ambient RNA Decontamination Set Up Standard QC ##################################
#Setting things up for quality control mitochondrial genes, ribosomal content, doublets----
NPC_87[["percent.mt"]] <- PercentageFeatureSet(NPC_87, pattern = "^mt-") 
NPC_87[["percent.rb"]] <- PercentageFeatureSet(NPC_87, pattern = "Rp[sl]")
NPC_88[["percent.mt"]] <- PercentageFeatureSet(NPC_88, pattern = "^mt-") 
NPC_88[["percent.rb"]] <- PercentageFeatureSet(NPC_88, pattern = "Rp[sl]")
NPC_91[["percent.mt"]] <- PercentageFeatureSet(NPC_91, pattern = "^mt-") 
NPC_91[["percent.rb"]] <- PercentageFeatureSet(NPC_91, pattern = "Rp[sl]")
NPC_92[["percent.mt"]] <- PercentageFeatureSet(NPC_92, pattern = "^mt-") 
NPC_92[["percent.rb"]] <- PercentageFeatureSet(NPC_92, pattern = "Rp[sl]")


#Add Conclusions from meta data as QC6 columns QC6 from PMID: 36901774; UMI 500<UMI<40000, Genes 500<GENES<6000, MT MT<25 ----
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nFeature_RNA > 6000,'Doublet','Pass')
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 300 & NPC_87@meta.data$QC6 == 'Pass' , paste('Low_nFeature',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 300 & NPC_87@meta.data$QC6 != 'Pass' & NPC_87@meta.data$QC6 != 'Low_nFeature',paste('Low_nFeature',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$percent.mt > 25 & NPC_87@meta.data$QC6 == 'Pass','High_MT',NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nCount_RNA > 40000 & NPC_87@meta.data$QC6 == 'Pass','High_UMI',NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nCount_RNA < 500 & NPC_87@meta.data$QC6 == 'Pass','Low_UMI',NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 300 & NPC_87@meta.data$QC6 != 'Pass'& NPC_87@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(log10(NPC_87@meta.data$nFeature_RNA) / log10(NPC_87@meta.data$nCount_RNA)>= 0.8 & NPC_87@meta.data$QC6 != 'Pass'& NPC_87@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(log10(NPC_87@meta.data$nFeature_RNA) / log10(NPC_87@meta.data$nCount_RNA)< 0.8 & NPC_87@meta.data$QC6 == 'Pass'& NPC_87@meta.data$QC6 != 'High_MT',paste('low complex',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
table(NPC_87[['QC6']])

NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nFeature_RNA > 6000,'Doublet','Pass')
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 300 & NPC_88@meta.data$QC6 == 'Pass' , paste('Low_nFeature',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 300 & NPC_88@meta.data$QC6 != 'Pass' & NPC_88@meta.data$QC6 != 'Low_nFeature',paste('Low_nFeature',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$percent.mt > 25 & NPC_88@meta.data$QC6 == 'Pass','High_MT',NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nCount_RNA > 40000 & NPC_88@meta.data$QC6 == 'Pass','High_UMI',NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nCount_RNA < 500 & NPC_88@meta.data$QC6 == 'Pass','Low_UMI',NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 300 & NPC_88@meta.data$QC6 != 'Pass'& NPC_88@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(log10(NPC_88@meta.data$nFeature_RNA) / log10(NPC_88@meta.data$nCount_RNA)>= 0.8 & NPC_88@meta.data$QC6 != 'Pass'& NPC_88@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(log10(NPC_88@meta.data$nFeature_RNA) / log10(NPC_88@meta.data$nCount_RNA)< 0.8 & NPC_88@meta.data$QC6 == 'Pass'& NPC_88@meta.data$QC6 != 'High_MT',paste('low complex',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
table(NPC_88[['QC6']])

NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nFeature_RNA > 6000,'Doublet','Pass')
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 300 & NPC_91@meta.data$QC6 == 'Pass' , paste('Low_nFeature',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 300 & NPC_91@meta.data$QC6 != 'Pass' & NPC_91@meta.data$QC6 != 'Low_nFeature',paste('Low_nFeature',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$percent.mt > 25 & NPC_91@meta.data$QC6 == 'Pass','High_MT',NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nCount_RNA > 40000 & NPC_91@meta.data$QC6 == 'Pass','High_UMI',NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nCount_RNA < 500 & NPC_91@meta.data$QC6 == 'Pass','Low_UMI',NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 300 & NPC_91@meta.data$QC6 != 'Pass'& NPC_91@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(log10(NPC_91@meta.data$nFeature_RNA) / log10(NPC_91@meta.data$nCount_RNA)>= 0.8 & NPC_91@meta.data$QC6 != 'Pass'& NPC_91@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(log10(NPC_91@meta.data$nFeature_RNA) / log10(NPC_91@meta.data$nCount_RNA)< 0.8 & NPC_91@meta.data$QC6 == 'Pass'& NPC_91@meta.data$QC6 != 'High_MT',paste('low complex',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
table(NPC_91[['QC6']])

NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nFeature_RNA > 6000,'Doublet','Pass')
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 300 & NPC_92@meta.data$QC6 == 'Pass' , paste('Low_nFeature',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 300 & NPC_92@meta.data$QC6 != 'Pass' & NPC_92@meta.data$QC6 != 'Low_nFeature',paste('Low_nFeature',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$percent.mt > 25 & NPC_92@meta.data$QC6 == 'Pass','High_MT',NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nCount_RNA > 40000 & NPC_92@meta.data$QC6 == 'Pass','High_UMI',NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nCount_RNA < 500 & NPC_92@meta.data$QC6 == 'Pass','Low_UMI',NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 300 & NPC_92@meta.data$QC6 != 'Pass'& NPC_92@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(log10(NPC_92@meta.data$nFeature_RNA) / log10(NPC_92@meta.data$nCount_RNA)>= 0.8 & NPC_92@meta.data$QC6 != 'Pass'& NPC_92@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(log10(NPC_92@meta.data$nFeature_RNA) / log10(NPC_92@meta.data$nCount_RNA)< 0.8 & NPC_92@meta.data$QC6 == 'Pass'& NPC_92@meta.data$QC6 != 'High_MT',paste('low complex',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
table(NPC_92[['QC6']])

#######################  Generate the Object "metadata" to plot Things more easily ############################################################################################
#https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html


metadata <- rbind(NPC_87@meta.data,NPC_88@meta.data, NPC_91@meta.data, NPC_92@meta.data)
metadata$log10GenesPerUMI <- log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)

metadata$cells <- rownames(metadata)
metadata <- metadata %>%  dplyr::rename(nUMI = nCount_RNA,nGene = nFeature_RNA)
metadata_QC6 <-subset(metadata, QC6 == "Pass")
########################### Visualizations of the QC Parameters on Data wo QC  #############################################################################
# Visualize number of cell per sample----
png(paste0("./03_plots/QC_noQC_Cells_per_sample_bar",i,".png"))
x <-metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells")
print(x)
dev.off()

# Visualize number of UMIS per sample violin----
png(paste0("./03_plots/QC_noQC_UMIs_per_sample_violin",i,".png"))
x <-metadata %>% 
  ggplot(aes(x=sample, fill=sample, y =nUMI)) + 
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 500)+
  geom_hline(yintercept = 40000)+
  ggtitle("nUMI")
print(x)
dev.off()

# Visualize number of UMIS per sample density plot----
png(paste0("./03_plots/QC_noQC_UMIs_per_sample_density",i,".png"))  
x <- metadata%>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+
  geom_vline(xintercept = 40000)
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample violin----
png(paste0("./03_plots/QC_noQC_Gene(log10)_per_sample_violin",i,".png"))  
x <- metadata %>% 
  ggplot(aes(x=sample, fill=sample, y =log10(nGene))) + 
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  geom_hline(yintercept = log10(6000))+
  geom_hline(yintercept = log10(300))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Log10 GENE")
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample BoxPlot----
png(paste0(("./03_plots/QC_noQC_Gene(log10)_per_sample_box",i,".png"))  
x <- metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = log10(6000))+
  geom_hline(yintercept = log10(300))+
  ggtitle("log10 NGenes")
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample DensityPlot----
png(paste0(("./03_plots/QC_noQC_Gene(log10)_per_sample_density",i,".png"))  
x <- metadata %>% 
  ggplot(aes(color=sample, x=log10(nGene), fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()+  
  geom_vline(xintercept = log10(6000))+
  geom_vline(xintercept = log10(300))+
  ggtitle("log10 NGenes")
print(x)
dev.off()

# Visualize number of Genes per sample ViolinPlot----
png(paste0(("./03_plots/QC_noQC_Gene_per_sample_violin",i,".png"))  
x <- metadata %>% 
  ggplot(aes(x=sample, fill=sample, y =nGene)) + 
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 6000)+
  geom_hline(yintercept = 300)+
  ggtitle("nGENE")
print(x)
dev.off()

# Visualize number of Genes per sample BoxPlot----
png(paste0(("./03_plots/QC_noQC_Gene_per_sample_box",i,".png"))  
x <- metadata %>% 
  ggplot(aes(x=sample, y=nGene, fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 6000)+
  geom_hline(yintercept = 300)+
  ggtitle("nGenes")
print(x)
dev.off()

# Visualize number of Genes per sample DensityPlot----
png(paste0(("./03_plots/QC_noQC_Gene_per_sample_density",i,".png"))  
x <- metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  ggtitle("nGenes")+
  geom_vline(xintercept = 6000)+
  geom_vline(xintercept = 300)
print(x)
dev.off()


# Visualize number of percent.met per sample Violin plot----
png(paste0(("./03_plots/QC_noQC_Mito_per_sample_violin",i,".png"))  
x <- metadata %>% 
  ggplot(aes(x=sample, fill=sample, y =percent.mt)) + 
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 25)+
  ggtitle("Percent mt")
print(x)
dev.off()

# Visualize number of percent.mt per sample Density plot----
png(paste0(("./03_plots/QC_noQC_Mito_per_sample_density",i,".png"))  
x <- metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 25)
print(x)
dev.off()

# Visualize number of percent.rb per sample Violin plot----
png(paste0(("./03_plots/QC_noQC_Rb_per_sample_violin",i,".png"))  
x <- metadata %>% 
  ggplot(aes(x=sample, fill=sample, y =percent.rb)) + 
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Percent rb")
print(x)
dev.off()

# Visualize number of Complexity of cells/sample as nGene vs nUMI----
png(paste0(("./03_plots/QC_noQC_Gene_vs_UMI_persample",i,".png"))  
x <- metadata%>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 40000) +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 6000) +
  facet_wrap(~sample)
print(x)
dev.off()

# Visualize number of Complexity of cells/sample as log10Genes/UMI----
png(paste0(("./03_plots/QC_noQC_Gene perUMI_density",i,".png"))  
x <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
print(x)
dev.off()


# ########################### Visualizations of the QC Parameters on Data with QC6############################################################################

# Visualize number of cell per sample----
png(paste0(("./03_plots/QC_QC6_Cells_per_sample_bar",i,".png"))  
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells")
print(x)
dev.off()

# Visualize number of UMIS per sample violin----
png(paste0(("./03_plots/QC_QC6_UMIs_per_sample_violin",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =nUMI)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 500)+
  geom_hline(yintercept = 40000)+
  ggtitle("nUMI")
print(x)
dev.off()

# Visualize number of UMIS per sample density plot----
png(paste0(("./03_plots/QC_QC6_UMIs_per_sample_density",i,".png")) 
x <- metadata_QC6%>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+
  geom_vline(xintercept = 40000)
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample violin----
png(paste0(("./03_plots/QC_QC6_Gene(log10)_per_sample_violin",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =log10(nGene))) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  geom_hline(yintercept = log10(6000))+
  geom_hline(yintercept = log10(300))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Log10 GENE")
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample BoxPlot----
png(paste0(("./03_plots/QC_QC6_Gene(log10)_per_sample_box",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = log10(6000))+
  geom_hline(yintercept = log10(300))+
  ggtitle("log10 NGenes")
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample DensityPlot----
png(paste0(("./03_plots/QC_QC6_Gene(log10)_per_sample_density",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(color=sample, x=log10(nGene), fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10()+
  geom_vline(xintercept = log10(6000))+
  geom_vline(xintercept = log10(300))+
  ggtitle("log10 NGenes")
print(x)
dev.off()


# Visualize number of Genes per sample ViolinPlot----
png(paste0(("./03_plots/QC_QC6_Gene_per_sample_violin",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =nGene)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 6000)+
  geom_hline(yintercept = 300)+
  ggtitle("nGENE")
print(x)
dev.off()

# Visualize number of Genes per sample BoxPlot----
png(paste0(("./03_plots/QC_QC6_Gene_per_sample_box",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, y=nGene, fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 6000)+
  geom_hline(yintercept = 300)+
  ggtitle("nGenes")
print(x)
dev.off()

# Visualize number of Genes per sample DensityPlot----
png(paste0(("./03_plots/QC_QC6_Gene_per_sample_density",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ggtitle("nGenes")+
  geom_vline(xintercept = 6000)+
  geom_vline(xintercept = 300)
print(x)
dev.off()


# Visualize number of percent.met per sample Violin plot----
png(paste0(("./03_plots/QC_QC6_Mito_per_sample_violin",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =percent.mt)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 25)+
  ggtitle("Percent mt")
print(x)
dev.off()

# Visualize number of percent.mt per sample Density plot----
png(paste0(("./03_plots/QC_QC6_Mito_per_sample_density",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 25)
print(x)
dev.off()

# Visualize number of percent.rb per sample Violin plot----
png(paste0(("./03_plots/QC_QC6_Rb_per_sample_violin",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =percent.rb)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Percent rb")
print(x)
dev.off()

# Visualize number of Complexity of cells/sample as nGene vs nUMI----
png(paste0(("./03_plots/QC_QC6_Gene_vs_UMI_persample",i,".png")) 
x <- metadata_QC6%>%
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 40000) +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 6000) +
  facet_wrap(~sample)
print(x)
dev.off()

# Visualize number of Complexity of cells/sample as log10Genes/UMI----
png(paste0(("./03_plots/QC_QC6_Gene perUMI_density",i,".png")) 
x <- metadata_QC6 %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
print(x)
dev.off()

#################### Normalize Data and Find Variable Features Data wo QC ############################################################
NPC_87 <- NormalizeData(NPC_87)
NPC_87 <- FindVariableFeatures(NPC_87, selection.method = "vst", nfeatures = 2000)
NPC_88 <- NormalizeData(NPC_88)
NPC_88 <- FindVariableFeatures(NPC_88, selection.method = "vst", nfeatures = 2000)
NPC_91 <- NormalizeData(NPC_91)
NPC_91 <- FindVariableFeatures(NPC_91, selection.method = "vst", nfeatures = 2000)
NPC_92 <- NormalizeData(NPC_92)
NPC_92 <- FindVariableFeatures(NPC_92, selection.method = "vst", nfeatures = 2000)

#################### Normalize Data and Find Variable Features Data with QC6 ############################################################

#normalize data set to account for sequencing depth, default scale to 10 000 and log2-transform
NPC_87_QC6 <- NormalizeData(subset(NPC_87, subset = QC6 == 'Pass'))
NPC_87_QC6 <- FindVariableFeatures(NPC_87_QC6, selection.method = "vst", nfeatures = 2000)
NPC_88_QC6 <- NormalizeData(subset(NPC_88, subset = QC6 == 'Pass'))
NPC_88_QC6 <- FindVariableFeatures(NPC_88_QC6, selection.method = "vst", nfeatures = 2000)
NPC_91_QC6 <- NormalizeData(subset(NPC_91, subset = QC6 == 'Pass'))
NPC_91_QC6 <- FindVariableFeatures(NPC_91_QC6, selection.method = "vst", nfeatures = 2000)
NPC_92_QC6 <- NormalizeData(subset(NPC_92, subset = QC6 == 'Pass'))
NPC_92_QC6 <- FindVariableFeatures(NPC_92_QC6, selection.method = "vst", nfeatures = 2000)
#Remove unused data from memory to save ram
rm(NPC_87.data,NPC_88.data, NPC_91.data, NPC_92.data, metadata_QC5, metadata_QC6, metadata)
saveRDS(NPC_87_QC6, paste0("./01_tidy_data/QC_1_QC6_NPC_87",i,".rds"))
saveRDS(NPC_88_QC6, paste0("./01_tidy_data/QC_1_QC6_NPC_88",i,".rds"))
saveRDS(NPC_91_QC6, paste0("./01_tidy_data/QC_1_QC6_NPC_91",i,".rds"))
saveRDS(NPC_92_QC6, paste0("./01_tidy_data/QC_1_QC6_NPC_92",i,".rds"))
rm(NPC_87_QC6,NPC_88_QC6,NPC_88_QC6,NPC_91_QC6,NPC_92_QC6)
#### for loop ends ####
}

########################## Define Anchors for Integration and Integrate Different Data Sets ####
## fill in correct FastCar CutOFF!
#NPC_87 <-readRDS("./01_tidy_data/QC_1_QC6_NPC_87",i".rds"))
#NPC_88 <-readRDS("./01_tidy_data/QC_1_QC6_NPC_88",i".rds"))
#NPC_91 <-readRDS("./01_tidy_data/QC_1_QC6_NPC_91",i".rds"))
#NPC_92 <-readRDS("./01_tidy_data/QC_1_QC6_NPC_92",i".rds"))
NPC.anchors <- FindIntegrationAnchors(object.list = list(NPC_87, NPC_88, NPC_91, NPC_92), dims = 1:20)
NPC.combined<- IntegrateData(anchorset = NPC.anchors, dims = 1:20)
NPC_QC6.anchors <- FindIntegrationAnchors(object.list = list(NPC_87_QC6, NPC_88_QC6, NPC_91_QC6, NPC_92_QC6), dims = 1:20)
NPC_QC6.combined<- IntegrateData(anchorset = NPC_QC6.anchors, dims = 1:20)
rm(NPC_87,NPC_88, NPC_91, NPC_92,NPC.anchors, NPC_QC5.anchors, NPC_QC6.anchors )
######################## Save Integrated Data Sets with QC (QC5 and QC6 and wo QC) #########
saveRDS(NPC.combined, "./01_tidy_data/QC_noQC_NPC.combined")
saveRDS(NPC_QC6.combined, "./01_tidy_data/QC_QC6_NPC_QC6.combined")
rm(, NPC_87_QC6, NPC_87, NPC_88, , NPC_88_QC6, NPC_91, , NPC_91_QC6, NPC_92, , NPC_92_QC6 )
rm(NPC.combined,  NPC_QC6.combined)
########################################################################
########################################################################

