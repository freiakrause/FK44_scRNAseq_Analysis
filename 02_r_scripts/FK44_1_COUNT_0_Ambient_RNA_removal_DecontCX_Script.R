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
library(gprofiler2)
library(decontX)
library(celda)
library(scater)
library(stringr)
library(Matrix)
library(gridExtra)
library(qlcMatrix)
library(pheatmap)
set.seed(42)

################### Load Output from CellRanger and generate decontaminated Counts with DeconX ############
#deconX= ambient RNA is floating around and got incorporated into GEMs and now cells which should not expres certain RNA seem to express it.
#e.g all cell types show Saa upregulation or whatever. but this is from dying heps that spread their RNA everywhere
#So deconX kinda filters
animals <-c("87","88","91","92")
interesting_genes <-c("Cxcl1","Ncam1","Fcgr3","Cd3e","Cd3d","Cd4","Cd8a","Alb","Saa1","Saa2","Hpd","C1qc","Ptprb", "Cd79a","Skap1", "Ccl5", "Tyrobp", "Spp1","Dapp1","S100a9", "Dcn", "Adgb", "Rbms3", "Snca" , "Cadm1")
markers <- list(Hep_Markers=c("Alb","Saa1","Saa2"),
                T_Markers =c("Cd3e","Cd3d","Cd4","Cd8a"),
                NK_Markers =c("Ncam1","Fcgr3"),
                B_Markers = c("Cd19", "Cd79a"),
                M_Markers = c("C1qc"),
                SAA = c("Saa1","Saa2"))
for (i in animals){
  filtered <- Read10X(paste0("./00_raw_data/Liver_NPC_iAL",i,"_transcriptome/filtered_feature_bc_matrix"))
  print(paste0(i))
  raw <-Read10X(paste0("./00_raw_data/Liver_NPC_iAL",i,"_transcriptome/raw_feature_bc_matrix/")) #CellRangerinitializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
  filtered <- SingleCellExperiment(list(counts = filtered))
  raw <- SingleCellExperiment(list(counts = raw))
  decon<- decontX(filtered, background = raw)
  saveRDS(decon, paste0("./01_tidy_data/NPC_",i,"_DecontX.rds"))
  umap <- reducedDim(decon, "decontX_UMAP")
  png(paste0("./03_plots/QC_0_Ambient_DecontX_",i,"_cluster.png"))
  p <-plotDimReduceCluster(x = decon$decontX_clusters,dim1 = umap[, 1], dim2 = umap[, 2])
  print(p)
  dev.off()
  png(paste0("./03_plots/QC_0_Ambient_DecontX_",i,"_cluster_decon.png"))
  p <-plotDecontXContamination(decon)
  print(p)
  dev.off()
  decon <- logNormCounts(decon)
  
  
  png((paste0("./03_plots/QC_0_Ambient_DecontX_",i,"_MarkerExpression_all.png"))) 
  p <-plotDecontXMarkerPercentage(decon,markers = markers,#groupClusters = cellTypeMappings,
                              assayName = c("counts", "decontXcounts"))
  print(p)
  dev.off()
  
  for (m in interesting_genes){
    png((paste0("./03_plots/QC_0_Ambient_DecontX_",i,"_MarkerExpression_",m,".png"))) 
    p <-plotDecontXMarkerExpression(decon,
                              markers = m,
                              #groupClusters = list(Tcells = 2, Bcells = 5, Hepatocytes = 1, NKcells = 6),
                              ncol = 4)
    print(p)
    dev.off()
  }
  #to look at normalized expression
  decon <- logNormCounts(decon,exprs_values = "decontXcounts",name = "decontXlogcounts")
  for (f in interesting_genes){
    png((paste0("./03_plots/QC_0_Ambient_DecontX_",i,"MarkerExpression_Normalized_",f,".png"))) 
    p<-plotDecontXMarkerExpression(decon,markers =f,#groupClusters = cellTypeMappings,
                                   ncol = 4,         assayName = c("logcounts", "decontXlogcounts"))
    print(p)
    dev.off()
  }
}
NPC_87 <-readRDS("./01_tidy_data/NPC_87_DecontX.rds")

###################################### Visualize the effect of DeconX #############################################################
#cellTypeMappings <- list(Tcells = 2, Bcells = 5, Hepatocytes = 1, NKcells = 6)
#The percetage of cells within a cluster that have detectable expression of marker genes
 


#Some helpful hints when using plotDecontXMarkerPercentage:
#
#  Cell clusters can be renamed and re-grouped using the groupCluster parameter, which also needs to be a named list.
#If groupCluster is used, cell clusters not included in the list will be excluded in the barplot. For example, if we wanted to group T-cells and NK-cells together,
#we could set cellTypeMappings <- list(NK_Tcells = c(2,6), Bcells = 5, Monocytes = 1)
#The level a gene that needs to be expressed to be considered detected in a cell can be adjusted using the threshold parameter.
#If you are not using a SingleCellExperiment, then you will need to supply the original counts matrix or the decontaminated counts matrix as the first argument to generate the barplots.


