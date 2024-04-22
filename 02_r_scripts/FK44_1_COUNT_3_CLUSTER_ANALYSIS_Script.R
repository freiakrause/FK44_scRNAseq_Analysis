#This is the script to analyse FK44.1 scRNAseq COUNT data provided by BSF; Data from 1_QC SCript is loaded
######################## Installed Packages #############################################
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
#renv::install("rlang")
#renv::install("patchwork")
#renv::install("igraph") # für das scheiß igraph (benötigt für seurat) nach tausend jahren troublehsooting gefunden: ich brauche: sudo apt install build-essential gfortran UND sudo apt install build-essential gfortran, dann gehts
#renv::install("Seurat")
#renv::install("HGNChelper")
#renv::install("openxlsx")
#if (!require("BiocManager", quietly = TRUE))
#  renv::install("BiocManager")
#renv::install("gprofiler2")
#BiocManager::install(version = "3.18")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install('glmGamPoi')
#remotes::install_github("Moonerss/scrubletR")
#remotes::install_github('immunogenomics/presto')
#remotes::install_github('mojaveazure/ggseurat')
#renv::install("svglite")
#################### Libraries #########################
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(scrubletR)
library(gprofiler2)
library(rlang)
library(ggseurat)
library(patchwork)
library(svglite)
library(HGNChelper)
library(openxlsx)
set.seed(42)

################################# Load Input Data  ##################################################
#Basic QC is done. QC5 and QC6 contain filtered data: UMI count, Gene Count, MT% and Complexity
#We did SCTransformation with regression in MT percentage and CellCycleScore
#runUMAP and FindeNeighbours was run with dim 1:30 each, findClusters was run with reoslution 0.1
NPC_CLUSTER <- readRDS("~/FK44_NPC_scRNAseq_Analysis/01_tidy_data/NPC_ALL_TRANSFORM.rds")
###################### Function Create_Vplots ####################################
#Create multiple ViolinPlots with Seurat Object and vector of displayed features as input as save the plots as png
Create_Vplots <- function(x,y){
  for (i in y){
    print(i)
    a <- VlnPlot(x, features = i)
    png(filename = paste0("./03_plots/Clustermarker_",i,".png"))
    print(a)
    dev.off()
    
    'svglite(filename = paste0("Clustermarker_",i,".svg"),
          width = 10,
          height = 10,
          bg="transparent",
          pointsize=10,
          fix_text_size=F)
  print(a)
  dev.off()'
  }
}
####### Create Violin Plots of all Clusters with shown potential marker genes ############
#From my TOP10 Clustermarker table
features <- c("Hpd","C1qc","Ptprb", "Cd79a","Skap1", "Ccl5", "Tyrobp", "Spp1","Dapp1","S100a9", "Dcn", "Adgb", "Rbms3", "Snca" , "Cadm1") 
#From Internet/head/FACS
features <- c("Alb","Clec4f","Ptprb", "Cd19","Cd3e", "Nkg7", "Cd14", "Lyz1","Lyz2","Timp1", "Col1a2", "Krt19", "Ly6g", "C1qa" , "Itgax", "Cd86") 
Create_Vplots(NPC_CLUSTER,features)
#FROM Scott Paper DEGs in mouse cells
VlnPlot(NPC_CLUSTER, features = "Mmrn2")
VlnPlot(NPC_CLUSTER, features = "Ntm")
VlnPlot(NPC_CLUSTER, features = "Fabp1")
VlnPlot(NPC_CLUSTER, features = "Apoa2")
VlnPlot(NPC_CLUSTER, features = "Ddit4l")
VlnPlot(NPC_CLUSTER, features = "Cd209a")
VlnPlot(NPC_CLUSTER, features = "Nudt17")
VlnPlot(NPC_CLUSTER, features = "Ly6i")

############################## Dimension Plot of Clusters ###########################################
DimPlot(NPC_CLUSTER, split.by = "orig.ident")
DimPlot(NPC_CLUSTER, split.by = "stim")
DimPlot(NPC_CLUSTER, split.by = "sex")
############################### Create Subsets by Cluster #########################################
NPC_0 <- subset(NPC_CLUSTER, seurat_clusters== "0")
NPC_1 <- subset(NPC_CLUSTER, seurat_clusters== "1")
NPC_2 <- subset(NPC_CLUSTER, seurat_clusters== "2")
NPC_3 <- subset(NPC_CLUSTER, seurat_clusters== "3")
NPC_4 <- subset(NPC_CLUSTER, seurat_clusters== "4")
NPC_5 <- subset(NPC_CLUSTER, seurat_clusters== "5")
NPC_6 <- subset(NPC_CLUSTER, seurat_clusters== "6")
NPC_7 <- subset(NPC_CLUSTER, seurat_clusters== "7")
NPC_8 <- subset(NPC_CLUSTER, seurat_clusters== "8")
NPC_9 <- subset(NPC_CLUSTER, seurat_clusters== "9")
NPC_10<- subset(NPC_CLUSTER, seurat_clusters== "10")
NPC_11<- subset(NPC_CLUSTER, seurat_clusters== "11")
NPC_12<- subset(NPC_CLUSTER, seurat_clusters== "12")

DimPlot(NPC_0, split.by = "stim")
DimPlot(NPC_1, split.by = "stim")
DimPlot(NPC_2, split.by = "stim")
DimPlot(NPC_3, split.by = "stim")
DimPlot(NPC_4, split.by = "stim")
DimPlot(NPC_5, split.by = "stim")
DimPlot(NPC_6, split.by = "stim")
DimPlot(NPC_7, split.by = "stim")
DimPlot(NPC_8, split.by = "stim")
DimPlot(NPC_9, split.by = "stim")
DimPlot(NPC_10, split.by = "stim")
DimPlot(NPC_11, split.by = "stim")
DimPlot(NPC_12, split.by = "stim")
VlnPlot(NPC_0, features= "Alb", split.by = "stim")
VlnPlot(NPC_0, features= "Saa1", split.by = "stim")
VlnPlot(NPC_1, features= "Ptprb", split.by = "stim")
VlnPlot(NPC_2, features= "Cd3e", split.by = "stim")
VlnPlot(NPC_2, features= "Cd4", split.by = "stim")
VlnPlot(NPC_2, features= "Cd8a", split.by = "stim")
VlnPlot(NPC_2, features= "Gata3", split.by = "stim")
VlnPlot(NPC_2, features= "Foxp3", split.by = "stim")
VlnPlot(NPC_2, features= "Pdcd1", split.by = "stim")
VlnPlot(NPC_2, features= "Il10", split.by = "stim")
VlnPlot(NPC_2, features= "Il6", split.by = "stim")
VlnPlot(NPC_4, features= "Il6", split.by = "stim")
VlnPlot(NPC_4, features= "Il10", split.by = "stim")
VlnPlot(NPC_4, features= "Clec4f", split.by = "stim")
#########################################################################
############# Try subclustering ###################
NPC_0 <- SCTransform(NPC_0, verbose = F,vars.to.regress = c("sex","stim"))
NPC_0 <- RunPCA(NPC_0, verbose = F)
NPC_0 <- RunUMAP(NPC_0, dims = 1:30, verbose = F)
NPC_0 <- FindNeighbors(NPC_0, dims = 1:30, verbose = F)
NPC_0 <- FindClusters(NPC_0, verbose = F,resolution = 0.2, save.SNN = TRUE)
NPC_0 <- SetIdent(NPC_0, value = "orig.ident")
NPC_0 <- SetIdent(NPC_0, value = "stim")
NPC_0 <- SetIdent(NPC_0, value = "sex")
NPC_0 <- SetIdent(NPC_0, value = "seurat_clusters")
DimPlot(NPC_0)

all.markers <- FindAllMarkers(NPC_0, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 0.5) %>%  slice_head(n = 10) %>%  ungroup() -> top10__MOUSEMAIN_CLUSTERMARKER


###########################################################################
#############################################################################
####### Try subclustering###################
NPC_2 <- SCTransform(NPC_2, verbose = F,vars.to.regress = c("sex","stim"))
NPC_2 <- SCTransform(NPC_2, verbose = F)
NPC_2 <- RunPCA(NPC_2, verbose = F)
NPC_2 <- RunUMAP(NPC_2, dims = 1:30, verbose = F)
NPC_2 <- FindNeighbors(NPC_2, dims = 1:30, verbose = F)
NPC_2 <- FindClusters(NPC_2, verbose = F,resolution = 0.2, save.SNN = TRUE)
NPC_2 <- SetIdent(NPC_2, value = "orig.ident")
NPC_2 <- SetIdent(NPC_2, value = "stim")
NPC_2 <- SetIdent(NPC_2, value = "sex")
NPC_2 <- SetIdent(NPC_2, value = "seurat_clusters")
DimPlot(NPC_2)

all.markers <- FindAllMarkers(NPC_2, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 0.5) %>%  slice_head(n = 10) %>%  ungroup() -> top10__MOUSEMAIN_CLUSTERMARKER
VlnPlot(NPC_2, features= "Ptprc", split.by = "stim")
VlnPlot(NPC_2, features= "Cd3e", split.by = "stim")
VlnPlot(NPC_2, features= "Cd8a", split.by = "stim")
VlnPlot(NPC_2, features= "Gzma", split.by = "stim")
VlnPlot(NPC_2, features= "Nkg2d", split.by = "stim")
VlnPlot(NPC_2, features= "Cd4", split.by = "stim")
VlnPlot(NPC_2, features= "Tbx21", split.by = "stim")
VlnPlot(NPC_2, features= "Gata3", split.by = "stim")
VlnPlot(NPC_2, features= "Rorc", split.by = "stim")

######## Potentiall Interesting link for subclustering or DE ######################
#####https://satijalab.org/seurat/articles/de_vignette
