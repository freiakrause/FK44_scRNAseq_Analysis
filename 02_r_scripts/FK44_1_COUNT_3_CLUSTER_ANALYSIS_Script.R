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
#renv::install("sccustomize")

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
#remotes::install_github('satijalab/seurat-data')
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




# ############################### Create Subsets by Cluster #########################################
# Cluster <-  unique(as.list(NPC_CLUSTER@meta.data$mouseRNA.main))
# for (c in Cluster){
#   x<-subset(NPC_CLUSTER, mouseRNA.main== c)
#   l<- length(x@meta.data$sex)
#   if  (l>10)
#     {
#     print(paste0("Cluster of ",c," has ",l," cells."))
#     p<-DimPlot(x, split.by = "stim")
#     print(p)
#     x <- SCTransform(x, verbose = F,vars.to.regress = c("S.Score","G2M.Score"))
#     x <- RunPCA(x, verbose = F)
#     x <- RunUMAP(x, dims = 1:30, verbose = F)
#     x <- FindNeighbors(x, dims = 1:30, verbose = F)
#     x <- FindClusters(x, verbose = F,resolution = 0.8, save.SNN = TRUE)
#     x <- SetIdent(x, value = "orig.ident")
#     DimPlot(x)
#     x <- SetIdent(x, value = "stim")
#     DimPlot(x)
#     x <- SetIdent(x, value = "sex")
#     DimPlot(x)
#     x <- SetIdent(x, value = "seurat_clusters")
#     DimPlot(x)
#     }
#     else{print(paste0("Cluster of ",x," is annoying and has less than 10 cells"))}
# }
# x<-subset(NPC_CLUSTER, mouseRNA.main== "Hepatocytes")
# l<- length(x@meta.data$sex)
# length(subset(NPC_CLUSTER, mouseRNA.main== "Hepatocytes")@meta.data$sex)
# length(subset(NPC_CLUSTER, mouseRNA.main== "Epithelial cells")@meta.data$sex)
# length(NPC_CLUSTER@meta.data$nFeature_RNA)
# length(NPC_CLUSTER@meta.data$sex)
# #########################################################################
