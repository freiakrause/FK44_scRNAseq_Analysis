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
#runUMAP and FindeNeighbours was run with dim 1:30 each, findClusters was run with resolution 0.3
NPC_CLUSTER <- readRDS("./01_tidy_data/4_NPC_ALL_TRANSFORM.rds")
###################### Function Create_Vplots ####################################
#Create multiple ViolinPlots with Seurat Object and vector of displayed features as input as save the plots as png
Create_Vplots <- function(DataSet,feature_list){
  for (i in feature_list){
    print(i)
    a <- VlnPlot(DataSet, features = i)
    #png(filename = paste0("./03_plots/Clustermarker_",i,".png"))
    print(a)
    #dev.off()
    
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
features_TOP10 <- c("Hpd","C1qc","Ptprb", "Cd79a","Skap1", "Ccl5", "Tyrobp", "Spp1","Dapp1","S100a9", "Dcn", "Adgb", "Rbms3", "Snca" , "Cadm1") 
#From Internet/head/FACS
features_FACS <- c("Alb","Clec4f", "Cd19","Cd3e", "Nkg7", "Cd14", "Lyz1","Lyz2","Timp1", "Col1a2", "Krt19", "Ly6g", "C1qa" , "Itgax", "Cd86") 
#FROM Scott Paper DEGs in mouse cells
features_DEG_Scott <-c("Mmrn2","Ntm","Fabp1","Apoa2","Ddit4l","Cd209a", "Nudt17","Ly6i")
features_Cytokine <-c("Il6","Il10","Il1ß","Tnfa")

NPC_CLUSTER<-SetIdent(NPC_CLUSTER,value = "mouseRNA.main")
Create_Vplots(NPC_CLUSTER,features_TOP10)
Create_Vplots(NPC_CLUSTER,features_FACS)
Create_Vplots(NPC_CLUSTER,features_DEG_Scott)


############################## Dimension Plot of Clusters ###########################################
DimPlot(NPC_CLUSTER, split.by = "orig.ident")
DimPlot(NPC_CLUSTER, split.by = "stim")
DimPlot(NPC_CLUSTER, split.by = "sex")
############################### Create Subsets by Cluster #########################################
Cluster <-  unique(as.list(NPC_CLUSTER@meta.data$mouseRNA.main))
for (c in Cluster){
  x<-subset(NPC_CLUSTER, mouseRNA.main== c)
  l<- length(x@meta.data$sex)
  if  (l>10)
    {
    print(paste0("Cluster of ",c," has ",l," cells."))
    p<-DimPlot(x, split.by = "stim")
    print(p)
    x <- SCTransform(x, verbose = F,vars.to.regress = c("S.Score","G2M.Score"))
    x <- RunPCA(x, verbose = F)
    x <- RunUMAP(x, dims = 1:30, verbose = F)
    x <- FindNeighbors(x, dims = 1:30, verbose = F)
    x <- FindClusters(x, verbose = F,resolution = 0.8, save.SNN = TRUE)
    x <- SetIdent(x, value = "orig.ident")
    DimPlot(x)
    x <- SetIdent(x, value = "stim")
    DimPlot(x)
    x <- SetIdent(x, value = "sex")
    DimPlot(x)
    x <- SetIdent(x, value = "seurat_clusters")
    DimPlot(x)
    }
    else{print(paste0("Cluster of ",x," is annoying and has less than 10 cells"))}
}
x<-subset(NPC_CLUSTER, mouseRNA.main== "Hepatocytes")
l<- length(x@meta.data$sex)
length(subset(NPC_CLUSTER, mouseRNA.main== "Hepatocytes")@meta.data$sex)
length(subset(NPC_CLUSTER, mouseRNA.main== "Epithelial cells")@meta.data$sex)
length(NPC_CLUSTER@meta.data$nFeature_RNA)
length(NPC_CLUSTER@meta.data$sex)
#########################################################################
############# Try subclustering ###################


all.markers <- FindAllMarkers(NPC_0, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 0.5) %>%  slice_head(n = 10) %>%  ungroup() -> top10__MOUSEMAIN_CLUSTERMARKER


###########################################################################



######## Potentiall Interesting link for subclustering or DE ######################
#####https://satijalab.org/seurat/articles/de_vignette
