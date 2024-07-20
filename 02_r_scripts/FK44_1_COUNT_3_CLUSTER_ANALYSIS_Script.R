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


#### Load Input Data #####
NPC_CLUSTER <-readRDS(file = "./01_tidy_data/5_NPC_ALL_TRANSFORM_Markers_on_integrated.rds")
Clustering_3_top10_mouseMAIN_CLUSTERMARKER<-read.csv("./99_other/Clustering_3_top10_mouseMAIN_CLUSTERMARKER_integrated.csv")

###################### Function Create_Vplots ####################################
#Create multiple ViolinPlots with Seurat Object and vector of displayed features as input as save the plots as png
Create_Vplots <- function(DataSet,feature_list){
  for (i in feature_list){
    print(i)
    a <- VlnPlot(DataSet, features = i)
    png(filename = paste0("./03_plots/2_Clustering/Clustermarker_",i,".png"))
    print(a)
    dev.off()
    
    'svglite(filename = paste0("Clustermarker_integrated_",i,".svg"),
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
features_TOP10 <- unique(Clustering_3_top10_mouseMAIN_CLUSTERMARKER$gene)
#From Internet/head/FACS
features_FACS <- c("Alb","Clec4f", "Cd19","Cd3e", "Nkg7", "Cd14", "Lyz1","Lyz2","Timp1", "Col1a2", "Krt19", "Ly6g", "C1qa" , "Itgax", "Cd86") 
#FROM Scott Paper DEGs in mouse cells
features_DEG_Scott <-c("Mmrn2","Ntm","Fabp1","Apoa2","Ddit4l","Cd209a", "Nudt17","Ly6i")
features_Cytokine <-c("Il6","Il10","Il1ß","Tnfa")
NPC_CLUSTER<-SetIdent(NPC_CLUSTER,value = "mouseRNA.main")
DefaultAssay(object = NPC_CLUSTER)<-"RNA"
Create_Vplots(NPC_CLUSTER,"Malat1")
Create_Vplots(NPC_CLUSTER,features_TOP10)

Create_Vplots(NPC_CLUSTER,features_FACS)
Create_Vplots(NPC_CLUSTER,features_DEG_Scott)
##### Bis hier alles seit dem neuen Integrieren durchgespielt. Nur DotPLot kram macht fehler 19.07.24
# https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/
#https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1223471/full#supplementary-material
Idents(NPC_CLUSTER) <-NPC_CLUSTER$mouseRNA.main
NPC_CLUSTER <- subset(NPC_CLUSTER, mouseRNA.main %in% na.omit(NPC_CLUSTER$mouseRNA.main)) #removing cells which have a NA in mouseRNAmain, it gave error in DotPlot
 Idents(NPC_CLUSTER) <-factor(Idents(NPC_CLUSTER),
  levels=c("Hepatocytes","Macrophages","Endothelial cells","B cells","T cells", "NK cells",   "Monocytes", "Microglia","Granulocytes", "Fibroblasts","Erythrocytes","Dendritic cells", "Epithelial cells","Adipocytes"))
#non immune: MUC5A KRT5 SFTPD EPCAM CDH5 COL1A2 ACTA2 PECAM1 COL3A1 TMSB10 CALD1 FTH1 COL6A2 FKBP1A
#immune PTPRC CD45 CSF3R FCGR3B KLRD1 FPR1 CD8A CD1B TNFRSF17 BANK1 FCRL2 PNOC CR2 FCN1 GNLY KIR2DL1 KIR3DL1 KIR3DL2 CD1E CD1A CD163 AIF1 CD79A JCHAIN
markers.to.plot <-c( "Alb","Hpd","Fabp1","Apoa2","Saa1", "Saa2",
                     "Cd14","Clec4f","Lyz1","Cd68","C1qa","Adgre1","C1qc","Ly6i","Adgb", "Cd36" , "Cadm1",
                     "Ptprb","Nrp1","Pecam1", 
                     "Ptprc","Cd52","Tyrobp","Cd1e",
                     "Cd19","Cd79a","Cd79b",
                     "Cd3d","Cd3e","Cd3g","Cd2","Cd7","Il7r","Ccr7","Cd4","Cd8a",
                     "Skap1","Nkg7",
                     "Cd86","Lyz2","Fcer1g","S100a9",
                     "Ly6g","Dapp1",
                   "Col1a2","Rbms3","Dcn",
                   "Spp1","Ddit4l",  "Krt19","Timp1",
                   "Mmrn2")
DotPlot(NPC_CLUSTER, features = features_TOP10)+RotatedAxis()
DotPlot(NPC_CLUSTER, features = markers.to.plot)+RotatedAxis()
DoHeatmap(subset(NPC_CLUSTER, downsample = 100),features = features_TOP10, size = 3)

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
