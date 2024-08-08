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
BiocManager::install("ComplexHeatmap")
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

NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/4_NPC_ALL_TRANSFORM_Markers.rds")
############################### Create Subsets by Cluster #########################################
Idents(NPC_ALL_TRANSFORMED) <-"mouseRNA.main"
Tcells <-subset(NPC_ALL_TRANSFORMED,idents = c("T cells", "NK cells"))
Idents(Tcells)<-"sample"
animals = c("87","88","91","92")

T_list<-list()
for(a in animals){
  print(a)
  NPC <-subset(Tcells,idents = paste0("iAL",a))
  all.genes <- rownames(NPC)
  #generated clusters to check cellcycle scoring
  #in older scrit I did this and then plotted results here i only do it and dont plot results
  NPC <- ScaleData(NPC, features = all.genes, verbose = F)
  NPC <- RunPCA(NPC, features = VariableFeatures(object = NPC), npcs = 20, verbose = F)
  NPC <- RunUMAP(NPC, dims = 1:20, verbose = F, reduction = "pca")
  NPC <- FindNeighbors(NPC, dims = 1:20, k.param = 10, verbose = F)
  NPC <- FindClusters(NPC, resolution = 0.1, verbose = F)
  #scaling, norm, UMAP and Clustering per ,mouse
  # SCTransfrom might be beneficial bc it gices better signal to noise ratio. regression is performed with Mt5 and cell cylce Scores bc they introduce unwanted variation
  NPC <- SCTransform(NPC,  vst.flavor= "v2",method = "glmGamPoi",  verbose = F, vars.to.regress = c("percent.rb","percent.mt","S.Score","G2M.Score"))
  print(paste0(" I did SCTransformation for NPC_",a,"."))
  T_list<-append(T_list,NPC)
}
features <-SelectIntegrationFeatures(object.list = T_list, nfeatures = 3000)
T_list <- PrepSCTIntegration(object.list = T_list, anchor.features = features)
T.anchors <- FindIntegrationAnchors(object.list = T_list, dims = 1:20, normalization.method = "SCT", anchor.features = features)
Tcells_integrated <- IntegrateData(anchorset = T.anchors, dims = 1:20, normalization.method = "SCT")
Tcells_integrated <- RunPCA(Tcells_integrated, verbose = F)
Tcells_integrated <- RunUMAP(Tcells_integrated, dims = 1:30, verbose = F)
Tcells_integrated <- FindNeighbors(Tcells_integrated, dims = 1:30, verbose = F)
Tcells_integrated <- FindClusters(Tcells_integrated, verbose = F,resolution =2) 
#Muss ich irgendiwe machen ,weil sonst verschiedne layer da sind? und dasnn kann zb Heatmap nicht arbeiten.
DefaultAssay(Tcells_integrated) <- "RNA"
Tcells_integrated <- NormalizeData(Tcells_integrated)
Tcells_integrated <- FindVariableFeatures(Tcells_integrated, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tcells_integrated)
Tcells_integrated <- ScaleData(Tcells_integrated, features = all.genes)
Tcells_integrated <-JoinLayers(Tcells_integrated)
#rm(NPC, NPC.anchors, NPC_list)
#saveRDS(NPC_ALL_TRANSFORMED, "./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
#rm(NPC_ALL_TRANSFORMED)
#NPC_QC<-readRDS("./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
DefaultAssay(Tcells_integrated) <- "integrated"
DimPlot(Tcells_integrated, group.by = "seurat_clusters")
DimPlot(Tcells_integrated, group.by = "sample")
DimPlot(Tcells_integrated, group.by = "stim")
DimPlot(Tcells_integrated, group.by = "sex")
# #########################################################################

DefaultAssay(Tcells) <-"SCT"
Tcells_integrated <- PrepSCTFindMarkers(Tcells_integrated)
CD8Tcells1 <-c("Cd3e","Cd8a","Cd8b1","Tnf","Ifng","Il2","Cxcr3","Tbx21")
CD8Tcells2 <-c("Cd3e","Cd8a","Cd8b1","Cd28","Il4","Il5","Ccr4","Gata3")
CD8Tcells9 <-c("Cd3e","Cd8a","Cd8b1","Cd28","Il9","Il10","Irf4")
CD8Tcells17 <-c("Cd3e","Cd8a","Cd8b1","Cd28","Ccr6","Klrb1","Il17a","Irf4","Rorc")
#myTcellsSorting <-c("0","1","2","3","4","5","6","7")
#Idents(Tcells_integrated) <-factor(Idents(Tcells_integrated),levels=myTcellsSorting)
### Cd4 also very low in Scott Atlas
TcellGenes <-c("Cd3e","Cd3d","Cd28","Cd4","Cd8a","Cd8b1","Foxp3","Il2ra","Stat5a","Ctla4","Il10","Tgfb1","Il6","Il10","Il17a","Il2")
#Tcells2 from https://www.nature.com/articles/ncomms15081 might be interesting
TcellGenes2<-c("Cst7","Gzma","Gzmb","Ifng","Nkg7","Prf1","Tnfsf10","Btla","Ctla4","Havcr2","Lag3","Pdcd1","Tigit","Il2ra","Il4ra","Il7","Tgfb1","Tgfb2","Tgfb3","Tgfbi","Ccr7","Lef1","Sell","Tcf7","Cd226","Icos","Slamf1","Tnfrsf14","Tnfrsf25","Tnfrsf9")
ApoptosisGenes <-c("Bax","Fas","Bid","Bcl2","Psrc1","Ccng1","Sesn2","Eda2r","Pmaip1","Bbc3","Cdkn1a","Mdm2","Tnfrsf10b","Trp53","Rarg","Eme2","Jag2","Zmat3","Traf3","Phlda3","Sfn","Lhx3","Ifit3")
p<-DoHeatmap(Tcells_integrated, assay = "RNA", slot = "scale.data", features = unique(c(TcellGenes2)),
             draw.lines = T,lines.width = NULL,
             label = F, group.bar =T)+
  scale_fill_viridis_c()+
  theme(axis.text.y.left = element_text(size = 5),
        legend.justification = "top", 
        legend.title = element_text(size=6),
        legend.key.height= unit(2, 'mm'), 
        legend.key.width= unit(1, 'mm'),
        legend.text = element_text(size=6))
print(p)
ggsave( filename = paste0("./03_plots/2_Clustering/Clustermarker5_HeatMap.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )

TConservedMarkers20<-data.frame()

b<-unique(Tcells_integrated@meta.data$seurat_clusters)
for (c in b){
  markers <- FindConservedMarkers(Tcells_integrated, assay = "RNA", ident.1 = c, grouping.var = "sample",verbose = T)
  write.csv(markers,paste0("./99_other/2_Clustering/Subset_Tcells_ConservedMarkers_",c,".csv"))
  markers20<-data.frame(cluster=c,genes=rownames(markers)[1:5])
  TConservedMarkers20<-rbind(TConservedMarkers20,markers20)
  print(paste0("I saved conserved Markers for ",c,"."))
}
TConservedMarkers20 <-TConservedMarkers20%>%arrange(cluster)
p<-DoHeatmap(Tcells_integrated, assay = "RNA", slot = "scale.data", features = unique(c(TConservedMarkers20$genes,ApoptosisGenes)),
             draw.lines = T,lines.width = NULL,
             label = F, group.bar =T)+
  scale_fill_viridis_c()+
  theme(axis.text.y.left = element_text(size = 5),
        legend.justification = "top", 
        legend.title = element_text(size=6),
        legend.key.height= unit(2, 'mm'), 
        legend.key.width= unit(1, 'mm'),
        legend.text = element_text(size=6))
print(p)
Tcells_integrated@meta.data
###Check if percent rb gives variability to clusters. if so, regress it out
x <- Tcells_integrated@meta.data %>%
  ggplot(aes(x=seurat_clusters, fill=seurat_clusters, y =percent.rb)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Percent Rb")
print(x)