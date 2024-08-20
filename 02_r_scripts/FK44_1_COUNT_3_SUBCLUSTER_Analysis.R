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
#BiocManager::install("ComplexHeatmap")
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
library(scran)
library(enrichR)
set.seed(42)

# NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/4_NPC_ALL_TRANSFORM_Markers.rds")
# ############################### Create Subsets by Cluster #########################################
# ####https://www.nature.com/articles/s41598-020-76972-9 : Article states that Celltype classification accuracy is low with low cell numbers. With T cell dataset they calim, you need ~2500 cells for good accuracy
# Idents(NPC_ALL_TRANSFORMED) <-"mouseRNA.main"
# Tcells <-subset(NPC_ALL_TRANSFORMED,idents = c("T cells", "NK cells"))
# Idents(Tcells)<-"sample"
# animals = c("87","88","91","92")
# DefaultAssay(Tcells)<-"RNA"
# T_list<-list()
# for(a in animals){
#   print(a)
#   NPC <-subset(Tcells,idents = paste0("iAL",a))
#   all.genes <- rownames(NPC)
#   #generated clusters to check cellcycle scoring
#   NPC <- ScaleData(NPC, features = all.genes, verbose = T)
#   NPC <- RunPCA(NPC, features = VariableFeatures(object = NPC), npcs = 50, verbose = F)
#   p <-DimHeatmap(NPC, dims = 1:10, nfeatures = 20, cells = 500, balanced = T)
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_DimHeatMap_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   p <-ElbowPlot(NPC)&labs(title = paste0("iAL",a," ElbowPLot PCA in T cell Subclusters"))
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_ElbowPlot_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   p <-DimPlot(NPC, reduction = "pca")&labs(title = paste0("iAL",a," DimPlot PCA T cell Subclusters"))
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_DimPlot_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   NPC <- RunUMAP(NPC, dims = 1:10, verbose = F, reduction = "pca")
#   NPC <- FindNeighbors(NPC, dims = 1:10, k.param = 10, verbose = T)
#   NPC <- FindClusters(NPC, resolution = 1, verbose =T)
#   p <-FeaturePlot(NPC,features = "percent.mt",label.size = 4,repel = T,label = T) +   theme(plot.title = element_text(size=5))
#   p <-VlnPlot(NPC,features = "percent.mt") & theme(plot.title = element_text(size=10))
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_MtPerc_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   p <-VlnPlot(NPC,features = c("nFeature_RNA")) & theme(plot.title = element_text(size=10))
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_FeatureRNA_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   p <-VlnPlot(NPC,features = c("nCount_RNA")) & theme(plot.title = element_text(size=10))
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_CountRNA_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   p <-VlnPlot(NPC,features = "G2M.Score") &   theme(plot.title = element_text(size=10))
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_G2Mscore_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   p <-VlnPlot(NPC,features = "S.Score") &   theme(plot.title = element_text(size=10))
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_Sscore_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   p <-VlnPlot(NPC,features = "percent.rb") &   theme(plot.title = element_text(size=10))
#   print(p)
#   ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_Rbperc_iAL",a,".png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
#   #scaling, norm, UMAP and Clustering per ,mouse
#   # SCTransfrom might be beneficial bc it gices better signal to noise ratio. regression is performed with Mt5 and cell cylce Scores bc they introduce unwanted variation
#   NPC <- SCTransform(NPC,  vst.flavor= "v2",method = "glmGamPoi",  verbose = F, vars.to.regress = c("percent.rb","percent.mt","S.Score","G2M.Score"))
#   print(paste0(" I did SCTransformation for NPC_",a,"."))
#   T_list<-append(T_list,NPC)
# }
# rm(NPC)
# features <-SelectIntegrationFeatures(object.list = T_list, nfeatures = 3000)
# T_list <- PrepSCTIntegration(object.list = T_list, anchor.features = features)
# T.anchors <- FindIntegrationAnchors(object.list = T_list, dims = 1:10, normalization.method = "SCT", anchor.features = features)
# Tcells_integrated <- IntegrateData(anchorset = T.anchors, dims = 1:10, normalization.method = "SCT")
# Tcells_integrated <- RunPCA(Tcells_integrated, verbose = F)
# p <-ElbowPlot(Tcells_integrated)&labs(title = paste0("ElbowPLot PCA in T cell Subclusters"))
# print(p)
# ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_ElbowPlot_FeatureRNA.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
# 
# Tcells_integrated <- RunUMAP(Tcells_integrated, dims = 1:30, verbose = F)
# Tcells_integrated <- FindNeighbors(Tcells_integrated, dims = 1:30, verbose = F)
# Tcells_integrated <- FindClusters(Tcells_integrated, verbose = F,resolution =1) 
# #Muss ich irgendiwe machen ,weil sonst verschiedne layer da sind? und dasnn kann zb Heatmap nicht arbeiten.
# DefaultAssay(Tcells_integrated) <- "RNA"
# Tcells_integrated <- NormalizeData(Tcells_integrated)
# Tcells_integrated <- FindVariableFeatures(Tcells_integrated, selection.method = "vst", nfeatures = 4000)
# all.genes <- rownames(Tcells_integrated)
# Tcells_integrated <- ScaleData(Tcells_integrated, features = all.genes)
# Tcells_integrated <-JoinLayers(Tcells_integrated)
# Tcells_integrated$clusters.stim <- paste(Tcells_integrated$seurat_clusters, Tcells_integrated$stim, sep = "_")
# DefaultAssay(Tcells_integrated) <-"SCT"
# Tcells_integrated <- PrepSCTFindMarkers(Tcells_integrated)
# rm(Tcells, T.anchors, T_list,NPC_ALL_TRANSFORMED)
#saveRDS(Tcells_integrated, "./01_tidy_data/3_Subclustering_Tcells.rds")
#rm(Tcells_integrated)
Tcells_integrated<-readRDS("./01_tidy_data/3_Subclustering_Tcells.rds")
#### Vizuals
p <-VlnPlot(Tcells_integrated,features = "percent.mt") & theme(plot.title = element_text(size=10))
print(p)
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_MtPerc.png"),  p,  width = 3.25,  height = 3.25,  dpi = 600,  bg="transparent"  )
p <-VlnPlot(Tcells_integrated,features = c("nFeature_RNA")) & theme(plot.title = element_text(size=10))
print(p)
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_FeatureRNA.png"),  p,  width = 3.25,  height = 3.25,  dpi = 600,  bg="transparent"  )
p <-VlnPlot(Tcells_integrated,features = c("nCount_RNA")) & theme(plot.title = element_text(size=10))
print(p)
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_CountRNA.png"),  p,  width = 3.25,  height = 3.25,  dpi = 600,  bg="transparent"  )
p <-VlnPlot(Tcells_integrated,features = "G2M.Score") &   theme(plot.title = element_text(size=10))
print(p)
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_G2Mscore.png"),  p,  width = 3.25,  height = 3.25,  dpi = 600,  bg="transparent"  )
p <-VlnPlot(Tcells_integrated,features = "S.Score") &   theme(plot.title = element_text(size=10))
print(p)
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_Sscore.png"),  p,  width = 3.25,  height = 3.25,  dpi = 600,  bg="transparent"  )
p <-VlnPlot(Tcells_integrated,features = "percent.rb") &   theme(plot.title = element_text(size=10))
print(p)
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_VlnPlot_Rbperc_.png"),  p,  width = 3.25,  height = 3.25,  dpi = 600,  bg="transparent"  )
#scaling, norm, UMAP and Clustering per ,mouse

DefaultAssay(Tcells_integrated) <- "integrated"
DimPlot(Tcells_integrated, group.by = "seurat_clusters")
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_DimPlot_SeuratClusters_.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
DimPlot(Tcells_integrated, group.by = "sample")
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_DimPlot_Sample_.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
DimPlot(Tcells_integrated, group.by = "stim")
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_DimPlot_Stim_.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )
DimPlot(Tcells_integrated, group.by = "sex")
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_DimPlot_Sex_.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )

# #########################################################################

#### Visualize manually assigned potential amrkers to identify tcells clusters
#### Define Cluster Sorting

#### Define some interesting geens (manually)
Cytokines_and_Stuff <-c("Lyve1","Flt4","Efnb2","Ephb4","Icam1","Selp","F3",
                        "Itgax","Cd163","Msr1","Mrc1","Vegfa","Maf","Cxcl9","Cxcl10","Cxcl11","Stat6","Socs1",
                        "Cxcl1","Il6","Il10","Tgfb1","Ifng","Cxcr6","Il2","Saa1","Saa2","Cd40",
                        "Cd28","Cd86","Stat3","Socs3","Gzma","Gzmb","Prf1","Il4","Cxcl15","Ccl2","Tnf",
                        "Runx3","Cd8a","Cd4","Cd3e","Il21","Il23a","Il17a","Il17f","Il22","Rorc","Rora","Tbx21","Gata3","Foxp3","Il2ra","Eomes","Il1b","Ifng","Il12a")
CD8Tcells1 <-c("Cd3e","Cd8a","Cd8b1","Tnf","Ifng","Il2","Cxcr3","Tbx21")
CD8Tcells2 <-c("Cd3e","Cd8a","Cd8b1","Cd28","Il4","Il5","Ccr4","Gata3")
CD8Tcells9 <-c("Cd3e","Cd8a","Cd8b1","Cd28","Il9","Il10","Irf4")
CD8Tcells17 <-c("Cd3e","Cd8a","Cd8b1","Cd28","Ccr6","Klrb1","Il17a","Irf4","Rorc")
myTcellsSorting <-c("0","1","2","3","4","5","6","7")
myTcellsSorting2 <-c("0_EtOH","0_TAM","1_EtOH","1_TAM","2_EtOH","2_TAM","3_EtOH","3_TAM","4_EtOH","4_TAM","5_EtOH","5_TAM","6_EtOH","6_TAM","7_EtOH","7_TAM")

#Idents(Tcells_integrated) <-factor(Idents(Tcells_integrated),levels=myTcellsSorting)
### 
TcellGenes <-c("Cd3e","Cd3d","Cd28","Cd4","Cd8a","Cd8b1","Foxp3","Il2ra","Stat5a","Ctla4","Il10","Tgfb1","Il6","Il10","Il17a","Il2")
#Tcells2 from https://www.nature.com/articles/ncomms15081 might be interesting
TcellGenes2<-c("Cst7","Gzma","Gzmb","Ifng","Nkg7","Prf1","Tnfsf10","Btla","Ctla4","Havcr2","Lag3","Pdcd1","Tigit","Il2ra","Il4ra","Il7","Tgfb1","Tgfb2","Tgfb3","Tgfbi","Ccr7","Lef1","Sell","Tcf7","Cd226","Icos","Slamf1","Tnfrsf14","Tnfrsf25","Tnfrsf9")
ApoptosisGenes <-c("Bax","Fas","Bid","Bcl2","Psrc1","Ccng1","Sesn2","Eda2r","Pmaip1","Bbc3","Cdkn1a","Mdm2","Tnfrsf10b","Trp53","Rarg","Eme2","Jag2","Zmat3","Traf3","Phlda3","Sfn","Lhx3","Ifit3")
Idents(Tcells_integrated)<-"clusters.stim"
Idents(Tcells_integrated)<-factor(Idents(Tcells_integrated),levels=myTcellsSorting2)

p<-DoHeatmap(Tcells_integrated, assay = "RNA", slot = "scale.data", features = unique(Cytokines_and_Stuff),
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
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_HeatMap_ClusterMarker_Tcell2Genes.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )

p<-DotPlot(Tcells_integrated,  assay = "RNA", scale= T, features =unique(Cytokines_and_Stuff))+
  RotatedAxis()+
  scale_size(breaks = c(0, 25, 50, 75, 100),range(0,14))+
  scale_colour_distiller(palette="Blues", trans="reverse")+
  coord_cartesian( ylim=c(0,16),clip = "off")+
 # annotate("text", y = 9.75, x = 0.5+length(unique(c(Leukocyte_Marker)))/2, label = "L", size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker)))-length(T_Marker)/2), label = "T",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker)))-length(NK_Marker)/2), label = "NK",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker)))-length(B_Marker)/2), label = "B",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker)))-length(myeloid_Marker)/2), label = "Mye",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker)))-length(Macro_Marker)/2), label = "Macro",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker)))-length(KC_Marker)/2), label = "KC",size=4)+
 # annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker, Mono_Marker)))-length(Mono_Marker)/2), label = "Mono",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker)))-length(Neutro_Marker)/2), label = "N",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker)))-(length(HSC_Marker)+length(Fibro_Marker))/2), label = "Fibro & HSC",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker,Endo_Marker)))-length(Endo_Marker)/2), label = "Endo",size=4)+
#  annotate("text", y = 9.75, x = 0.5+(length(unique(Canonical_ClusterMarker))-length(Hep_Marker)/2), label = "Hep",size=4)+
  guides(colour = guide_colourbar(reverse = TRUE))+
  theme(panel.background = element_rect(fill = "gray95", linewidth = 0.8,color = "black"),
        axis.line.y.left =element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y.left = element_text(size = 9),
        axis.text.x.bottom = element_text(size = 8, angle = 90,vjust = 0.5),
        legend.justification = "top", 
        legend.key.height= unit(0.4, 'cm'), 
        legend.key.width= unit(0.2, 'cm'),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        axis.line.x.bottom =element_blank(),
        axis.text.x.top =element_text(angle =0,vjust = 0.5))+
  xlab("Canonical marker genes")+
  ylab("Cell type")
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker)))+0.5)+
  #geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker,Endo_Marker)))+0.5)+
  #geom_hline(yintercept=9.4)
print(p)


#### Caclucalte Conseverd Markers
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
  theme(axis.text.y.left = element_text(size = 4),
        legend.justification = "top", 
        legend.title = element_text(size=6),
        legend.key.height= unit(2, 'mm'), 
        legend.key.width= unit(1, 'mm'),
        legend.text = element_text(size=6))
print(p)
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_HeatMap_ClusterMarker_Conserved_Apopotosis.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )

###Check if percent rb gives variability to clusters. if so, regress it out
x <- Tcells_integrated@meta.data %>%
  ggplot(aes(x=seurat_clusters, fill=seurat_clusters, y =percent.mt)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Percent Rb")
print(x)
#### Caclulate AllMarkers to potentially get Subcluster Marker
AllMarkersW <-FindAllMarkers(Tcells_integrated,test.use = "wilcox")
write.csv(AllMarkersW,paste0("./99_other/3_Subclustering/Subcluster_T_AllMarkers_Wilcox.csv"))
AllMarkersL <-FindAllMarkers(Tcells_integrated,test.use = "wilcox_limma")
write.csv(AllMarkersL,paste0("./99_other/3_Subclustering/Subcluster_T_AllMarkersLimma.csv"))
AllMarkersB <-FindAllMarkers(Tcells_integrated,test.use = "bimod")
write.csv(AllMarkersB,paste0("./99_other/3_Subclustering/Subcluster_T_AllMarkers_bimod.csv"))
AllMarkersM <-FindAllMarkers(Tcells_integrated,test.use = "MAST")
write.csv(AllMarkersM,paste0("./99_other//3_Subclustering/Subcluster_T_AllMarkers_MAST.csv"))
AllMarkersM <-AllMarkersM%>%group_by(cluster)%>%slice_head(n=20)
#AllMarkersB <-read.csv(paste0("./99_other/3_Subclustering/Subcluster_T_AllMarkers_bimod.csv"))
#AllMarkersM <-read.csv(paste0("./99_other//3_Subclustering/Subcluster_T_AllMarkers_MAST.csv"))
#AllMarkersL <-read.csv(paste0("./99_other/3_Subclustering/Subcluster_T_AllMarkersLimma.csv"))
#AllMarkersW <-read.csv(paste0("./99_other/3_Subclustering/Subcluster_T_AllMarkers_Wilcox.csv"))

p<-DoHeatmap(Tcells_integrated, assay = "RNA", slot = "scale.data", features = unique(c("Ptprc","Cd3e","Cd3d","Cd8a","Cd4","Cd28",AllMarkersM)),
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
ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_HeatMap_ClusterMarker_AllMarkersB.png"),  p,  width = 3.25,  height = 3.25,  dpi = 1200,  bg="transparent"  )

####try to fin enriched gensets to determine t cell subclusters
dbs <- listEnrichrDbs()  
dbs_selection <- unique(dbs$libraryName)[grep("Mouse",unique(dbs$libraryName))]

for(c in Tcells_integrated$seurat_clusters){
  for (i in 1:length(dbs_selection[1:3])) {
    database<-dbs_selection[i]
    print(database)
    p<-DEenrichRPlot(Tcells_integrated,ident.1= c,test.use = "MAST",enrich.database = database, max.genes = 1000, verbose =T)+
      print(p)
    ggsave( filename = paste0("./03_plots/3_Subclustering/Subcluster_T_DEnrichPlot_Cluster_",c,"_",database,".png"),  p,  width = 10,  height = 5,  dpi = 300,  bg="white"  )
      }
}
########
#####https://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html#assigning-cluster-labels-from-markers
