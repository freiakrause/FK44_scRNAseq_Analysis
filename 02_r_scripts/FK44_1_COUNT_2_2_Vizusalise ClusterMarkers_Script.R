#This is the script to analyse FK44.1 scRNAseq COUNT data provided by BSF; Data from 1_QC SCript is loaded

##################### Install Packages ############################################
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
#renv::install('metap')
#if (!require("BiocManager", quietly = TRUE))
#  renv::install("BiocManager")
#renv::install("gprofiler2")
#renv::install("DESeq2")
#renv::install("MAST")

#BiocManager::install(version = "3.18")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install('glmGamPoi')
#BiocManager::install("EnhancedVolcano")
#BiocManager::install('multtest')

#remotes::install_github("Moonerss/scrubletR")
#remotes::install_github('immunogenomics/presto')
#remotes::install_github("ianmoran11/mmtable2")

################################# Load Libraries #########################################################
library(Seurat)
library(ggplot2)
library(SingleR)
library(tidyr)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(scrubletR)
library(gprofiler2)
library(patchwork)
library(EnhancedVolcano)
library(readr)
library(mmtable2)
library(gt)
source("02_r_scripts/malat1_function.R")
source("02_r_scripts/VlnPlot_Function.R") 
set.seed(42)
#### Load Input Data ####
### Results from FK44.1_COUNT_2_0_TRANSFORM and CLUSTERING
### SoupX, QC, SCT, Integration, Clustering, Annotation, Removal of Clusters with low cell numbers
#NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced.rds)
#y<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())
#DefaultAssay(NPC_ALL_TRANSFORMED) <-"SCT"
#NPC_ALL_TRANSFORMED <- PrepSCTFindMarkers(NPC_ALL_TRANSFORMED)
#saveRDS(NPC_ALL_TRANSFORMED, file = "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced.rds"")
NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced.rds")

#### Define Cluster Sorting
myClusterSorting <-c("T cells","NK cells","B cells","Macrophages","Monocytes","Granulocytes","Fibroblasts","Endothelial cells","Hepatocytes")
myClusterSorting2 <-c("T cells_EtOH","T cells_TAM","NK cells_EtOH","NK cells_TAM","B cells_EtOH","B cells_TAM","Macrophages_EtOH","Macrophages_TAM",
                      "Monocytes_EtOH","Monocytes_TAM","Granulocytes_EtOH","Granulocytes_TAM",
                      "Fibroblasts_EtOH","Fibroblasts_TAM","Endothelial cells_EtOH","Endothelial cells_TAM","Hepatocytes_EtOH","Hepatocytes_TAM")
#### Define some interesting geens (manually)
Cytokines_and_Stuff <-c("Lyve1","Flt4","Efnb2","Ephb4","Icam1","Selp","F3",
                        "Itgax","Cd163","Msr1","Mrc1","Vegfa","Maf","Cxcl9","Cxcl10","Cxcl11","Stat6","Socs1",
                        "Cxcl1","Il6","Il10","Tgfb1","Ifng","Cxcr6","Il2","Saa1","Saa2","Cd40",
                        "Cd28","Cd86","Stat3","Socs3","Gzma","Gzmb","Prf1","Il4","Cxcl15","Ccl2","Tnf",
                        "Runx3","Cd8a","Cd4","Cd3e","Il21","Il23a","Il17a","Il17f","Il22","Rorc","Rora","Tbx21","Gata3","Foxp3","Il2ra","Eomes","Il1b","Ifng","Il12a",
                        "Col1a2","Col3a1","Fgg","Fga","Fgb","Fbln5","Apoa1","Fabp1","Gnmt","Selenbp2")
Leukocyte_Marker <-c("Ptprc","Cd52")#alle noch nicht überprüft
T_Marker <-c("Cd4","Cd8a","Il7r","Cd3d","Cd3e","Cd3g","Lat", "Lck")
NK_Marker <-c("Ccl5","Nkg7","Gzmb","Gzma","Ncr1","Prf1","Ccl4","Ccl3")
B_Marker <-c("Cd19","Cd79a","Ighm","Ighd","Fcmr","Ly6d","Ebf1","Ms4a1")#,"Cd86"
myeloid_Marker <-c("Itgax","Cd14","Fcgr3","Itgam")
Macro_Marker <-c("Cd74","Cyth4","Lyz2","Csf1r","Cd68","Aif1","Cybb","Ccl6")
KC_Marker <-c("Clec4f","C1qc","C1qa","C1qb","Cd163","Clec1b","Adgre1")
Mono_Marker <-c("Ly6c1","Ccr2","Itga2","Sell","Cx3cr1","Clec7a","Ifitm2","Cxcr4")
Neutro_Marker <-c("Mmp8","Mmp9","Ly6g","Hdc","Il1r2","Ccr1") #,"Hp","Lcn2"
Fibro_Marker <-c("Gsn","Egr1","Col6a2","Fstl1","Dcn","Mmp2","Lum")
HSC_Marker <-c("Col1a1","Col1a2","Col3a1","Igfbp3","Igfbp7","Bgn")
Endo_Marker <-c("Id3","Ptprb","Pecam1","Egfl7","Gng11","Flt1","Cldn5","Adgrf5","Eng")
Hep_Marker <-c("Ass1","Orm1","Apoa1","Apoa2","Alb","Aldh6a1","Ambp")
Canonical_ClusterMarker <-unique(c(Leukocyte_Marker,T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker, KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker,Endo_Marker,Hep_Marker))
#### Do HeatMap of manually assigned Canonical Markers ####
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
Idents(NPC_ALL_TRANSFORMED) <-factor(Idents(NPC_ALL_TRANSFORMED),levels=myClusterSorting)
p<-DoHeatmap(NPC_ALL_TRANSFORMED, assay = "RNA",slot = "scale.data", features = Canonical_ClusterMarker,
             draw.lines = T,lines.width = NULL,
             label = F, group.bar =T)+
  scale_fill_viridis_c()+
  theme(axis.text.y.left = element_text(size = 2.8),
        legend.justification = "top", 
        legend.title = element_text(size=6),
        legend.key.height= unit(0.2, 'cm'), 
        legend.key.width= unit(0.1, 'cm'),
        legend.text = element_text(size=4))
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_Canonnical_HeatMap.png"), p,width = 3.25, height = 3.25,dpi = 1200, bg="transparent")
#### Do DotPlot of manually assgined Canonical Markers
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
Idents(NPC_ALL_TRANSFORMED) <-factor(Idents(NPC_ALL_TRANSFORMED),levels=myClusterSorting)
p<-DotPlot(NPC_ALL_TRANSFORMED,  assay = "RNA", scale= T, features =unique(Canonical_ClusterMarker))+
  RotatedAxis()+
  scale_size(breaks = c(0, 25, 50, 75, 100),range(0,10))+
  scale_colour_distiller(palette="Blues", trans="reverse")+
  coord_cartesian( ylim=c(0,9.5),clip = "off")+
  annotate("text", y = 9.75, x = 0.5+length(unique(c(Leukocyte_Marker)))/2, label = "L", size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker)))-length(T_Marker)/2), label = "T",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker)))-length(NK_Marker)/2), label = "NK",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker)))-length(B_Marker)/2), label = "B",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker)))-length(myeloid_Marker)/2), label = "Mye",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker)))-length(Macro_Marker)/2), label = "Macro",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker)))-length(KC_Marker)/2), label = "KC",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker, Mono_Marker)))-length(Mono_Marker)/2), label = "Mono",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker)))-length(Neutro_Marker)/2), label = "N",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker)))-(length(HSC_Marker)+length(Fibro_Marker))/2), label = "Fibro & HSC",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker,Endo_Marker)))-length(Endo_Marker)/2), label = "Endo",size=4)+
  annotate("text", y = 9.75, x = 0.5+(length(unique(Canonical_ClusterMarker))-length(Hep_Marker)/2), label = "Hep",size=4)+
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
  ylab("Cell type")+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker)))+0.5)+
  geom_vline(linetype= "dotted",xintercept=length(unique(c(Leukocyte_Marker, T_Marker,NK_Marker,B_Marker,myeloid_Marker,Macro_Marker,KC_Marker,Mono_Marker,Neutro_Marker,Fibro_Marker,HSC_Marker,Endo_Marker)))+0.5)+
  geom_hline(yintercept=9.4)
  
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_DotPlot_Canonical_Markers.png"), p,width = 18, height = 4, dpi = 600,bg="transparent")





#### Test some potential interesting marker stuff #####
Idents(NPC_ALL_TRANSFORMED) <- "celltype.stim"
Idents(NPC_ALL_TRANSFORMED) <-factor(Idents(NPC_ALL_TRANSFORMED),levels=myClusterSorting2)
VlnPlot(NPC_ALL_TRANSFORMED,features = c("Col1a2"),assay = "RNA")

p<-DotPlot(NPC_ALL_TRANSFORMED,  assay = "RNA", features =unique(Cytokines_and_Stuff),cols=c("pink","green"))+
  RotatedAxis()+
  scale_size(breaks = c(0, 25, 50, 75, 100),range(0,10))+
  scale_colour_distiller(palette="Blues", trans="reverse")+
  guides(colour = guide_colourbar(reverse = TRUE))+
  theme(panel.background = element_rect(fill = "gray95",colour="black", linewidth = 1),
        axis.line.y.left =element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y.left = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8),
        legend.justification = "top", 
        legend.key.height= unit(0.4, 'cm'), 
        legend.key.width= unit(0.2, 'cm'),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        axis.line.x.bottom =element_blank(),
        axis.text.x =element_text(angle = 90,vjust = 0.5))+
  xlab("Marker genes")+
  ylab("Cell Type")
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_DotPlot_Cytokines Test.png"), p,width = 12, height = 3, dpi = 800,bg="transparent")

######################## Find Conserved Markers in Clusters across Stimulation ########
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
Idents(NPC_ALL_TRANSFORMED) <-factor(Idents(NPC_ALL_TRANSFORMED),levels=myClusterSorting)

ConservedMarkers5<-data.frame()
ConservedMarkers10<-data.frame()
ConservedMarkers20<-data.frame()

b<-unique(NPC_ALL_TRANSFORMED@meta.data$mouseRNA.main)
for (c in b){
  markers <- FindConservedMarkers(NPC_ALL_TRANSFORMED, assay = "RNA", ident.1 = c, grouping.var = "sample",verbose = T)
  write.csv(markers,paste0("./99_other/2_Clustering/ConservedMarkers_",c,".csv"))
  markers5<-data.frame(cluster=c,genes=rownames(markers)[1:5])
  ConservedMarkers5<-rbind(ConservedMarkers5,markers5)
  markers10<-data.frame(cluster=c,genes=rownames(markers)[1:10])
  ConservedMarkers10<-rbind(ConservedMarkers10,markers10)
  markers20<-data.frame(cluster=c,genes=rownames(markers)[1:20])
  ConservedMarkers20<-rbind(ConservedMarkers20,markers20)
  print(paste0("I saved conserved Markers for ",c,"."))
}
write.csv(ConservedMarkers5,paste0("./99_other/2_Clustering/ConservedMarkers_Top5.csv"))
write.csv(ConservedMarkers10,paste0("./99_other/2_Clustering/ConservedMarkers_Top10.csv"))
write.csv(ConservedMarkers20,paste0("./99_other/2_Clustering/ConservedMarkers_Top20.csv"))

ConservedMarkers5<- read.csv(paste0("./99_other/2_Clustering/ConservedMarkers_Top5.csv"))
ConservedMarkers10 <-read.csv(paste0("./99_other/2_Clustering/ConservedMarkers_Top10.csv"))
ConservedMarkers20 <-read.csv(paste0("./99_other/2_Clustering/ConservedMarkers_Top20.csv"))

ConservedMarkers1<-c("Cd3e","Gzma","Cd79a","Clec4f","Ank2","Havcr2","S100a8","Col3a1","Ptprb","Alb","Saa1")
ConservedMarkers5<-ConservedMarkers5[order(factor(ConservedMarkers5$cluster, levels=myClusterSorting)),]                        
ConservedMarkers10<-ConservedMarkers10[order(factor(ConservedMarkers10$cluster, levels=myClusterSorting)),]
ConservedMarkers20<-ConservedMarkers20[order(factor(ConservedMarkers20$cluster, levels=myClusterSorting)),]

#### Do HeatMap of Marker 5 Marker genes per Clusters ####
p<-DoHeatmap(NPC_ALL_TRANSFORMED, assay = "RNA", slot = "scale.data", features = c("Ptprc",ConservedMarkers5$gene),
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

#### Do HeatMap of Marker 10 Marker genes per Clusters ####
p<-DoHeatmap(NPC_ALL_TRANSFORMED, assay = "RNA", slot = "scale.data",features = c("Ptprc",ConservedMarkers10$gene),
             draw.lines = T,lines.width = NULL, label = F, group.bar =T)+
            scale_fill_viridis_c()+
            theme(axis.text.y.left = element_text(size = 2),
            legend.justification = "top", 
            legend.title = element_text(size=6),
            legend.key.height= unit(0.2, 'cm'), 
            legend.key.width= unit(0.1, 'cm'),
          legend.text = element_text(size=4))
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker10_HeatMap.png"),p,width = 3.25, height = 3.25,dpi = 1200, bg="transparent"  )

#### Do HeatMap of Marker 20 Marker genes per Clusters ####
p<-DoHeatmap(NPC_ALL_TRANSFORMED, assay = "RNA",slot = "scale.data", features = c(unique(ConservedMarkers20$gene)),
             draw.lines = T,lines.width = NULL,
             label = F, group.bar =T)+
  scale_fill_viridis_c()+
  theme(axis.text.y.left = element_text(size = 1),
        legend.justification = "top", 
        legend.title = element_text(size=6),
        legend.key.height= unit(0.2, 'cm'), 
        legend.key.width= unit(0.1, 'cm'),
        legend.text = element_text(size=4))
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker20_HeatMap.png"), p,width = 3.25, height = 3.25,dpi = 1200, bg="transparent")

#### Do DotPLot(BubblePlot) of 5 Marker genes per Clusters ####
p<-DotPlot(NPC_ALL_TRANSFORMED,  assay = "RNA",features =unique(ConservedMarkers5$gene))+
  RotatedAxis()+
  scale_size(breaks = c(0, 25, 50, 75, 100),range(0,10))+
  scale_colour_distiller(palette="Blues", trans="reverse")+
  guides(colour = guide_colourbar(reverse = TRUE))+
  theme(panel.background = element_rect(fill = "gray95",colour="black", linewidth = 1),
        axis.line.y.left =element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y.left = element_text(size = 8),
        axis.text.x.bottom = element_text(size = 8),
        legend.justification = "top", 
        legend.key.height= unit(0.4, 'cm'), 
        legend.key.width= unit(0.2, 'cm'),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        axis.line.x.bottom =element_blank(),
        axis.text.x =element_text(angle = 90,vjust = 0.5))+
        xlab("Marker genes")+
        ylab("Cell Type")
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_DotPlot5.png"), p,width = 12, height = 3, dpi = 800,bg="transparent")
#### Do DotPLot(BubblePlot) of 1 Marker genes per Clusters ####
p<-DotPlot(NPC_ALL_TRANSFORMED,  assay = "RNA",features =ConservedMarkers1)+
  RotatedAxis()+
  scale_size(breaks = c(0, 25, 50, 75, 100),range(0,10))+
  scale_colour_distiller(palette="Blues", trans="reverse")+
  guides(colour = guide_colourbar(reverse = TRUE))+
  theme(panel.background = element_rect(fill = "gray95",colour="black", linewidth = 1),
        axis.line.y.left =element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y.left = element_text(size = 10),
        axis.text.x.bottom = element_text(size = 10),
        legend.justification = "top", 
        legend.key.height= unit(0.6, 'cm'), 
        legend.key.width= unit(0.3, 'cm'),
        legend.title = element_text(size=9),
        legend.text = element_text(size=7),
        axis.line.x.bottom =element_blank(),
        axis.text.x =element_text(angle = 90,vjust = 0.5))+
  xlab("Marker genes")+
  ylab("Cell Type")
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_DotPlot1.png"), p,width = 8, height = 4, dpi = 800,bg="transparent")
#### Do Stacked VlnPlot of 1 Marker genes per Clusters ####
DefaultAssay(NPC_ALL_TRANSFORMED) <-"RNA"
p<-VlnPlot(NPC_ALL_TRANSFORMED, features = ConservedMarkers1, assay= "RNA", layer= "scale.data",log = T, stack = T,flip = F, fill.by = "ident")+
    theme_classic()+NoLegend()+theme(axis.title = element_blank())+
   labs(y="Cell Type")
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_VlnPlot1.png"), p,width = 8, height = 4, dpi = 800,bg="transparent")

#### Do Stacked VlnPlot ofSome special Fatty Genes####
Idents(NPC_ALL_TRANSFORMED) <- "celltype.stim"
Idents(NPC_ALL_TRANSFORMED) <-factor(Idents(NPC_ALL_TRANSFORMED),levels=myClusterSorting2)
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
Idents(NPC_ALL_TRANSFORMED) <-factor(Idents(NPC_ALL_TRANSFORMED),levels=myClusterSorting)

FattyGenes <-c("Saa1","Saa2","Cxcl1","Fabp1","Alb","Apoa1","Serpina1e","Fgg")
DefaultAssay(NPC_ALL_TRANSFORMED) <-"RNA"
p<-VlnPlot(NPC_ALL_TRANSFORMED, features = FattyGenes, assay= "RNA", layer= "scale.data",log = F, stack = T,flip = F, fill.by = "ident")+
  theme_classic()+NoLegend()+theme(axis.title = element_blank())+
  labs(y="Cell Type")
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_VlnPlot1.png"), p,width = 8, height = 4, dpi = 800,bg="transparent")

#### Do Stacked VlnPlot of 5 Marker per Cluster
DefaultAssay(NPC_ALL_TRANSFORMED) <-"RNA"
for (c in unique(NPC_ALL_TRANSFORMED$mouseRNA.main)){
  mm <- ConservedMarkers5%>%filter(cluster == c)
  p<-VlnPlot(NPC_ALL_TRANSFORMED, features = mm$genes, assay= "RNA", layer= "scale.data",log = T, stack = T,flip = F, fill.by = "ident")+
  theme_classic()+NoLegend()+theme(axis.title = element_blank())+
  labs(y="Cell Type")
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_VlnPlot",c,".png"), p,width = 8, height = 4, dpi = 800,bg="transparent")
}

#### Try to get Markers from Panglao DB to have independet clustermakers to proove my clusters
#Code from interent to get human geens to mouse gene conversion
mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
# separate human and mouse 
mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
# remove some columns
mouse <- mouse[,c(1,4)]
human <- human[,c(1,4)]
# merge the 2 dataset  (note that the human list is longer than the mouse one)
mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 
mh_data<-mh_data%>%arrange(desc(Symbol.y))
head(mh_data)
#Load Panglao Markers
PanglaoMarkers <- read.table("./99_other/PanglaoDB_markers_27_Mar_2020.tsv",h=T,sep="\t",quote="")
#Shrink Dataset to contain right species, necessary cell types and organs and some selection in sensitivity/canonical marker
PanglaoMarkers <-PanglaoMarkers%>%
  filter(species == "Mm" | species == "Mm Hs" & sensitivity_mouse > 0.01,cell.type !="Adipocytes",cell.type !="Chondrocytes",cell.type !="Plateletes",
                        !is.na(sensitivity_mouse ), !is.na(specificity_mouse), !is.na(canonical.marker),canonical.marker==1,
                        organ == "Liver" |organ == "Immune system" | organ == "Blood" |organ == "Connective tissue" | organ=="Epithelium"|
                        organ == "Vasculature" |organ == "Other")
#Get the human genes names fromPanlga transformed to mouseGene names
PanglaoMarkers[,"mouseGene"] <- NA
for (g in mh_data$Symbol.y){
  if (g %in% PanglaoMarkers$official.gene.symbol){
    pos_in_Pang<-grep(g,PanglaoMarkers$official.gene.symbol)
    pos_in_mh <-match(g,mh_data$Symbol.y)
    for(p in pos_in_Pang){
      PanglaoMarkers$mouseGene[pos_in_Pang]<-(mh_data$Symbol.x [pos_in_mh])
    }
    }
  }
#reorder columns
PanglaoMarkers%>%select(cell.type,mouseGene,official.gene.symbol,species,organ,sensitivity_mouse,specificity_mouse,ubiquitousness.index,nicknames,product.description,gene.type,germ.layer) 
#check if there a unassigned mosueGenes
nanana<-PanglaoMarkers%>%filter(is.na(mouseGene) )
#### need to find mouse gene names for the nananana s ###
write.csv(nanana,paste0("./99_other/2_Clustering/PanglaoDBMarkes_woMouseName.csv"))
#### added mouse gene names by hand #####
# convert the NA mosuegenes with the hand filled datatable to mouseGene names
Panglao_noNA <-read.csv("99_other/2_Clustering/PanglaoDBMarkers_manually_assigned.csv", sep = ";")
for (g in Panglao_noNA$official.gene.symbol){
  if (g %in% PanglaoMarkers$official.gene.symbol){
    pos_in_Pang<-grep(g,PanglaoMarkers$official.gene.symbol)
    pos_in_noNA <-match(g,Panglao_noNA$official.gene.symbol)
    for(p in pos_in_Pang){
      PanglaoMarkers$mouseGene[pos_in_Pang]<-(mh_data$Symbol.x [pos_in_noNA])
    }
  }
}
#Are there NAs left?
nonono<-PanglaoMarkers%>%filter(is.na(mouseGene) )
#Sort genes in Panlga bei sensitivity in mouse
PanglaoMarkers <-PanglaoMarkers%>%group_by(cell.type)%>%arrange((desc(sensitivity_mouse)),.by_group = T)

#PanglaoMarkers<-PanglaoMarkers%>%
#  mutate(mouseGene = ifelse(mouseGene %in% NPC_ALL_TRANSFORMED@assays$RNA@meta.data$var.features, mh_data$Symbol.x, NA))%>%
#  filter(!is.na(mouseGene))
for(c in unique(PanglaoMarkers$cell.type)){
  CP<-subset(PanglaoMarkers,PanglaoMarkers$cell.type==c)
  p<-DoHeatmap(NPC_ALL_TRANSFORMED, assay = "RNA",slot = "scale.data", features = c(unique(CP$mouseGene)),
               draw.lines = T,lines.width = NULL,
               label = F, group.bar =T)+
    scale_fill_viridis_c()+
    theme(axis.text.y.left = element_text(size = 4),
          legend.justification = "top", 
          legend.title = element_text(size=6),
          legend.key.height= unit(0.2, 'cm'), 
          legend.key.width= unit(0.1, 'cm'),
          legend.text = element_text(size=4))
  print(p)
  ggsave(filename = paste0("./03_plots/2_Clustering/ClustermarkerPanglao",c,".png"), p,width = 3.25, height = 3.25,dpi = 1200, bg="transparent")
  
}
##### Panglao Markers kinda specific but also not. need to have good metrics to select genes to represent "Canonical Cluster Markers"#
#### I looked at PanglaoMarkers and on top of this script I defined manually marker genes per cluster. Here I want to get the product descriptions for these defined markers
CanoMarkers_Table<-PanglaoMarkers%>%ungroup()%>%filter(mouseGene %in% Canonical_ClusterMarker)%>%
  select(mouseGene, official.gene.symbol, product.description,cell.type)%>%
  arrange(factor(mouseGene,levels=Canonical_ClusterMarker))%>%
filter(official.gene.symbol!="CD44" & official.gene.symbol!="CD40" & official.gene.symbol!="CD47" & official.gene.symbol!="CD48")%>%
select(mouseGene,product.description,official.gene.symbol)%>%mutate(Reference= "1")%>%
  rename( humanGene =official.gene.symbol,Protein = product.description)
CanoMM<-list(Leukocytes=Leukocyte_Marker,
             `T Cells`=T_Marker,
             `NK Cells`=NK_Marker,
             `B Cells`=B_Marker,
             `myeloid Cells`=myeloid_Marker,
             Macropahges=Macro_Marker,
             `Kupffer Cells`=KC_Marker,
             Monocytes=Mono_Marker,
             Neutrophils=Neutro_Marker,
             Fibroblasts=Fibro_Marker,
             `Hepatic Stellate Cells`= HSC_Marker,
             `Endothelial Cells`=Endo_Marker,
             Hepatocytes= Hep_Marker)
MarkerFrame <-data.frame()

for (e in 1:length(CanoMM)){
  print(e)
 print(CanoMM[[e]])
print(names(CanoMM[e])[1])
x<-data.frame(mouseGene=CanoMM[[e]], Marker.for =names(CanoMM[e])[1])
MarkerFrame <-rbind(MarkerFrame,x)
}
CanoMarkers_Table<-merge(MarkerFrame,CanoMarkers_Table, by = "mouseGene", all.x =T)%>%arrange(factor(mouseGene,levels=Canonical_ClusterMarker))%>%distinct(mouseGene,.keep_all = T)
#### Generate nice table as output for canonical Markers
CanoMarkers_Table%>%rename(`Murine Gene`=mouseGene,`Human Gene`=humanGene,`Marker for`=Marker.for)%>%
  gt()%>%
  tab_header(title = "Canonical Cluster Markers")%>%
  tab_style(style = cell_text(color = "black", weight = "bold", align = "center"),locations = cells_title("title"))%>%
  tab_style(style= cell_text(weight= "bold"),locations = cells_column_labels(everything()))%>%
   opt_table_lines()%>%
  tab_options(table.font.size="small")%>%
  gtsave("Canonical_Markers_Table.docx",path="./03_plots/2_Clustering/")



##################### Plotting Things ####
#NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "mouseRNA.fine")
Idents(NPC_ALL_TRANSFORMED) <-"mouseRNA.main"
DefaultAssay(NPC_ALL_TRANSFORMED) <- "RNA"
Hep_genes <-c("Saa1", "Saa2","Alb","mt-Cytb","mt-Co1", "2cdnaFLAG-LGP130", "1cdnaKANA","3cdnaZSGREEN","mt-Atp6")
T_genes <-c("Cd3e","Cd3d", "Cd4", "Cd8a") 
B_genes <-c("Cd19")
M_genes <-c("Clec4f")
Gene_List <-list(Hep_genes,T_genes,B_genes,M_genes)
for (GL in Gene_List){
  for (g in GL){
    p<-FeaturePlot(NPC_ALL_TRANSFORMED, features = g, split.by = "stim", max.cutoff = 5,cols = c("grey", "red"))
    png(paste0("./03_plots/2_Clustering/FeaturePlot_",g,".png"))
    print(p)
    print(paste0("I did FeaturePlot of ",g,". Now I will do VlnPlot of ",g,"."))
    dev.off()
    p<-VlnPlot(NPC_ALL_TRANSFORMED, features = g, split.by = "stim",group.by = "mouseRNA.main", pt.size = 0, combine = FALSE)
    png(paste0("./03_plots/2_Clustering/VlnPlot_",g,".png"))
    print(p)
    dev.off()
  }
}
