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

################################# Load Libraries #########################################################
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(scrubletR)
library(gprofiler2)
library(patchwork)
library(EnhancedVolcano)
source("02_r_scripts/malat1_function.R")
set.seed(42)

#### Load Input Data ####
### Results from FK44.1_COUNT_2_0
### SoupX, QC, SCT, Integration, Clustering
NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
y<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())
#Remove Clusters with low cell numbers
NPC_ALL_TRANSFORMED<-subset(NPC_ALL_TRANSFORMED, subset=mouseRNA.main!= "Adipocytes"&mouseRNA.main!= "Epithelial cells"&mouseRNA.main!= "NA"&mouseRNA.main!= "Erythrocytes"&mouseRNA.main!= "Dendritic cells")
z<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())
#### Save and vizualizes TOP10 ClusterMarker #####
images <-list()
DefaultAssay(NPC_ALL_TRANSFORMED) <-"SCT"
NPC_ALL_TRANSFORMED <- PrepSCTFindMarkers(NPC_ALL_TRANSFORMED)
#Order Clusters in my way
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
Idents(NPC_ALL_TRANSFORMED) <-factor(Idents(NPC_ALL_TRANSFORMED),
                                     levels=c("T cells", "NK cells","B cells","Macrophages","Microglia","Monocytes","Granulocytes", 
                                              "Fibroblasts","Endothelial cells","Hepatocytes"))
######################## Find Conserved Markers inClusters across Stimulation ########

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
ConservedMarkers1<-c("Cd3e","Gzma","Cd79a","Clec4f","Ank2","Havcr2","S100a8","Col3a1","Ptprb","Spp1")
ConservedMarkers5<-ConservedMarkers5[order(factor(ConservedMarkers5$cluster, levels=c("T cells", "NK cells","B cells","Macrophages", "Microglia","Monocytes","Granulocytes", 
                                                                                      "Fibroblasts","Endothelial cells","Hepatocytes"))),]                        
ConservedMarkers10<-ConservedMarkers10[order(factor(ConservedMarkers10$cluster, levels=c("T cells", "NK cells","B cells","Macrophages", "Microglia","Monocytes","Granulocytes", 
                                                                                         "Fibroblasts","Endothelial cells","Hepatocytes"))),]
ConservedMarkers20<-ConservedMarkers20[order(factor(ConservedMarkers20$cluster, levels=c("T cells", "NK cells","B cells","Macrophages", "Microglia","Monocytes","Granulocytes", 
                                                                                         "Fibroblasts","Endothelial cells","Hepatocytes"))),]

#### Do HeatMap of Marker 5 Marker genes per Clusters ####

p<-DoHeatmap(NPC_ALL_TRANSFORMED, assay = "RNA", 
             slot = "scale.data",
             features = c("Ptprc",ConservedMarkers5$gene),
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
p<-DoHeatmap(NPC_ALL_TRANSFORMED, assay = "RNA", 
             slot = "scale.data",
             features = c(unique(ConservedMarkers20$gene)),
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
#### Do DotPLot(BubblePlot) of Marker Cleanest? Marker genes per Clusters ####
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


DefaultAssay(NPC_ALL_TRANSFORMED) <-"RNA"
p<-VlnPlot(NPC_ALL_TRANSFORMED, features = ConservedMarkers1, assay= "RNA", layer= "scale.data",log = T, stack = T,flip = F, fill.by = "ident")+
    theme_classic()+NoLegend()+theme(axis.title = element_blank())+
   labs(y="Cell Type")
print(p)
ggsave(filename = paste0("./03_plots/2_Clustering/Clustermarker_VlnPlot1.png"), p,width = 8, height = 4, dpi = 800,bg="transparent")


Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.fine"
all.markers_fine <- FindAllMarkers(NPC_ALL_TRANSFORMED, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
all.markers_fine %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 5) %>%  ungroup() -> top5__MOUSEFINE_CLUSTERMARKER
write.csv(top5__MOUSEMAIN_CLUSTERMARKER, "./99_other/Clustering_3_top5_mouseMAIN_CLUSTERMARKER.csv", row.names=FALSE)
write.csv(top5__MOUSEFINE_CLUSTERMARKER, "./99_other/Clustering_3_top5_mouseFINE_CLUSTERMARKER.csv", row.names=FALSE)
saveRDS(NPC_ALL_TRANSFORMED, file = "./01_tidy_data/4_NPC_ALL_TRANSFORM_Markers.rds")


#### Load Input Data #####
NPC_CLUSTER <-readRDS(file =  "./01_tidy_data/4_NPC_ALL_TRANSFORM_Markers.rds")
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





#### Load RDS with Cluster Markers found in SCT assay  #####
NPC_ALL_TRANSFORMED <- readRDS("./01_tidy_data/4_NPC_ALL_TRANSFORM_Markers.rds")
########################## Identify differentially expressed Genes in Clusters across Conditions Enhanced Volcano ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
NPC_ALL_TRANSFORMED<-subset(NPC_ALL_TRANSFORMED, subset=mouseRNA.main!= "Adipocytes"&mouseRNA.main!= "Epithelial cells"&mouseRNA.main!= "NA"&mouseRNA.main!= "Erythrocytes")
NPC_ALL_TRANSFORMED$celltype.stim <- paste(NPC_ALL_TRANSFORMED$mouseRNA.main, NPC_ALL_TRANSFORMED$stim, sep = "_")
NPC_ALL_TRANSFORMED$celltype.sex <- paste(NPC_ALL_TRANSFORMED$mouseRNA.main, NPC_ALL_TRANSFORMED$sex, sep = "_")
Idents(NPC_ALL_TRANSFORMED) <- "celltype.stim"

NPC_ALL_TRANSFORMED <- PrepSCTFindMarkers(NPC_ALL_TRANSFORMED)
a<-unique(NPC_ALL_TRANSFORMED@meta.data$mouseRNA.main)
y<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())
for (c in a){
  single_l.de <-FindMarkers(NPC_ALL_TRANSFORMED, assay = NULL, ident.2 = paste0(c,"_EtOH"), ident.1 = paste0(c,"_TAM"), 
                            verbose = T, recorrect_umi = FALSE, min.cells.feature = 3, min.pct= 0.2,
                            test.use="wilcox_limma")
  single_M.de <-FindMarkers(NPC_ALL_TRANSFORMED, assay = NULL, ident.2 = paste0(c,"_EtOH"), ident.1 = paste0(c,"_TAM"), 
                            verbose = T,recorrect_umi = FALSE, min.cells.feature = 3, min.pct= 0.2,
                            test.use="MAST")
  write.csv(single_l.de,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_single_limma",c,".csv"))
  write.csv(single_M.de,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_single_MAST",c,".csv"))
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/EnhancedVolcano_Single_Limma",c,".png"))
  p<-  EnhancedVolcano(single_l.de,lab = rownames(single_l.de), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj", pCutoff = 1e-05, FCcutoff = 1.0 ,
                       title = paste0("DE ",c," TAM vs EtOH"), 
                       caption = 'FC cutoff, 1.0; p-value cutoff: p_val_adj<1e-05')
  print(p)
  dev.off()
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/EnhancedVolcano_Single_MAST",c,".png"))
  p<-  EnhancedVolcano(single_M.de,lab = rownames(single_M.de), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj", pCutoff = 1e-05, FCcutoff = 1.0 ,
                       title = paste0("DE ",c," TAM vs EtOH"), 
                       caption = 'FC cutoff, 1.0; p-value cutoff: p_val_adj<1e-05')
  print(p)
  dev.off()
  print(paste0("I just saved Enhanced Volcano of singeDE Analysis of ",c,"."))
  
  names(single_l.de) <- paste0(names(single_l.de), "_l.sc")
  single_l.de$gene <- rownames(single_l.de)
  names(single_M.de) <- paste0(names(single_M.de), "_M.sc")
  single_M.de$gene <- rownames(single_M.de)
  merge_dat <- merge(single_M.de, single_l.de,by = "gene")
  merge_dat <- merge_dat[order(merge_dat$p_val_M.sc), ]
  # Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
  common <- merge_dat$gene[which(merge_dat$p_val_l.sc < 0.05&
                                   merge_dat$p_val_M.sc < 0.05)]
  
  only_sc_l <- merge_dat$gene[which(merge_dat$p_val_M.sc > 0.05 &
                                      merge_dat$p_val_l.sc < 0.05)]
  
  only_sc_M <- merge_dat$gene[which(  merge_dat$p_val_l.sc > 0.05 &
                                        merge_dat$p_val_M.sc < 0.05)]
  print(paste0('# Common: ',length(common)))
  print(merge_dat[merge_dat$gene%in%common[1:10],c('gene','p_val_M.sc','p_val_l.sc')])
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Vln_Common_Limma_MAST",c,".png"))
  p<-VlnPlot(NPC_ALL_TRANSFORMED, features = common[1:16], idents = c(paste0(c,"_EtOH"), paste0(c,"_TAM")), group.by = "stim")
  print(p)
  dev.off()
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Vln_Only_Limma",c,".png"))
  p<-VlnPlot(NPC_ALL_TRANSFORMED, features = only_sc_M[1:16], idents = c(paste0(c,"_EtOH"), paste0(c,"_TAM")), group.by = "stim")
  print(p)
  dev.off()
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Vln_Only_MAST",c,".png"))
  p<-VlnPlot(NPC_ALL_TRANSFORMED, features = only_sc_l[1:16], idents = c(paste0(c,"_EtOH"), paste0(c,"_TAM")), group.by = "stim")
  print(p)
  dev.off()
}

########################## Identify differentially expressed Genes in Clusters across in Females Enhanced Volcano ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
NPC_female <-subset(NPC_ALL_TRANSFORMED, subset = sex == "female")
Idents(NPC_ALL_TRANSFORMED) <- "celltype.stim"
NPC_female <- PrepSCTFindMarkers(NPC_female)
a<-unique(NPC_female@meta.data$mouseRNA.main)
NPC_female@meta.data%>%group_by(mouseRNA.main)%>%summarise(n=n())
b<- list()
for (c in a){
  x <-subset(NPC_female, idents = c(paste0(c,"_EtOH"), paste0(c,"_TAM")))
  x <-FindMarkers(x, assay = "SCT", ident.2 = paste0(c,"_EtOH"), ident.1 = paste0(c,"_TAM"), verbose = FALSE, recorrect_umi = FALSE)
  write.csv(x,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_female",c,".csv"))
  p<-  EnhancedVolcano(x,lab = rownames(x), 
                       x = "avg_log2FC", y = "p_val", 
                       pCutoffCol = "p_val_adj" ,FCcutoff = 1.0, pCutoff = 1e-05,
                       title = paste0("DE ",c," TAM vs EtOH only females"), 
                       caption = 'FC cutoff: 1.0; p-value cutoff: p-val adj <1e-05')
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Volcano_Main_female",c,".png"))
  print(p)
  dev.off()
  print(paste0("I just saved Enhanced Volcano of ",c," for female samples."))
  
}

########################## Identify differentially expressed Genes in Clusters across in male Enhanced Volcano ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
NPC_male <-subset(NPC_ALL_TRANSFORMED, subset = sex == "male")
Idents(NPC_ALL_TRANSFORMED) <- "celltype.stim"
NPC_male <- PrepSCTFindMarkers(NPC_male)
a<-unique(NPC_male@meta.data$mouseRNA.main)
NPC_male@meta.data%>%group_by(mouseRNA.main)%>%summarise(n=n())
a<-subset(a, subset=a!= "Dendritic cells")
b<- list()
for (c in a){
  x <-subset(NPC_male, idents = c(paste0(c,"_EtOH"), paste0(c,"_TAM")))
  x <-FindMarkers(x, assay = "SCT", ident.2 = paste0(c,"_EtOH"), ident.1 = paste0(c,"_TAM"), verbose = FALSE, recorrect_umi = FALSE)
  write.csv(x,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_male",c,".csv"))
  p<-  EnhancedVolcano(x,lab = rownames(x), x = "avg_log2FC", y = "p_val",
                       pCutoffCol = "p_val_adj" ,FCcutoff = 1.0 ,pCutoff = 1e-05,
                       title = paste0("DE ",c," TAM vs EtOH only males"), 
                       caption = 'FC cutoff: 1.5; p-value cutoff: p-val adj<1e-05')
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Volcano_Main_male",c,".png"))
  print(p)
  dev.off()
  print(paste0("I just saved Enhanced Volcano of ",c," for male samples."))
  
}

########################## Identify differentially expressed Genes in Clusters across in 87 vs 91 Enhanced Volcano ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
NPC_TAM <-subset(NPC_ALL_TRANSFORMED, subset = stim == "TAM")
Idents(NPC_TAM) <- "celltype.sex"
NPC_TAM <- PrepSCTFindMarkers(NPC_TAM)
a<-unique(NPC_TAM@meta.data$mouseRNA.main)
NPC_TAM@meta.data%>%group_by(mouseRNA.main)%>%summarise(n=n())
a<-subset(a, subset=a!= "Dendritic cells")
b<- list()
for (c in a){
  x <-subset(NPC_TAM, idents = c(paste0(c,"_female"), paste0(c,"_male")))
  x <-FindMarkers(x, assay = "SCT", ident.2 = paste0(c,"_female"), ident.1 = paste0(c,"_male"), verbose = FALSE, recorrect_umi = FALSE)
  write.csv(x,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_87vs91",c,".csv"))
  p<-  EnhancedVolcano(x,lab = rownames(x), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,FCcutoff = 1.0,pCutoff = 1e-05,
                       title = paste0("DE ",c," 91 vs 87"), caption = 'FC cutoff, 1.5; p-value cutoff: p-val_adj <1e-05')
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Volcano_Main_87vs91",c,".png"))
  print(p)
  dev.off()
  print(paste0("I just saved Enhanced Volcano of ",c," for TAM samples."))
  
}


##################### Find TOP Markers that define Clusters ########

######################## Find Conserved Markers in Clusters across Stimulation on Integrated Assay ########
ConservedMarkers<-list()
b<-unique(NPC_ALL_TRANSFORMED@meta.data$mouseRNA.main)
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
for (c in b){
  markers <- FindConservedMarkers(NPC_ALL_TRANSFORMED, assay = "integrated", ident.1 = c, grouping.var = "stim",verbose = FALSE)
  write.csv(markers,paste0("./99_other/2_Clustering/ConservedMarkers_integrated_",c,".csv"))
  ConservedMarkers<-append(ConservedMarkers,markers)
  print(paste0("I saved conserved Markers for ",c,"."))
}
##################### Find TOP Markers that define Clusters on Integrated Assay########
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
DefaultAssay(NPC_ALL_TRANSFORMED) <-"integrated"
all.markers <- FindAllMarkers(NPC_ALL_TRANSFORMED, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 10) %>%  ungroup() -> top10__MOUSEMAIN_CLUSTERMARKER
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.fine"
all.markers <- FindAllMarkers(NPC_ALL_TRANSFORMED, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 10) %>%  ungroup() -> top10__MOUSEFINE_CLUSTERMARKER
dim(all.markers)
write.csv(top10__MOUSEMAIN_CLUSTERMARKER, "./99_other/Clustering_3_top10_mouseMAIN_CLUSTERMARKER_integrated.csv", row.names=FALSE)
write.csv(top10__MOUSEFINE_CLUSTERMARKER, "./99_other/Clustering_3_top10_mouseFINE_CLUSTERMARKER_integrated.csv", row.names=FALSE)

#### Save RDS with Cluster Markers found in Integrated Dataset #####
saveRDS(NPC_ALL_TRANSFORMED, file = "./01_tidy_data/5_NPC_ALL_TRANSFORM_Markers_on_integrated.rds")

####################################################################################################################
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
