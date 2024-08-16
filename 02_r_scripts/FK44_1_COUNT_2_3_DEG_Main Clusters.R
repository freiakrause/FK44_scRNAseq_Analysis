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
library(readr)
source("02_r_scripts/malat1_function.R")
source("02_r_scripts/VlnPlot_Function.R") 
set.seed(42)
#### Load Input Data ####
### Results from FK44.1_COUNT_2_0
### SoupX, QC, SCT, Integration, Clustering
NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/4_NPC_ALL_TRANSFORM_Markers.rds")

########################## Identify differentially expressed Genes in Clusters across Conditions Enhanced Volcano ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
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
  names(single_l.de) <- paste0(names(single_l.de), "_l.sc")
  single_l.de$gene <- rownames(single_l.de)
  names(single_M.de) <- paste0(names(single_M.de), "_M.sc")
  single_M.de$gene <- rownames(single_M.de)
  merge_dat <- merge(single_M.de, single_l.de,by = "gene")
  merge_dat <- merge_dat[order(merge_dat$p_val_M.sc), ]
  # Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
  common <- merge_dat$gene[which(merge_dat$p_val_l.sc < 0.05&
                                   merge_dat$p_val_M.sc < 0.05)]
  write.csv(common,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_Common",c,".csv"))
  only_sc_l <- merge_dat$gene[which(merge_dat$p_val_M.sc > 0.05 &
                                      merge_dat$p_val_l.sc < 0.05)]
  write.csv(only_sc_l,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_Only_Limma",c,".csv"))
  
  only_sc_M <- merge_dat$gene[which(  merge_dat$p_val_l.sc > 0.05 &
                                        merge_dat$p_val_M.sc < 0.05)]
  write.csv(only_sc_M,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_Only_MAST",c,".csv"))
  
}

MAST<-read.csv("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_single_MASTT cells.csv")
MAST<-MAST%>%mutate(pct.diff=(pct.1-pct.2))%>%arrange(desc(avg_log2FC),desc(pct.1))
p<-DoHeatmap(subset(NPC_ALL_TRANSFORMED, mouseRNA.main=="T cells"), assay = "RNA",slot = "scale.data", features = c(unique(MAST$X[10:100])),
             draw.lines = T,lines.width = NULL,
             label = F, group.bar =T)+
  scale_fill_viridis_c()+
  theme(axis.text.y.left = element_text(size = 6),
        legend.justification = "top", 
        legend.title = element_text(size=6),
        legend.key.height= unit(0.2, 'cm'), 
        legend.key.width= unit(0.1, 'cm'),
        legend.text = element_text(size=4))
print(p)



p<-DotPlot(subset(NPC_ALL_TRANSFORMED, mouseRNA.main=="T cells"),  assay = "RNA",scale=F,features =c(unique(MAST$X[1:50])))+
  RotatedAxis()+
  scale_size(breaks = c(0, 25, 50, 75, 100),range(0,10))+
  scale_colour_distiller(palette="Blues", trans="reverse")+
  guides(colour = guide_colourbar(reverse = TRUE))+
  theme(panel.background = element_rect(fill = "red",colour="black", linewidth = 1),
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









png(paste0("./03_plots/3_DEG_Analysis_MainCluster/EnhancedVolcano_Single_Limma",c,".png"))
p<-  EnhancedVolcano(single_l.de,lab = rownames(single_l.de), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj", pCutoff = 1e-05, FCcutoff = 1.0 ,
                     title = paste0("DE ",c," TAM vs EtOH"), 
                     caption = 'FC cutoff, 1.0; p-value cutoff: p_val_adj<1e-05')
print(p)
dev.off()
png(paste0("./03_plots/3_DEG_Analysis_MainCluster/EnhancedVolcano_Single_MAST",c,".png"))
p<-  EnhancedVolcano(MAST,lab = MAST$X, x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj", pCutoff = 0.05, FCcutoff = 0.5 ,
                     title = paste0("DE ",c," TAM vs EtOH"), 
                     caption = 'FC cutoff, 0.5; p-value cutoff: p_val_adj<1e-05')
print(p)
dev.off()
print(paste0("I just saved Enhanced Volcano of singeDE Analysis of ",c,"."))