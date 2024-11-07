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
#renv::install("paletteer")
#renv::install("ggsci")
#renv::install("harmony")
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
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

library(Seurat)
library(tidyverse)
library(ggplot2)
source("02_r_scripts/Function_Malat1.R")
library(SingleR)
library(dplyr)
library(stringr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(gt)
library(paletteer)
library(ggsci)
library(SoupX)
set.seed(42)

############################################# Load Input Data ###############################################
#Basic QC is done. QC6 contains filtered data: UMI count, Gene Count, MT% and Complexity
#Since we have performed extensive QC with doublet and empty cell removal, we can now test
#if  cell cycle score or MT percentage change very dramatically between clusters
#so that we know, if we remove them by regression during SCTransform,we will not remove biological signal,
#but only some unwanted variation.

#To get cell cycle to need to convert cell cycle gene list which is for human into mouse
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

animals<-c("87","88","91","92")
#I do this here per animal/sample bc before I integrated samples after QC and then die SCT and clustering.
#but i saw visual differences per cluster per animal and in SCTransform i can not regress out sex,stim or sample.When I look at integraedd data w/o SCTransform after integration,
#then clusters would not be split by sample. But after sctransform yes
#Online i found that some people to normalization etc per sample and then integration after that might resolve per animal differences
#https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette
for(a in animals){
  print(a)
  NPC <-readRDS(paste0("./01_tidy_data/1_QC_1_QC6_NPC_",a,".rds"))
  all.genes <- rownames(NPC)
  #generated clusters to check cellcycle scoring
  #in older scrit I did this and then plotted results here i only do it and dont plot results
  NPC <- ScaleData(NPC, features = all.genes, verbose = F)
  NPC <- RunPCA(NPC, features = VariableFeatures(object = NPC), npcs = 20, verbose = F)
  NPC <- RunUMAP(NPC, dims = 1:20, verbose = F, reduction = "pca")
  NPC <- FindNeighbors(NPC, dims = 1:20, k.param = 10, verbose = F)
  NPC <- FindClusters(NPC, resolution = 0.1, verbose = F)
  NPC <-CellCycleScoring(NPC, s.features = mmus_s, g2m.features = mmus_g2m)
#### Vizualize PCA results ----
  png(paste0("./03_plots/1_QC/QC_2_PCADIMENSIONS_PCA1_iAL",a,".png"))
  p<-VizDimLoadings(NPC, dims = 1:9, reduction = "pca") &  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_PCADIMENSIONS_HEAT_iAL",a,".png"))
  p<-DimHeatmap(NPC, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_PCADIMENSIONS_PCA2_iAL",a,".png"))
  p<-DimPlot(NPC, reduction = "pca")
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_PCADIMENSIONS_ELBOW_iAL",a,".png"))
  p<-ElbowPlot(NPC) #It’s often good to find how many PCs can be used without much information loss. Look were the big drop happens
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_CLUSTER for REGRESSION Test_iAL",a,".png"))
  p<-DimPlot(NPC,label.size = 4,repel = T,label = T)
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_CLUSTER for REGRESSION MT_iAL",a,".png"))
  p<-FeaturePlot(NPC,features = "percent.mt",label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_MT_iAL",a,".png"))
  p<-VlnPlot(NPC,features = "percent.mt") & theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_CLUSTER for REGRESSION RB_iAL",a,".png"))
  p<-FeaturePlot(NPC,features = "percent.rb",label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_GENES_iAL",a,".png"))
  p<-VlnPlot(NPC,features = c("nFeature_RNA")) & theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_UMIs_iAL",a,".png"))
  p<-VlnPlot(NPC,features = c("nCount_RNA")) & theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_CLUSTER for REGRESSION CELLCYCLESCORE_SG2M_iAL",a,".png"))
  p<-FeaturePlot(NPC,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_CELLCYCLESCORE_iAL",a,".png"))
  p<-VlnPlot(NPC,features = c("S.Score","G2M.Score")) &   theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_CELLCYCLESCORE_S_iAL",a,".png"))
  p<-VlnPlot(NPC,features = "S.Score") &   theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(paste0("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_CELLCYCLESCORE_G2M_iAL",a,".png"))
  p<-VlnPlot(NPC,features = "G2M.Score") &   theme(plot.title = element_text(size=10))
  print(p)
  dev.off()
  png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_iAL",a,".png"))
  p <-RidgePlot(NPC, features = "Malat1")
  print(p)
  dev.off()
  norm_counts <- NPC@assays$RNA$data["Malat1",]
  threshold <- define_malat1_threshold(norm_counts)
  malat1_threshold <- norm_counts > threshold
  NPC$malat1_threshold <- malat1_threshold
  NPC$malat1_threshold <- factor(NPC$malat1_threshold, levels = c("TRUE","FALSE"))
  png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter2_iAL",a,".png"))
  p <-DimPlot(NPC, group.by = "malat1_threshold")
  print(p)
  dev.off()
  saveRDS(NPC, file = paste0("./01_tidy_data/2_1_2_NPC_",a,"_afterCellCycleScoring.rds"))
  print(paste0(" I saved the RDS with CellCycleScoring for NPC_",a,"."))
}
rm(mmus_g2m, mmus_s, NPC)
NPC_87<-readRDS("./01_tidy_data/2_1_2_NPC_87_afterCellCycleScoring.rds")
NPC_88<-readRDS("./01_tidy_data/2_1_2_NPC_88_afterCellCycleScoring.rds")
NPC_91<-readRDS("./01_tidy_data/2_1_2_NPC_91_afterCellCycleScoring.rds")
NPC_92<-readRDS("./01_tidy_data/2_1_2_NPC_92_afterCellCycleScoring.rds")

####Due to change in Seurat integration I had to adjust code here (integrate anchors changed to Integrate Layers ~16.10.24. 
NPC <-merge(NPC_87,y=c(NPC_88,NPC_91,NPC_92), add.cell.ids=c("87","88","91","92"))
NPC <-SCTransform(NPC,  vst.flavor= "v2",method = "glmGamPoi",verbose = F, 
                  vars.to.regress = c("percent.mt","S.Score","G2M.Score"))
NPC <-RunPCA(NPC, npcs=30, verbose = F)
options(future.globals.maxSize = 3e+10)
NPC_ALL_TRANSFORMED <-IntegrateLayers(object= NPC, method = HarmonyIntegration, 
                                      verbose = FALSE, normalization.method = "SCT",
                                      orig.reduction = "pca", new.reduction = 'integrated.harmony')

NPC_ALL_TRANSFORMED <- FindNeighbors(NPC_ALL_TRANSFORMED, dims = 1:30, verbose = F,reduction = "integrated.harmony")
NPC_ALL_TRANSFORMED <- FindClusters(NPC_ALL_TRANSFORMED, verbose = F,resolution = 0.2)
NPC_ALL_TRANSFORMED <- RunUMAP(NPC_ALL_TRANSFORMED, dims = 1:30, verbose = F,reduction = "integrated.harmony")

########### Don't understand what I did here and why it was necessary #####################
###addition 10.10.24: maybe need to do this to have scaled data layer in RNA assay. In HeatMap might be nice to use for visualization
DefaultAssay(NPC_ALL_TRANSFORMED) <- "RNA"
NPC_ALL_TRANSFORMED <- NormalizeData(NPC_ALL_TRANSFORMED)
NPC_ALL_TRANSFORMED <- FindVariableFeatures(NPC_ALL_TRANSFORMED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_ALL_TRANSFORMED)
NPC_ALL_TRANSFORMED <- ScaleData(NPC_ALL_TRANSFORMED, features = all.genes)
NPC_ALL_TRANSFORMED <-JoinLayers(NPC_ALL_TRANSFORMED)
###Test expression of example genes to see if souo removal etc is beter
APP <-c("Alb","Saa1","Hp","Fgg")
Recruitment <-c("Vegfa","Cxcl1","Csf1","Vcam1","Cxcl2","Ccr2") #"Cxcl12",
ECM <-c("Serpina1e","Tgfb1","Col1a1","Mmp2","Fn1","Timp1","Mmp9")#,"Des","Acta2"
Fatty <-c("Fabp1","Apoa1","Apoa2")
DefaultAssay(NPC_ALL_TRANSFORMED) <-"RNA"
p<-VlnPlot(NPC_ALL_TRANSFORMED, 
           features = APP, 
           assay= "RNA", 
           layer= "scale.data",
           log = T, 
           stack = T,
           flip = F, 
           fill.by = "ident",
           split.by="stim",
           combine=T)+
  coord_cartesian( ylim=c(1,16),clip = "off")+
  labs(title= "Acute Phase Response Genes",
       x= "Expression Level", y="Cell Type")+scale_fill_manual(values=c("#90bff9","#99cc99"))+
  theme_classic()+
  theme(plot.title = element_text(size=9.5,hjust=0.5),
        axis.title=element_text(size=9))
print(p)
ggsave(filename = paste0("./03_plots/240711_Test_Clustering before_woSoup.png"),
       p,width =(1+(length(unique(APP))*1.5)),
       height = 5, dpi = 600,bg="transparent")
####
saveRDS(NPC_ALL_TRANSFORMED, "./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
rm(NPC, NPC_87, NPC_88, NPC_91, NPC_92)
#rm(NPC_ALL_TRANSFORMED)
#NPC_ALL_TRANSFORMED<-readRDS("./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")

png(filename ="./03_plots/1_QC/QC_Malat1-Filter_RidgePlot_ALLTransformed.png")
RidgePlot(NPC_ALL_TRANSFORMED, features = "Malat1")
dev.off()

png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_DimPLot_ALLTranformed.png"))
DimPlot(NPC_ALL_TRANSFORMED, group.by = "malat1_threshold")
dev.off()

######## Annotate Clusters using cellDex #################################
sce <- as.SingleCellExperiment(DietSeurat(NPC_ALL_TRANSFORMED))
mouseRNA.ref <- celldex::MouseRNAseqData()
mouseRNA.main <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.main)
mouseRNA.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.fine)
NPC_ALL_TRANSFORMED@meta.data$mouseRNA.main <- mouseRNA.main$pruned.labels
NPC_ALL_TRANSFORMED@meta.data$mouseRNA.fine <- mouseRNA.fine$pruned.labels
x<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main)%>%summarise(n=n())
rm(mouseRNA.fine, mouseRNA.main, mouseRNA.ref, sce)
################### Save the SCT Transformed Seurat Object with Annotations############
saveRDS(NPC_ALL_TRANSFORMED, file = "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated.rds")
################## Load SC Transformed Seurat Object with Annotations ###############
#NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated.rds")
y<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())

##### Throw out Cells/Cell Clusters, that are very low in number and cells that do not pass Malat1 Filter
NPC_ALL_TRANSFORMED<-subset(NPC_ALL_TRANSFORMED, 
                            subset=mouseRNA.main!= "Adipocytes"& #bc only 2 0 vs 2
                              mouseRNA.main!= "Epithelial cells"& #bc only 4 4vs0
                              mouseRNA.main!= "Dendritic cells"& #bc only 17 9vs 8
                              mouseRNA.main!= "Erythrocytes")#& #bc only 16 4 vs 12
                              #(malat1_threshold=="TRUE" | mouseRNA.main=="Hepatocytes"))         #substracts empty droplets/cells wo nucleus, but not for heps bc heps always look like sh*?$%$t
NPC_ALL_TRANSFORMED$mouseRNA.main[grepl("Microglia", NPC_ALL_TRANSFORMED$mouseRNA.main)] <- "Macrophages"       # Mikroglia a brain marcos. i checked some of Conserved micro genes and they also fit kupffer cells, so i just merge these clusters here; I read that HSC also express some mikroglia genes, so maybe my microglia here are HSC?                      
z<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())

myClusterSorting <-c("Hepatocytes","Endothelial cells","B cells",
                     "T cells","Macrophages","Granulocytes","NK cells",
                     "Monocytes","Fibroblasts")
z2<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim,sex)%>%summarise(n=n())%>%ungroup()%>%
  bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~'Total')))%>%
  group_by(mouseRNA.main)%>%
  bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~'Total')))%>%
  ungroup()%>%
  arrange(factor(mouseRNA.main, levels = myClusterSorting))


z2%>%
  gt()%>%
  tab_header(title = "Number of Cells per Cluster and Condition")%>%
  tab_style(style = cell_text(color = "black", weight = "bold", align = "center"),locations = cells_title("title"))%>%
  tab_style(style= cell_text(weight= "bold"),locations = cells_column_labels(everything()))%>%
  opt_table_lines()%>%
  tab_options(table.font.size="small")%>%
  gtsave("Cell_Numbers.docx",path="./03_plots/2_Clustering/")

#Add colums for celltype_stim and celltype_sex and sex_stim for potential future analysis
NPC_ALL_TRANSFORMED$celltype.stim <- paste(NPC_ALL_TRANSFORMED$mouseRNA.main, NPC_ALL_TRANSFORMED$stim, sep = "_")
NPC_ALL_TRANSFORMED$celltype.sex <- paste(NPC_ALL_TRANSFORMED$mouseRNA.main, NPC_ALL_TRANSFORMED$sex, sep = "_")
NPC_ALL_TRANSFORMED$sex.stim <- paste(NPC_ALL_TRANSFORMED$sex, NPC_ALL_TRANSFORMED$stim, sep = "_")


saveRDS(NPC_ALL_TRANSFORMED, file = "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced_woMALAT_Filter.rds")
######
###########################Load Data and Decontaminate with SoupX ##################
#### Soup Decont ----
ConservedMarkers_Top20_woMALAT_Filter <- read_csv("99_other/2_Clustering_bad soup_25.10.24/ConservedMarkers_Top20_woMALAT_Filter.csv")
Hep_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "Hepatocytes")[3]%>%pull(genes)
T_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "T cells")[3]%>%pull(genes)
B_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "B cells")[3]%>%pull(genes)
M_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "Macrophages")[3]%>%pull(genes)
Mono_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "Monocytes")[3]%>%pull(genes)
N_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "Granulocytes")[3]%>%pull(genes)
E_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "Endothelial cells")[3]%>%pull(genes)
F_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "Fibroblasts")[3]%>%pull(genes)
NK_genes <-subset(ConservedMarkers_Top20_woMALAT_Filter,subset=ConservedMarkers_Top20_woMALAT_Filter$cluster== "NK cells")[3]%>%pull(genes)


Gene_List <-list(Hep_genes,T_genes,B_genes,M_genes, Mono_genes,N_genes,E_genes,F_genes,NK_genes)
animals <-c("87","88","91","92")


for (i in animals){ 
  i="87"
  NPC_metadata <-NPC_ALL_TRANSFORMED@meta.data
  NPC_metadata<-filter(NPC_metadata,startsWith(rownames(NPC_metadata),"87"))
  NPC_metadata<-mutate(NPC_metadata,ID = str_extract(string = rownames(NPC_metadata),pattern = "[AGTC]+-1$"))%>%
    remove_rownames %>%
    column_to_rownames(var="ID")
  Clusters<-as.character(pull(NPC_metadata,seurat_clusters))
  names(Clusters)<-rownames(NPC_metadata)
  sc = load10X(paste0("./00_raw_data/biomedical-sequencing.at/projects/BSA_0873_FK44_1_LiverMet_A_1_1_51ddbfd228ec40b096e110101b219cb0/COUNT/Liver_NPC_iAL",i,"_transcriptome"))
  sc = setClusters(sc, Clusters)
  Soup <-sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ]
  write.csv(Soup,paste0("./99_other/0_Decont_SoupX/0_Decont_SoupX_Soup_Genes_iAL",i,".csv"))
  Soup_Genes <-head(rownames(Soup), n=10)
  #sc = autoEstCont(sc)
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG= Hep_genes,T_genes,B_genes,M_genes, Mono_genes,N_genes,E_genes,F_genes,NK_genes))
  sc = calculateContaminationFraction(sc, list(IG= Hep_genes,T_genes,B_genes,M_genes, Mono_genes,N_genes,E_genes,F_genes,NK_genes), useToEst = useToEst,forceAccept=TRUE)
  out = adjustCounts(sc)
  srat <-CreateSeuratObject(out)
  saveRDS(srat, paste0("./01_tidy_data/0_iAL",i,"_SoupX.rds"))
  rm(srat)
  
  ##Vizuals----
  dd = sc$metaData[colnames(sc$toc), ]
  mids = aggregate(cbind(tSNE1, tSNE2) ~ clustersFine, data = dd, FUN = mean)
  gg = ggplot(dd, aes(tSNE1, tSNE2))+
    geom_point(aes(colour = clustersFine), size = 0.2) +
    geom_label(data = mids, aes(label = clustersFine)) + ggtitle(paste0(i))+
    guides(colour = guide_legend(override.aes = list(size = 1)))
  png(paste0("./03_plots/0_Ambient_RNA_removal_SoupX/QC_0_SoupX_",i,"_ClustersFine.png"))
  print(plot(gg))
  dev.off()
  
  for (GL in Gene_List){     
    for (g in GL){
      dd$val = sc$toc[g, ]
      gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
      png(paste0("./03_plots/0_Ambient_RNA_removal_SoupX/QC_0_SoupX_",i,"_Clusters_",g,".png"))
      print(plot(gg))
      dev.off()
      gg = plotMarkerMap(sc, g)
      png(paste0("./03_plots/0_Ambient_RNA_removal_SoupX/QC_0_SoupX_",i,"_Clusters_",g,"_SIG.png"))
      print(plot(gg))
      dev.off()
      gg <-plotChangeMap(sc, out, g)
      png(paste0("./03_plots/0_Ambient_RNA_removal_SoupX/QC_0_SoupX_",i,"_Clusters_",g,"_CHANGE.png"))
      print(gg)
      dev.off()
    }
  }
  for (g in Soup_Genes){
    dd$val = sc$toc[g, ]
    gg = ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = val > 0))
    png(paste0("./03_plots/0_Ambient_RNA_removal_SoupX/QC_0_SoupX_",i,"_Clusters_",g,".png"))
    print(plot(gg))
    dev.off()
    gg = plotMarkerMap(sc, g)
    png(paste0("./03_plots/0_Ambient_RNA_removal_SoupX/QC_0_SoupX_",i,"_Clusters_",g,"_SIG.png"))
    print(plot(gg))
    dev.off()
    gg <-plotChangeMap(sc, out, g)
    png(paste0("./03_plots/0_Ambient_RNA_removal_SoupX/QC_0_SoupX_",i,"_Clusters_",g,"_CHANGE.png"))
    print(gg)
    dev.off()
  }
}
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
####
#NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced_woMALAT_Filter.rds.rds")
#### Vizuals Malat1 Filter ###

Idents(NPC_ALL_TRANSFORMED)<-"mouseRNA.main"
png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_RidgePlot_ALLTransformed_mouseRNAMain_woMALAT_Filter.png"))
RidgePlot(NPC_ALL_TRANSFORMED, features = "Malat1")
dev.off()

png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_DimPLot_ALLTranformed_mouseRNAMain_woMALAT_Filter.png"))
DimPlot(NPC_ALL_TRANSFORMED, group.by = "malat1_threshold")
dev.off()

png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_VlnPLot_ALLTranformed_mouseRNAMain_woMALAT_Filter.png"))
VlnPlot(NPC_ALL_TRANSFORMED,features = "Malat1")
dev.off()

png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_RidgePlot_ALLTransformed_mouseRNAMain_after_MALAT1_and_Number_woMALAT_Filter.png"))
RidgePlot(NPC_ALL_TRANSFORMED, features = "Malat1")
dev.off()

png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_DimPLot_ALLTranformed_mouseRNAMain_MALAT1_and_Number_woMALAT_Filter.png"))
DimPlot(NPC_ALL_TRANSFORMED, group.by = "malat1_threshold")
dev.off()

png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_VlnPLot_ALLTranformed_mouseRNAMain_MALAT1_and_Number_woMALAT_Filter.png"))
VlnPlot(NPC_ALL_TRANSFORMED,features = "Malat1")
dev.off()

############## Vizualise annotated Clusters ###################################
#Vizualise Cluster No ----
png("./03_plots/2_Clustering/Clustering_1_ClusterNo_woLegend_woMALAT_Filter.png")
DimPlot(NPC_ALL_TRANSFORMED, label = T, repel = T , group.by = "seurat_clusters") + ggtitle("Unsupervised clustering")+ NoLegend()
dev.off()

png("./03_plots/2_Clustering/Clustering_1_ClusterNo_wLegend_woMALAT_Filter.png")
NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "seurat_clusters")
DimPlot(NPC_ALL_TRANSFORMED, label = T , repel = T, label.size = 3)+ggtitle("Unsupervised clustering")
dev.off()

#Vizualise Clusters with Fine Annotation ----
png("./03_plots/2_Clustering/Clustering_1_ClusterMouseFine_woLegend_woMALAT_Filter.png")
DimPlot(NPC_ALL_TRANSFORMED, label = T, repel = T, group.by = "mouseRNA.fine") + ggtitle("Annotation Fine")+NoLegend()
dev.off()

png("./03_plots/2_Clustering/Clustering_1_ClusterMouseFine_wLegend_woMALAT_Filter.png")
NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "mouseRNA.fine")
DimPlot(NPC_ALL_TRANSFORMED, label = F , repel = T, label.size = 3)
dev.off()
#Vizualise Cluster with Main Annoation----
p<-DimPlot(NPC_ALL_TRANSFORMED, 
        label = T, 
        label.size = 5,
        repel = T, 
        group.by = "mouseRNA.main",
        pt.size = 1)+
  NoLegend()+
  scale_color_observable()
ggsave(filename = paste0("./03_plots/2_Clustering/Clustering_1_ClusterMouseMain_woLegend_woMALAT_Filter.png"), p,width = 10, height = 10, dpi = 600)


NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "mouseRNA.main")
png("./03_plots/2_Clustering/Clustering_1_ClusterMouseMain_wLegend_woMALAT_Filter.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F , repel = T,group.by = "mouseRNA.main", label.size = 3)+
  scale_color_paletteer_d("peRReo::planb")
dev.off()

#Vizualise Clusters with Annotation of Sex ----
png("./03_plots/2_Clustering/Clustering_1_ClusterSex_woMALAT_Filter.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F, repel = T, group.by = "sex") + ggtitle("Sex")
dev.off()

#Vizualise Cluster with Annotation of Animal ----
png("./03_plots/2_Clustering/Clustering_1_ClusterAnimal_woMALAT_Filter.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F, repel = T, group.by = "sample") + ggtitle("Animal")
dev.off()

NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "stim")
png("./03_plots/2_Clustering/Clustering_3_QC6_CR1_Stim_woMALAT_Filter.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F , group.by = "stim",repel = T, label.size = 3)+ ggtitle("Stimulation")
dev.off()

