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

############################################# Load Input Data ###############################################
#Basic QC is done. QC6 contains filtered data: UMI count, Gene Count, MT% and Complexity
#Since we have performed extensive QC with doublet and empty cell removal, we can now test
#if  cell cycle score or MT percentage change very dramatically between clusters
#so that we know, if we remove them by regression during SCTransform,we will not remove biological signal,
#but only some unwanted variation.

#To get cell cycle to need to convert cell cycle gene list which is for human into mouse
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
NPC_list<-list()
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
  saveRDS(NPC, file = paste0("./01_tidy_data/2_1_2_NPC_",a,"_afterCellCycleScoring.rds"))
  print(paste0(" I saved the RDS with CellCycleScoring for NPC_",a,"."))
  #scaling, norm, UMAP and Clustering per ,mouse
  # SCTransfrom might be beneficial bc it gices better signal to noise ratio. regression is performed with Mt5 and cell cylce Scores bc they introduce unwanted variation
  NPC <- SCTransform(NPC,  vst.flavor= "v2",method = "glmGamPoi",  verbose = F, vars.to.regress = c("percent.mt","S.Score","G2M.Score"))
  print(paste0(" I did SCTransformation for NPC_",a,"."))
  png(filename = paste0("./03_plots/1_QC/Malat1-Filter_iAL",a,".png"))
  RidgePlot(NPC, features = "Malat1")
  dev.off()
  norm_counts <- NPC@assays$SCT@data["Malat1",]
  threshold <- define_malat1_threshold(norm_counts)
  malat1_threshold <- norm_counts > threshold
  NPC$malat1_threshold <- malat1_threshold
  NPC$malat1_threshold <- factor(NPC$malat1_threshold, levels = c("TRUE","FALSE"))
  png(filename = paste0("./03_plots/1_QC/Malat1-Filter2_iAL",a,".png"))
  DimPlot(NPC, group.by = "malat1_threshold")
  dev.off()
  NPC <- subset(NPC, malat1_threshold == "TRUE")
  print(paste0(" I removed cells expressing low Malat1 from NPC_",a,"."))
  NPC_list<-append(NPC_list,NPC)
}



rm(mmus_g2m, mmus_s)
features <-SelectIntegrationFeatures(object.list = NPC_list, nfeatures = 3000)
NPC_list <- PrepSCTIntegration(object.list = NPC_list, anchor.features = features)
NPC.anchors <- FindIntegrationAnchors(object.list = NPC_list, dims = 1:20, normalization.method = "SCT", anchor.features = features)
NPC_ALL_TRANSFORMED <- IntegrateData(anchorset = NPC.anchors, dims = 1:20, normalization.method = "SCT")
NPC_ALL_TRANSFORMED <- RunPCA(NPC_ALL_TRANSFORMED, verbose = F)
NPC_ALL_TRANSFORMED <- RunUMAP(NPC_ALL_TRANSFORMED, dims = 1:30, verbose = F)
NPC_ALL_TRANSFORMED <- FindNeighbors(NPC_ALL_TRANSFORMED, dims = 1:30, verbose = F)
NPC_ALL_TRANSFORMED <- FindClusters(NPC_ALL_TRANSFORMED, verbose = F,resolution = 0.2, save.SNN = TRUE) 
#rm(NPC, NPC.anchors, NPC_list)
#saveRDS(NPC_ALL_TRANSFORMED, "./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
#rm(NPC_ALL_TRANSFORMED)
#NPC_QC<-readRDS("./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
DefaultAssay(NPC_ALL_TRANSFORMED) <- "integrated"
p1 <- DimPlot(NPC_ALL_TRANSFORMED, group.by = "sample")
plot(p1)
###### These image files are from the time, when I did integration first and then cellcycle scoring etc. The code is not yet adjusted to all 4 animals in single object etc. 
  # #### Vizualize PCA results ----
  # png("./03_plots/1_QC/QC_2_QC6_PCADIMENSIONS_PCA1.png")
  # VizDimLoadings(NPC_QC, dims = 1:9, reduction = "pca") &  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
  # dev.off()
  # 
  # png("./03_plots/1_QC/QC_2_QC6_PCADIMENSIONS_HEAT.png")
  # DimHeatmap(NPC_QC, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
  # dev.off()
  # 
  # png("./03_plots/1_QC/QC_2_QC6_PCADIMENSIONS_PCA2.png")
  # DimPlot(NPC_QC, reduction = "pca")
  # dev.off()
  # 
  # png("./03_plots/1_QC/QC_2_QC6_PCADIMENSIONS_ELBOW.png")
  # ElbowPlot(NPC_QC) #It’s often good to find how many PCs can be used without much information loss. Look were the big drop happens
  # dev.off()
  # 
  # 
  # 
  # png("./03_plots/1_QC/QC_2_REGRESSION_CLUSTER for REGRESSION Test.png")
  # DimPlot(NPC_QC,label.size = 4,repel = T,label = T)
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_CLUSTER for REGRESSION MT.png")
  # FeaturePlot(NPC_QC,features = "percent.mt",label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_MT.png")
  # VlnPlot(NPC_QC,features = "percent.mt") & theme(plot.title = element_text(size=10))
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_CLUSTER for REGRESSION RB.png")
  # FeaturePlot(NPC_QC,features = "percent.rb",label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_GENES.png")
  # VlnPlot(NPC_QC,features = c("nFeature_RNA")) & theme(plot.title = element_text(size=10))
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_UMIs.png")
  # VlnPlot(NPC_QC,features = c("nCount_RNA")) & theme(plot.title = element_text(size=10))
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_CLUSTER for REGRESSION CELLCYCLESCORE_SG2M.png")
  # FeaturePlot(NPC_QC,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_CELLCYCLESCORE.png")
  # VlnPlot(NPC_QC,features = c("S.Score","G2M.Score")) &   theme(plot.title = element_text(size=10))
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_CELLCYCLESCORE_S.png")
  # VlnPlot(NPC_QC,features = "S.Score") &   theme(plot.title = element_text(size=10))
  # dev.off()
  # png("./03_plots/1_QC/QC_2_REGRESSION_VIOLIN_CELLCYCLESCORE_G2M.png")
  # VlnPlot(NPC_QC,features = "G2M.Score") &   theme(plot.title = element_text(size=10))
  # dev.off()

######################### Perform Transformation of Data Set ###############################################


########### Don't understand what I did here and why it was necessary #####################
DefaultAssay(NPC_ALL_TRANSFORMED) <- "RNA"
NPC_ALL_TRANSFORMED <- NormalizeData(NPC_ALL_TRANSFORMED)
NPC_ALL_TRANSFORMED <- FindVariableFeatures(NPC_ALL_TRANSFORMED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_ALL_TRANSFORMED)
NPC_ALL_TRANSFORMED <- ScaleData(NPC_ALL_TRANSFORMED, features = all.genes)
NPC_ALL_TRANSFORMED <-JoinLayers(NPC_ALL_TRANSFORMED)
######## Annotate Clusters using cellDex #################################
sce <- as.SingleCellExperiment(DietSeurat(NPC_ALL_TRANSFORMED))
#sce
mouseRNA.ref <- celldex::MouseRNAseqData()
mouseRNA.main <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.main)
mouseRNA.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.fine)
NPC_ALL_TRANSFORMED@meta.data$mouseRNA.main <- mouseRNA.main$pruned.labels
NPC_ALL_TRANSFORMED@meta.data$mouseRNA.fine <- mouseRNA.fine$pruned.labels
rm(mouseRNA.fine, mouseRNA.main, mouseRNA.ref, sce)
x<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main)%>%summarise(n=n())
################### Save the SCT Transformed Seurat Object with Annotations############
saveRDS(NPC_ALL_TRANSFORMED, file = "./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
################## Load SC Transformed Seurat Object with Annotations ###############
NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
y<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())
NPC_ALL_TRANSFORMED<-subset(NPC_ALL_TRANSFORMED, subset=mouseRNA.main!= "Adipocytes"&mouseRNA.main!= "Epithelial cells"&mouseRNA.main!= "NA"&mouseRNA.main!= "Erythrocytes")
z<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())
#### Save and vizualizes TOP10 ClusterMarker #####
DefaultAssay(NPC_ALL_TRANSFORMED) <-"SCT"
NPC_ALL_TRANSFORMED <- PrepSCTFindMarkers(NPC_ALL_TRANSFORMED)
Idents(NPC_ALL_TRANSFORMED) <-factor(Idents(NPC_ALL_TRANSFORMED),
                                     levels=c("T cells", "NK cells","B cells","Macrophages","Monocytes", "Microglia","Granulocytes", "Fibroblasts","Endothelial cells","Erythrocytes","Dendritic cells","Hepatocytes", "Epithelial cells","Adipocytes"))

#Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"

all.markers_main <- FindAllMarkers(NPC_ALL_TRANSFORMED, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
all.markers_main %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 5) %>%  ungroup() -> top5__MOUSEMAIN_CLUSTERMARKER
p<-DoHeatmap(NPC_ALL_TRANSFORMED, assay = "RNA", slot = "data",features = c("Ptprc","Cd3e",top5__MOUSEMAIN_CLUSTERMARKER$gene),
             label = F, group.bar =T)
print(p)
a<-unique(NPC_ALL_TRANSFORMED@meta.data$mouseRNA.main)
DefaultAssay(NPC_ALL_TRANSFORMED) <-"RNA"

for (c in a){
t<-subset(top5__MOUSEMAIN_CLUSTERMARKER,cluster==c)
p<-VlnPlot(NPC_ALL_TRANSFORMED, features = t$gene, log = TRUE, ncol=5)
print(p)
}
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
################################# Load Input Data  ##################################################
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
######################## Find Conserved Markers inClusters across Stimulation ########
ConservedMarkers<-list()
b<-unique(NPC_ALL_TRANSFORMED@meta.data$mouseRNA.main)
Idents(NPC_ALL_TRANSFORMED) <- "mouseRNA.main"
for (c in b){
  markers <- FindConservedMarkers(NPC_ALL_TRANSFORMED, assay = "SCT", ident.1 = c, grouping.var = "stim",verbose = FALSE)
  write.csv(markers,paste0("./99_other/2_Clustering/ConservedMarkers_",c,".csv"))
  ConservedMarkers<-append(ConservedMarkers,markers)
  print(paste0("I saved conserved Markers for ",c,"."))
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

############## Vizualise annotated Clusters ###################################
png("./03_plots/1_QC/QC_3_CLUSTERING.png")
DimPlot(NPC_ALL_TRANSFORMED, label = T , repel = T, label.size = 3,reduction = "umap", group.by = "stim")
dev.off()
#Vizualise Cluster No ----
png("./03_plots/2_Clustering/Clustering_1_ClusterNo.png")
DimPlot(NPC_ALL_TRANSFORMED, label = T, repel = T ) + ggtitle("Unsupervised clustering")+ NoLegend()
dev.off()
#Vizualise Clusters with Fine Annotation ----
png("./03_plots/2_Clustering/Clustering_1_ClusterMouseFine.png")
DimPlot(NPC_ALL_TRANSFORMED, label = T, repel = T, group.by = "mouseRNA.fine") + ggtitle("Annotation Fine")+NoLegend()
dev.off()
#Vizualise Cluster with Main Annoation----
png("./03_plots/2_Clustering/Clustering_1_ClusterMouseMain.png")
DimPlot(NPC_ALL_TRANSFORMED, label = T, repel = T, group.by = "mouseRNA.main") + ggtitle("Annotation Main")+NoLegend()
dev.off()
#Vizualise Clusters with Annotation of Sex ----
png("./03_plots/2_Clustering/Clustering_1_ClusterSex.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F, repel = T, group.by = "sex") + ggtitle("Sex")
dev.off()
#Vizualise Cluster with Annotation of Stimulus ----
png("./03_plots/2_Clustering/Clustering_1_ClusterStim.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F, repel = T, group.by = "stim") + ggtitle("Stimulation")
dev.off()
#Vizualise Cluster with Annotation of Animal ----
png("./03_plots/2_Clustering/Clustering_1_ClusterAnimal.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F, repel = T, group.by = "sample") + ggtitle("Animal")
dev.off()

png("./03_plots/Clustering_3_QC6_CR1_mouseRNA_Anno_fine.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F , repel = T, label.size = 3)
dev.off()

NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "mouseRNA.main")
png("./03_plots/Clustering_3_QC6_CR1_mouseRNA_Anno_main.png")
DimPlot(NPC_ALL_TRANSFORMED, label = T , repel = T, label.size = 3)
dev.off()

NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "sex")
png("./03_plots/Clustering_3_QC6_CR1_Sex.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F , repel = T, label.size = 3)
dev.off()

NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "stim")
png("./03_plots/Clustering_3_QC6_CR1_Stim.png")
DimPlot(NPC_ALL_TRANSFORMED, label = F , repel = T, label.size = 3)
dev.off()

png("./03_plots/Clustering_3_QC6_CR1_ClusterNo.png")
NPC_ALL_TRANSFORMED <- SetIdent(NPC_ALL_TRANSFORMED, value = "seurat_clusters")
DimPlot(NPC_ALL_TRANSFORMED, label = T , repel = T, label.size = 3)
dev.off()

