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
#if (!require("BiocManager", quietly = TRUE))
#  renv::install("BiocManager")
#renv::install("gprofiler2")
#BiocManager::install(version = "3.18")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install('glmGamPoi')
#BiocManager::install("EnhancedVolcano")
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
for(a in animals){
  print(a)
  NPC <-readRDS(paste0("./01_tidy_data/1_QC_1_QC6_NPC_",a,".rds"))
  all.genes <- rownames(NPC)
  NPC <- ScaleData(NPC, features = all.genes, verbose = F)
  NPC <- RunPCA(NPC, features = VariableFeatures(object = NPC), npcs = 20, verbose = F)
  NPC <- RunUMAP(NPC, dims = 1:20, verbose = F, reduction = "pca")
  NPC <- FindNeighbors(NPC, dims = 1:20, k.param = 10, verbose = F)
  NPC <- FindClusters(NPC, resolution = 0.1, verbose = F)
  NPC <-CellCycleScoring(NPC, s.features = mmus_s, g2m.features = mmus_g2m)
  saveRDS(NPC, file = paste0("./01_tidy_data/2_1_2_NPC_",a,"_afterCellCycleScoring.rds"))
  NPC <- SCTransform(NPC,  vst.flavor= "v2",method = "glmGamPoi",  verbose = F, vars.to.regress = c("percent.mt","S.Score","G2M.Score"))
  NPC <- RunPCA(NPC, verbose = F)
  NPC <- RunUMAP(NPC, dims = 1:30, verbose = F)
  NPC <- FindNeighbors(NPC, dims = 1:30, verbose = F)
  NPC <- FindClusters(NPC, verbose = F,resolution = 0.2, save.SNN = TRUE)
  p1 <- DimPlot(NPC, reduction = "umap")
  plot(p1)
  NPC_list<-append(NPC_list,NPC)
}

rm(mmus_g2m, mmus_s)
NPC.anchors <- FindIntegrationAnchors(object.list = NPC_list, dims = 1:20)
NPC_ALL_TRANSFORMED<- IntegrateData(anchorset = NPC.anchors, dims = 1:20)
rm(NPC, NPC.anchors, NPC_list)
saveRDS(NPC_ALL_TRANSFORMED, "./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
#rm(NPC_ALL_TRANSFORMED)
#NPC_QC<-readRDS("./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
DefaultAssay(NPC_ALL_TRANSFORMED) <- "integrated"
p1 <- DimPlot(NPC_ALL_TRANSFORMED, reduction = "umap", group.by = "sample")
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
# SCTransfrom might be beneficial bc it gices better signal to noise ratio. regression is performed with Mt5 and cell cylce Scores bc they introduce unwanted variation
# Should I also regress for Sex and Stimulation? Sex variation is also unwanted at the moment
## regression mit vst.flavors= "v2" klappt nicht, wenn mit sex oder stim anwenden Fehler in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : Kontraste können nur auf Faktoren mit 2 oder mehr Stufen angewendet werden 


png("./03_plots/1_QC/QC_3_CLUSTERING.png")
DimPlot(NPC_ALL_TRANSFORM, label = T , repel = T, label.size = 3,reduction = "umap", group.by = "stim")
dev.off()

########### Don't understand what I did here and why it was necessary #####################
DefaultAssay(NPC_ALL_TRANSFORM) <- "RNA"
NPC_ALL_TRANSFORM <- NormalizeData(NPC_ALL_TRANSFORM)
NPC_ALL_TRANSFORM <- FindVariableFeatures(NPC_ALL_TRANSFORM, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_ALL_TRANSFORM)
NPC_ALL_TRANSFORM <- ScaleData(NPC_ALL_TRANSFORM, features = all.genes)
NPC_ALL_TRANSFORM <-JoinLayers(NPC_ALL_TRANSFORM)
######## Annotate Clusters using cellDex #################################
sce <- as.SingleCellExperiment(DietSeurat(NPC_ALL_TRANSFORM))
#sce
mouseRNA.ref <- celldex::MouseRNAseqData()
mouseRNA.main <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.main)
mouseRNA.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.fine)

table(mouseRNA.main$pruned.labels)
table(mouseRNA.fine$pruned.labels)
NPC_ALL_TRANSFORM@meta.data$mouseRNA.main <- mouseRNA.main$pruned.labels
NPC_ALL_TRANSFORM@meta.data$mouseRNA.fine <- mouseRNA.fine$pruned.labels
rm(mouseRNA.fine, mouseRNA.main, mouseRNA.ref, sce)
#"MouseRNAseqData set fits better than MonacoImmuneData and ImmGenData set # So deleted code for these datasets.

# #### Try to Use Scott "Spatial Protegenomics Macro niches as Ref data set for Annotations-------
# expression_matrix_Scott <- ReadMtx(  feature.column = 1, "./99_other/rawData_mouseStSt/countTable_mouseStSt/matrix.mtx.gz", features = "./99_other/rawData_mouseStSt/countTable_mouseStSt/features.tsv.gz",
#                                      cells = "./99_other/rawData_mouseStSt/countTable_mouseStSt/barcodes.tsv.gz")
# Matrix::readMM("./99_other/rawData_mouseStSt/countTable_mouseStSt/matrix.mtx.gz")
# 
# Scott.main <- SingleR(test = sce,assay.type.test = 1,ref = expression_matrix_Scott,labels = mouseRNA.ref$label.main)
# 
# 
# 
# 

# 
############## Vizualise annotated Clusters ###################################
#Vizualise Cluster No ----
png("./03_plots/2_Clustering/Clustering_1_ClusterNo.png")
DimPlot(NPC_ALL_TRANSFORM, label = T, repel = T ) + ggtitle("Unsupervised clustering")+ NoLegend()
dev.off()
#Vizualise Clusters with Fine Annotation ----
png("./03_plots/2_Clustering/Clustering_1_ClusterMouseFine.png")
DimPlot(NPC_ALL_TRANSFORM, label = T, repel = T, group.by = "mouseRNA.fine") + ggtitle("Annotation Fine")+NoLegend()
dev.off()
#Vizualise Cluster with Main Annoation----
png("./03_plots/2_Clustering/Clustering_1_ClusterMouseMain.png")
DimPlot(NPC_ALL_TRANSFORM, label = T, repel = T, group.by = "mouseRNA.main") + ggtitle("Annotation Main")+NoLegend()
dev.off()
#Vizualise Clusters with Annotation of Sex ----
png("./03_plots/2_Clustering/Clustering_1_ClusterSex.png")
DimPlot(NPC_ALL_TRANSFORM, label = F, repel = T, group.by = "sex") + ggtitle("Sex")
dev.off()
#Vizualise Cluster with Annotation of Stimulus ----
png("./03_plots/2_Clustering/Clustering_1_ClusterStim.png")
DimPlot(NPC_ALL_TRANSFORM, label = F, repel = T, group.by = "stim") + ggtitle("Stimulation")
dev.off()
#Vizualise Cluster with Annotation of Animal ----
png("./03_plots/2_Clustering/Clustering_1_ClusterAnimal.png")
DimPlot(NPC_ALL_TRANSFORM, label = F, repel = T, group.by = "sample") + ggtitle("Animal")
dev.off()
################### Save the SCT Transformed Seurat Object with Annotations############
saveRDS(NPC_ALL_TRANSFORM, file = "./01_tidy_data/3_NPC_ALL_TRANSFORM.rds")


################## Load SC Transformed Seurat Object with Annotations ###############
NPC_ALL_TRANSFORM <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORM.rds")

########################## Identify differentially expressed Genes in Clusters across Conditions ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
NPC_ALL_TRANSFORM$celltype.stim <- paste(NPC_ALL_TRANSFORM$mouseRNA.main, NPC_ALL_TRANSFORM$stim, sep = "_")
NPC_ALL_TRANSFORM$celltype.sex <- paste(NPC_ALL_TRANSFORM$mouseRNA.main, NPC_ALL_TRANSFORM$sex, sep = "_")
Idents(NPC_ALL_TRANSFORM) <- "celltype.stim"
NPC_ALL_TRANSFORM <- PrepSCTFindMarkers(NPC_ALL_TRANSFORM)
a<-unique(NPC_ALL_TRANSFORM@meta.data$mouseRNA.main)
a<-subset(a, subset=a!= "Adipocytes"&a!= "Epithelial cells")
print(p)
b<- list()
for (c in a){
  x <-subset(NPC_ALL_TRANSFORM, idents = c(paste0(c,"_EtOH"), paste0(c,"_TAM")))
  x <-FindMarkers(x, assay = "SCT", ident.2 = paste0(c,"_EtOH"), ident.1 = paste0(c,"_TAM"), verbose = FALSE, recorrect_umi = FALSE)
  write.csv(x,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_",c,".csv"))
  p<-  EnhancedVolcano(x,lab = rownames(x), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val",FCcutoff = 1.5 ,title = paste0("DE ",c," TAM vs EtOH"), caption = 'FC cutoff, 1.5; p-value cutoff: p-val')
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Volcano_Main_ALL",c,".png"))
  print(p)
  dev.off()
}

########################## Identify differentially expressed Genes in Clusters across in Females ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
NPC_female <-subset(NPC_ALL_TRANSFORM, subset = sex == "female")
Idents(NPC_ALL_TRANSFORM) <- "celltype.stim"
NPC_female <- PrepSCTFindMarkers(NPC_female)
a<-unique(NPC_female@meta.data$mouseRNA.main)
a<-subset(a, subset=a!= "Adipocytes"&a!= "Epithelial cells")
b<- list()
for (c in a){
  x <-subset(NPC_female, idents = c(paste0(c,"_EtOH"), paste0(c,"_TAM")))
  x <-FindMarkers(x, assay = "SCT", ident.2 = paste0(c,"_EtOH"), ident.1 = paste0(c,"_TAM"), verbose = FALSE, recorrect_umi = FALSE)
  write.csv(x,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_female",c,".csv"))
  p<-  EnhancedVolcano(x,lab = rownames(x), 
                       x = "avg_log2FC", y = "p_val", 
                       pCutoffCol = "p_val" ,FCcutoff = 1.5, 
                       title = paste0("DE ",c," TAM vs EtOH only females"), 
                       caption = 'FC cutoff: 1.5; p-value cutoff: p-val')
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Volcano_Main_female",c,".png"))
  print(p)
  dev.off()
}

########################## Identify differentially expressed Genes in Clusters across in male ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
NPC_male <-subset(NPC_ALL_TRANSFORM, subset = sex == "male")
Idents(NPC_ALL_TRANSFORM) <- "celltype.stim"
NPC_male <- PrepSCTFindMarkers(NPC_male)
a<-unique(NPC_male@meta.data$mouseRNA.main)
a<-subset(a, subset=a!= "Adipocytes"&a!= "Epithelial cells")
b<- list()
for (c in a){
  x <-subset(NPC_male, idents = c(paste0(c,"_EtOH"), paste0(c,"_TAM")))
  x <-FindMarkers(x, assay = "SCT", ident.2 = paste0(c,"_EtOH"), ident.1 = paste0(c,"_TAM"), verbose = FALSE, recorrect_umi = FALSE)
  write.csv(x,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_male",c,".csv"))
  p<-  EnhancedVolcano(x,lab = rownames(x), x = "avg_log2FC", y = "p_val",
                       pCutoffCol = "p_val" ,FCcutoff = 1.5 ,
                       title = paste0("DE ",c," TAM vs EtOH only males"), 
                       caption = 'FC cutoff: 1.5; p-value cutoff: p-val')
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Volcano_Main_male",c,".png"))
  print(p)
  dev.off()
}

########################## Identify differentially expressed Genes in Clusters across in 87 vs 91 ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
NPC_TAM <-subset(NPC_ALL_TRANSFORM, subset = stim == "TAM")
Idents(NPC_TAM) <- "celltype.sex"
NPC_TAM <- PrepSCTFindMarkers(NPC_TAM)
a<-unique(NPC_TAM@meta.data$mouseRNA.main)
a<-subset(a, subset=a!= "Adipocytes"&a!= "Epithelial cells")
b<- list()
for (c in a){
  x <-subset(NPC_TAM, idents = c(paste0(c,"_female"), paste0(c,"_male")))
  x <-FindMarkers(x, assay = "SCT", ident.2 = paste0(c,"_female"), ident.1 = paste0(c,"_male"), verbose = FALSE, recorrect_umi = FALSE)
  write.csv(x,paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_87vs91",c,".csv"))
  p<-  EnhancedVolcano(x,lab = rownames(x), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val" ,FCcutoff = 1.5,title = paste0("DE ",c," 91 vs 87"), caption = 'FC cutoff, 1.5; p-value cutoff: p-val')
  png(paste0("./03_plots/3_DEG_Analysis_MainCluster/Volcano_Main_87vs91",c,".png"))
  print(p)
  dev.off()
}

####################################################################################################################
Idents(NPC_ALL_TRANSFORM) <-"mouseRNA.main"
DefaultAssay(NPC_ALL_TRANSFORM) <- "SCT"
Hep_genes <-c("Saa1", "Saa2","Alb","mt-Cytb","mt-Co1", "2cdnaFLAG-LGP130", "1cdnaKANA","3cdnaZSGREEN","mt-Atp6")
T_genes <-c("Cd3e","Cd3d", "Cd4", "Cd8a") 
B_genes <-c("Cd19")
M_genes <-c("Clec4f")
Gene_List <-list(Hep_genes,T_genes,B_genes,M_genes)
for (GL in Gene_List){
  for (g in GL){
  p<-FeaturePlot(NPC_ALL_TRANSFORM, features = g, split.by = "stim", max.cutoff = 5,cols = c("grey", "red"))
  print(p)
  p<-VlnPlot(NPC_ALL_TRANSFORM, features = g, split.by = "stim",group.by = "mouseRNA.main", pt.size = 0, combine = FALSE)
  print(p)
  }
}

all.markers <- FindAllMarkers(NPC_ALL_TRANSFORM, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 10) %>%  ungroup() -> top10_CLUSTERMARKER
all.markers %>%  group_by(cluster) %>%  dplyr::filter(cluster == "Hepatocytes") %>%  slice_head(n = 100) %>%  ungroup() -> Hep_CLUSTERMARKER

NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "mouseRNA.fine")
png("./03_plots/Clustering_3_QC6_CR1_mouseRNA_Anno_fine.png")
DimPlot(NPC_ALL_TRANSFORM, label = F , repel = T, label.size = 3)
dev.off()

all.markers <- FindAllMarkers(NPC_ALL_TRANSFORM, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 10) %>%  ungroup() -> top10__MOUSEFINE_CLUSTERMARKER
all.markers %>%  group_by(cluster) %>%  dplyr::filter(cluster == "Hepatocytes") %>%  slice_head(n = 100) %>%  ungroup() -> Hep_CLUSTERMARKER


NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "mouseRNA.main")
png("./03_plots/Clustering_3_QC6_CR1_mouseRNA_Anno_main.png")
DimPlot(NPC_ALL_TRANSFORM, label = T , repel = T, label.size = 3)
dev.off()

all.markers <- FindAllMarkers(NPC_ALL_TRANSFORM, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 10) %>%  ungroup() -> top10__MOUSEMAIN_CLUSTERMARKER
write.csv(top10__MOUSEMAIN_CLUSTERMARKER, "./99_other/Clustering_3_top10_mouseMAIN_CLUSTERMARKER.csv", row.names=FALSE)
write.csv(top10__MOUSEFINE_CLUSTERMARKER, "./99_other/Clustering_3_op10_mouseFINE_CLUSTERMARKER.csv", row.names=FALSE)

NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "sex")
png("./03_plots/Clustering_3_QC6_CR1_Sex.png")
DimPlot(NPC_ALL_TRANSFORM, label = F , repel = T, label.size = 3)
dev.off()

NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "stim")
png("./03_plots/Clustering_3_QC6_CR1_Stim.png")
DimPlot(NPC_ALL_TRANSFORM, label = F , repel = T, label.size = 3)
dev.off()

png("./03_plots/Clustering_3_QC6_CR1_ClusterNo.png")
NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "seurat_clusters")
DimPlot(NPC_ALL_TRANSFORM, label = T , repel = T, label.size = 3)
dev.off()
saveRDS(NPC_ALL_TRANSFORM, file = "./01_tidy_data/4_NPC_ALL_TRANSFORM_Markers.rds")

