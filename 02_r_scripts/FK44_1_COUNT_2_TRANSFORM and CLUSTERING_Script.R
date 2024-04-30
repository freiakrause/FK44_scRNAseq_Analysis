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
#Basic QC is done. QC5 and QC6 contain filtered data: UMI count, Gene Count, MT% and Complexity
#Since we have performed extensive QC with doublet and empty cell removal, we can now test
#if  cell cycle score or MT percentage change very dramatically between clusters
#so that we know, if we remove them by regression during SCTransform,we will not remove biological signal,
#but only some unwanted variation.
#NPC_QC<-readRDS("./01_tidy_data/QC_noQC_NPC.combined")
#NPC_QC<-readRDS("./01_tidy_data/QC_QC5_NPC_QC5.combined")
NPC_QC<-readRDS("./01_tidy_data/QC_QC6_NPC_QC6.combined")

########################## Check if Regression of MT% and CellCycle Score is appropriate ##################################
all.genes <- rownames(NPC_QC)
NPC_QC <- ScaleData(NPC_QC, features = all.genes)
NPC_QC <- RunPCA(NPC_QC, features = VariableFeatures(object = NPC_QC))

png("./03_plots/QC_2_QC6_PCADIMENSIONS_PCA1.png")
VizDimLoadings(NPC_QC, dims = 1:9, reduction = "pca") &  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
dev.off()

png("./03_plots/QC_2_QC6_PCADIMENSIONS_HEAT.png")
DimHeatmap(NPC_QC, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
dev.off()

png("./03_plots/QC_2_QC6_PCADIMENSIONS_PCA2.png")
DimPlot(NPC_QC, reduction = "pca")
dev.off()

png("./03_plots/QC_2_QC6_PCADIMENSIONS_ELBOW.png")
ElbowPlot(NPC_QC) #It’s often good to find how many PCs can be used without much information loss. Look were the big drop happens
dev.off()
## Möchte nun schauen,
##- wie viele PCs ich brauch zum rechnen der Untscheidlichkeit.
##-ob cell cycle score wichtig ist für clustering oder nur noch mehr variation reinbring
## ob mt und rb score wichtig sind für clustering, oder noch mehr variation reinbringen

#To get cell cycle to need to convert cell cycle gene list which is for human into mouse
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
NPC_QC <-CellCycleScoring(NPC_QC, s.features = mmus_s, g2m.features = mmus_g2m)
#print(NPC_QC@meta.data)
rm(mmus_g2m, mmus_s)
#now do clustering, resolution matters, default is 0.8 but play around with it to get best resolution
NPC_QC <- FindNeighbors(NPC_QC, dims = 1:10)
NPC_QC_R1 <- FindClusters(NPC_QC, resolution = 0.1)
# NPC_QC_R8 <- FindClusters(NPC_QC, resolution = 0.8)
# NPC_QC_R5 <- FindClusters(NPC_QC, resolution = 0.5)
# NPC_QC_R10 <- FindClusters(NPC_QC, resolution = 1)
NPC_QC_R1 <- RunUMAP(NPC_QC_R1, dims = 1:10, verbose = F)
# NPC_QC_R5 <- RunUMAP(NPC_QC_R5, dims = 1:10, verbose = F)
# NPC_QC_R8 <- RunUMAP(NPC_QC_R8, dims = 1:10, verbose = F)
# NPC_QC_R10 <- RunUMAP(NPC_QC_R10, dims = 1:10, verbose = F)
#table(NPC_QC@meta.data$seurat_clusters)
png("./03_plots/QC_2_REGRESSION_CLUSTERR1 for REGRESSION Test.png")
DimPlot(NPC_QC_R1,label.size = 4,repel = T,label = T)
dev.off()
# DimPlot(NPC_QC_R5,label.size = 4,repel = T,label = T)
# DimPlot(NPC_QC_R8,label.size = 4,repel = T,label = T)
# DimPlot(NPC_QC_R10,label.size = 4,repel = T,label = T)
#FeaturePlot(NPC_QC_R10,features = "percent.mt",label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
#FeaturePlot(NPC_QC,features = "percent.rb",label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))

#table(NPC_QC[[]]$Phase)
png("./03_plots/QC_2_REGRESSION_GENES_UMIs.png")
VlnPlot(NPC_QC_R1,features = c("nCount_RNA","nFeature_RNA")) & theme(plot.title = element_text(size=10))
dev.off()
#VlnPlot(NPC_QC_R5,features = c("nCount_RNA","nFeature_RNA")) & theme(plot.title = element_text(size=10))
#VlnPlot(NPC_QC_R8,features = c("nCount_RNA","nFeature_RNA")) & theme(plot.title = element_text(size=10))
#VlnPlot(NPC_QC_R10,features = c("nCount_RNA","nFeature_RNA")) & theme(plot.title = element_text(size=10))

#FeaturePlot(NPC_QC,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
png("./03_plots/QC_2_REGRESSION_CELLCYCLESCORE.png")
VlnPlot(NPC_QC_R1,features = c("S.Score","G2M.Score")) &   theme(plot.title = element_text(size=10))
dev.off()
#VlnPlot(NPC_QC_R5,features = c("S.Score","G2M.Score")) &   theme(plot.title = element_text(size=10))
#VlnPlot(NPC_QC_R8,features = c("S.Score","G2M.Score")) &   theme(plot.title = element_text(size=10))
#VlnPlot(NPC_QC_R10,features = c("S.Score","G2M.Score")) &   theme(plot.title = element_text(size=10))

rm(NPC_QC_R1, NPC_QC_R10, NPC_QC_R5, NPC_QC_R8)

######################### Perform Transformation of Data Set ###############################################
# SCTransfrom might be beneficial bc it gices better signal to noise ratio. regression is performed with Mt5 and cell cylce Scores bc they introduce unwanted variation
# Should I also regress for Sex and Stimulation? Sex variation is also unwanted at the moment
## regression mit vst.flavors= "v2" klappt nicht, wenn mit sec oder stim anwenden Fehler in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : Kontraste können nur auf Faktoren mit 2 oder mehr Stufen angewendet werden 
NPC_ALL_TRANSFORM <- SCTransform(NPC_QC,  vst.flavor= "v2",method = "glmGamPoi", vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F) 
NPC_ALL_TRANSFORM <- RunPCA(NPC_ALL_TRANSFORM, verbose = F)
NPC_ALL_TRANSFORM <- RunUMAP(NPC_ALL_TRANSFORM, dims = 1:30, verbose = F)
NPC_ALL_TRANSFORM <- FindNeighbors(NPC_ALL_TRANSFORM, dims = 1:30, verbose = F)
NPC_ALL_TRANSFORM <- FindClusters(NPC_ALL_TRANSFORM, verbose = F,resolution = 0.1, save.SNN = TRUE)
#NPC_ALL_TRANSFORM_R5 <- FindClusters(NPC_ALL_TRANSFORM, verbose = F,resolution = 0.5, save.SNN = TRUE)
#NPC_ALL_TRANSFORM_R8 <- FindClusters(NPC_ALL_TRANSFORM, verbose = F,resolution = 0.8, save.SNN = TRUE)
#NPC_ALL_TRANSFORM_R10 <- FindClusters(NPC_ALL_TRANSFORM, verbose = F,resolution = 1, save.SNN = TRUE)
#PrintFindClustersParams(object = NPC_ALL_TRANSFORM)
#table(NPC_ALL_TRANSFORM[[]]$seurat_clusters)

png("./03_plots/QC_3_CLUSTERING_R1.png")
DimPlot(NPC_ALL_TRANSFORM, label = T , repel = T, label.size = 3)
dev.off()
# png("./03_plots/QC_3_CLUSTERING_R5.png")
# DimPlot(NPC_ALL_TRANSFORM_R5, label = T , repel = T, label.size = 3)
# dev.off()
# png("./03_plots/QC_3_CLUSTERING_R8.png")
# DimPlot(NPC_ALL_TRANSFORM_R8, label = T , repel = T, label.size = 3)
# dev.off()
# png("./03_plots/QC_3_CLUSTERING_R10.png")
# DimPlot(NPC_ALL_TRANSFORM_R10, label = T , repel = T, label.size = 3)
# dev.off()
###choose R1 bc it gave me the basic cell types and cluster numbers for rough clusters. Plan is to later go deeper into each rough cluster and do sub clustering

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
# rm(mouseRNA.fine, mouseRNA.main, mouseRNA.ref, sce)
# #"MouseRNAseqData set fits better than MonacoImmuneData and ImmGenData set # So deleted code for these datasets.

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
png("./03_plots/Clustering_1_ClusterNo.png")
DimPlot(NPC_ALL_TRANSFORM, label = T, repel = T ) + ggtitle("Unsupervised clustering")+ NoLegend()
dev.off()
#Vizualise Clusters with Fine Annotation ----
png("./03_plots/Clustering_1_ClusterMouseFine.png")
DimPlot(NPC_ALL_TRANSFORM, label = T, repel = T, group.by = "mouseRNA.fine") + ggtitle("Annotation Fine")+NoLegend()
dev.off()
#Vizualise Cluster with Main Annoation----
png("./03_plots/Clustering_1_ClusterMouseMain.png")
DimPlot(NPC_ALL_TRANSFORM, label = T, repel = T, group.by = "mouseRNA.main") + ggtitle("Annotation Main")+NoLegend()
dev.off()
#Vizualise Clusters with Annotation of Sex ----
png("./03_plots/Clustering_1_ClusterSex.png")
DimPlot(NPC_ALL_TRANSFORM, label = F, repel = T, group.by = "sex") + ggtitle("Sex")
dev.off()
#Vizualise Cluster with Annotation of Stimulus ----
png("./03_plots/Clustering_1_ClusterStim.png")
DimPlot(NPC_ALL_TRANSFORM, label = F, repel = T, group.by = "stim") + ggtitle("Stimulation")
dev.off()
#Vizualise Cluster with Annotation of Animal ----
png("./03_plots/Clustering_1_ClusterAnimal.png")
DimPlot(NPC_ALL_TRANSFORM, label = F, repel = T, group.by = "sample") + ggtitle("Animal")

################### Save the SCT Transformed Seurat Object with Annotations############
#saveRDS(NPC_ALL_TRANSFORM, file = "./01_tidy_data/NPC_ALL_TRANSFORM.rds")
#saveRDS(NPC_ALL_TRANSFORM, file = "./01_tidy_data/NPC_ALL_TRANSFORM_QC5.rds")
#saveRDS(NPC_ALL_TRANSFORM, file = "./01_tidy_data/NPC_ALL_TRANSFORM_QC6.rds")

################## Load SC Transformed Seurat Object with Annotations ###############
#NPC_ALL_TRANSFORM <- readRDS( "./01_tidy_data/NPC_ALL_TRANSFORM.rds")
#NPC_ALL_TRANSFORM <- readRDS( "./01_tidy_data/NPC_ALL_TRANSFORM_QC5.rds")
#NPC_ALL_TRANSFORM <- readRDS( "./01_tidy_data/NPC_ALL_TRANSFORM_QC6.rds")


########################## Identify differentially expressed Genes in Clusters across Conditions ############################
#Add "celltype.stim" to meta data and PrepSCT FIndMarkers----
NPC_ALL_TRANSFORM$celltype.stim <- paste(NPC_ALL_TRANSFORM$mouseRNA.main, NPC_ALL_TRANSFORM$stim, sep = "_")
Idents(NPC_ALL_TRANSFORM) <- "celltype.stim"
NPC_ALL_TRANSFORM <- PrepSCTFindMarkers(NPC_ALL_TRANSFORM)

#Find Markers T cells ----
Tcell_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("T cells_EtOH", "T cells_TAM"))
Tcell_TAM_response <- FindMarkers(Tcell_TAM_response, assay = "SCT", ident.2 = "T cells_EtOH", ident.1 = "T cells_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers Hepatocytes ----
Hep_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("Hepatocytes_EtOH", "Hepatocytes_TAM"))
Hep_TAM_response <- FindMarkers(Hep_TAM_response, assay = "SCT", ident.2 = "Hepatocytes_EtOH", ident.1 = "Hepatocytes_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers NK cells ----
NKcell_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("NK cells_EtOH", "NK cells_TAM"))
NKcell_TAM_response <- FindMarkers(NKcell_TAM_response, assay = "SCT", ident.2 = "NK cells_EtOH", ident.1 = "NK cells_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers Granulocytes ----
Granulo_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("Granulocytes_EtOH", "Granulocytes_TAM"))
Granulo_TAM_response <- FindMarkers(Granulo_TAM_response, assay = "SCT", ident.2 = "Granulocytes_EtOH", ident.1 = "Granulocytes_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers Endothelial Cells ----
Endo_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("Endothelial cells_EtOH", "Endothelial cells_TAM"))
Endo_TAM_response <- FindMarkers(Endo_TAM_response, assay = "SCT", ident.2 = "Endothelial cells_EtOH", ident.1 = "Endothelial cells_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers Fibroblasts/ HSC ----
HSC_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("Fibroblasts_EtOH", "Fibroblasts_TAM"))
HSC_TAM_response <- FindMarkers(HSC_TAM_response, assay = "SCT", ident.2 = "Fibroblasts_EtOH", ident.1 = "Fibroblasts_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers Macrophages ----
Macro_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("Macrophages_EtOH", "Macrophages_TAM"))
Macro_TAM_response <- FindMarkers(Macro_TAM_response, assay = "SCT", ident.2 = "Macrophages_EtOH", ident.1 = "Macrophages_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers Monocytes  ----
Mono_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("Monocytes_EtOH","Monocytes_TAM"))
Mono_TAM_response <- FindMarkers(Mono_TAM_response, assay = "SCT", ident.2 = "Monocytes_EtOH", ident.1 = "Monocytes_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Marker Dendritic Cells ----
DC_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("Dendritic cells_EtOH", "Dendritic cells_TAM"))
DC_TAM_response <- FindMarkers(DC_TAM_response, assay = "SCT", ident.2 = "Dendritic cells_EtOH", ident.1 = "Dendritic cells_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers B cells ----
Bcell_TAM_response <- subset(NPC_ALL_TRANSFORM, idents = c("B cells_EtOH", "B cells_TAM"))
Bcell_TAM_response <- FindMarkers(Bcell_TAM_response, assay = "SCT", ident.2 = "B cells_EtOH", ident.1 = "B cells_TAM", verbose = FALSE, recorrect_umi = FALSE)
#Find Markers Adipocytes ----
#Find Markers Erythrocytes ----
#Find Markers Microglia ----


########################### Volcano Plots DE Genes between clusters ###################################################
png("./03_plots/Clustering_2_Volcano_Hep.png")
EnhancedVolcano(Hep_TAM_response,lab = rownames(Hep_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE Hepatocytes TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_T.png")
EnhancedVolcano(Tcell_TAM_response,lab = rownames(Tcell_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE T cells TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_B.png")
EnhancedVolcano(Bcell_TAM_response,lab = rownames(Bcell_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE B cells TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_NK.png")
EnhancedVolcano(NKcell_TAM_response,lab = rownames(NKcell_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE NK cells TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_Granulo.png")
EnhancedVolcano(Granulo_TAM_response,lab = rownames(Granulo_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE Granulocytes TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_Endo.png")
EnhancedVolcano(Endo_TAM_response,lab = rownames(Endo_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE Endothelial cells TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_HSC.png")
EnhancedVolcano(HSC_TAM_response,lab = rownames(HSC_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE Hepatic Stellate Cells TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_Macros.png")
EnhancedVolcano(Macro_TAM_response,lab = rownames(Macro_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE Macrophages TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_Monos.png")
EnhancedVolcano(Mono_TAM_response,lab = rownames(Mono_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE Monocytes TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
png("./03_plots/Clustering_2_Volcano_DC.png")
EnhancedVolcano(DC_TAM_response,lab = rownames(DC_TAM_response), x = "avg_log2FC", y = "p_val", pCutoffCol = "p_val_adj" ,title = "DE Dendritic Cells TAM vs EtOH", caption = 'FC cutoff, 1; p-value cutoff, 10e-5')
dev.off()
####################################################################################################################
Idents(NPC_ALL_TRANSFORM) <-"mouseRNA.main"
DefaultAssay(NPC_ALL_TRANSFORM) <- "SCT"
FeaturePlot(NPC_ALL_TRANSFORM, features = c("Saa2", "Saa1", "Alb"), split.by = "stim", max.cutoff = 5,cols = c("grey", "red"))
VlnPlot(NPC_ALL_TRANSFORM, features = c("Cd3e", "Cd19","Alb","Saa1","Clec4f"), split.by = "stim",group.by = "mouseRNA.main", pt.size = 0, combine = FALSE)

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

