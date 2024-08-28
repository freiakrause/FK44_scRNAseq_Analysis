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
library(viridis)
library(paletteer)
library(ggsci)

source("02_r_scripts/malat1_function.R")
source("02_r_scripts/VlnPlot_Function.R")
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
  ###### These image files are from the time, when I did integration first and then cellcycle scoring etc. The code is not yet adjusted to all 4 animals in single object etc. 
  # #### Vizualize PCA results ----
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
  ##### SCTransfrom, (scaling, norm) UMAP and Clustering per ,mouse, MALAT1 treshold
  NPC <- SCTransform(NPC,  vst.flavor= "v2",method = "glmGamPoi",  verbose = F, vars.to.regress = c("percent.mt","S.Score","G2M.Score"))
  print(paste0(" I did SCTransformation for NPC_",a,"."))
  png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_iAL",a,".png"))
  p <-RidgePlot(NPC, features = "Malat1")
  print(p)
  dev.off()
  norm_counts <- NPC@assays$SCT@data["Malat1",]
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

png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_RidgePlot_ALLTransformed.png"))
RidgePlot(NPC_ALL_TRANSFORMED, features = "Malat1")
dev.off()

png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_DimPLot_ALLTranformed",a,".png"))
DimPlot(NPC, group.by = "malat1_threshold")
dev.off()

# NPC <- subset(NPC, malat1_threshold == "TRUE"|)
# print(paste0(" I removed cells expressing low Malat1 from NPC_",a,"."))


#rm(NPC, NPC.anchors, NPC_list)
#saveRDS(NPC_ALL_TRANSFORMED, "./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")
#rm(NPC_ALL_TRANSFORMED)
NPC_ALL_TRANSFORMED<-readRDS("./01_tidy_data/3_NPC_ALL_TRANSFORMED.rds")


DefaultAssay(NPC_ALL_TRANSFORMED) <- "integrated"
p1 <- DimPlot(NPC_ALL_TRANSFORMED, group.by = "sample")
plot(p1)

########### Don't understand what I did here and why it was necessary #####################
DefaultAssay(NPC_ALL_TRANSFORMED) <- "RNA"
NPC_ALL_TRANSFORMED <- NormalizeData(NPC_ALL_TRANSFORMED)
NPC_ALL_TRANSFORMED <- FindVariableFeatures(NPC_ALL_TRANSFORMED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_ALL_TRANSFORMED)
NPC_ALL_TRANSFORMED <- ScaleData(NPC_ALL_TRANSFORMED, features = all.genes)
NPC_ALL_TRANSFORMED <-JoinLayers(NPC_ALL_TRANSFORMED)
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
NPC_ALL_TRANSFORMED$mouseRNA.main[grepl("Microglia", NPC_ALL_TRANSFORMED$mouseRNA.main)] <- "Macrophages or HSC"       # Mikroglia a brain marcos. i checked some of Conserved micro genes and they also fit kupffer cells, so i just merge these clusters here; I read that HSC also express some mikroglia genes, so maybe my microglia here are HSC?                      
z<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())
myClusterSorting <-c("Hepatocytes","Endothelial cells","B cells","T cells","Macrophages","Granulocytes","NK cells","Monocytes","Fibroblasts","Macrophages or HSC")
z<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim,sex)%>%summarise(n=n())%>%ungroup()%>%
  group_by(mouseRNA.main)%>%
  bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~'Total')))%>%
  ungroup()%>%
  bind_rows(summarise(., across(where(is.numeric), sum),across(where(is.character), ~'Total')))%>% 
  arrange(factor(mouseRNA.main, levels = myClusterSorting))


z%>%
  gt()%>%
  tab_header(title = "Number of Cells per Cluster and Condition")%>%
  tab_style(style = cell_text(color = "black", weight = "bold", align = "center"),locations = cells_title("title"))%>%
  tab_style(style= cell_text(weight= "bold"),locations = cells_column_labels(everything()))%>%
  opt_table_lines()%>%
  tab_options(table.font.size="small")%>%
  gtsave("Cell_Numbers.docx",path="./03_plots/2_Clustering/")

#Add colums for celltype_stim and celltype_sex for potential future analysis

NPC_ALL_TRANSFORMED$celltype.stim <- paste(NPC_ALL_TRANSFORMED$mouseRNA.main, NPC_ALL_TRANSFORMED$stim, sep = "_")
NPC_ALL_TRANSFORMED$celltype.sex <- paste(NPC_ALL_TRANSFORMED$mouseRNA.main, NPC_ALL_TRANSFORMED$sex, sep = "_")
#saveRDS(NPC_ALL_TRANSFORMED, file = "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced.rds")
saveRDS(NPC_ALL_TRANSFORMED, file = "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced_woMALAT_Filter.rds")

#NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced.rds")

#### Vizuals Malat1 Filter ###
Idents(NPC_ALL_TRANSFORMED)<-"mouseRNA.main"
#png(filename = paste0("./03_plots/1_QC/QC_Malat1-Filter_RidgePlot_ALLTransformed_mouseRNAMain.png"))
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

