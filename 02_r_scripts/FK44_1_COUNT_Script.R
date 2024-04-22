#This is the script to analyse FK44.1 scRNAseq COUNT data provided by BSF
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
#remotes::install_github("Moonerss/scrubletR")

#bioc manager installieren hab es nicht geachaftt die BiocManager Packages über renv zu installieren aber während installation wurde revn folder aufgerufen. also vielleicht trotzdem gepasst
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(scrubletR)
library(gprofiler2)
set.seed(42)

#Loading and Combining output (non normalized) fromCellRanger

NPC_87.data <- Read10X(data.dir = "./00_raw_data/iAL87") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
NPC_87 <- CreateSeuratObject(count = NPC_87.data, project = "FK44_NPC_87", min.cells = 3, min.features = 200)
NPC_88.data <- Read10X(data.dir = "./00_raw_data/iAL88") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
NPC_88 <- CreateSeuratObject(count = NPC_88.data, project = "FK44_NPC_88", min.cells = 3, min.features = 200)
NPC_91.data <- Read10X(data.dir = "./00_raw_data/iAL91") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
NPC_91 <- CreateSeuratObject(count = NPC_91.data, project = "FK44_NPC_91", min.cells = 3, min.features = 200)
NPC_92.data <- Read10X(data.dir = "./00_raw_data/iAL92") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
NPC_92 <- CreateSeuratObject(count = NPC_92.data, project = "FK44_NPC_92", min.cells = 3, min.features = 200)
#Kategorien stimlulation und sex hinzufügen, evtl noch age?
NPC_87$stim <- "TAM"
NPC_88$stim <- "EtOH"
NPC_91$stim <- "TAM"
NPC_92$stim <- "EtOH"
NPC_87$sex <- "female"
NPC_88$sex <- "female"
NPC_91$sex <- "male"
NPC_92$sex <- "male"


#Setting things up for qualitycontrol mitochondrial genes, ribosomal content, doublets
NPC_87[["percent.mt"]] <- PercentageFeatureSet(NPC_87, pattern = "^mt-") #
NPC_87[["percent.rb"]] <- PercentageFeatureSet(NPC_87, pattern = "Rp[sl]")
NPC_88[["percent.mt"]] <- PercentageFeatureSet(NPC_88, pattern = "^mt-") #
NPC_88[["percent.rb"]] <- PercentageFeatureSet(NPC_88, pattern = "Rp[sl]")
NPC_91[["percent.mt"]] <- PercentageFeatureSet(NPC_91, pattern = "^mt-") #
NPC_91[["percent.rb"]] <- PercentageFeatureSet(NPC_91, pattern = "Rp[sl]")
NPC_92[["percent.mt"]] <- PercentageFeatureSet(NPC_92, pattern = "^mt-") #
NPC_92[["percent.rb"]] <- PercentageFeatureSet(NPC_92, pattern = "Rp[sl]")

#Add Conclusions from meta data as QC columns
#NPC_combined[['QC']] <- ifelse(NPC_combined@meta.data$Is_doublet == 'True','Doublet','Pass')
NPC_87[['QC']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 500 , 'Low_nFeature','Pass')
NPC_87[['QC']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 500 & NPC_87@meta.data$QC != 'Pass' & NPC_87@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',NPC_87@meta.data$QC,sep = ','),NPC_87@meta.data$QC)
NPC_87[['QC']] <- ifelse(NPC_87@meta.data$percent.mt > 15 & NPC_87@meta.data$QC == 'Pass','High_MT',NPC_87@meta.data$QC)
NPC_87[['QC']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 500 & NPC_87@meta.data$QC != 'Pass'& NPC_87@meta.data$QC != 'High_MT',paste('High_MT',NPC_87@meta.data$QC,sep = ','),NPC_87@meta.data$QC)
table(NPC_87[['QC']])

NPC_88[['QC']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 500 , 'Low_nFeature','Pass')
NPC_88[['QC']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 500 & NPC_88@meta.data$QC != 'Pass' & NPC_88@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',NPC_88@meta.data$QC,sep = ','),NPC_88@meta.data$QC)
NPC_88[['QC']] <- ifelse(NPC_88@meta.data$percent.mt > 15 & NPC_88@meta.data$QC == 'Pass','High_MT',NPC_88@meta.data$QC)
NPC_88[['QC']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 500 & NPC_88@meta.data$QC != 'Pass'& NPC_88@meta.data$QC != 'High_MT',paste('High_MT',NPC_88@meta.data$QC,sep = ','),NPC_88@meta.data$QC)
table(NPC_88[['QC']])

NPC_91[['QC']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 500 , 'Low_nFeature','Pass')
NPC_91[['QC']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 500 & NPC_91@meta.data$QC != 'Pass' & NPC_91@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',NPC_91@meta.data$QC,sep = ','),NPC_91@meta.data$QC)
NPC_91[['QC']] <- ifelse(NPC_91@meta.data$percent.mt > 15 & NPC_91@meta.data$QC == 'Pass','High_MT',NPC_91@meta.data$QC)
NPC_91[['QC']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 500 & NPC_91@meta.data$QC != 'Pass'& NPC_91@meta.data$QC != 'High_MT',paste('High_MT',NPC_91@meta.data$QC,sep = ','),NPC_91@meta.data$QC)
table(NPC_91[['QC']])

NPC_92[['QC']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 500 , 'Low_nFeature','Pass')
NPC_92[['QC']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 500 & NPC_92@meta.data$QC != 'Pass' & NPC_92@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',NPC_92@meta.data$QC,sep = ','),NPC_92@meta.data$QC)
NPC_92[['QC']] <- ifelse(NPC_92@meta.data$percent.mt > 15 & NPC_92@meta.data$QC == 'Pass','High_MT',NPC_92@meta.data$QC)
NPC_92[['QC']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 500 & NPC_92@meta.data$QC != 'Pass'& NPC_92@meta.data$QC != 'High_MT',paste('High_MT',NPC_92@meta.data$QC,sep = ','),NPC_92@meta.data$QC)
table(NPC_92[['QC']])

#Plot only cells that pass QC
png("./output/QC_Features_AFTER_500_15_iAL88.png")
VlnPlot(subset(NPC_88, subset = QC == 'Pass'), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10))
dev.off()
png("./output/QC_Features_AFTER_500_15_iAL87.png")
VlnPlot(subset(NPC_87, subset = QC == 'Pass'), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10))
dev.off()

png("./output/QC_Features_AFTER_500_15_iAL91.png")
VlnPlot(subset(NPC_91, subset = QC == 'Pass'), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10))
dev.off()

png("./output/QC_Features_AFTER_500_15_iAL92.png")
VlnPlot(subset(NPC_92, subset = QC == 'Pass'), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10))
dev.off()


#Plot ALL cells that
png("./output/QC_Features_woQC_iAL88.png")
VlnPlot(subset(NPC_88), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10))
dev.off()
png("./output/QC_Features_woQC_iAL87.png")
VlnPlot(subset(NPC_87), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10))
dev.off()

png("./output/QC_Features_woQC_iAL91.png")
VlnPlot(subset(NPC_91), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10))
dev.off()

png("./output/QC_Features_woQC_iAL92.png")
VlnPlot(subset(NPC_92), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10))
dev.off()

#normalize data set to account for sequencing depth, default scale to 10 000 and log2-transform
NPC_87 <- NormalizeData(NPC_87)
NPC_87 <- FindVariableFeatures(NPC_87, selection.method = "vst", nfeatures = 2000)
NPC_88 <- NormalizeData(NPC_88)
NPC_88 <- FindVariableFeatures(NPC_88, selection.method = "vst", nfeatures = 2000)
NPC_91 <- NormalizeData(NPC_91)
NPC_91 <- FindVariableFeatures(NPC_91, selection.method = "vst", nfeatures = 2000)
NPC_92 <- NormalizeData(NPC_92)
NPC_92 <- FindVariableFeatures(NPC_92, selection.method = "vst", nfeatures = 2000)
NPC.anchors <- FindIntegrationAnchors(object.list = list(NPC_87, NPC_88, NPC_91, NPC_92), dims = 1:20)
NPC.combined<- IntegrateData(anchorset = NPC.anchors, dims = 1:20)
#Remove unused data from memory to save ram
rm(NPC_87.data, NPC_88.data, NPC_91.data, NPC_92.data,NPC_87, NPC_88, NPC_91, NPC_92)

#Scale data
all.genes <- rownames(NPC.combined)
NPC_ALL_Scaled <- ScaleData(NPC.combined, features = all.genes)
#PCA 
NPC_ALL_Scaled <- RunPCA(NPC_ALL_Scaled, features = VariableFeatures(object = NPC_ALL_Scaled))

VizDimLoadings(NPC_ALL_Scaled, dims = 1:9, reduction = "pca") &  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
DimHeatmap(NPC_ALL_Scaled, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(NPC_ALL_Scaled, reduction = "pca")
ElbowPlot(NPC_ALL_Scaled) #It’s often good to find how many PCs can be used without much information loss. Look were the big drop happens

#now do clustering, resolution matters, default is 0.8 but play around with it to get best resolution
NPC_ALL_Scaled <- FindNeighbors(NPC_ALL_Scaled, dims = 1:10)
NPC_ALL_Scaled <- FindClusters(NPC_ALL_Scaled, resolution = 0.5)
NPC_ALL_Scaled <- RunUMAP(NPC_ALL_Scaled, dims = 1:10, verbose = F)
table(NPC_ALL_Scaled@meta.data$seurat_clusters)
DimPlot(NPC_ALL_Scaled,label.size = 4,repel = T,label = T)


NPC_ALL_Scaled_QC <- subset(NPC_ALL_Scaled, subset = QC == 'Pass')
rm(NPC_ALL_Scaled)

DimPlot(NPC_ALL_Scaled_QC,label.size = 4,repel = T,label = T)
FeaturePlot(NPC_ALL_Scaled_QC,features = "percent.mt",label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
FeaturePlot(NPC_ALL_Scaled_QC,features = "percent.rb",label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))

#Also want to get cell cycle to nedd to convert cell cycle gene list which is for human into mouse
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
NPC_ALL_Scaled_QC <- CellCycleScoring(NPC_ALL_Scaled_QC, s.features = mmus_s, g2m.features = mmus_g2m)
table(NPC_ALL_Scaled_QC[[]]$Phase)



VlnPlot(NPC_ALL_Scaled_QC,features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(NPC_ALL_Scaled_QC,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
VlnPlot(NPC_ALL_Scaled_QC,features = "percent.rb") & theme(plot.title = element_text(size=10))

VlnPlot(NPC_ALL_Scaled_QC,features = c("nCount_RNA","nFeature_RNA")) & theme(plot.title = element_text(size=10))

FeaturePlot(NPC_ALL_Scaled_QC,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) &   theme(plot.title = element_text(size=10))
VlnPlot(NPC_ALL_Scaled_QC,features = c("S.Score","G2M.Score")) &   theme(plot.title = element_text(size=10))

#Since we have performed extensive QC with doublet and empty cell removal, we can now apply SCTransform normalization, 
#that was shown to be beneficial for finding rare cell populations by improving signal/noise ratio
#Single SCTransform command replaces NormalizeData, ScaleData, and FindVariableFeatures
#We will also correct for % MT genes and cell cycle scores using vars.to.regress variables
#our previous exploration has shown that neither cell cycle score nor MT percentage change very dramatically between clusters,
#so we will not remove biological signal, but only some unwanted variation.
NPC_ALL_TRANSFORM <- SCTransform(NPC_ALL_Scaled_QC, method = "glmGamPoi", ncells = 8824, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)

NPC_ALL_TRANSFORM <- RunPCA(NPC_ALL_TRANSFORM, verbose = F)
NPC_ALL_TRANSFORM <- RunUMAP(NPC_ALL_TRANSFORM, dims = 1:30, verbose = F)
NPC_ALL_TRANSFORM <- FindNeighbors(NPC_ALL_TRANSFORM, dims = 1:30, verbose = F)
NPC_ALL_TRANSFORM <- FindClusters(NPC_ALL_TRANSFORM, verbose = F, save.SNN = TRUE)
PrintFindClustersParams(object = NPC_ALL_TRANSFORM)
table(NPC_ALL_TRANSFORM[[]]$seurat_clusters)
png("./output/UMAP_Clusters.png")
DimPlot(NPC_ALL_TRANSFORM,label = T)
dev.off()
#Schauen, ob wir Cluster noch correct habeb, durch expression von marker genen von kleinen Populationen platelets und DCs hier (einfüllen,w as ich will) außerdem wurde Farbschmea geändert

FeaturePlot(NPC_ALL_TRANSFORM,"Ptprc") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(NPC_ALL_TRANSFORM,"Cd8a") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(NPC_ALL_TRANSFORM,"Lilrb4a") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(NPC_ALL_TRANSFORM,"Adgre1") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(NPC_ALL_TRANSFORM,"Saa1") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(NPC_ALL_TRANSFORM,"Saa2") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(NPC_ALL_TRANSFORM,"Alb") &  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

DefaultAssay(NPC_ALL_TRANSFORM) <- "RNA"
NPC_ALL_TRANSFORM <- NormalizeData(NPC_ALL_TRANSFORM)
NPC_ALL_TRANSFORM <- FindVariableFeatures(NPC_ALL_TRANSFORM, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_ALL_TRANSFORM)
NPC_ALL_TRANSFORM <- ScaleData(NPC_ALL_TRANSFORM, features = all.genes)
NPC_ALL_TRANSFORM <-JoinLayers(NPC_ALL_TRANSFORM)
all.markers <- FindAllMarkers(NPC_ALL_TRANSFORM, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
all.markers %>%  group_by(cluster) %>%  dplyr::filter(avg_log2FC > 1) %>%  slice_head(n = 10) %>%  ungroup() -> top10_CLUSTERMARKER
png("./output/UMAP_MARKER__Cxcl1.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Cxcl1"))
dev.off()
png("./output/UMAP_MARKER_Saa1.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Saa1"))
dev.off()
VlnPlot(NPC_ALL_TRANSFORM, features = c("Saa2"))
png("./output/UMAP_MARKER_Alb.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Alb"))
dev.off()
VlnPlot(NPC_ALL_TRANSFORM, features = c("Gfap"))
png("./output/UMAP_MARKER__Nkg7.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Saa1"), split.by = "stim")
dev.off()
png("./output/UMAP_MARKER__Ptprc.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Ptprc"),split.by = "stim")
dev.off()
png("./output/UMAP_MARKER_Cxcl1.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Cxcl1"),split.by = "stim")
dev.off()
#HSC Markers
#Extensive studies using the CCl4 model of liver fibrosis in mice have identified 3 distinct HSC phenotypes; 
#quiescent (defined by genes such as Ngfr, Hspa1a, and Hspa1b), 
#activated (defined by genes such as Col1a1, Acta2, and Timp1), and 
#inactivated (defined by genes such as Smoc2, Gabra3, and Gsn) (3,4,33). PMID: 33550587
#Acta 2, Cryab, Spp1, Prnp, and Pai-1 PMID: 17072980
#Gpc3, Lox, and Mgp PMID: 17072980
png("./output/UMAP_MARKER_HSC_Cryab.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Cryab"))
dev.off()
png("./output/UMAP_MARKER_HSC_Prnp.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Prnp"))
dev.off()
png("./output/UMAP_MARKER_HSC_Spp1.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Spp1"))
dev.off()
png("./output/UMAP_MARKER_HSC_Gpc3.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Gpc3"))
dev.off()
png("./output/UMAP_MARKER_HSC_Lox.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Lox"))
dev.off()
png("./output/UMAP_MARKER_HSC_Mgp.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Mgp"))
dev.off()
#quiescent not present
VlnPlot(NPC_ALL_TRANSFORM, features = c("Ngfr"))
dev.off()
VlnPlot(NPC_ALL_TRANSFORM, features = c("Hspa1a"))
dev.off()
VlnPlot(NPC_ALL_TRANSFORM, features = c("Hspa1b"))
dev.off()
#activated looks like it
png("./output/UMAP_MARKER_HSC_activated_Col1a1.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Col1a1"))
dev.off()
png("./output/UMAP_MARKER_HSC_activated_Acta2.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Acta2"))
dev.off()
png("./output/UMAP_MARKER_HSC_activated_Timp1.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Timp1"))
dev.off()
#maybe also some inactivated
png("./output/UMAP_MARKER_HSC_inactivated_Smoc2.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Smoc2"))
dev.off()
png("./output/UMAP_MARKER_HSC_inactivated_Gabra3.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Gabra3"))
dev.off()
png("./output/UMAP_MARKER_HSC_inactivated_Gsn.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Gsn"))
dev.off()
#Kupffer Cell Lyz, Gzmb, and Il1b,
VlnPlot(NPC_ALL_TRANSFORM, features = c("Lyz1"))
VlnPlot(NPC_ALL_TRANSFORM, features = c("Gzmb"))
VlnPlot(NPC_ALL_TRANSFORM, features = c("Il1b"))
png("./output/UMAP_MARKER_Il6st.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Il6st"))
dev.off()
png("./output/UMAP_MARKER_Il6ra.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Il6ra"))
dev.off()
png("./output/UMAP_MARKER_Il6.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Il6"))
dev.off()
png("./output/UMAP_MARKER_Il10.png")
VlnPlot(NPC_ALL_TRANSFORM, features = c("Il10"))
dev.off()
sce <- as.SingleCellExperiment(DietSeurat(NPC_ALL_TRANSFORM))
sce
mouseRNA.ref <- celldex::MouseRNAseqData()
mouseRNA.main <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.main)
mouseRNA.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.fine)
table(mouseRNA.main$pruned.labels)
table(mouseRNA.fine$pruned.labels)
NPC_ALL_TRANSFORM@meta.data$mouseRNA.main <- mouseRNA.main$pruned.labels
NPC_ALL_TRANSFORM@meta.data$mouseRNA.fine <- mouseRNA.fine$pruned.labels
monaco.ref <- celldex::MonacoImmuneData()
monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
table(monaco.main$pruned.labels)
table(monaco.fine$pruned.labels)
NPC_ALL_TRANSFORM@meta.data$monaco.main <- monaco.main$pruned.labels
NPC_ALL_TRANSFORM@meta.data$monaco.fine <- monaco.fine$pruned.labels
ImmGenData.ref <- celldex::ImmGenData()
ImmGenData.main <- SingleR(test = sce,assay.type.test = 1,ref = ImmGenData.ref,labels = ImmGenData.ref$label.main)
ImmGenData.fine <- SingleR(test = sce,assay.type.test = 1,ref = ImmGenData.ref,labels = ImmGenData.ref$label.fine)
table(ImmGenData.main$pruned.labels)
table(ImmGenData.fine$pruned.labels)
NPC_ALL_TRANSFORM@meta.data$ImmGenData.main <- ImmGenData.main$pruned.labels
NPC_ALL_TRANSFORM@meta.data$ImmGenData.fine <- ImmGenData.fine$pruned.labels

saveRDS(NPC_ALL_TRANSFORM, file = "./output/NPC_ALL_TRANSFORM.rds")
rm(monaco.fine, monaco.main,monaco.ref, mouseRNA.fine,mouseRNA.main,mouseRNA.ref,ImmGenData.fine,ImmGenData.ref,ImmGenData.main)

NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "mouseRNA.fine")

png("./output/UMAP_Clusters_MouseRNASeq_FINE_Annotated.png")
DimPlot(NPC_ALL_TRANSFORM, label = F , repel = T, label.size = 3)
dev.off()
NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "mouseRNA.main")
png("./output/UMAP_Clusters_MouseRNASeq_main_Annotated.png")
DimPlot(NPC_ALL_TRANSFORM, label = F , repel = T, label.size = 3)
dev.off()


NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "monaco.fine")
png("./output/UMAP_Clusters_MonacoImmune_Annotated.png")
DimPlot(NPC_ALL_TRANSFORM, label = T , repel = T, label.size = 3) + NoLegend()
dev.off()


NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "ImmGenData.main")
png("./output/UMAP_Clusters_ImmGenData_Annotated.png")
DimPlot(NPC_ALL_TRANSFORM, label = F , repel = F, label.size = 3) 
dev.off()

NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "sex")
png("./output/UMAP_Clusters_sex.png")
DimPlot(NPC_ALL_TRANSFORM, label = F , repel = T, label.size = 3)
dev.off()

NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "stim")
png("./output/UMAP_Clusters_stim.png")
DimPlot(NPC_ALL_TRANSFORM, label = F , repel = T, label.size = 3)
dev.off()
png("./output/UMAP_Clusters.png")
NPC_ALL_TRANSFORM <- SetIdent(NPC_ALL_TRANSFORM, value = "seurat_clusters")
DimPlot(NPC_ALL_TRANSFORM, label = T , repel = T, label.size = 3)
dev.off()

