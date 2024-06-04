#This is the script to perform QC on FK44.1 scRNAseq COUNT data provided by BSF

############################ Install Packages #################################
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
#renv::install("SoupX")
#if (!require("BiocManager", quietly = TRUE))
#  renv::install("BiocManager")
#renv::install("gprofiler2")
#BiocManager::install(version = "3.18")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install('glmGamPoi')
#BiocManager::install('decontX')
#remotes::install_github("Moonerss/scrubletR")
#remotes::install_github('immunogenomics/presto')
#remotes::install_github("cysouw/qlcMatrix")
#remotes::install_github("Nawijn-Group-Bioinformatics/FastCAR")
############################# Load Libraries ###################################

library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(scrubletR)
library(gprofiler2)
#library(decontX)
library(SoupX)
library(celda)
library(scater)
library(stringr)
library(Matrix)
library(gridExtra)
library(qlcMatrix)
library(FastCAR)
library(pheatmap)
set.seed(42)
#####Load Output from CellRanger and generated Decon with FastCAR #######
#### FastCar iAL87----
cellExpressionFolder  = paste0("./00_raw_data/Liver_NPC_iAL87_transcriptome/filtered_feature_bc_matrix/")
fullMatrixFolder      = paste0("./00_raw_data/Liver_NPC_iAL87_transcriptome/raw_feature_bc_matrix/")
cellMatrix     = Read10X(cellExpressionFolder)
fullMatrix     = Read10X(fullMatrixFolder)
ambProfile = describe.ambient.RNA.sequence(fullMatrix = fullMatrix,  start = 10, stop = 1000,  by = 10, contaminationChanceCutoff = 0.05)
correctionEffectProfile_87 = describe.correction.effect(fullMatrix, cellMatrix, 50, 800, 100, 0.05)
x <-subset(correctionEffectProfile_87, subset=correctionEffectProfile_87$ctsScores>0)
annoying_genes <-rownames(x) 
emptyDropletCutoff = 2 ### smaller than 2 is not possible. it gives error
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_002     = remove.background(cellMatrix, ambientProfile)
NPC_87_002= CreateSeuratObject(cellMatrix_002)
NPC_87_002$sample <- "iAL87_002"
emptyDropletCutoff = 150
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_150     = remove.background(cellMatrix, ambientProfile)
NPC_87_150= CreateSeuratObject(cellMatrix_150)
NPC_87_150$sample <- "iAL87_150"
emptyDropletCutoff = 250
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_250     = remove.background(cellMatrix, ambientProfile)
NPC_87_250= CreateSeuratObject(cellMatrix_250)
NPC_87_250$sample <- "iAL87_250"
emptyDropletCutoff = 350
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_350     = remove.background(cellMatrix, ambientProfile)
NPC_87_350= CreateSeuratObject(cellMatrix_350)
emptyDropletCutoff = 450
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_450     = remove.background(cellMatrix, ambientProfile)
NPC_87_450= CreateSeuratObject(cellMatrix_450)
NPC_87_450$sample <- "iAL87_450"
emptyDropletCutoff = 550
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_550     = remove.background(cellMatrix, ambientProfile)
NPC_87_550= CreateSeuratObject(cellMatrix_550)
NPC_87_550$sample <- "iAL87_550"
#### Prepare Results from FastCAR 87 for Expression Analysis ####
NPC_87_CLEANED <- merge(NPC_87_002, y= c(NPC_87_150, NPC_87_250, NPC_87_350,NPC_87_450, NPC_87_550), add.cell.ids= c("002","150", "250","350","450","550"))
saveRDS(NPC_87_002, file = "./01_tidy_data/iAL87_FastCar_002.rds")
saveRDS(NPC_87_150, file = "./01_tidy_data/iAL87_FastCar_150.rds")
saveRDS(NPC_87_250, file = "./01_tidy_data/iAL87_FastCar_250.rds")
saveRDS(NPC_87_350, file = "./01_tidy_data/iAL87_FastCar_350.rds")
saveRDS(NPC_87_450, file = "./01_tidy_data/iAL87_FastCar_450.rds")
saveRDS(NPC_87_550, file = "./01_tidy_data/iAL87_FastCar_550.rds")
rm(NPC_87_002, NPC_87_150, NPC_87_250, NPC_87_350, NPC_87_450, NPC_87_550)
rm(cellMatrix_002,cellMatrix_150,cellMatrix_250,cellMatrix_350,cellMatrix_450,cellMatrix_550)

gc()
NPC_87_CLEANED$stim <- "TAM"
NPC_87_CLEANED <- NormalizeData(subset(NPC_87_CLEANED))
NPC_87_CLEANED <- FindVariableFeatures(NPC_87_CLEANED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_87_CLEANED)
NPC_87_CLEANED <- ScaleData(NPC_87_CLEANED, features = all.genes)
NPC_87_CLEANED <- RunPCA(NPC_87_CLEANED, features = VariableFeatures(object = NPC_87_CLEANED))
NPC_87_CLEANED <- FindNeighbors(NPC_87_CLEANED, dims = 1:11)
NPC_87_CLEANED <- FindClusters(NPC_87_CLEANED, resolution = 0.1)
NPC_87_CLEANED <- RunUMAP(NPC_87_CLEANED, dims = 1:11, verbose = F)
#### Perform Transformation of Data Set 87 ----
NPC_87_CLEANED <- SCTransform(NPC_87_CLEANED,  vst.flavor= "v2",method = "glmGamPoi", verbose = F) 
NPC_87_CLEANED <- RunPCA(NPC_87_CLEANED, verbose = F)
NPC_87_CLEANED <- RunUMAP(NPC_87_CLEANED, dims = 1:30, verbose = F)
NPC_87_CLEANED <- FindNeighbors(NPC_87_CLEANED, dims = 1:30, verbose = F)
NPC_87_CLEANED <- FindClusters(NPC_87_CLEANED, verbose = F,resolution = 0.1, save.SNN = TRUE)
##### Don't understand what I did here and why it was necessary 87 ----
DefaultAssay(NPC_87_CLEANED) <- "RNA"
NPC_87_CLEANED <- NormalizeData(NPC_87_CLEANED)
NPC_87_CLEANED <- FindVariableFeatures(NPC_87_CLEANED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_87_CLEANED)
NPC_87_CLEANED <- ScaleData(NPC_87_CLEANED, features = all.genes)
NPC_87_CLEANED <-JoinLayers(NPC_87_CLEANED)
#### Annotate Clusters using cellDex 87 ####
sce <- as.SingleCellExperiment(DietSeurat(NPC_87_CLEANED))
mouseRNA.ref <- celldex::MouseRNAseqData()
mouseRNA.main <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.main)
mouseRNA.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.fine)
NPC_87_CLEANED@meta.data$mouseRNA.main <- mouseRNA.main$pruned.labels
NPC_87_CLEANED@meta.data$mouseRNA.fine <- mouseRNA.fine$pruned.labels
NPC_87_CLEANED$celltype.stim <- paste(NPC_87_CLEANED$mouseRNA.main, NPC_87_CLEANED$stim, sep = "_")
NPC_87_CLEANED <- SetIdent(NPC_87_CLEANED, value = NPC_87_CLEANED@meta.data$mouseRNA.main)
png("./03_plots/QC_0_Ambient_FeaturePlot87_Saa2.png")
FeaturePlot(NPC_87_CLEANED, features = c("Saa2"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
png("./03_plots/QC_0_Ambient_FeaturePlot87_Saa1.png")
FeaturePlot(NPC_87_CLEANED, features = c("Saa1"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
png("./03_plots/QC_0_Ambient_FeaturePlot87_Alb.png")
FeaturePlot(NPC_87_CLEANED, features = c("Alb"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
for (i in annoying_genes){
  png(paste0("./03_plots/QC_0_Ambient_VlnPlot_87_",i,".png"))
  x=VlnPlot(NPC_87_CLEANED, features = paste0(i), split.by = "sample",group.by = "mouseRNA.main", pt.size = 1, combine = FALSE)
  print(x)
  dev.off()
}
interesting_genes=  c("Cd19","Cd3e","Clec4f","Cxcl1","Ptprb","Saa1","Alb","mt-Cytb")
for (i in interesting_genes){
  png(paste0("./03_plots/QC_0_Ambient_VlnPlot_87_",i,".png"))
  x=VlnPlot(NPC_87_CLEANED, features = paste0(i), split.by = "sample",group.by = "mouseRNA.main", pt.size = 1, combine = FALSE)
  print(x)
  dev.off()
  }

#### FastCar iAL88----
cellExpressionFolder  = paste0("./00_raw_data/Liver_NPC_iAL88_transcriptome/filtered_feature_bc_matrix/")
fullMatrixFolder      = paste0("./00_raw_data/Liver_NPC_iAL88_transcriptome/raw_feature_bc_matrix/")
cellMatrix     = Read10X(cellExpressionFolder)
fullMatrix     = Read10X(fullMatrixFolder)
ambProfile = describe.ambient.RNA.sequence(fullMatrix = fullMatrix,  start = 10, stop = 1000,  by = 10, contaminationChanceCutoff = 0.05)
correctionEffectProfile_88 = describe.correction.effect(fullMatrix, cellMatrix, 50, 800, 100, 0.05)
x <-subset(correctionEffectProfile_88, subset=correctionEffectProfile_88$ctsScores>0)
annoying_genes <-rownames(x) 
emptyDropletCutoff = 2 ### smaller than 2 is not possible. it gives error
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_002     = remove.background(cellMatrix, ambientProfile)
NPC_88_002= CreateSeuratObject(cellMatrix_002)
NPC_88_002$sample <- "iAL88_002"
emptyDropletCutoff = 150
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_150     = remove.background(cellMatrix, ambientProfile)
NPC_88_150= CreateSeuratObject(cellMatrix_150)
NPC_88_150$sample <- "iAL88_150"
emptyDropletCutoff = 250
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_250     = remove.background(cellMatrix, ambientProfile)
NPC_88_250= CreateSeuratObject(cellMatrix_250)
NPC_88_250$sample <- "iAL88_250"
emptyDropletCutoff = 350
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_350     = remove.background(cellMatrix, ambientProfile)
NPC_88_350= CreateSeuratObject(cellMatrix_350)
NPC_88_350$sample <- "iAL88_350"
emptyDropletCutoff = 450
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_450     = remove.background(cellMatrix, ambientProfile)
NPC_88_450= CreateSeuratObject(cellMatrix_450)
NPC_88_450$sample <- "iAL88_450"
emptyDropletCutoff = 550
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_550     = remove.background(cellMatrix, ambientProfile)
NPC_88_550= CreateSeuratObject(cellMatrix_550)
NPC_88_550$sample <- "iAL88_550"
#### Prepare Results from FastCAR 88 for Expression Analysis ####
NPC_88_CLEANED <- merge(NPC_88_002, y= c(NPC_88_150, NPC_88_250, NPC_88_350, NPC_88_450, NPC_88_550), add.cell.ids= c("2","150", "250","350","450","550"))

saveRDS(NPC_88_002, file = "./01_tidy_data/iAL88_FastCar_002.rds")
saveRDS(NPC_88_150, file = "./01_tidy_data/iAL88_FastCar_150.rds")
saveRDS(NPC_88_250, file = "./01_tidy_data/iAL88_FastCar_250.rds")
saveRDS(NPC_88_350, file = "./01_tidy_data/iAL88_FastCar_350.rds")
saveRDS(NPC_88_450, file = "./01_tidy_data/iAL88_FastCar_450.rds")
saveRDS(NPC_88_550, file = "./01_tidy_data/iAL88_FastCar_550.rds")
rm(cellMatrix_002,cellMatrix_150,cellMatrix_250,cellMatrix_350,cellMatrix_450,cellMatrix_550)

gc()
rm(NPC_88_002,NPC_88_150,NPC_88_250,NPC_88_350,NPC_88_450,NPC_88_550)
NPC_88_CLEANED$stim <- "TAM"
NPC_88_CLEANED <- NormalizeData(subset(NPC_88_CLEANED))
NPC_88_CLEANED <- FindVariableFeatures(NPC_88_CLEANED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_88_CLEANED)
NPC_88_CLEANED <- ScaleData(NPC_88_CLEANED, features = all.genes)
NPC_88_CLEANED <- RunPCA(NPC_88_CLEANED, features = VariableFeatures(object = NPC_88_CLEANED))
NPC_88_CLEANED <- FindNeighbors(NPC_88_CLEANED, dims = 1:11)
NPC_88_CLEANED <- FindClusters(NPC_88_CLEANED, resolution = 0.1)
NPC_88_CLEANED <- RunUMAP(NPC_88_CLEANED, dims = 1:11, verbose = F)
#### Perform Transformation of Data Set 88 ----
NPC_88_CLEANED <- SCTransform(NPC_88_CLEANED,  vst.flavor= "v2",method = "glmGamPoi", verbose = F) 
NPC_88_CLEANED <- RunPCA(NPC_88_CLEANED, verbose = F)
NPC_88_CLEANED <- RunUMAP(NPC_88_CLEANED, dims = 1:30, verbose = F)
NPC_88_CLEANED <- FindNeighbors(NPC_88_CLEANED, dims = 1:30, verbose = F)
NPC_88_CLEANED <- FindClusters(NPC_88_CLEANED, verbose = F,resolution = 0.1, save.SNN = TRUE)
##### Don't understand what I did here and why it was necessary 88 ----
DefaultAssay(NPC_88_CLEANED) <- "RNA"
NPC_88_CLEANED <- NormalizeData(NPC_88_CLEANED)
NPC_88_CLEANED <- FindVariableFeatures(NPC_88_CLEANED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_88_CLEANED)
NPC_88_CLEANED <- ScaleData(NPC_88_CLEANED, features = all.genes)
NPC_88_CLEANED <-JoinLayers(NPC_88_CLEANED)
#### Annotate Clusters using cellDex 88 ####
sce <- as.SingleCellExperiment(DietSeurat(NPC_88_CLEANED))
mouseRNA.ref <- celldex::MouseRNAseqData()
mouseRNA.main <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.main)
mouseRNA.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.fine)
NPC_88_CLEANED@meta.data$mouseRNA.main <- mouseRNA.main$pruned.labels
NPC_88_CLEANED@meta.data$mouseRNA.fine <- mouseRNA.fine$pruned.labels
NPC_88_CLEANED$celltype.stim <- paste(NPC_88_CLEANED$mouseRNA.main, NPC_88_CLEANED$stim, sep = "_")
NPC_88_CLEANED <- SetIdent(NPC_88_CLEANED, value = NPC_88_CLEANED@meta.data$mouseRNA.main)
png("./03_plots/QC_0_Ambient_FeaturePlot88_Saa2.png")
FeaturePlot(NPC_88_CLEANED, features = c("Saa2"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
png("./03_plots/QC_0_Ambient_FeaturePlot88_Saa1.png")
FeaturePlot(NPC_88_CLEANED, features = c("Saa1"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
png("./03_plots/QC_0_Ambient_FeaturePlot88_Alb.png")
FeaturePlot(NPC_88_CLEANED, features = c("Alb"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
for (i in annoying_genes){
  png(paste0("./03_plots/QC_0_Ambient_VlnPlot_88_",i,".png"))
  x=VlnPlot(NPC_88_CLEANED, features = paste0(i), split.by = "sample",group.by = "mouseRNA.main", pt.size = 1, combine = FALSE)
  print(x)
  dev.off()
}
interesting_genes=  c("Cd19","Cd3e","Clec4f","Cxcl1","Ptprb","Saa1","Alb","mt-Cytb")
for (i in interesting_genes){
  png(paste0("./03_plots/QC_0_Ambient_VlnPlot_88_",i,".png"))
  x=VlnPlot(NPC_88_CLEANED, features = paste0(i), split.by = "sample",group.by = "mouseRNA.main", pt.size = 1, combine = FALSE)
  print(x)
  dev.off()
}

#### FastCar iAL91----
cellExpressionFolder  = paste0("./00_raw_data/Liver_NPC_iAL91_transcriptome/filtered_feature_bc_matrix/")
fullMatrixFolder      = paste0("./00_raw_data/Liver_NPC_iAL91_transcriptome/raw_feature_bc_matrix/")
cellMatrix     = Read10X(cellExpressionFolder)
fullMatrix     = Read10X(fullMatrixFolder)
ambProfile = describe.ambient.RNA.sequence(fullMatrix = fullMatrix,  start = 10, stop = 1000,  by = 10, contaminationChanceCutoff = 0.05)
correctionEffectProfile_91 = describe.correction.effect(fullMatrix, cellMatrix, 50, 800, 100, 0.05)
x <-subset(correctionEffectProfile_91, subset=correctionEffectProfile_91$ctsScores>0)
annoying_genes <-rownames(x) 
emptyDropletCutoff = 2 ### smaller than 2 is not possible. it gives error
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_002     = remove.background(cellMatrix, ambientProfile)
NPC_91_002= CreateSeuratObject(cellMatrix_002)
NPC_91_002$sample <- "iAL91_002"
emptyDropletCutoff = 150
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_150     = remove.background(cellMatrix, ambientProfile)
NPC_91_150= CreateSeuratObject(cellMatrix_150)
NPC_91_150$sample <- "iAL91_150"
emptyDropletCutoff = 250
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_250     = remove.background(cellMatrix, ambientProfile)
NPC_91_250= CreateSeuratObject(cellMatrix_250)
NPC_91_250$sample <- "iAL91_250"
emptyDropletCutoff = 350
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_350     = remove.background(cellMatrix, ambientProfile)
NPC_91_350= CreateSeuratObject(cellMatrix_350)
NPC_91_350$sample <- "iAL91_350"
emptyDropletCutoff = 450
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_450     = remove.background(cellMatrix, ambientProfile)
NPC_91_450= CreateSeuratObject(cellMatrix_450)
NPC_91_450$sample <- "iAL91_450"
emptyDropletCutoff = 550
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_550     = remove.background(cellMatrix, ambientProfile)
NPC_91_550= CreateSeuratObject(cellMatrix_550)
NPC_91_550$sample <- "iAL91_550"
#### Prepare Results from FastCAR 91 for Expression Analysis ####


saveRDS(NPC_91_002, file = "./01_tidy_data/iAL91_FastCar_002.rds")
saveRDS(NPC_91_150, file = "./01_tidy_data/iAL91_FastCar_150.rds")
saveRDS(NPC_91_250, file = "./01_tidy_data/iAL91_FastCar_250.rds")
saveRDS(NPC_91_350, file = "./01_tidy_data/iAL91_FastCar_350.rds")
saveRDS(NPC_91_450, file = "./01_tidy_data/iAL91_FastCar_450.rds")
saveRDS(NPC_91_550, file = "./01_tidy_data/iAL91_FastCar_550.rds")
NPC_91_CLEANED <- merge(NPC_91_002, y= c(NPC_91_150, NPC_91_250, NPC_91_350,NPC_91_450, NPC_91_550), add.cell.ids= c("2","150","250","350","450","550"))

rm(NPC_91_002, NPC_91_150,NPC_91_250,NPC_91_350,NPC_91_450,NPC_91_550)
# NPC_91_002 <- readRDS("./01_tidy_data/iAL91_FastCar_002.rds")
# NPC_91_150 <- readRDS("./01_tidy_data/iAL91_FastCar_150.rds")
# NPC_91_250 <- readRDS("./01_tidy_data/iAL91_FastCar_250.rds")
# NPC_91_350 <- readRDS("./01_tidy_data/iAL91_FastCar_350.rds")
# NPC_91_450 <- readRDS("./01_tidy_data/iAL91_FastCar_450.rds")
# NPC_91_550 <- readRDS("./01_tidy_data/iAL91_FastCar_550.rds")
rm(cellMatrix_002,cellMatrix_150,cellMatrix_250,cellMatrix_350,cellMatrix_450,cellMatrix_550)

gc()
NPC_91_CLEANED$stim <- "TAM"
NPC_91_CLEANED <- NormalizeData(subset(NPC_91_CLEANED))
NPC_91_CLEANED <- FindVariableFeatures(NPC_91_CLEANED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_91_CLEANED)
NPC_91_CLEANED <- ScaleData(NPC_91_CLEANED, features = all.genes)
NPC_91_CLEANED <- RunPCA(NPC_91_CLEANED, features = VariableFeatures(object = NPC_91_CLEANED))
NPC_91_CLEANED <- FindNeighbors(NPC_91_CLEANED, dims = 1:11)
NPC_91_CLEANED <- FindClusters(NPC_91_CLEANED, resolution = 0.1)
NPC_91_CLEANED <- RunUMAP(NPC_91_CLEANED, dims = 1:11, verbose = F)
#### Perform Transformation of Data Set 91 ----
NPC_91_CLEANED <- SCTransform(NPC_91_CLEANED,  vst.flavor= "v2",method = "glmGamPoi", verbose = F) 
NPC_91_CLEANED <- RunPCA(NPC_91_CLEANED, verbose = F)
NPC_91_CLEANED <- RunUMAP(NPC_91_CLEANED, dims = 1:30, verbose = F)
NPC_91_CLEANED <- FindNeighbors(NPC_91_CLEANED, dims = 1:30, verbose = F)
NPC_91_CLEANED <- FindClusters(NPC_91_CLEANED, verbose = F,resolution = 0.1, save.SNN = TRUE)
##### Don't understand what I did here and why it was necessary 91 ----
DefaultAssay(NPC_91_CLEANED) <- "RNA"
NPC_91_CLEANED <- NormalizeData(NPC_91_CLEANED)
NPC_91_CLEANED <- FindVariableFeatures(NPC_91_CLEANED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_91_CLEANED)
NPC_91_CLEANED <- ScaleData(NPC_91_CLEANED, features = all.genes)
NPC_91_CLEANED <-JoinLayers(NPC_91_CLEANED)
saveRDS(NPC_91_CLEANED, "./01_tidy_data/NPC_91_CLEANED.rds")
#NPC_91_CLEANED <-readRDS( "./01_tidy_data/NPC_91_CLEANED.rds")
#### Annotate Clusters using cellDex 91 ####
sce <- as.SingleCellExperiment(DietSeurat(NPC_91_CLEANED))
mouseRNA.ref <- celldex::MouseRNAseqData()
mouseRNA.main <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.main)
mouseRNA.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.fine)
NPC_91_CLEANED@meta.data$mouseRNA.main <- mouseRNA.main$pruned.labels
NPC_91_CLEANED@meta.data$mouseRNA.fine <- mouseRNA.fine$pruned.labels
NPC_91_CLEANED$celltype.stim <- paste(NPC_91_CLEANED$mouseRNA.main, NPC_91_CLEANED$stim, sep = "_")
NPC_91_CLEANED <- SetIdent(NPC_91_CLEANED, value = NPC_91_CLEANED@meta.data$mouseRNA.main)
png("./03_plots/QC_0_Ambient_FeaturePlot91_Saa2.png")
FeaturePlot(NPC_91_CLEANED, features = c("Saa2"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
png("./03_plots/QC_0_Ambient_FeaturePlot91_Saa1.png")
FeaturePlot(NPC_91_CLEANED, features = c("Saa1"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
png("./03_plots/QC_0_Ambient_FeaturePlot91_Alb.png")
FeaturePlot(NPC_91_CLEANED, features = c("Alb"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
for (i in annoying_genes){
  png(paste0("./03_plots/QC_0_Ambient_VlnPlot_91_",i,".png"))
  x=VlnPlot(NPC_91_CLEANED, features = paste0(i), split.by = "sample",group.by = "mouseRNA.main", pt.size = 1, combine = FALSE)
  print(x)
  dev.off()
}
interesting_genes=  c("Cd19","Cd3e","Clec4f","Cxcl1","Ptprb","Saa1","Alb","mt-Cytb")
for (i in interesting_genes){
  png(paste0("./03_plots/QC_0_Ambient_VlnPlot_91_",i,".png"))
  x=VlnPlot(NPC_91_CLEANED, features = paste0(i), split.by = "sample",group.by = "mouseRNA.main", pt.size = 1, combine = FALSE)
  print(x)
  dev.off()
}
rm(NPC_91_CLEANED)
#### FastCar iAL92----
cellExpressionFolder  = paste0("./00_raw_data/Liver_NPC_iAL92_transcriptome/filtered_feature_bc_matrix/")
fullMatrixFolder      = paste0("./00_raw_data/Liver_NPC_iAL92_transcriptome/raw_feature_bc_matrix/")
cellMatrix     = Read10X(cellExpressionFolder)
fullMatrix     = Read10X(fullMatrixFolder)
ambProfile = describe.ambient.RNA.sequence(fullMatrix = fullMatrix,  start = 10, stop = 1000,  by = 10, contaminationChanceCutoff = 0.05)
correctionEffectProfile_92 = describe.correction.effect(fullMatrix, cellMatrix, 50, 800, 100, 0.05)
x <-subset(correctionEffectProfile_92, subset=correctionEffectProfile_92$ctsScores>0)
annoying_genes <-rownames(x) 
emptyDropletCutoff = 2 ### smaller than 2 is not possible. it gives error
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_002     = remove.background(cellMatrix, ambientProfile)
NPC_92_002= CreateSeuratObject(cellMatrix_002)
NPC_92_002$sample <- "iAL92_002"
emptyDropletCutoff = 150
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_150     = remove.background(cellMatrix, ambientProfile)
NPC_92_150= CreateSeuratObject(cellMatrix_150)
NPC_92_150$sample <- "iAL92_150"
emptyDropletCutoff = 250
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_250     = remove.background(cellMatrix, ambientProfile)
NPC_92_250= CreateSeuratObject(cellMatrix_250)
NPC_92_250$sample <- "iAL92_250"
emptyDropletCutoff = 350
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_350     = remove.background(cellMatrix, ambientProfile)
NPC_92_350= CreateSeuratObject(cellMatrix_350)
NPC_92_350$sample <- "iAL92_350"
emptyDropletCutoff = 450
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_450     = remove.background(cellMatrix, ambientProfile)
NPC_92_450= CreateSeuratObject(cellMatrix_450)
NPC_92_450$sample <- "iAL92_450"
emptyDropletCutoff = 550
contaminationChanceCutoff = 0.05
ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix_550     = remove.background(cellMatrix, ambientProfile)
NPC_92_550= CreateSeuratObject(cellMatrix_550)
NPC_92_550$sample <- "iAL92_550"
#### Prepare Results from FastCAR 92 for Expression Analysis ####
NPC_92_CLEANED <- merge(NPC_92_002, y= c(NPC_92_150, NPC_92_250, NPC_92_350, NPC_92_450,NPC_92_550), add.cell.ids= c("2","150", "250","350","450","550"))

saveRDS(NPC_92_002, file = "./01_tidy_data/iAL92_FastCar_002.rds")
saveRDS(NPC_92_150, file = "./01_tidy_data/iAL92_FastCar_150.rds")
saveRDS(NPC_92_250, file = "./01_tidy_data/iAL92_FastCar_250.rds")
saveRDS(NPC_92_350, file = "./01_tidy_data/iAL92_FastCar_350.rds")
saveRDS(NPC_92_450, file = "./01_tidy_data/iAL92_FastCar_450.rds")
saveRDS(NPC_92_550, file = "./01_tidy_data/iAL92_FastCar_550.rds")
rm(NPC_92_002,NPC_92_150,NPC_92_250,NPC_92_350,NPC_92_450,NPC_92_550)
rm(cellMatrix_002,cellMatrix_150,cellMatrix_250,cellMatrix_350,cellMatrix_450,cellMatrix_550)

gc()
NPC_92_CLEANED$stim <- "TAM"
NPC_92_CLEANED <- NormalizeData(subset(NPC_92_CLEANED))
NPC_92_CLEANED <- FindVariableFeatures(NPC_92_CLEANED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_92_CLEANED)
NPC_92_CLEANED <- ScaleData(NPC_92_CLEANED, features = all.genes)
NPC_92_CLEANED <- RunPCA(NPC_92_CLEANED, features = VariableFeatures(object = NPC_92_CLEANED))
NPC_92_CLEANED <- FindNeighbors(NPC_92_CLEANED, dims = 1:11)
NPC_92_CLEANED <- FindClusters(NPC_92_CLEANED, resolution = 0.1)
NPC_92_CLEANED <- RunUMAP(NPC_92_CLEANED, dims = 1:11, verbose = F)
#### Perform Transformation of Data Set 92 ----
NPC_92_CLEANED <- SCTransform(NPC_92_CLEANED,  vst.flavor= "v2",method = "glmGamPoi", verbose = F) 
NPC_92_CLEANED <- RunPCA(NPC_92_CLEANED, verbose = F)
NPC_92_CLEANED <- RunUMAP(NPC_92_CLEANED, dims = 1:30, verbose = F)
NPC_92_CLEANED <- FindNeighbors(NPC_92_CLEANED, dims = 1:30, verbose = F)
NPC_92_CLEANED <- FindClusters(NPC_92_CLEANED, verbose = F,resolution = 0.1, save.SNN = TRUE)
##### Don't understand what I did here and why it was necessary 92 ----
DefaultAssay(NPC_92_CLEANED) <- "RNA"
NPC_92_CLEANED <- NormalizeData(NPC_92_CLEANED)
NPC_92_CLEANED <- FindVariableFeatures(NPC_92_CLEANED, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NPC_92_CLEANED)
NPC_92_CLEANED <- ScaleData(NPC_92_CLEANED, features = all.genes)
NPC_92_CLEANED <-JoinLayers(NPC_92_CLEANED)
#### Annotate Clusters using cellDex 92 ####
sce <- as.SingleCellExperiment(DietSeurat(NPC_92_CLEANED))
mouseRNA.ref <- celldex::MouseRNAseqData()
mouseRNA.main <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.main)
mouseRNA.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouseRNA.ref,labels = mouseRNA.ref$label.fine)
NPC_92_CLEANED@meta.data$mouseRNA.main <- mouseRNA.main$pruned.labels
NPC_92_CLEANED@meta.data$mouseRNA.fine <- mouseRNA.fine$pruned.labels
NPC_92_CLEANED$celltype.stim <- paste(NPC_92_CLEANED$mouseRNA.main, NPC_92_CLEANED$stim, sep = "_")
NPC_92_CLEANED <- SetIdent(NPC_92_CLEANED, value = NPC_92_CLEANED@meta.data$mouseRNA.main)
png("./03_plots/QC_0_Ambient_FeaturePlot92_Saa2.png")
FeaturePlot(NPC_92_CLEANED, features = c("Saa2"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
png("./03_plots/QC_0_Ambient_FeaturePlot92_Saa1.png")
FeaturePlot(NPC_92_CLEANED, features = c("Saa1"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
png("./03_plots/QC_0_Ambient_FeaturePlot92_Alb.png")
FeaturePlot(NPC_92_CLEANED, features = c("Alb"), split.by = "sample", max.cutoff = 5,cols = c("grey", "red"))
dev.off()
for (i in annoying_genes){
  png(paste0("./03_plots/QC_0_Ambient_VlnPlot_92_",i,".png"))
  x=VlnPlot(NPC_92_CLEANED, features = paste0(i), split.by = "sample",group.by = "mouseRNA.main", pt.size = 1, combine = FALSE)
  print(x)
  dev.off()
}
interesting_genes=  c("Cd19","Cd3e","Clec4f","Cxcl1","Ptprb","Saa1","Alb","mt-Cytb")
for (i in interesting_genes){
  png(paste0("./03_plots/QC_0_Ambient_VlnPlot_92_",i,".png"))
  x=VlnPlot(NPC_92_CLEANED, features = paste0(i), split.by = "sample",group.by = "mouseRNA.main", pt.size = 1, combine = FALSE)
  print(x)
  dev.off()
}
rm(NPC_92_CLEANED)






