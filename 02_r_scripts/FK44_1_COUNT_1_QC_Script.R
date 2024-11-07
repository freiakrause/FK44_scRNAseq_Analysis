#This is the script to perform QC on FK44.1 scRNAseq COUNT data provided by BSF

########################### Install Packages #################################
# install.packages("renv")
# renv::init()    #creates project library,  lockfile and Rprofile
# renv::snapshot() # updates lockfile abou current-used packages-> can be shared and reproduced by others when restore() is used
# renv::restore() # installs the exact same versions of packages determined in lockfile
# renv::update()) #
# renv::history() #
# renv::install("usethis")
# usethis::create_github_token()
#renv::install("gitcreds")
#gitcreds::gitcreds_set()
#renv::install("remotes")
#renv::install("ggplot2", "dplyr" ,"RColorBrewer") # package "SingleR" required by tutorial but not found when trying isntallation
#renv::install("igraph","Seurat","SoupX") # für das scheiß igraph (benötigt für seurat) nach tausend jahren troublehsooting gefunden: ich brauche: sudo apt install build-essential gfortran UND sudo apt install build-essential gfortran, dann gehts
# if (!require("BiocManager", quietly = TRUE))
#  renv::install("BiocManager")
# renv::install("gprofiler2")
#BiocManager::install(version = "3.18")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install('glmGamPoi')
#remotes::install_github("Moonerss/scrubletR")
#remotes::install_github('immunogenomics/presto')
#remotes::install_github("cysouw/qlcMatrix")
############################# Load Libraries ###################################

rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(gprofiler2)
source("02_r_scripts/Function_Malat1.R")
set.seed(42)


########################## Kategorien Stimlulation und Sex hinzufügen, evtl noch age? ####################

NPC_87 <-Read10X(data.dir = "./00_raw_data/biomedical-sequencing.at/projects/BSA_0873_FK44_1_LiverMet_A_1_1_51ddbfd228ec40b096e110101b219cb0/COUNT/Liver_NPC_iAL87_transcriptome/filtered_feature_bc_matrix/")%>%
CreateSeuratObject(project="FK44.1",min.cells=10, min.features=200) 
NPC_88 <-Read10X(data.dir = "./00_raw_data/biomedical-sequencing.at/projects/BSA_0873_FK44_1_LiverMet_A_1_1_51ddbfd228ec40b096e110101b219cb0/COUNT/Liver_NPC_iAL88_transcriptome/filtered_feature_bc_matrix/")%>%
  CreateSeuratObject(project="FK44.1",min.cells=10, min.features=200) 
NPC_91 <-Read10X(data.dir = "./00_raw_data/biomedical-sequencing.at/projects/BSA_0873_FK44_1_LiverMet_A_1_1_51ddbfd228ec40b096e110101b219cb0/COUNT/Liver_NPC_iAL91_transcriptome/filtered_feature_bc_matrix/")%>%
  CreateSeuratObject(project="FK44.1",min.cells=10, min.features=200) 
NPC_92 <-Read10X(data.dir = "./00_raw_data/biomedical-sequencing.at/projects/BSA_0873_FK44_1_LiverMet_A_1_1_51ddbfd228ec40b096e110101b219cb0/COUNT/Liver_NPC_iAL92_transcriptome/filtered_feature_bc_matrix/")%>%
  CreateSeuratObject(project="FK44.1",min.cells=10, min.features=200) 
NPC_87$stim <- "TAM"
NPC_88$stim <- "EtOH"
NPC_91$stim <- "TAM"
NPC_92$stim <- "EtOH"
NPC_87$sex <- "female"
NPC_88$sex <- "female"
NPC_91$sex <- "male"
NPC_92$sex <- "male"
NPC_87$sample <- "iAL87"
NPC_88$sample <- "iAL88"
NPC_91$sample <- "iAL91"
NPC_92$sample <- "iAL92"

############################### After Applying ambient RNA Decontamination Set Up Standard QC ##################################
NPC_87[["percent.mt"]] <- PercentageFeatureSet(NPC_87, pattern = "^mt-") 
NPC_87[["percent.rb"]] <- PercentageFeatureSet(NPC_87, pattern = "Rp[sl]")
NPC_88[["percent.mt"]] <- PercentageFeatureSet(NPC_88, pattern = "^mt-") 
NPC_88[["percent.rb"]] <- PercentageFeatureSet(NPC_88, pattern = "Rp[sl]")
NPC_91[["percent.mt"]] <- PercentageFeatureSet(NPC_91, pattern = "^mt-") 
NPC_91[["percent.rb"]] <- PercentageFeatureSet(NPC_91, pattern = "Rp[sl]")
NPC_92[["percent.mt"]] <- PercentageFeatureSet(NPC_92, pattern = "^mt-") 
NPC_92[["percent.rb"]] <- PercentageFeatureSet(NPC_92, pattern = "Rp[sl]")


#Add Conclusions from meta data as QC6 columns QC6 from PMID: 36901774; UMI 500<UMI<40000, Genes 500<GENES<6000, MT MT<25 ----
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nFeature_RNA > 6000,'Doublet','Pass')
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 300 & NPC_87@meta.data$QC6 == 'Pass' , paste('Low_nFeature',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 300 & NPC_87@meta.data$QC6 != 'Pass' & NPC_87@meta.data$QC6 != 'Low_nFeature',paste('Low_nFeature',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$percent.mt > 10 & NPC_87@meta.data$QC6 == 'Pass','High_MT',NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nCount_RNA > 40000 & NPC_87@meta.data$QC6 == 'Pass','High_UMI',NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nCount_RNA < 500 & NPC_87@meta.data$QC6 == 'Pass','Low_UMI',NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(NPC_87@meta.data$nFeature_RNA < 300 & NPC_87@meta.data$QC6 != 'Pass'& NPC_87@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(log10(NPC_87@meta.data$nFeature_RNA) / log10(NPC_87@meta.data$nCount_RNA)>= 0.8 & NPC_87@meta.data$QC6 != 'Pass'& NPC_87@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
NPC_87[['QC6']] <- ifelse(log10(NPC_87@meta.data$nFeature_RNA) / log10(NPC_87@meta.data$nCount_RNA)< 0.8 & NPC_87@meta.data$QC6 == 'Pass'& NPC_87@meta.data$QC6 != 'High_MT',paste('low complex',NPC_87@meta.data$QC6,sep = ','),NPC_87@meta.data$QC6)
table(NPC_87[['QC6']])

NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nFeature_RNA > 6000,'Doublet','Pass')
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 300 & NPC_88@meta.data$QC6 == 'Pass' , paste('Low_nFeature',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 300 & NPC_88@meta.data$QC6 != 'Pass' & NPC_88@meta.data$QC6 != 'Low_nFeature',paste('Low_nFeature',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$percent.mt > 10 & NPC_88@meta.data$QC6 == 'Pass','High_MT',NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nCount_RNA > 40000 & NPC_88@meta.data$QC6 == 'Pass','High_UMI',NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nCount_RNA < 500 & NPC_88@meta.data$QC6 == 'Pass','Low_UMI',NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(NPC_88@meta.data$nFeature_RNA < 300 & NPC_88@meta.data$QC6 != 'Pass'& NPC_88@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(log10(NPC_88@meta.data$nFeature_RNA) / log10(NPC_88@meta.data$nCount_RNA)>= 0.8 & NPC_88@meta.data$QC6 != 'Pass'& NPC_88@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
NPC_88[['QC6']] <- ifelse(log10(NPC_88@meta.data$nFeature_RNA) / log10(NPC_88@meta.data$nCount_RNA)< 0.8 & NPC_88@meta.data$QC6 == 'Pass'& NPC_88@meta.data$QC6 != 'High_MT',paste('low complex',NPC_88@meta.data$QC6,sep = ','),NPC_88@meta.data$QC6)
table(NPC_88[['QC6']])

NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nFeature_RNA > 6000,'Doublet','Pass')
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 300 & NPC_91@meta.data$QC6 == 'Pass' , paste('Low_nFeature',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 300 & NPC_91@meta.data$QC6 != 'Pass' & NPC_91@meta.data$QC6 != 'Low_nFeature',paste('Low_nFeature',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$percent.mt > 10 & NPC_91@meta.data$QC6 == 'Pass','High_MT',NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nCount_RNA > 40000 & NPC_91@meta.data$QC6 == 'Pass','High_UMI',NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nCount_RNA < 500 & NPC_91@meta.data$QC6 == 'Pass','Low_UMI',NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(NPC_91@meta.data$nFeature_RNA < 300 & NPC_91@meta.data$QC6 != 'Pass'& NPC_91@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(log10(NPC_91@meta.data$nFeature_RNA) / log10(NPC_91@meta.data$nCount_RNA)>= 0.8 & NPC_91@meta.data$QC6 != 'Pass'& NPC_91@meta.data$QC6 != 'High_MT',paste('High_MT',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
NPC_91[['QC6']] <- ifelse(log10(NPC_91@meta.data$nFeature_RNA) / log10(NPC_91@meta.data$nCount_RNA)< 0.8 & NPC_91@meta.data$QC6 == 'Pass'& NPC_91@meta.data$QC6 != 'High_MT',paste('low complex',NPC_91@meta.data$QC6,sep = ','),NPC_91@meta.data$QC6)
table(NPC_91[['QC6']])

NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nFeature_RNA > 6000,'Doublet','Pass')
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 300 & NPC_92@meta.data$QC6 == 'Pass' , 
                          paste('Low_nFeature',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 300 & NPC_92@meta.data$QC6 != 'Pass' & NPC_92@meta.data$QC6 != 'Low_nFeature',
                          paste('Low_nFeature',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$percent.mt > 10 & NPC_92@meta.data$QC6 == 'Pass','High_MT',NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nCount_RNA > 40000 & NPC_92@meta.data$QC6 == 'Pass','High_UMI',NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nCount_RNA < 500 & NPC_92@meta.data$QC6 == 'Pass','Low_UMI',NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(NPC_92@meta.data$nFeature_RNA < 300 & NPC_92@meta.data$QC6 != 'Pass'& NPC_92@meta.data$QC6 != 'High_MT',
                          paste('High_MT',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(log10(NPC_92@meta.data$nFeature_RNA) / log10(NPC_92@meta.data$nCount_RNA)>= 0.8 &
                            NPC_92@meta.data$QC6 != 'Pass'&
                            NPC_92@meta.data$QC6 != 'High_MT',
                          paste('High_MT',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
NPC_92[['QC6']] <- ifelse(log10(NPC_92@meta.data$nFeature_RNA) / log10(NPC_92@meta.data$nCount_RNA)< 0.8 &
                            NPC_92@meta.data$QC6 == 'Pass'&
                            NPC_92@meta.data$QC6 != 'High_MT',
                          paste('low complex',NPC_92@meta.data$QC6,sep = ','),NPC_92@meta.data$QC6)
table(NPC_92[['QC6']])

#######################  Generate the Object "metadata" to plot Things more easily ############################################################################################
#https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

NPC_87@meta.data
NPC_88@meta.data
metadata <- rbind(NPC_87@meta.data,NPC_88@meta.data, NPC_91@meta.data, NPC_92@meta.data)
metadata$log10GenesPerUMI <- log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)
metadata$cells <- rownames(metadata)
metadata <- metadata %>%  dplyr::rename(nUMI = nCount_RNA,nGene = nFeature_RNA)
metadata[["samples_QC"]]<-ifelse(metadata$QC6 == "Pass",paste0("PASSED"),paste0("FAILED"))
metadata_QC6 <-subset(metadata, QC6 == "Pass")

########################### Visualizations of the QC Parameters on Data wo QC  #############################################################################
# Visualize number of cell per sample----
png(paste0("./03_plots/1_QC/QC_1_noQC_Cells_per_sample_bar.png"))
x <-metadata %>%
  ggplot(aes(x=sample, fill=samples_QC)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells")
print(x)
dev.off()

# Visualize number of UMIS per sample violin----
png(paste0("./03_plots/1_QC/QC_1_noQC_UMIs_per_sample_violin.png"))
x <-metadata %>%
  ggplot(aes(x=sample, fill=sample, y =nUMI)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 500)+
  geom_hline(yintercept = 40000)+
  ggtitle("nUMI")
print(x)
dev.off()

# Visualize number of UMIS per sample density plot----
png(paste0("./03_plots/1_QC/QC_1_noQC_UMIs_per_sample_density.png"))
x <- metadata%>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+
  geom_vline(xintercept = 40000)
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample violin----
png(paste0("./03_plots/1_QC/QC_1_noQC_Gene(log10)_per_sample_violin.png"))
x <- metadata %>%
  ggplot(aes(x=sample, fill=sample, y =log10(nGene))) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  geom_hline(yintercept = log10(6000))+
  geom_hline(yintercept = log10(300))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Log10 GENE")
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample BoxPlot----
png(paste0("./03_plots/1_QC/QC_1_noQC_Gene(log10)_per_sample_box.png"))
x <- metadata %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = log10(6000))+
  geom_hline(yintercept = log10(300))+
  ggtitle("log10 NGenes")
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample DensityPlot----
png(paste0("./03_plots/1_QC/QC_1_noQC_Gene(log10)_per_sample_density.png"))
x <- metadata %>%
  ggplot(aes(color=sample, x=log10(nGene), fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10()+
  geom_vline(xintercept = log10(6000))+
  geom_vline(xintercept = log10(300))+
  ggtitle("log10 NGenes")
print(x)
dev.off()

# Visualize number of Genes per sample ViolinPlot----
png(paste0("./03_plots/1_QC/QC_1_noQC_Gene_per_sample_violin.png"))
x <- metadata %>%
  ggplot(aes(x=sample, fill=sample, y =nGene)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 6000)+
  geom_hline(yintercept = 300)+
  ggtitle("nGENE")
print(x)
dev.off()

# Visualize number of Genes per sample BoxPlot----
png(paste0("./03_plots/1_QC/QC_1_noQC_Gene_per_sample_box.png"))
x <- metadata %>%
  ggplot(aes(x=sample, y=nGene, fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 6000)+
  geom_hline(yintercept = 300)+
  ggtitle("nGenes")
print(x)
dev.off()

# Visualize number of Genes per sample DensityPlot----
png(paste0("./03_plots/1_QC/QC_1_noQC_Gene_per_sample_density.png"))
x <- metadata %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ggtitle("nGenes")+
  geom_vline(xintercept = 6000)+
  geom_vline(xintercept = 300)
print(x)
dev.off()


# Visualize number of percent.met per sample Violin plot----
png(paste0("./03_plots/1_QC/QC_1_noQC_Mito_per_sample_violin.png"))
x <- metadata %>%
  ggplot(aes(x=sample, fill=sample, y =percent.mt)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 10)+
  ggtitle("Percent mt")
print(x)
dev.off()

# Visualize number of percent.mt per sample Density plot----
png(paste0("./03_plots/1_QC/QC_1_noQC_Mito_per_sample_density.png"))
x <- metadata %>%
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 10)
print(x)
dev.off()

# Visualize number of percent.rb per sample Violin plot----
png(paste0("./03_plots/1_QC/QC_1_noQC_Rb_per_sample_violin.png"))
x <- metadata %>%
  ggplot(aes(x=sample, fill=sample, y =percent.rb)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Percent rb")
print(x)
dev.off()

# Visualize number of Complexity of cells/sample as nGene vs nUMI----
png(paste0("./03_plots/1_QC/QC_1_noQC_Gene_vs_UMI_persample.png"))
x <- metadata%>%
  ggplot(aes(x=nUMI, y=nGene)) +
  geom_point(aes(colour=percent.mt)) +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 40000) +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 6000) +
  facet_wrap(~sample)
print(x)
dev.off()

# Visualize number of Complexity of cells/sample as log10Genes/UMI----
png(paste0("./03_plots/1_QC/QC_1_noQC_Gene perUMI_density.png"))
x <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
print(x)
dev.off()

# ########################### Visualizations of the QC Parameters on Data with QC6############################################################################
# Visualize number of cell per sample----
png(paste0("./03_plots/1_QC/QC_1_QC6_Cells_per_sample_bar.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells")
print(x)
dev.off()

# Visualize number of UMIS per sample violin----
png(paste0("./03_plots/1_QC/QC_1_QC6_UMIs_per_sample_violin.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =nUMI)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 500)+
  geom_hline(yintercept = 40000)+
  ggtitle("nUMI")
print(x)
dev.off()

# Visualize number of UMIS per sample density plot----
png(paste0("./03_plots/1_QC/QC_1_QC6_UMIs_per_sample_density.png"))
x <- metadata_QC6%>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+
  geom_vline(xintercept = 40000)
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample violin----
png(paste0("./03_plots/1_QC/QC_1_QC6_Gene(log10)_per_sample_violin.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =log10(nGene))) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  geom_hline(yintercept = log10(6000))+
  geom_hline(yintercept = log10(300))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Log10 GENE")
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample BoxPlot----
png(paste0("./03_plots/1_QC/QC_1_QC6_Gene(log10)_per_sample_box.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = log10(6000))+
  geom_hline(yintercept = log10(300))+
  ggtitle("log10 NGenes")
print(x)
dev.off()

# Visualize number of Log10(Genes) per sample DensityPlot----
png(paste0("./03_plots/1_QC/QC_1_QC6_Gene(log10)_per_sample_density.png"))
x <- metadata_QC6 %>%
  ggplot(aes(color=sample, x=log10(nGene), fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10()+
  geom_vline(xintercept = log10(6000))+
  geom_vline(xintercept = log10(300))+
  ggtitle("log10 NGenes")
print(x)
dev.off()


# Visualize number of Genes per sample ViolinPlot----
png(paste0("./03_plots/1_QC/QC_1_QC6_Gene_per_sample_violin.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =nGene)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 6000)+
  geom_hline(yintercept = 300)+
  ggtitle("nGENE")
print(x)
dev.off()

# Visualize number of Genes per sample BoxPlot----
png(paste0("./03_plots/1_QC/QC_1_QC6_Gene_per_sample_box.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, y=nGene, fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 6000)+
  geom_hline(yintercept = 300)+
  ggtitle("nGenes")
print(x)
dev.off()

# Visualize number of Genes per sample DensityPlot----
png(paste0("./03_plots/1_QC/QC_1_QC6_Gene_per_sample_density.png"))
x <- metadata_QC6 %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ggtitle("nGenes")+
  geom_vline(xintercept = 6000)+
  geom_vline(xintercept = 300)
print(x)
dev.off()


# Visualize number of percent.met per sample Violin plot----
png(paste0("./03_plots/1_QC/QC_1_QC6_Mito_per_sample_violin.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =percent.mt)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 10)+
  ggtitle("Percent mt")
print(x)
dev.off()

# Visualize number of percent.mt per sample Density plot----
png(paste0("./03_plots/1_QC/QC_1_QC6_Mito_per_sample_density.png"))
x <- metadata_QC6 %>%
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 10)
print(x)
dev.off()

# Visualize number of percent.rb per sample Violin plot----
png(paste0("./03_plots/1_QC/QC_1_QC6_Rb_per_sample_violin.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=sample, fill=sample, y =percent.rb)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.1) +
  geom_violin(alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Percent rb")
print(x)
dev.off()

# Visualize number of Complexity of cells/sample as nGene vs nUMI----
png(paste0("./03_plots/1_QC/QC_1_QC6_Gene_vs_UMI_persample.png"))
x <- metadata_QC6%>%
  ggplot(aes(x=nUMI, y=nGene)) +
  geom_point(aes(color=percent.mt)) +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 40000) +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 6000) +
  facet_wrap(~sample)
print(x)
dev.off()

# Visualize number of Complexity of cells/sample as log10Genes/UMI----
png(paste0("./03_plots/1_QC/QC_1_QC6_Gene perUMI_density.png"))
x <- metadata_QC6 %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
print(x)
dev.off()


#################### Normalize Data and Find Variable Features Data with QC6 ############################################################

#normalize data set to account for sequencing depth, default scale to 10 000 and log2-transform
NPC_87 <-subset(NPC_87, subset = QC6 == 'Pass')%>%
  NormalizeData(verbose=F)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F)
saveRDS(NPC_87, "./01_tidy_data/1_QC_1_QC6_NPC_87.rds")

NPC_88 <-subset(NPC_88, subset = QC6 == 'Pass')%>%NormalizeData(verbose=F)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F)
saveRDS(NPC_88, "./01_tidy_data/1_QC_1_QC6_NPC_88.rds")

NPC_91 <-subset(NPC_91, subset = QC6 == 'Pass')%>%NormalizeData(verbose=F)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F)
saveRDS(NPC_91, "./01_tidy_data/1_QC_1_QC6_NPC_91.rds")

NPC_92 <-subset(NPC_92, subset = QC6 == 'Pass')%>%NormalizeData(verbose=F)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F)
saveRDS(NPC_92, "./01_tidy_data/1_QC_1_QC6_NPC_92.rds")

#Remove unused data from memory to save ram
rm( NPC_87,NPC_88,NPC_91,NPC_92, metadata_QC6, metadata)
