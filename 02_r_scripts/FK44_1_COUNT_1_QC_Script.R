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
library(ggplot2)
library(Seurat)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(scrubletR)
library(gprofiler2)
library(SoupX)
library(celda)
library(scater)
library(stringr)
library(Matrix)
library(gridExtra)
library(qlcMatrix)
library(pheatmap)
source("02_r_scripts/malat1_function.R")
set.seed(42)
# output_dir <- file.path("./", "99_other/0_Decont_SoupX")
# 
# if (!dir.exists(output_dir)){
#   dir.create(output_dir)
# } else {
#   print("Dir already exists!")
# }

# ###########################Load Data and Decontaminate with SoupX ##################
#  #### Soup Decont ----
Hep_genes <-c("Saa1", "Saa2","Alb","Tat")
T_genes <-c("Cd3e","Cd3d", "Cd4", "Cd8a") 
B_genes <-c("Cd19")
M_genes <-c("Clec4f")
Gene_List <-list(Hep_genes,T_genes,B_genes,M_genes)
animals <-c("87","88","91","92")
for (i in animals){  
  sc = load10X(paste0("./00_raw_data/biomedical-sequencing.at/projects/BSA_0873_FK44_1_LiverMet_A_1_1_51ddbfd228ec40b096e110101b219cb0/COUNT/Liver_NPC_iAL",i,"_transcriptome"))
  sc = setClusters(sc, sc$metaData$clustersFine)
  Soup <-sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ]
  write.csv(Soup,paste0("./99_other/0_Decont_SoupX/0_Decont_SoupX_Soup_Genes_iAL",i,".csv"))
  Soup_Genes <-head(rownames(Soup), n=10)
  #  sc = autoEstCont(sc)
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = Hep_genes, T_genes,B_genes,M_genes))
  sc = calculateContaminationFraction(sc, list(IG = Hep_genes,T_genes,B_genes,M_genes), useToEst = useToEst,forceAccept=TRUE)
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


 
########################## Kategorien Stimlulation und Sex hinzufügen, evtl noch age? ####################
animals <-c("_")#"87","88","91","92")

NPC_87 <-readRDS(paste0("./01_tidy_data/0_iAL87_SoupX.rds")) 
NPC_88 <- readRDS(paste0("./01_tidy_data/0_iAL88_SoupX.rds")) 
NPC_91 <- readRDS(paste0("./01_tidy_data/0_iAL91_SoupX.rds")) 
NPC_92 <- readRDS(paste0("./01_tidy_data/0_iAL92_SoupX.rds")) 
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
#Setting things up for quality control mitochondrial genes, ribosomal content, doublets----
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
  scale_x_log10() +
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
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) +
  geom_point() +
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
  scale_x_log10() +
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
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) +
  geom_point() +
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
rm(NPC_87.data,NPC_88.data, NPC_91.data, NPC_92.data, metadata_QC5, metadata_QC6, metadata)

########################## Define Anchors for Integration and Integrate Different Data Sets ####

#NPC_87 <-readRDS("./01_tidy_data/1_QC_1_QC6_NPC_87.rds")
#NPC_88 <-readRDS("./01_tidy_data/1_QC_1_QC6_NPC_88.rds")
#NPC_91 <-readRDS("./01_tidy_data/1_QC_1_QC6_NPC_91.rds")
#NPC_92 <-readRDS("./01_tidy_data/1_QC_1_QC6_NPC_92.rds")

########################################################################
########################################################################