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
#renv::install("ggpubr")
#renv::install("tidyverse")
#renv::install("hrbrthemes")
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

#### Correlate Soup Genes per Sample with Clusters per smaples and check wich cluster contriubtes to soup most #####
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(stringr)
library(ggpubr)
library(hrbrthemes)
library(tidyverse)
library(rlang)
library(reshape2)
#BiocManager::install("GeneExpressionSignature")

set.seed(42)


animals <-c("87","88","91","92")
# Import SCtransformed, integrated, clustered, annotated,  dataset from FK44_1_COUNT_2_TRANSFORM and CLUSTERING_Script.R
NPC_CLUSTER <- readRDS("./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced_woMALAT_Filter.rds")
# Do Tutorial from https://divingintogeneticsandgenomics.com/post/how-to-do-gene-correlation-for-single-cell-rnaseq-data-part-1/
#Functions they define in the tutorial
matrix_to_expression_df<- function(x, obj){
  df<- x %>%
    as.matrix() %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var= "gene") %>%
    tidyr::pivot_longer(cols = -1, names_to = "cell", values_to = "expression") %>%
    tidyr::pivot_wider(names_from = "gene", values_from = expression) %>%
    left_join(obj@meta.data %>% 
                tibble::rownames_to_column(var = "cell"))
  return(df)
}
get_expression_data<- function(obj, assay = "RNA", slot = "data", genes = NULL, cells = NULL){
  if (is.null(genes) & !is.null(cells)){
    df<- GetAssayData(obj, assay = assay, slot = slot)[, cells, drop = FALSE] %>%
      matrix_to_expression_df(obj = obj)
  } else if (!is.null(genes) & is.null(cells)){
    df <- GetAssayData(obj, assay = assay, slot = slot)[genes, , drop = FALSE] %>%
      matrix_to_expression_df(obj = obj)
  } else if (is.null(genes & is.null(cells))){
    df <- GetAssayData(obj, assay = assay, slot = slot)[, , drop = FALSE] %>%
      matrix_to_expression_df(obj = obj)
  } else {
    df<- GetAssayData(obj, assay = assay, slot = slot)[genes, cells, drop = FALSE] %>%
      matrix_to_expression_df(obj = obj)
  }
  return(df)
}
Idents(NPC_CLUSTER)<- NPC_CLUSTER$mouseRNA.main
DimPlot(NPC_CLUSTER, reduction = "umap", label = TRUE)
for(a in animals){
# Import Soup Profiles generated in Script FK44_1_COUNT_1_QC_Script.R by SoupX
  Soup <- read_csv(paste0("99_other/0_Decont_SoupX/0_Decont_SoupX_Soup_Genes_iAL",a,".csv"))
# Import Seurat Object generated in FK44_1_COUNT_2_TRANSFORM and CLUSTERING.R 
  dataset<-subset(NPC_CLUSTER, sample== paste0("iAL",a))
  umapplot<-DimPlot(dataset, reduction = "umap", label = TRUE)
  png(paste0("./03_plots/2_1_Soup_Correlation_with_Clusters/Soup_Correlation_UMAP_iAL",a,".png"))
  print(umapplot)
  dev.off()
#Get the Genes from Soup file  
  Genes <-as.data.frame(Soup)%>%filter(counts>10) #Counts für Gene aus der Soup >10
  row.names(Genes)<-Genes[,1]
  vGenes <-as.vector(unlist(Genes[,1]))
#Get Mean expression per gene per cluster inSeurat object 
  expression <-get_expression_data(dataset, genes = vGenes)
  Mean_SoupGenes_Clusters<-expression%>%group_by(mouseRNA.main)%>%summarize_at(vGenes,mean)
# Transpose Clusters Dataframe and merge Cleaned cluster df with soup df
  Mean_SoupGenes_Clusters <- as.data.frame(Mean_SoupGenes_Clusters)
  Clusters_SoupGenes <- as.data.frame(t(Mean_SoupGenes_Clusters))
  names(Clusters_SoupGenes) <-lapply(Clusters_SoupGenes[1,], as.character) # set first row as column names
  Clusters_SoupGenes <-Clusters_SoupGenes[-1,] # delete first row
  Genes <- as.data.frame(Genes)
#merge gene expression from soup and cleaned up dataset/clusters
  CellySoup <-merge(Clusters_SoupGenes,Genes, by = "row.names", all = T)
  CellySoup <-data.frame(lapply(CellySoup,as.numeric))
#Make Scatter Plots with ggplot and calculate correlation of expression in soup vs clusters
  plot_list<-list()
  
## Now the Plots of different clusters vs Soup follow. was not able to put clusters in a loop to loop over cluster vector  
  Hepatocytes<-ggplot(CellySoup,aes(x=Hepatocytes, y= log2(counts)))+
      geom_point(size = 0.8) +
      geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
      theme_ipsum()+
      ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
      ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
      annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
      annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
    plot_list$Hepatocytes<-Hepatocytes
  
  Tcells<-ggplot(CellySoup,aes(x=T.cells, y= log2(counts)))+
    geom_point(size = 0.8) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_ipsum()+
  ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
    ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
    annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
    annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
  plot_list$Tcells<-Tcells
  
  Bcells<-ggplot(CellySoup,aes(x=B.cells, y= log2(counts)))+
    geom_point(size = 0.8) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_ipsum()+
  ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
    ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
    annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
    annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
  plot_list$Bcells<-Bcells
  
  EndothelialCells<-ggplot(CellySoup,aes(x=Endothelial.cells, y= log2(counts)))+
    geom_point(size = 0.8) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_ipsum()+
  ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
    ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
    annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
    annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
  plot_list$EndothelialCells<-EndothelialCells
  
  Fibroblasts<-ggplot(CellySoup,aes(x=Fibroblasts, y= log2(counts)))+
    geom_point(size = 0.8) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_ipsum()+
  ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
    ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
    annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
    annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
  plot_list$Fibroblasts<-Fibroblasts
  
  Macrophages<-ggplot(CellySoup,aes(x=Macrophages, y= log2(counts)))+
    geom_point(size = 0.8) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_ipsum()+
  ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
    ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
    annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
    annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
  plot_list$Macrophages<-Macrophages
  
  NKcells<-ggplot(CellySoup,aes(x=NK.cells, y= log2(counts)))+
    geom_point(size = 0.8) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_ipsum()+
  ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
    ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
    annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
    annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
  plot_list$NKcells<-NKcells
  
  Monocytes<-ggplot(CellySoup,aes(x=Monocytes, y= log2(counts)))+
    geom_point(size = 0.8) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_ipsum()+
  ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
    ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
    annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
    annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
  plot_list$Monocytes<-Monocytes
  
  Granulocytes<-ggplot(CellySoup,aes(x=Granulocytes, y= log2(counts)))+
    geom_point(size = 0.8) +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    theme_ipsum()+
  ggpubr::stat_cor(method = "pearson",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=1)+
    ggpubr::stat_cor(method = "spearman",label.x = -Inf, label.y=Inf,hjust=-0.7,vjust=2.5)+
    annotate(geom="text",-Inf,Inf, label="Pearson",color="black", hjust=0,vjust=1)+
    annotate(geom="text",-Inf,Inf, label="Spearman",color="black", hjust=0,vjust=3)
  plot_list$Granulocytes<-Granulocytes
  
#Print and save plots  
  for(i in 1:length(names(plot_list))){
    cluster<-names(plot_list[i])
    png(paste0("./03_plots/2_1_Soup_Correlation_with_Clusters/Soup_Correlation_",cluster,"_iAL",a,".png"))
    print(plot_list[[i]])
    dev.off()
    }
  }

