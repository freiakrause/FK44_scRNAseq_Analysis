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
#BiocManager::install('fgsea')

#remotes::install_github("Moonerss/scrubletR")
#remotes::install_github('immunogenomics/presto')
#remotes::install_github("ianmoran11/mmtable2")

################################# Load Libraries #########################################################
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
library(Seurat)
library(fgsea)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(readr)

set.seed(42)
#everything from
#https://biostatsquid.com/fgsea-tutorial-gsea/
in_path <- "./01_tidy_data/"
GS_path<-"./99_other/5_GeneSets_for_GSEA"
out_path <- "./03_plots/3_GSEA_MainCluster/"
#### Load Input Data - latest NPC dataset
NPC_ALL_TRANSFORMED <- readRDS( "./01_tidy_data/3_NPC_ALL_TRANSFORMED_Annotated_Reduced_woMALAT_Filter_PrepSCT.rds")
y<-NPC_ALL_TRANSFORMED@meta.data%>%group_by(mouseRNA.main,stim)%>%summarise(n=n())
# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
## Function: prepare_gmt --------------------------------------
library(fgsea)
prepare_gmt <- function(gmt_file, my_genes, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(my_genes, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}
# Analysis ====================================================
## 1. Read in data -----------------------------------------------------------
clusters<-unique(NPC_ALL_TRANSFORMED@meta.data$mouseRNA.main)

DEG_Genes<-sapply(X=clusters,
       FUN= function(c){c <- read.csv(paste0("./99_other/3_DEG_Analysis_MainCluster/1_DEG_Analysis_single_MAST",c,".csv"), row.names = 1)},
       simplify = F)
  

## 2. Prepare background genes -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files 
my_genes <- rownames(NPC_ALL_TRANSFORMED@assays$RNA$counts)
gmt_files <- list.files(path = GS_path, pattern = '.gmt', full.names = TRUE)
gmt_files

  #### Prepare Ranked Gene list
  for (c in clusters){
    for(gmt in gmt_files){
    bg_genes <- prepare_gmt(gmt,my_genes)
    a<-substring(paste0(gmt[[1]]),32,33)
    DEG1<-DEG_Genes[[c]]
    rankings <- sign(DEG1$avg_log2FC)*(-log10(DEG1$p_val)) # we will use the signed p values from spatial DGE as ranking
    names(rankings) <- rownames(DEG1) # genes as names#
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
  # Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
    max_ranking <- max(rankings[is.finite(rankings)])
    min_ranking <- min(rankings[is.finite(rankings)])
    rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
    rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
    ###GSEA needs named vector
    #### Run GSEA
    GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                     stats = rankings,
                     eps=0.0,
                     scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                     minSize = 10,
                     maxSize = 500,
                     nproc = 1) # for parallelisation
    
    sum(GSEAres[, padj < 0.01])
    sum(GSEAres[, pval < 0.01])
    saveRDS(GSEAres, file = paste0("./01_tidy_data/3_GSEA_results",c,"_",a,".RDS"))
    
    # Select only independent pathways, removing redundancies/similar pathways
    collapsedPathways <- collapsePathways(GSEAres[order(pval)][pval < 0.05], bg_genes, rankings)
    mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
 
    t<-plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)+
    theme_classic()+
    theme(title=element_text())+
    labs(title=paste0("Top 5 Up and Down Pathways in ",c," against List ",a),
         subtitle = paste0("Number of Genesets enriched with p_val<0.01:",sum(GSEAres[, pval < 0.01])))
    ggsave(filename =paste0("./03_plots/3_GSEA_MainCluster/GSEA_Top30pathways_",c,"_List_",a,"_selectedPathways.png"),
         t, width = 23, height = 15, dpi = 150, bg= "transparent")
  
    #### Plot top pathways
    topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = sum(GSEAres[, pval < 0.01])), pathway]
    topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = sum(GSEAres[, pval < 0.01])), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    ### Plot Top Enriched Pathways
    t<-plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)+
      theme_classic()+
      theme(title=element_text())+
      labs(title=paste0("Top 5 Up and Down Pathways in ",c," against List ",a),
           subtitle = paste0("Number of Genesets enriched with p_val<0.01:",sum(GSEAres[, pval < 0.01])))
    ggsave(filename =paste0("./03_plots/3_GSEA_MainCluster/GSEA_Top30pathways_",c,"_List_",a,".png"),
           t, width = 23, height = 15, dpi = 150, bg= "transparent")
    #### Plot Nerichement scores for Top pathway
      p <- plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
                   rankings) + 
        labs(title = head(GSEAres[order(padj), ], 1)$pathway) + 
        theme_classic() +
        theme(text = element_text(size=5),
          title = element_text(size=5),
              axis.title = element_text(size =5))+
        scale_x_continuous('Rank') +
        scale_y_continuous(paste0('Enrichment score (ES) in ',c))
        
    ggsave(filename = paste0("./03_plots/3_GSEA_MainCluster/GSEA_Enrichment_",c,"_List_",a,".png"),
           p,width = 6, height = 3,dpi = 300, bg="transparent")
    }
    ### Other viszualizaions were easier for me with dataframe
    rankings<-as.data.frame(rankings)
    rankings$position<-c(1:length(rankings$rankings))
    rankings$gene<-(c(rownames(rankings)))
    p<-ggplot(rankings,aes(x=position,y=rankings))+
      geom_point()
    print(p)
    ggsave(filename = paste0("./03_plots/3_GSEA_MainCluster/Ranking_",c,"_Raw.png"),
           p,width = 6, height = 3,dpi = 300, bg="transparent")
    p<-ggplot(rankings[1:20,], aes(x=gene, rankings)) + 
    geom_point() +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(p)
    ggsave(filename = paste0("./03_plots/3_GSEA_MainCluster/Ranking_",c,"_Top20.png"),
           p,width = 6, height = 3,dpi = 300, bg="transparent")
    p<-ggplot(rankings[(nrow(rankings)-19):nrow(rankings),], aes(x=gene, rankings)) + 
      geom_point() +
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(p)
    ggsave(filename = paste0("./03_plots/3_GSEA_MainCluster/Ranking_",c,"_TopMinus20.png"),
           p,width = 6, height = 3,dpi = 300, bg="transparent")
}

