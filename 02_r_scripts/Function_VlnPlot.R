###################### Function Create_Vplots ####################################
#Create multiple ViolinPlots with Seurat Object and vector of features , also give folder for saving as
#input to save the plots as png
Create_Vplots <- function(DataSet,feature_list,folder){
  for (i in feature_list){
    print(i)
    a <- VlnPlot(DataSet, features = i)
    png(filename = paste0("./03_plots/",folder,"/VlnPlot_",i,".png"))
    print(a)
    dev.off()
  }
}
