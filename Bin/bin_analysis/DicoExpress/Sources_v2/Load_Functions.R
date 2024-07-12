Load_Functions<-function(source_path)
{
    library(edgeR)
    library(ggplot2)
    library(FactoMineR)
    library(gplots)
    library(reshape2)
    library(ggpubr)
    library(plyr)
    library(dplyr)
    library(RColorBrewer)
    library(coseq)
    library(data.table)
    library(UpSetR)
    library(pdftools)
    
  
    # source_path is the name of the folder containing the Source code.
    # If you want to use version 1, you must modify source_path in the Template_script.
    Sources_Directory=paste0("./",source_path)
    source(paste0(Sources_Directory,"/Load_Data_Files.R"))
    source(paste0(Sources_Directory,"/Quality_Control.R"))
    source(paste0(Sources_Directory,"/GLM_Contrasts.R"))
    source(paste0(Sources_Directory,"/DiffAnalysis_edgeR.R"))
    source(paste0(Sources_Directory,"/Contrast_Comparison.R"))
    source(paste0(Sources_Directory,"/Coexpression_coseq.R"))
    source(paste0(Sources_Directory,"/Enrichment.R"))
    source(paste0(Sources_Directory,"/Save_Parameters.R"))
    source(paste0(Sources_Directory,"/Salmon_to_CSV.R"))
    source(paste0(Sources_Directory,"/xlx_to_CSV.R"))
    
}
