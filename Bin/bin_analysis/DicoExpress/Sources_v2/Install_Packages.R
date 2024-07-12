
## install CRAN packages
packages <- c("ggplot2","FactoMineR","gplots","reshape2","ggpubr","plyr","dplyr","RColorBrewer","data.table","UpSetR","pdftools")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

## install Bioconductor packages
BCpackages <- c("edgeR","coseq")
if (length(setdiff(BCpackages, rownames(installed.packages()))) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(setdiff(BCpackages, rownames(installed.packages())))
}

print(sapply(c(packages, BCpackages), require, character.only=T))

