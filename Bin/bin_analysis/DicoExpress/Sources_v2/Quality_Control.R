###### Filtering and Normalization ######
######### Data Quality Control ##########

Quality_Control <- function(Data_Directory, Results_Directory, Project_Name, 
                            Target, Raw_Counts, Color_Group=NULL, Min.count=15,
                            Normalization_Method="TMM",marginal_sample, Subset=NULL, ...){
    ## Subdirectory
    if (!I(Project_Name %in% dir(Results_Directory))) dir.create(paste0(Results_Directory,"/",Project_Name), showWarnings=FALSE)

    subdirectory <- "Quality_Control"
    path=paste0(Results_Directory,"/",Project_Name,"/",subdirectory)
    
    if (!I(subdirectory %in% dir(paste0(Results_Directory,"/",Project_Name)))) dir.create(path, showWarnings=FALSE)
    
    ## Nb of biological factors, genes and samples 
    NbBioFactors <- ncol(Target)-1
    NbGenes <- nrow(Raw_Counts)
    NbSamples <- ncol(Raw_Counts)
    

    ## Name of the conditions and Nb od replicate per conidtion
    if(NbBioFactors >= 2) ## Modification Marie de NbBioFactors == 2 en >=2
    {
        BioFactors <- unique(apply(Target[,1:NbBioFactors],1,function(x) paste(x,collapse="_")))
        NbReplicate <- table(apply(Target[,1:NbBioFactors],1,function(x) paste(x,collapse="_")))[BioFactors]
    }
    else
    {
        BioFactors <- unique(Target[,1])
        NbReplicate <- table(Target[,1])
    }
    ## Nb of conditions 
    NbConditions = length(BioFactors)
    if(length(Color_Group)!=NbConditions)
      {
        Color_Group=NULL
        cat("The color number is diffrent from the condition number. DiCoExpress redefines a new vector \n")   
      }
    ## Vector of colors for biological factors
    if(is.null(Color_Group))
    {
        if (NbConditions<=12 && NbConditions>=3)
            Color_Group <- brewer.pal(NbConditions, name = ifelse(nlevels(Target[,1])==2,'Paired','Set3'))
        else
            Color_Group <-colors()[seq(from=3,length=NbConditions,by=4)]
    }
    
    ## Filtering
    
    ## Creation group vector
    group = c()
    for (i in 1:length(NbReplicate)){
      group = c(group,rep(0+i,as.vector(NbReplicate)[i]))
    }
    Counts <- DGEList(counts = Raw_Counts, genes = rownames(Raw_Counts),group = group)

    ## Filtering strategy
    ## filterByExpr function of edgeR package
    if(Min.count >= 0)
    {
      keep <- filterByExpr(Counts, min.count=Min.count,... )
    }
    # Error message if Min.count < 0
    else{
      cat("Min.count can't be negative\n")
    }

    ## Filtering
    fcounts <- Counts[keep,]
    ## File.txt for low counts genes
    id.filteredGenes <- rownames(Counts)[!keep]
    NbFilteredGenes <- length(id.filteredGenes)
    
    
    CPM_Cutoff <- Min.count/median(fcounts$samples$lib.size)*1e6

    
    fileout=paste0(path,"/","Low_count_genes.txt")
    write.table(id.filteredGenes, file=fileout, sep="\t",quote=F,row.names=F,col.names=F)
  
    ## Normalization
    fcounts$samples$lib.size=colSums(fcounts$counts)
    normCounts <- calcNormFactors(fcounts, method=Normalization_Method)
    tmm <- normCounts$samples$norm.factors
    N <- normCounts$samples$lib.size
    f <- tmm * N/mean(tmm * N) 
    ## normalized read counts are obtained by dividing raw read counts by these re-scaled normalization factors
    NormCounts <- scale(normCounts$counts, center=FALSE, scale=f)
    
    write.csv2(NormCounts,paste0(path,"/","NormCounts.csv"),quote=F)
    ## Log2 + 1 transformation of NormCounts
    NormCounts_log2 <- log2(NormCounts+1)
    
    ## Write result file
    Output <- file(paste0(path,"/","Normalization_Results.txt"), open="wt")
    sink(Output)
    sink(Output, type = "message")
    
    cat("################################################\nFiltering and Normalization\n################################################\n\n")
    cat("#### Description of Raw counts table ####\n")
    cat("Number of samples:",NbSamples,"\n")
    cat("Number of genes:",NbGenes,"\n\n")
    cat("#### Filtering ####\n")
    cat("Number of genes discarded by the filtering:", NbFilteredGenes,"\n")
    cat("Number of genes analyzed after filtering:",nrow(fcounts$counts),"\n")
    cat("CPM_Cutoff:",CPM_Cutoff,"\n\n")
    
    cat("################################################\n Statistics on the normalization factors \n################################################\n\n")
    print(summary(normCounts$samples$norm.factors))
    cat("\n")
    sink(type = "message")
    sink()
    
    cat("################################################\nFiltering\n################################################\n\n")
    cat("#### Description of Raw counts table ####\n")
    cat("Number of samples:",NbSamples,"\n")
    cat("Number of genes:",NbGenes,"\n\n")
    cat("Number of genes discarded by the filtering:", NbFilteredGenes,"\n")
    cat("Number of genes analyzed after filtering:",nrow(fcounts$counts),"\n\n")
    cat("################################################\n Statistics on the normalization factors \n################################################\n\n")
    print(summary(normCounts$samples$norm.factors))
    cat("\n")
    
#### Data Quality Control
    ## Pdf file where plots are stored
    #pdf_args_names <- names(formals(pdf))
    #args <- list(...)
    #save_args <- args[names(args) %in% pdf_args_names]
    #args$min.total.count <- NULL
    pdf(file=paste0(path,"/","Data_Quality_Control.pdf"))
    
    pdf_args_names <- c("width", "height", "onefile", "family", "title", "fonts", "paper",
      "encoding", "pointsize", "bg", "fg", "pagecentre", "useDingbats","colormodel", "fillOddEven", "compress")
    i = 0
    arg <- c()
    for (e in names(list(...))){
      i = i+1
      if (e %in% pdf_args_names){
        arg <- c(arg,e,...elt(i))
      }
    }
    #for (i in 1:length(arg)){}
    
    colorBioFactors <- rep(Color_Group,NbReplicate) 
    
        
#### Barplot ####
    ## Barplot of raw counts
    df_Total<-as.data.frame(as.table(colSums(Raw_Counts)))
    names(df_Total)[1:2] <- c("Samples","Library_Size")
    title <- "Total raw counts per sample"
    ylab <- "Library size (millions)"
    xlab <- "Samples"
    p <- ggplot(data=df_Total,aes(x=Samples , y=Library_Size/1e+06, fill=Samples))
    p <- p + geom_bar(stat="identity")
    p <- p + labs(title = title, x=xlab, y = ylab)
    p <- p + scale_fill_manual(values = colorBioFactors)
    p <- p + theme(title = element_text(size=15),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10),axis.text.y=element_text(size=10))
    p <- p + theme(legend.position='none')
    print(p)
    
    ## Barplot of normalized count library sizes (millions)
    df_Total <- as.data.frame(as.table(colSums(NormCounts)))
    names(df_Total)[1:2] <- c("Samples","Library_Size")
    title <- "Total normalized counts per sample"
    ylab <- "Library size (millions)"
    xlab <- "Samples"
    p <- ggplot(data=df_Total,aes(x=Samples , y=Library_Size/1e+06, fill=Samples))
    p <- p + geom_bar(stat="identity")
    p <- p + labs(title = title, x=xlab, y = ylab)
    p <- p + scale_fill_manual(values = colorBioFactors)
    p <- p + theme(title = element_text(size=15),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10),axis.text.y=element_text(size=10))
    p <- p + theme(legend.position='none')
    print(p)

#### Boxplot ####
    ## Boxplot of filtered counts
    log2fcounts <- as.data.frame(log2(fcounts$counts+1))
    df_fcount = reshape2::melt(log2fcounts, variable.name = "Samples", value.name = "Counts")
    title <- "Boxplot of raw counts"
    ylab <- "Gene expression (log2)"
    xlab <- "Samples"
    p <- ggplot(df_fcount, aes(x = Samples, y = Counts, fill = Samples)) 
    p <- p + geom_boxplot()
    p <- p + labs(title = title, x=xlab, y = ylab)
    p <- p + scale_fill_manual(values = colorBioFactors)
    p <- p + theme(title = element_text(size=15),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10),axis.text.y=element_text(size=10))
    p <- p + theme(legend.position='none')
    print(p)

    ## Boxplot of normalized counts
    title <- "Boxplot of normalized counts"
    log2norm <- as.data.frame(NormCounts_log2)
    df_norm = reshape2::melt(log2norm, variable.name = "Samples", value.name = "Counts")
    p <- ggplot(df_norm, aes(x = Samples, y = Counts, fill = Samples))
    p <- p + geom_boxplot()
    p <- p + labs(title = title, x=xlab, y = ylab)
    p <- p + scale_fill_manual(values = colorBioFactors)
    p <- p + theme(title = element_text(size=15),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10),axis.text.y=element_text(size=10))
    p <- p + theme(legend.position='none')
    print(p)

    
#### Clustering of the sample-to-sample distance ####
    ## Raw counts
    mat.dist <- as.matrix(dist(t(Raw_Counts)))
    par(cex.main=.9)
    heatmap.2(mat.dist,
              hclustfun = function(x) hclust(x, method = "ward.D"),
              margins = c(9,9),
              trace="none",
              cexRow = 0.25,  # Adjust the value as needed to make the labels fit
              cexCol = 0.25)
    title(main = "Heatmap of raw counts [Euclidean dist, Ward link]",
          cex.main = 1,
          adj = 1,
          line = 0)
    
    ## Normalized counts
    mat.dist <- as.matrix(dist(t(NormCounts)))
    par(cex.main=.9)
    output_file <- paste0(path, "/", "mat_dist_output.txt")
    
    # Écrire la matrice de distance dans un fichier texte
    write.csv2(mat.dist, file = output_file, sep = "\t", row.names = FALSE)
    library(gplots)
    
    # Assuming 'mat.dist' is already defined as your distance matrix
    heatmap.2(mat.dist, 
              hclustfun = function(x) hclust(x, method = "ward.D"), 
              margins = c(9, 9),
              trace = "none",
              cexRow = 0.25,  # Adjust the value as needed to make the labels fit
              cexCol = 0.25)  
    
    title(main = "Heatmap of normalized counts [Euclidean dist, Ward link]", 
          cex.main = 1,  
          adj = 1, 
          line = 0)
    
    #get the sum,mean of distance  for each sample and global mean
    sum_distances <- rowSums(mat.dist)
    mean_distances <- rowMeans(mat.dist)
    global_mean_distance <- mean(mat.dist)
    
    #get 5 more distal sample
    max_sum_indices <- order(-sum_distances)[1:marginal_sample]
    
    #extract 
    max_sum_samples <- rownames(mat.dist)[max_sum_indices]
    max_sum_values <- sum_distances[max_sum_indices]
    max_mean_values <- mean_distances[max_sum_indices]

    
    # Créer un dataframe contenant les résultats
    marginal_data <- data.frame(
      Sample = max_sum_samples,
      Sum_Distance = max_sum_values,
      Mean_Distance = max_mean_values)
    
    marginal_data <- marginal_data %>% 
      mutate("deviation_from_avg" = Mean_Distance-global_mean_distance )
    
    # save data in txt
    output_file <- paste0(path, "/", "marginal_samples.txt")
    write.table(marginal_data, file = output_file, sep = "\t", row.names = FALSE)
    
#### PCA ####
    qualitativeFactors <- lapply(Target, factor)
    cnames <- names(qualitativeFactors)
    nb.batch = length(cnames)
    Biological_Factors <- rep(BioFactors,NbReplicate)
    ## raw counts
    raw.PCA <- cbind.data.frame(qualitativeFactors,sqrt(t(Raw_Counts)))
    res.raw = PCA(raw.PCA, quali.sup=1:nb.batch, ncp=2, graph=FALSE)
    PC1.raw <- res.raw$ind$coord[,1]
    PC2.raw <- res.raw$ind$coord[,2]
    labs.raw <- rownames(res.raw$ind$coord)
    PCs.raw <- data.frame(cbind(PC1.raw,PC2.raw))
    rownames(PCs.raw) <- labs.raw
    PCs.raw <- cbind(PCs.raw,qualitativeFactors[],Biological_Factors)
    PCs.raw$Biological_Factors <- factor(PCs.raw$Biological_Factors, levels = BioFactors)
    ## Normalized counts
    norm.PCA <- cbind.data.frame(qualitativeFactors,sqrt(t(NormCounts)))
    res.norm = PCA(norm.PCA, quali.sup=1:nb.batch, ncp=2, graph=F)
    PC1.norm <- res.norm$ind$coord[,1]
    PC2.norm <- res.norm$ind$coord[,2]
    labs.norm <- rownames(res.norm$ind$coord)
    PCs.norm <- data.frame(cbind(PC1.norm,PC2.norm))
    rownames(PCs.norm) <- labs.norm
    PCs.norm <- cbind(PCs.norm,qualitativeFactors[],Biological_Factors)
    PCs.norm$Biological_Factors <- factor(PCs.norm$Biological_Factors, levels = BioFactors)

    for(ii in 3:ncol(PCs.raw))
    {
        name.quali <- colnames(PCs.raw)[ii]
        PC1 <- PCs.raw$PC1
        PC2 <- PCs.raw$PC2
        p <- ggplot(data=PCs.raw,aes(x=PC1 , y=PC2))
        # p <- p + geom_point (data=PCs.raw, aes_string (x="PC1", y="PC2", color=name.quali))
        p <- p + geom_point(aes(color = !!sym(name.quali)))
        
        p <- p + geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_vline(aes(xintercept = 0), linetype = "dashed")
        concat = cbind.data.frame(PCs.raw[,ii],PC1,PC2)
        ellipse.coord = coord.ellipse(concat,bary=T)
        ellipse.coord <-as.data.frame(ellipse.coord)
        setDT(ellipse.coord)
        setnames(ellipse.coord,c("res.PCs.raw...ii.","res.PC1","res.PC2"), c("group", "x", "y"))
        p <- p + geom_path(data=ellipse.coord, aes(x=x, y=y, color = group),show.legend = FALSE) 
        p <- p + theme_bw()
        if (name.quali=="Biological_Factors")
        {
            p <- p + scale_color_manual(values=Color_Group)
        }
        p <- p + labs(title = "PCA on raw counts", x=paste((format(round((res.raw$eig[1,2]), 2), nsmall = 2)),"%"), y = paste((format(round((res.raw$eig[2,2]), 2), nsmall = 2)),"%"))
        p <- p + guides(colour = guide_legend(override.aes = list(size= 3)))
        p <- p + theme(legend.text = element_text(size = 15))
        p <- p + theme(title = element_text(size=15),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20), axis.title = element_text(family = "Helvetica", size=18))
        print(p)
    }
    
    for(ii in 3:ncol(PCs.norm))
    {
        name.quali <- colnames(PCs.norm)[ii]
        PC1 <- PCs.norm$PC1
        PC2 <- PCs.norm$PC2
        p <- ggplot(data=PCs.norm,aes(x=PC1 , y=PC2))
        # p <- p + geom_point (data=PCs.norm, aes_string (x="PC1", y="PC2", color=name.quali))
        p <- p + geom_point(aes(color = !!sym(name.quali)))
        
        p <- p + geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_vline(aes(xintercept = 0), linetype = "dashed")
        concat = cbind.data.frame(PCs.norm[,ii],PC1,PC2)
        ellipse.coord = coord.ellipse(concat,bary=T)
        ellipse.coord <-as.data.frame(ellipse.coord)
        setDT(ellipse.coord)
        setnames(ellipse.coord,c("res.PCs.norm...ii.","res.PC1","res.PC2"), c("group", "x", "y"))
        p <- p + geom_path(data=ellipse.coord, aes(x=x, y=y, color = group),show.legend = FALSE) 
        p <- p + theme_bw()
        if (name.quali=="Biological_Factors")
        {
            p <- p + scale_color_manual(values=Color_Group)
        }
        p <- p + labs(title = "PCA on normalized counts", x=paste((format(round((res.norm$eig[1,2]), 2), nsmall = 2)),"%"), y = paste((format(round((res.norm$eig[2,2]), 2), nsmall = 2)),"%"))
        p <- p + guides(colour = guide_legend(override.aes = list(size= 3)))
        p <- p + theme(legend.text = element_text(size = 15))
        p <- p + theme(title = element_text(size=15),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20), axis.title = element_text(family = "Helvetica", size=18))
        print(p)
    }
    
    ## subset of normalized counts
    if (!is.null(Subset))
    {
      NormCounts.f=NormCounts[is.element(rownames(NormCounts),Subset),]
      norm.PCA <- cbind.data.frame(qualitativeFactors,sqrt(t(NormCounts.f)))
      res.norm = PCA(norm.PCA, quali.sup=1:nb.batch, ncp=2, graph=F)
      PC1.norm <- res.norm$ind$coord[,1]
      PC2.norm <- res.norm$ind$coord[,2]
      labs.norm <- rownames(res.norm$ind$coord)
      PCs.norm <- data.frame(cbind(PC1.norm,PC2.norm))
      rownames(PCs.norm) <- labs.norm
      PCs.norm <- cbind(PCs.norm,qualitativeFactors[],Biological_Factors)
      PCs.norm$Biological_Factors <- factor(PCs.norm$Biological_Factors, levels = BioFactors)
      for(ii in 3:ncol(PCs.norm))
      {
        name.quali <- colnames(PCs.norm)[ii]
        PC1 <- PCs.norm$PC1
        PC2 <- PCs.norm$PC2
        p <- ggplot(data=PCs.norm,aes(x=PC1 , y=PC2))
        # p <- p + geom_point (data=PCs.norm, aes_string (x="PC1", y="PC2", color=name.quali))
        p <- p + geom_point(aes(color = !!sym(name.quali)))
        
        p <- p + geom_hline(aes(yintercept = 0), linetype = "dashed") + geom_vline(aes(xintercept = 0), linetype = "dashed")
        concat = cbind.data.frame(PCs.norm[,ii],PC1,PC2)
        ellipse.coord = coord.ellipse(concat,bary=T)
        ellipse.coord <-as.data.frame(ellipse.coord)
        setDT(ellipse.coord)
        setnames(ellipse.coord,c("res.PCs.norm...ii.","res.PC1","res.PC2"), c("group", "x", "y"))
        p <- p + geom_path(data=ellipse.coord, aes(x=x, y=y, color = group),show.legend = FALSE) 
        p <- p + theme_bw()
        if (name.quali=="Biological_Factors")
        {
          p <- p + scale_color_manual(values=Color_Group)
        }
        p <- p + labs(title = "PCA on some normalized counts", x=paste((format(round((res.norm$eig[1,2]), 2), nsmall = 2)),"%"), y = paste((format(round((res.norm$eig[2,2]), 2), nsmall = 2)),"%"))
        p <- p + guides(colour = guide_legend(override.aes = list(size= 3)))
        p <- p + theme(legend.text = element_text(size = 15))
        p <- p + theme(title = element_text(size=15),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20), axis.title = element_text(family = "Helvetica", size=18))
        print(p)
      } 
    }
    dev.off()
}
