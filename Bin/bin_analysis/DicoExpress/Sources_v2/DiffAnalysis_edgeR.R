###### Differential expression analysis using GLM edgeR ######

DiffAnalysis.edgeR <- function(Data_Directory, Results_Directory, Project_Name,
                               Target, Raw_Counts, GLM_Model, Contrasts,
                               Index_Contrast=(1:nrow(Contrasts)), 
                               Alpha_DiffAnalysis=0.05, NbGenes_Profiles=20, NbGenes_Clustering=50,
                               Min.count=15,Normalization_Method="TMM",... )
{
    
    ## subdirectory
    subdirectory <- "DiffAnalysis"
    path=paste0(Results_Directory,"/",Project_Name,"/",subdirectory)
    if (!I(subdirectory %in% dir(paste0(Results_Directory,"/",Project_Name))))
        dir.create(path, showWarnings=FALSE)
    
  ## Contrasts of interest
    if(nrow(Contrasts)==1){
        Contrasts.interest = Contrasts
    } else {Contrasts.interest <- Contrasts[Index_Contrast,]}
    ## Write Contrasts of interest in table file
    fileout <- paste0(path,"/","Contrasts_Interest_Matrix.txt")
    write.table(Contrasts.interest,fileout,
                row.names=FALSE,col.names=TRUE, sep="\t", quote = FALSE)
    
    ## Description of the dataset
    NbBioFactors <- ncol(Target)-1
    NbGenes <- nrow(Raw_Counts)
    NbSamples <- ncol(Raw_Counts)

    ## Name of the conditions and Nb od replicate per conidtion
    if(NbBioFactors == 2)
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

    ## Creation group vector
    group = c()
    for (i in 1:length(NbReplicate)){
      group = c(group,rep(0+i,as.vector(NbReplicate)[i]))
    }
    ## Building dge object
    dge <- DGEList(counts = Raw_Counts, genes = rownames(Raw_Counts),group = group)
    
    ## Filtering strategy
    ## filterByExpr function of edgeR package
    if(Min.count >= 0)
    {
      keep <- filterByExpr(dge, min.count=Min.count,...)
    }
    # Error message if Min.count < 0
    else{
      cat("Min.count can't be negative\n")
    }
    
    ## Pour savoir combien de gènes sont filtrés
    fcounts <- dge[keep,keep.lib.size=FALSE]
    id.filteredGenes <- rownames(dge)[!keep]
    NbFilteredGenes <- length(id.filteredGenes)
    fcounts$samples$lib.size=colSums(fcounts$counts)
    cat("Number of genes discarded by the filtering:", NbFilteredGenes,"\n")
    cat("Number of genes analyzed after filtering:",nrow(fcounts$counts),"\n\n")
    
    ## Filtering
    dge <- dge[keep,]
    dge$samples$lib.size=colSums(dge$counts)
    # pareil que dge <- dge[keep,keep.lib.size=FALSE]
    
    ## Normalization
    dge <- calcNormFactors(dge, method=Normalization_Method)
    tmm <- dge$samples$norm.factors 
    N <- dge$samples$lib.size
    f <- tmm * N/mean(tmm * N) 
    NormCounts <- scale(dge$counts, center=FALSE, scale=f)
    Gene_ID <- rownames(NormCounts)
    ## Log2 + 1 transformation of NormCounts
    NormCounts_log2 <- log2(NormCounts+1)
    NormCounts_log2 <- data.frame(Gene_ID,NormCounts_log2)    
    ## Log2 NormCounts Mean SD
    NormCounts_Log2_Mean <- matrix(NA,nrow(NormCounts_log2),NbConditions)
    NormCounts_Log2_SD <- matrix(NA,nrow(NormCounts_log2),NbConditions)
    cpt=1
    for (p in BioFactors)
    {
        NormCounts_tmp <- NormCounts[,grep(p,colnames(NormCounts))]
        NormCounts_Log2_tmp <- log2(NormCounts_tmp+1)
        NormCounts_Log2_Mean[,cpt] <- rowMeans(NormCounts_Log2_tmp,na.rm=TRUE)
        NormCounts_Log2_SD[,cpt] <- apply(NormCounts_Log2_tmp,1,sd,na.rm=TRUE)
        cpt=cpt+1
    }
    colnames(NormCounts_Log2_Mean) <-  paste0(BioFactors,"_Mean")
    colnames(NormCounts_Log2_SD) <- paste0(BioFactors, "_SD")
    NormCounts_log2_Mean_SD <- cbind(NormCounts_Log2_Mean,NormCounts_Log2_SD)
    NormCounts_log2_Mean_SD <- data.frame(Gene_ID,NormCounts_log2_Mean_SD)

    ## estimating dispersions
    rownames(GLM_Model) <- GLM_Model[,1]
    GLM_Model <- GLM_Model[,-1]
    dge <- estimateDisp(dge, GLM_Model)
    
    write.table(data.frame(Gene_ID=rownames(dge$counts),Tagwise_Dispersion=dge$tagwise.dispersion),file=paste0(path,"/","Estimated_Dispersion.txt"),row.names=FALSE,col.names=TRUE,sep="\t")
  ## inference: parameter estimation
    fit <- glmFit(dge, GLM_Model)
    fittedval <- as.data.frame(log2(fit$fitted.values))
    fittedval <- tibble::rownames_to_column(fittedval, "Gene_ID")
    write.table(fittedval,file=paste0(path,"/","Fitted_Values.txt"),row.names=FALSE,col.names=TRUE,sep="\t")
  ## Control quality of the differential analysis
    lrt<-list()
    pdf(paste0(path,"/","Raw_pvalues_histograms.pdf"))
    par(mfrow=c(2,2))
    rownames(Contrasts.interest) <- Contrasts.interest[,1]
    Contrasts.interest <- Contrasts.interest[,-1]
    GLM_Model <- as.matrix(GLM_Model)
    Contrasts.interest <- as.matrix(Contrasts.interest)
    for (i in 1:(nrow(Contrasts.interest)))
    {
        lrt[[i]] <- glmLRT(fit,contrast=Contrasts.interest[i,])
        hist(lrt[[i]]$table$PValue,100,main=rownames(Contrasts.interest)[i],xlab="Raw pvalue",cex.main=0.8)
    }
    dev.off()
    
    ## Results of differential expression analysis
    Compare = dge$genes 
    Down_Up_table <- data.frame("Contrast"=rep(rownames(Contrasts.interest),each=2),"Expression"=rep(c("Up","Down"),nrow(Contrasts.interest)),"Nb_DEG"=0)
    DEG_table<- data.frame("Contrast"=rownames(Contrasts.interest),"Nb_DEG"=0)
    
    for (i in 1:nrow(Contrasts.interest))
    {
        Compare[,i+1]=0
        ## subdirectory for results
        subdirectory_contrast <- rownames(Contrasts.interest)[i]
        if (!I(subdirectory_contrast %in% dir(path))) dir.create(paste0(path,"/",subdirectory_contrast), showWarnings=FALSE)
        ## Table of results
        res=topTags(lrt[[i]],n=nrow(lrt[[i]]$genes),adjust.method="BH",sort.by="none")$table
        
        ## Add Gene_Name
        if(!is.null(Annotation))
            res <- merge(Annotation, res, by.x="Gene_ID", by.y="genes",all.y=TRUE)
        ## Write table for results
        fileout <- paste0(path,"/",subdirectory_contrast,"/","LRT_BH.txt")
        write.table(res,fileout,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
        if(sum(res$FDR<=Alpha_DiffAnalysis)!=0)
        {
            ## Table of Differential expressed genes by BH methods (DE.BH)
            DE.BH <- res[which(res$FDR<=Alpha_DiffAnalysis),]
            DE.BH <- DE.BH[order(DE.BH$FDR),]
            ## All DE.BH id genes
            IdGenes <- as.character(DE.BH[,1])
            ## Write table for DE.BH
            fileout <- paste0(path,"/",subdirectory_contrast,"/","DEG.BH.txt")
            write.table(DE.BH,fileout,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
            ## Write ID of DE genes 
            fileout <- paste0(path,"/",subdirectory_contrast,"/","Id_DEG.txt")
            write.table(IdGenes,fileout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
            
            ## Gene division and Compare table
            #Compare[is.element(Compare[,1],IdGenes),i+1]=1
            DEG_table[i,2] <- nrow(DE.BH)
            if(sum(is.element(Compare[,1],DE.BH[which(sign(DE.BH$logFC)==(-1)),1]))!=0)
                Compare[is.element(Compare[,1],DE.BH[which(sign(DE.BH$logFC)==(-1)),1]),i+1]=-1
            if(sum(is.element(Compare[,1],DE.BH[which(sign(DE.BH$logFC)==(1)),1]))!=0)
                Compare[is.element(Compare[,1],DE.BH[which(sign(DE.BH$logFC)==1),1]),i+1]=1
             
            Down_Up_table[(2*i-1),3] <- sum(sign(DE.BH$logFC)==1)
            Down_Up_table[2*i,3] <- sum(sign(DE.BH$logFC)==(-1))
            ## NormCounts table for DEG 
            DEG_NormCounts <- NormCounts[IdGenes,]
            DEG_log2NormCounts <- log2(DEG_NormCounts+1)
            ## Add Gene_Name
            if(!is.null(Annotation))
            {
                if (length(IdGenes) == 1) 
                    Gene_ID <- IdGenes
                else 
                    Gene_ID <- rownames(DEG_NormCounts)
                
                DEG_log2NormCounts <- data.frame(Gene_ID,DEG_log2NormCounts)
                DEG_log2NormCounts <- merge(DEG_log2NormCounts, Annotation, by="Gene_ID", all.x = TRUE,sort=FALSE) 
                DEG_log2NormCounts <- DEG_log2NormCounts[,c(1,ncol(DEG_log2NormCounts),2:(ncol(DEG_log2NormCounts)-1))]
            }
            
            ## plotSmear
            pdf(paste0(path,"/",subdirectory_contrast,"/","plotSmear.pdf"))
            plotSmear(lrt[[i]], de.tags = IdGenes ,cex=0.5,ylim=c(-8,8))
            dev.off()
            
            ## Differential Expressed Gene profiles
            nbGenes<-pmin(NbGenes_Profiles,nrow(DE.BH))
            pdf(paste0(path,"/",subdirectory_contrast,"/","Top",nbGenes,"_Profile.pdf"))
            if(nbGenes==1)
            {
                GeneID <- IdGenes[1:nbGenes]
                title <- GeneID
                if(!is.null(Annotation))
                {
                    title <- paste0(GeneID,"\n",DE.BH$Gene_Name)
                    table <- data.frame(Target[,1:NbBioFactors],t(DEG_log2NormCounts[-(1:2)]))
                }
                else
                    table <- data.frame(Target[,1:NbBioFactors],DEG_log2NormCounts)
                colnames(table)[ncol(table)]="log2NormCounts"
                if (NbBioFactors == 1){
                    p <- ggline(table, x = colnames(table[1]), y = "log2NormCounts", add = "mean_se",title=title, xlab = "Condition")
                }
                if (NbBioFactors == 2){
                    p <- ggline(table, x = colnames(table[2]), y = "log2NormCounts", add = "mean_se", color = colnames(table[1]), xlab = "Condition",title=title)
                }
                print(p)
            }
            if(nbGenes!=1)
            {
                GeneID <- IdGenes[1:nbGenes]
                title <- GeneID
                if(!is.null(Annotation))
                    title <- paste0(GeneID,"\n",DE.BH$Gene_Name)
                
                ## Profile graphs
                for (j in 1:nbGenes)
                {
                    if(!is.null(Annotation))
                        table <- data.frame(Target[,1:NbBioFactors],t(DEG_log2NormCounts[j,-(1:2)]))
                    else
                        table <- data.frame(Target[,1:NbBioFactors],DEG_log2NormCounts[j,])
                    colnames(table)[ncol(table)]="log2NormCounts"
                    if (NbBioFactors == 1){
                        p <- ggline(table, x = colnames(table[1]), y = "log2NormCounts", add = "mean_se",title=title[j], xlab = "Condition")
                    }
                    if (NbBioFactors == 2){
                        p <- ggline(table, x = colnames(table[2]), y = "log2NormCounts", add = "mean_se", color = colnames(table[1]), xlab = "Condition",title=title[j])
                    }
                    print(p)
                }
            }
            dev.off()
            ## Gene clustering (euclidean distance / Ward linkage)
            nbGenes=pmin(NbGenes_Clustering,nrow(DE.BH))
            if(nbGenes>5)
            {
                ClustID <-  as.character(DE.BH[1:nbGenes,1])
                if(!is.null(Annotation)){
                    Annot <- paste0(DE.BH[,1],"-",DE.BH[,2])[1:nbGenes]
                } else {
                    Annot=DE.BH[1:nbGenes,1]
                }
                
                if(NbBioFactors == 2){
                    Samples <- apply(Target[,1:NbBioFactors],1,function(x) paste(x,collapse="_"))
                } else {
                    Samples <- Target[,1]
                }
                
                
                mycol <- colorpanel(1000,"blue","white","red")
                 
                if(!is.null(Annotation))
                {
                    rownames(DEG_log2NormCounts) <- DEG_log2NormCounts[,1]
                    DEG_log2NormCounts <- DEG_log2NormCounts[,-c(1:2)]
                }
                
                pdf(paste0(path,"/",subdirectory_contrast,"/","Top",nbGenes,"_Clustering.pdf"))
                heatmap.2(as.matrix(DEG_log2NormCounts[ClustID,]),
                          hclustfun = function(x) hclust(x, method = "ward.D"),
                          scale="row", labRow=Annot, labCol=Samples,
                          cexCol = 0.5,
                          col=mycol, trace="none", density.info="none",
                          dendrogram="both", key = FALSE, lhei=c(2,10),
                          margins = c(4,18))
                dev.off()
            }
        }
    }
    ## Comparisons of results
    colnames(Compare) <- c("Gene_ID",rownames(Contrasts.interest))
    
    ## Save results
    fileout <- paste0(path,"/","Compare_table.txt")
    write.table(Compare,fileout,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
    
    write.csv2(DEG_table, paste0(path,"/","DiffAnalysis_Comparisons.csv"))

    cat("################################################\n Results of the differential analysis \n################################################\n\n")
    cat("Number of DEGs for each contrast:\n")
    print(DEG_table)
    cat("\n")
    print(Down_Up_table)
    cat("################################################\n\n")
    

    
    pdf(paste0(path,"/","Down_Up_DEG.pdf"))
    title <- "Number of differentially genes expressed for each contrast"
    xlab <- "Contrasts"
    ylab <- "Number of DEG"
    Down_Up_table$z <- 1:nrow(Down_Up_table)
    Down_Up_table <- arrange(ddply(Down_Up_table, .(Contrast), transform, pos = cumsum(Nb_DEG) - (0.5 * Nb_DEG), z=z), z)
    p <- ggplot(data=Down_Up_table, aes(x=Contrast, y=Nb_DEG, fill=Expression)) 
    p <- p + geom_bar(stat="identity")
    p <- p + scale_fill_manual(values=c("mediumturquoise","indianred2"))
    p <- p + scale_x_discrete(limits=Down_Up_table$Contrast)
    p <- p + labs(title = title, x=xlab, y = ylab)
    p <- p + theme(axis.text.x=element_text(angle = +45, vjust = 1, hjust=1))
    p <- p + geom_text(data=Down_Up_table, aes(x = Contrast, y = pos, label = Nb_DEG), size=1)
    print(p)
    dev.off()
    
    ## Add Annotations
    if(!is.null(Annotation)){
        NormCounts_log2 <- merge(Annotation, NormCounts_log2, by="Gene_ID", all.y=TRUE)
        ## Mean SD
        NormCounts_log2_Mean_SD <- merge(Annotation, NormCounts_log2_Mean_SD, by="Gene_ID", all.y=TRUE)
    }
    
### Write All NormCounts tables

    ## Log2 NormCounts table
    fileout <- paste0(path,"/","NormCounts_log2.txt")
    write.table(NormCounts_log2,fileout,row.names=FALSE,col.names=TRUE, sep="\t", quote = FALSE)

    ## Mean/SD Log2 table
    fileout <- paste0(path,"/","NormCounts_log2_Mean_SD.txt")
    write.table(NormCounts_log2_Mean_SD,fileout,row.names=FALSE,col.names=TRUE, sep="\t", quote = FALSE)
}
