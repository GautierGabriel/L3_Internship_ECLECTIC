###### Co-expression analysis using Coseq ######
Coexpression_coseq <- function(Data_Directory, Results_Directory, Project_Name,
                               Title, Groups, Operation, Target, Raw_Counts, Color_Group=NULL,
                               A=5, B=40, K=c(2,seq(5,30,by=5)), ...)
{
    ## subdirectories
    if (!I(Project_Name %in% dir(Results_Directory)))
        dir.create(paste0(Results_Directory,"/",Project_Name), showWarnings=FALSE)
    subdirectory <- "Coexpression"
    path=paste0(Results_Directory,"/",Project_Name,"/",subdirectory)
    if (!I(subdirectory %in% dir(paste0(Results_Directory,"/",Project_Name))))
        dir.create(path, showWarnings=FALSE)
    if (!I(Title %in% dir(paste0(Results_Directory,"/",Project_Name,"/",subdirectory))))
        dir.create(paste0(path,"/",Title), showWarnings=FALSE)
    
    ## Load Compare table
      Compare <- read.csv2(paste0(Results_Directory,"/",Project_Name,"/DiffAnalysis/Compare_table.txt"), header=TRUE, sep="\t",check.names=FALSE)
    
    ## check if the groups are really studied contrasts
    if (sum(Groups%in%colnames(Compare[,Groups]))!=length(Groups))
    {
      cat("Groups is bad defined. The constrast(s) \n")
      cat(Groups[!Groups%in%colnames(Compare)])
      cat("do(es) not exist \n")
      break
    }
    
    ## Union List
    if (Operation=="Union"){
      newCol <- ncol(Compare)+1
      Compare[,newCol]=apply(Compare[,Groups],1,function(a) 1*(sum(abs(a))!=0)) 
      colnames(Compare)[newCol]<- paste0(Title,"_Union")
      IdGene <- Compare$Gene_ID[which(Compare[,newCol]==1)]
      }
    if (Operation=="Intersection"){
      newCol <- ncol(Compare)+1
      Compare[,newCol]=apply(Compare[,Groups],1,function(a) 1*(sum(abs(a))==length(Groups))) 
      colnames(Compare)[newCol]<- paste0(Title,"_Intersection")
      IdGene <- Compare$Gene_ID[which(Compare[,newCol]==1)]
    }
      fileout <- paste0(path,"/",Title,"/","Gene_List.txt")
      write.table(IdGene,fileout,row.names=FALSE,col.names=FALSE, sep="\t", quote = FALSE)
      
      ## Add a column to the compare table
    if (file.exists(paste0(path,"/","Coexpression_Input_table.txt")))
    {
      tmp=read.table(paste0(path,"/","Coexpression_Input_table.txt"),header=TRUE,sep="\t",check.names = FALSE)
      Compare=merge(Compare,tmp)    
    }     
    fileout <-  paste0(path,"/","Coexpression_Input_table.txt")
    write.table(Compare,fileout,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

  ## List of gene for the coexpression analysis = subset
    subset<-(1:nrow(Raw_Counts))[is.element(rownames(Raw_Counts),IdGene)]
  
  ###### First grid to select a region for the number of components (K)
    Results.logLike <- data.frame()
    Results.icl<- data.frame()
    Results.1 <- list()
    Results.1_min_icl <- list()
   
    for (a in 1:A){
        try(coseq.res <- coseq(Raw_Counts, K=K, model="Normal", transformation="arcsin",
                               normFactors="TMM", subset=subset, seed=6789, ...))
        metadata <- metadata(coseq.res)
        Results.logLike[1:length(K),a]<- data.frame(metadata$logLike)
        Results.icl[1:length(K),a]<- data.frame(metadata$ICL)
        Results.1[[a]]=coseq.res
        Results.1_min_icl[[a]] <- min(metadata(Results.1[[a]])$ICL,na.rm=TRUE)
    }
  ## Maximum logLike and minimum ICL 
    colnames(Results.logLike) <- paste0("iter_",1:A)
    colnames(Results.icl) <- paste0("iter_",1:A)
    icl=apply(Results.icl,1,min,na.rm=TRUE)
    logLike=apply(Results.logLike,1,max,na.rm=TRUE)
    
  ## Table of loglike and table of ICL 
    Output <- file(paste0(path,"/",Title,"/","Results_First_Loop.txt"), open="wt")
    sink(Output)
    sink(Output, type="message")
    cat("################################################\n")
    cat("Co-expression analysis \n")
    cat("################################################\n\n")
    cat(paste0("Number of genes for the co-expression analysis: ",length(IdGene),"\n\n"))
    cat("########## First Loop ##########\n")
    cat("\n\nLogLikelihood table:\n")
    print(Results.logLike)
    cat("\nMaximum LogLikelihood table:\n")
    print(logLike)
    cat("\n\nICL table:\n")
    print(Results.icl)
    cat("\nMinimum ICL table:\n")
    print(icl)
    sink(type="message")
    sink()
    
  ## LogLike and ICL curves (pdf plot)
    pdf(paste0(path,"/",Title,"/","Loop_1.pdf"))
    matplot(K,icl,pch=2,ylab="ICL",xlab="Number of clusters")
    dev.off()
  
  ## Second loop to perform coexpression analysis on the new set of K 
  ## find the min ICL in the first loop results = index
    index=which.min(icl)
## Selection of new set of K close to the index
    if (index!=1 && index!=length(K)){
        K <- (K[index-1]+1):(K[index+1]-1)
    }
    if (index==1){
        K <- 2:K[index+1]
    }
    
    if (index==length(K)){
        cat("The first research was not relevant.
 The research will be done between",  K[length(K)]," and ",K[length(K)]+20,"\n")
        K <- K[length(K)]:(K[length(K)]+20)
    }
  
  ## Perform coexpression analysis on new set of K
    Results.2 <- list()
    print("DEBUT")
    for (k in K)
    {
        try(ref<-coseq(Raw_Counts, K=k, model="Normal", transformation="arcsin", 
                       normFactors="TMM", subset=subset, seed=6798, ...))
        for (b in 1:B)
        {

            try(new <- coseq(Raw_Counts, K=k, model="Normal", transformation="arcsin", 
                            normFactors="TMM", subset=subset, seed=6789,...))

            if(ICL(ref)>ICL(new)){ref=new}
        }
        Results.2[[k]]=ref
        ## save RData for the second loop results
        save(Results.2, file = paste0(path,"/",Title,"/","coseq_loop_2.RData"))
    }

    
    ## Results (text and plot files)
    Output <- file(paste0(path,"/",Title,"/","Results_Second_Loop.txt"), open="wt")
    sink(Output)
    sink(Output, type="message")
    cat("################################################\n")
    cat("Co-expression analysis \n")
    cat("################################################\n\n")
    cat("Index:\n")
    print(index)
    cat("\nInterval of K for the second loop:\n")
    print(K)
    cat("\n\n########## Second Loop ##########\n")
    cat("\nMinimum ICL values for each K:\n")
    print(sapply(Results.2,ICL)[K])
    cat("\nMinimum ICL value:\n")
    print(which.min(sapply(Results.2,ICL)))
    sink(type = "message")
    sink()
    
    pdf(paste0(path,"/",Title,"/","Loop_2.pdf"))
    matplot(K,sapply(Results.2,ICL)[K],ylab="ICL",xlab="Number of clusters")
    dev.off()
  
### Final Results
    Output <- file(paste0(path,"/",Title,"/","Results_Final.txt"), open="wt")
    sink(Output)
    sink(Output, type="message")
    cat("################################################\n")
    cat("Co-expression analysis \n")
    cat("################################################\n\n")
    
    ## Find the min icl from all results (first and second loop)
    Results.1_final <- min(unlist(Results.1_min_icl))
    Results.2_final <- min(sapply(Results.2,ICL), na.rm=TRUE)
    if(Results.1_final < Results.2_final){
    final <- Results.1[[which.min(Results.1_min_icl)]]
  } else{ final <- Results.2[[which.min(sapply(Results.2,ICL))]]}
    
### Final.Rdata
    save(final, file = paste0(path,"/",Title,"/","coseq_final.RData"))
    
    ## Text file output with the summary of final results
    cat("########## Final Results ##########\n")
    cat("\nSummary coseq:\n")
    print(summary(final))
    cat("\n")
  
  ## Write ID Gene list for all clusters 
    Clusters <- clusters(final)
    Clusters <- data.frame(Clusters)
    for (b in 1:max(Clusters))
    {
        Id_Genes <- rownames(Clusters)[which(Clusters == b)]
        fileout <- paste0(path,"/",Title,"/","Cluster",b,"_GeneID.txt")
        write.table(Id_Genes,fileout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    }
    
    ## Write a summary of clusters results
    Gene_ID <- rownames(Clusters)
    rownames(Clusters) <- NULL
    Summary_coseq <- cbind(Gene_ID,Clusters)
    IdCoseq <- Summary_coseq[,1]
    IdGene <- as.factor(IdGene)
    diff <- setdiff(IdGene, IdCoseq)
    Nb_no_cluster <- length(diff)
    cat("\n Number of genes not included in the co-expression analysis:\n")
    print(Nb_no_cluster)
    if(Nb_no_cluster != 0){
        diff_table <- data.frame(Gene_ID=diff,Clusters=0)
        Summary_coseq <- rbind(Summary_coseq, diff_table)
    }
  
  # Add annotation
    if(!is.null(Annotation))
    {
        Summary_coseq <- merge(Annotation, Summary_coseq, by="Gene_ID", all.y=TRUE)
    }
    
    ## All maximum conditional probabilities 
    Proba_Cond <- assay(final)
    Max_Proba_Cond <- data.frame()
    for (i in 1:nrow(Proba_Cond)){
        Max_Proba_Cond[i,1] <- rownames(Proba_Cond)[i]
        Max_Proba_Cond[i,2] <- max(Proba_Cond[i,])
    }
    colnames(Max_Proba_Cond) <- c("Gene_ID","Proba_Cond")
    Summary_coseq <- merge(Summary_coseq, Max_Proba_Cond, by = "Gene_ID", all.x=TRUE)
    cat("\n")
    cat("Head of Summary_coseq table :\n")
    print(head(Summary_coseq),10)
    cat("\n")
    sink(type="message")
    sink() 
    
    fileout <- paste0(path,"/",Title,"/","AllClusters.txt")
    write.table(Summary_coseq,fileout,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
    
### graphical outputs for the exploration of clusters
    ## Nb of biological factors, genes and samples 
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
    
    ## Conditions
    Conditions=rep(BioFactors,NbReplicate)
    ## Vector of colors for biological factors
    if(is.null(Color_Group))
    {
        if (NbConditions<=12 && NbConditions>=3)
            Color_Group <- brewer.pal(NbConditions, name = ifelse(nlevels(Target[,1])==2,'Paired','Set3'))
        else
            Color_Group <-colors()[seq(from=3,length=NbConditions,by=4)]
    }
    
    ## all plots
    pdf(paste0(path,"/",Title,"/","Final_Coseq.pdf"))
    p <- plot(final, conds=Conditions, #collapse_reps = "average",
              graphs = c("probapost_histogram","probapost_boxplots","probapost_barplots","profiles"),add_lines = FALSE)
    print(p)
    p <- plot(final, conds=Conditions, graphs = "boxplots", add_lines = FALSE)
    p <- p$boxplots + scale_fill_manual(values = Color_Group) + theme(axis.text.x = element_blank())
    print(p)
    dev.off()
    
  ## boxplot of profiles for each cluster
    pdf(paste0(path,"/",Title,"/","Boxplot_profiles_Coseq.pdf"))
    p <- plot(final, conds=Conditions, collapse_reps = "average", graphs = c("boxplots"), add_lines=FALSE)
    p <- p$boxplots + scale_fill_manual(values = Color_Group) + theme(axis.text.x = element_blank())
    print(p)
    dev.off()
    
    Clusters_sizes <- as.data.frame(table(Clusters))
    colnames(Clusters_sizes) <- c("Cluster","Number of genes")
    cat("################################################\nNumber of genes for each clusters\n################################################\n\n")
    print(Clusters_sizes)
   
}



