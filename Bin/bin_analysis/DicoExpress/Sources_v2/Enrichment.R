#### Clusters Enrichment analysis
enrichment_analysis <- function(Reference,Gene_List,Alpha,Subset)
  {
    if(nrow(Gene_List)==0)
      return(NULL)
    ## success in the urn /  For each annotation term, number of annotated genes in the Reference file
    m=table(Reference[,2])
    ## failures in the urn / For each annotation term, number of not annotated genes in the Reference file
    n=length(unique(Reference[,1]))-m

    trial<-merge(Gene_List,Reference)
    ## number of genes in the gene list also observed in the Reference
    k=length(unique(trial[,1]))
    ## trial success /  For each annotation term, number of annotated genes in the gene list file
    x=table(factor(trial[,2],levels=rownames(m)))
    ## Result files
    res=NaN
    Term=rownames(m)
    m=as.numeric(m)
    n=as.numeric(n)
    x=as.numeric(x)
    res=data.frame(Subset,Term,Urn_Success=m,Urn_Failures=n,Trial_Success=x,Trial_effective=k,Urn_percentage_Success=signif(100*m/(m+n),3),Trial_percentage_Success=signif(100*x/k,3), Pvalue_over=phyper(x-1,m,n,k,lower.tail=FALSE))
    res=data.frame(res,Bonferroni_over=p.adjust(res$Pvalue_over,method="bonferroni"),BH_over=p.adjust(res$Pvalue_over,method="BH"),
                   Pvalue_under=phyper(x,m,n,k,lower.tail=TRUE))
    res=data.frame(res,Bonferroni_under=p.adjust(res$Pvalue_under,method="bonferroni"),BH_under=p.adjust(res$Pvalue_under,method="BH"))
    
    res_over_under <-NULL
    index=which(res$Bonferroni_over<Alpha)
    if(length(index)!=0)
    {
      res_over <- res[index,]
      res_over[,ncol(res_over)+1] <- "overrepresented"
      colnames(res_over)[ncol(res_over)] <- c("Bonferroni_Decision")
      res_over_under <- res_over
    }
    
    index=which(res$Bonferroni_under<Alpha)
    if(length(index)!=0)
    {
      res_under <- res[index,]
      res_under[,ncol(res_under)+1] <- "underrepresented"
      colnames(res_under)[ncol(res_under)] <- c("Bonferroni_Decision")
      res_over_under <- rbind(res_over_under,res_under)
    }
    Bonferroni_over_under = res_over_under
    res_over_under <-NULL
    index=which(res$BH_over<Alpha)
    if(length(index)!=0)
    {
      res_over <- res[index,]
      res_over[,ncol(res_over)+1] <- "overrepresented"
      colnames(res_over)[ncol(res_over)] <- c("BH_Decision")
      res_over_under <- res_over
    }
    
    index=which(res$BH_under<Alpha)
    if(length(index)!=0)
    {
      res_under <- res[index,]
      res_under[,ncol(res_under)+1] <- "underrepresented"
      colnames(res_under)[ncol(res_under)] <- c("BH_Decision")
      res_over_under <- rbind(res_over_under,res_under)
    }
    BH_over_under = res_over_under
    Results=list("All_Results" = res, "Significant_Results" = merge(Bonferroni_over_under,BH_over_under,all=TRUE))
    return(Results)
}





Enrichment <- function(Results_Directory, Project_Name, Title, Reference_Enrichment,
                       Alpha_Enrichment = 0.05, append=FALSE){

    if(is.null(Reference_Enrichment))
    {
        stop("\n Please put a file ", paste0(Project_Name,"_Enrichment.csv"), " in the Data directory and re-run the function Load_Data_Files before performing an enrichment test\n")
    }
  
    File=read.table(paste0(Results_Directory,"/",Project_Name,"/DiffAnalysis/","Compare_table.txt"),h=T,sep="\t",check.names = FALSE)
#   ## restriction of the reference file to the genes present int the compare table 
    Ref=merge(Reference_Enrichment,File,by="Gene_ID",all.y=TRUE)[,1:2]
    
    if(is.null(Title))
    {
        Contrast_Name=colnames(File)[-1]
        ## Summary results of all DEG lists
        Nb_Term=length(unique(Ref$Annotation))
        Resume_Bonferroni=data.frame(Term=unique(Ref$Annotation),matrix(0,nrow=Nb_Term,ncol=length(Contrast_Name)+1))
        colnames(Resume_Bonferroni)=c("Term",Contrast_Name,"Total")
        Resume_BH=Resume_Bonferroni
        for (cc in 1:length(Contrast_Name))
        {
              ## Perform the 3 enrichment analyses
                Gene_List <- data.frame(Gene_ID=File[which(File[,cc+1]==1),1])
                Results_Up<- enrichment_analysis(Ref,Gene_List,Alpha_Enrichment,Subset="Up")
                Gene_List <- data.frame(Gene_ID=File[which(File[,cc+1]==(-1)),1])
                Results_Down<- enrichment_analysis(Ref,Gene_List,Alpha_Enrichment,Subset="Down")
                Gene_List <- data.frame(Gene_ID=File[which(abs(File[,cc+1])==1),1])
                Results_DE<- enrichment_analysis(Ref,Gene_List,Alpha_Enrichment,Subset="DE")
                ## Save all results
                path=paste0(Results_Directory,"/",Project_Name,"/DiffAnalysis/",Contrast_Name[cc],"/")
                fileout=paste0(path,"All_Enrichment_Results.txt")
                col.names.boolean=(!file.exists(fileout)||!append)
                tmp=rbind(Results_Up$All_Results,Results_Down$All_Results,Results_DE$All_Results)
                if(!is.null(tmp))
                  write.table(tmp,fileout,row.names=FALSE,col.names=col.names.boolean, sep="\t", quote = FALSE, append=append)
                

                tmp <- rbind(Results_Up$Significant_Results,Results_Down$Significant_Results,Results_DE$Significant_Results)
                index.col<-is.element(colnames(Resume_Bonferroni),Contrast_Name[cc])
                if(!is.null(tmp))
                  {
                    fileout=paste0(path,"Significant_Enrichments.txt")
                    col.names.boolean=(!file.exists(fileout)||!append)
                    write.table(tmp,fileout,col.names=col.names.boolean,row.names=FALSE, sep="\t", quote = FALSE, append=append)
                    
                    res.character=rep(0,Nb_Term)
                    
                    ## UP ##
                    index.row <-(1:Nb_Term)[is.element(Resume_Bonferroni$Term,Results_Up$Significant_Results$Term[which(Results_Up$Significant_Results$Bonferroni_Decision=="overrepresented")])]
                    if (length(index.row) == 0) {
                      index.row <- 1:Nb_Term
                      res.character[index.row]=0                    }
                    else{
                      res.character[index.row]=1
                    }
                    
                    ## DOWN ##
                    index.row <-(1:Nb_Term)[is.element(Resume_Bonferroni$Term,Results_Down$Significant_Results$Term[which(Results_Down$Significant_Results$Bonferroni_Decision=="overrepresented")])]
                    if (length(index.row) == 0) {
                      index.row <- 1:Nb_Term
                      res.character[index.row]=paste0(res.character[index.row],"0")
                    }
                    else{
                      res.character[index.row]=paste0(res.character[index.row],"1")
                      res.character[-index.row]=paste0(res.character[-index.row],"0")
                    }
 
                    ## DE ##
                    
                    index.row <-(1:Nb_Term)[is.element(Resume_Bonferroni$Term,Results_DE$Significant_Results$Term[which(Results_DE$Significant_Results$Bonferroni_Decision=="overrepresented")])]
                    if (length(index.row) == 0) {
                      index.row <- 1:Nb_Term
                      res.character[index.row]=paste0(res.character[index.row],"0")
                    }
                    else{
                      res.character[index.row]=paste0(res.character[index.row],"1")
                      res.character[-index.row]=paste0(res.character[-index.row],"0")
                    }
                    
                    Resume_Bonferroni[,index.col] <- res.character
                }
                else
                  Resume_Bonferroni[,index.col] <- "000"
                
                
                tmp <- rbind(Results_Up$Significant_Results,Results_Down$Significant_Results,Results_DE$Significant_Results)
                index.col<-is.element(colnames(Resume_BH),Contrast_Name[cc])
                if(!is.null(tmp))
                {
                  res.character=rep(0,Nb_Term)
                  
                  ## UP ##
                  index.row <-(1:Nb_Term)[is.element(Resume_BH$Term,Results_Up$Significant_Results$Term[which(Results_Up$Significant_Results$BH_Decision=="overrepresented")])]
                  if (length(index.row) == 0) {
                    index.row <- 1:Nb_Term
                    res.character[index.row]=0
                  }
                  else{
                    res.character[index.row]=1
                  }
                  
                  ## DOWN ##
                  index.row <-(1:Nb_Term)[is.element(Resume_BH$Term,Results_Down$Significant_Results$Term[which(Results_Down$Significant_Results$BH_Decision=="overrepresented")])]
                  if (length(index.row) == 0) {
                    index.row <- 1:Nb_Term
                    res.character[index.row]=paste0(res.character[index.row],"0")
                  }
                  else{
                    res.character[index.row]=paste0(res.character[index.row],"1")
                    res.character[-index.row]=paste0(res.character[-index.row],"0")
                  }
                  
                  ## DE ##
                  index.row <-(1:Nb_Term)[is.element(Resume_BH$Term,Results_DE$Significant_Results$Term[which(Results_DE$Significant_Results$BH_Decision=="overrepresented")])]
                  if (length(index.row) == 0) {
                    index.row <- 1:Nb_Term
                    res.character[index.row]=paste0(res.character[index.row],"0")
                  }
                  else{
                    res.character[index.row]=paste0(res.character[index.row],"1")
                    res.character[-index.row]=paste0(res.character[-index.row],"0")
                  }
                  
                  Resume_BH[,index.col] <- res.character
                }
                else
                  Resume_BH[,index.col] <- "000"
                
        }
        if(ncol(Resume_Bonferroni)==3){
          
          Resume_Bonferroni$Total <- ifelse(Resume_Bonferroni[,2] != "000", 1, 0)
          Resume_Bonferroni=Resume_Bonferroni[which(Resume_Bonferroni$Total!=0),]
        }
        else
        {
          Resume_Bonferroni$Total <- apply(Resume_Bonferroni[, -c(1, ncol(Resume_Bonferroni))], 1, function(x) sum(x != "000"))
          Resume_Bonferroni=Resume_Bonferroni[which(Resume_Bonferroni$Total!=0),]
        }

        fileout=paste0(Results_Directory,"/",Project_Name,"/DiffAnalysis/","Summary_Bonferroni_Overrepresentation.txt")
        write.table(Resume_Bonferroni,file=fileout,row.names=FALSE,col.names=!append, sep="\t", quote = FALSE, append=append)
        
        
        if(ncol(Resume_BH)==3){
          Resume_BH$Total <- ifelse(Resume_BH[,2] != "000", 1, 0)
          Resume_BH=Resume_BH[which(Resume_BH$Total!=0),]
        }
        else
        {        
          Resume_BH$Total=apply(Resume_BH[,-c(1,ncol(Resume_BH))],1,function(x) sum(x!="000"))
          Resume_BH=Resume_BH[which(Resume_BH$Total!=0),]
        }
        
        
        
        fileout=paste0(Results_Directory,"/",Project_Name,"/DiffAnalysis/","Summary_BH_Overrepresentation.txt")
        write.table(Resume_BH,file=fileout,row.names=FALSE,col.names=!append, sep="\t", quote = FALSE, append=append)
    }
    
    ## Title != NULL
    else
    {
#       ## Load Cluster file (in Coexpression directory)
        path=paste0(Results_Directory,"/",Project_Name,"/Coexpression/",Title,"/")
        Clusters <- read.csv2(paste0(path,"AllClusters.txt"), header = TRUE, sep="\t",quote="")
        Clusters <- Clusters[,c("Gene_ID","Clusters")]
        ##  Summary results of all clusters
        Nb_Term=length(unique(Ref$Annotation))
        Resume_Bonferroni=data.frame(Term=unique(Ref$Annotation),matrix(0,nrow=Nb_Term,ncol=max(Clusters$Clusters)+2))
        colnames(Resume_Bonferroni)=c("Term",paste0("Cluster_",1:max(Clusters$Clusters)),"AllClusters","Nb_enriched_clusters")
        Resume_BH=Resume_Bonferroni
        ## Perform enrichment analysis
        Gene_List <-as.data.frame(Clusters[,1])
        colnames(Gene_List) <- c("Gene_ID")
        Results <- enrichment_analysis(Ref,Gene_List,Alpha_Enrichment,Subset="Coexpression")
 
         ## Save results
        fileout=paste0(path,"AllClusters_All_Enrichment_Results.txt")
        write.table(Results$All_Results,file=fileout,row.names=FALSE,col.names=!append, sep="\t", quote = FALSE,append=append)
        if(!is.null(Results$Significant_Results))
        {
            ## Save significant results (over and under enriched)
            fileout=paste0(path,"AllClusters_Significant_Enrichments.txt")
            col.names.boolean=(!file.exists(fileout)||!append)
            write.table(Results$Significant_Results,fileout,col.names=col.names.boolean,row.names=FALSE, sep="\t", quote = FALSE,append=append)
            ## Summary results 
           
            index.col=grep("AllClusters", colnames(Resume_Bonferroni))
            index.row <-is.element(Resume_Bonferroni$Term,Results$Significant_Results$Term[which(Results$Significant_Results$Bonferroni_Decision=="overrepresented")])
            Resume_Bonferroni[index.row,index.col] <- 1
            index.row <-is.element(Resume_BH$Term,Results$Significant_Results$Term[which(Results$Significant_Results$BH_Decision=="overrepresented")])
            Resume_BH[index.row,index.col] <- 1
            }
                   
        ## Enrichment for each cluster
        for (i in 1:max(Clusters$Clusters))
        {
          
#           # Keep only the genes from the co-expression of each cluster
            Reference_Enrichment = Ref[Ref$Gene_ID %in% Gene_List$Gene_ID, ]
            
            
            Gene_List <-as.data.frame(Clusters[which(Clusters$Clusters == i),1])
            colnames(Gene_List) <- c("Gene_ID")

            
            ## Perform enrichment analysis
            Results <- enrichment_analysis(Reference_Enrichment,Gene_List,Alpha_Enrichment,Subset=paste0("Cluster_",i))
                      
            ## Save results
            fileout=paste0(path,"Cluster_",i,"_All_Enrichment_Results.txt")
            write.table(Results$All_Results,fileout,row.names=FALSE,col.names=!append, sep="\t", quote = FALSE,append=append)
            ## Save significant results (over and under enriched)
            if(!is.null(Results$Significant_Results))
            {
                fileout=paste0(path,"Cluster_",i,"_Significant_Enrichments.txt")
                col.names.boolean=(!file.exists(fileout)||!append)
                write.table(Results$Significant_Results,fileout,col.names=col.names.boolean,row.names=FALSE, sep="\t", quote = FALSE,append=append)
                index.col=grep(paste0("Cluster_",i), colnames(Resume_Bonferroni))
                index.row <-is.element(Resume_Bonferroni$Term,Results$Significant_Results$Term[which(Results$Significant_Results$Bonferroni_Decision=="overrepresented")])
                Resume_Bonferroni[index.row,index.col] <- 1
                index.row <-is.element(Resume_BH$Term,Results$Significant_Results$Term[which(Results$Significant_Results$BH_Decision=="overrepresented")])
                Resume_BH[index.row,index.col] <- 1
            }
        }
        Resume_Bonferroni$Nb_enriched_clusters=rowSums(Resume_Bonferroni[-c(1,(ncol(Resume_Bonferroni)-1),ncol(Resume_Bonferroni))])
        Resume_Bonferroni=Resume_Bonferroni[which(Resume_Bonferroni$Nb_enriched_clusters!=0),]
        fileout=paste0(path,"Summary_Cluster_Bonferroni_Overrepresentation.txt")
        write.table(Resume_Bonferroni,fileout,row.names=FALSE,col.names=!append, sep="\t", quote = FALSE,append=append)
        
        Resume_BH$Nb_enriched_clusters=rowSums(Resume_BH[-c(1,(ncol(Resume_BH)-1),ncol(Resume_BH))])
        Resume_BH=Resume_BH[which(Resume_BH$Nb_enriched_clusters!=0),]
        fileout=paste0(path,"Summary_Cluster_BH_Overrepresentation.txt")
        write.table(Resume_BH,fileout,row.names=FALSE,col.names=!append, sep="\t", quote = FALSE,append=append)
    }
}


