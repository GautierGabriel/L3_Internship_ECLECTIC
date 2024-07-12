###### Create Venn diagram ######

Contrast_Comparison <- function(Data_Directory, Results_Directory, Project_Name, 
                                Title,  Groups, Operation)
{    
    ## subdirectory "Project"
    if (!I(Project_Name %in% dir(Results_Directory))){
        cat("This project does not exist\n")
        break
    }
    ## subdirectory
    subdirectory <- "Contrast_Comparison"
    path=paste0(Results_Directory,"/",Project_Name,"/",subdirectory)
    if (!I(subdirectory %in% dir(paste0(Results_Directory,"/",Project_Name)))) dir.create(path, showWarnings=FALSE)
    if (!I(Title %in% dir(paste0(Results_Directory,"/",Project_Name,"/",subdirectory)))) dir.create(paste0(path,"/",Title), showWarnings=FALSE)
    ## Load Compare table
    Compare <- read.csv2(paste0(Results_Directory,"/",Project_Name,"/DiffAnalysis/Compare_table.txt"), header=TRUE, sep="\t",check.names=FALSE)
 
    ## check if the groups are really studied contrasts
    if (sum(Groups%in%colnames(Compare[,Groups]))!=length(Groups))
    {
        cat("Groups is badly defined. The constrast(s) \n")
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
        fileout <- paste0(path,"/",Title,"/","Union_List.txt")
        write.table(IdGene,fileout,row.names=FALSE,col.names=FALSE, sep="\t", quote = FALSE)

        Compare_subset <- Compare[is.element(Compare$Gene_ID,IdGene),Groups]
        Summary_Union <- data.frame("Gene_ID"=IdGene,DE_Group=apply(Compare_subset,1,function(a) paste(LETTERS[which(a!=0)],collapse='')))
        
        # Determine Up or Down regulation for each gene in the Intersection list
        Regulation <- apply(Compare_subset, 1, function(a) {
          if (any(a > 0) && !any(a < 0)) {
            return("Up")
          } else if (any(a < 0) && !any(a > 0)) {
            return("Down")
          } else {
            return("Mixed")
          }
        })
        
        # Create a dataframe with Gene_ID and Regulation
        Regulation_DF <- data.frame(Gene_ID = IdGene, Regulation = Regulation)
        
        #merge
        Summary_Union <- merge(Regulation_DF, Summary_Union, by="Gene_ID",all.y=TRUE)
        ## Add Gene_Name
        if(!is.null(Annotation))
            Summary_Union <- merge(Annotation, Summary_Union, by="Gene_ID",all.y=TRUE)
        fileout <- paste0(path,"/",Title,"/","Union_Summary_Table.txt")
        write.table(Summary_Union,fileout,row.names=FALSE,col.names=TRUE, sep="\t", quote = FALSE)
    }

    ## Intersection List
    if (Operation == "Intersection")
    {
        newCol <- ncol(Compare)+1
        Compare[,newCol]=apply(Compare[,Groups],1,function(a) 1*(sum(abs(a))==length(Groups)))
        colnames(Compare)[newCol]<- paste0(Title,"_Intersection")
        IdGene <- Compare$Gene_ID[which(Compare[,newCol]==1)]
        ## Intersection List
        fileout <- paste0(path,"/",Title,"/","Intersection_List.txt")
        write.table(IdGene,fileout,row.names=FALSE,col.names=FALSE, sep="\t", quote = FALSE)

        Compare_subset <- Compare[is.element(Compare$Gene_ID,IdGene),Groups]
        Summary_Intersection <- data.frame("Gene_ID"=IdGene,DE_Group=apply(Compare_subset,1,function(a) paste(LETTERS[which(a!=0)],collapse='')))
        
        # Determine Up or Down regulation for each gene in the Intersection list
        Regulation <- apply(Compare_subset, 1, function(a) {
          if (any(a > 0) && !any(a < 0)) {
            return("Up")
          } else if (any(a < 0) && !any(a > 0)) {
            return("Down")
          } else {
            return("Mixed")
          }
        })
        
        # Create a dataframe with Gene_ID and Regulation
        Regulation_DF <- data.frame(Gene_ID = IdGene, Regulation = Regulation)
        
        #merge
        Summary_Intersection <- merge(Regulation_DF, Summary_Intersection, by="Gene_ID",all.y=TRUE)
        
        ## Add Gene_Name
        if(!is.null(Annotation))
          Summary_Intersection <- merge(Annotation, Summary_Intersection, by="Gene_ID",all.y=TRUE)
        fileout <- paste0(path,"/",Title,"/","Intersection_Summary_Table.txt")
        write.table(Summary_Intersection,fileout,row.names=FALSE,col.names=TRUE, sep="\t", quote = FALSE)
    }    
    
        ## Add a column to the compare table
    if (file.exists(paste0(path,"/","Contrast_Comparison_Output_table.txt")))
    {
      tmp=read.table(paste0(path,"/","Contrast_Comparison_Output_table.txt"),header=TRUE,sep="\t",check.names = FALSE)
      tmp=merge(Compare,tmp)    
    } 
    fileout <-  paste0(path,"/","Contrast_Comparison_Output_table.txt")
    write.table(Compare,fileout,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

    if(Operation != "Union" && Operation != "Intersection")
    {
        cat("The Operation is misinformed. The possible values for the Operation argument are Union or Intersection")
    }

    ## check of the number of groups that are compared
    if (length(Groups) >=6)
    {
        cat("Cannot plot Venn diagram for more than 5 sets \n")
    }
    else
    {
    ## Create Venn diagram

    ## DEG ##
    pdf(paste0(path,"/",Title,"/","UpSetR2.pdf"),width = 7,height = 5)
    Venn <- vennDiagram(Compare[,Groups], circle.col = brewer.pal(length(Groups),"Set2"),
                        cex=1, names=LETTERS[1:length(Groups)],main="DEG")
    for (a in 1:length(Groups))
      title(paste0(LETTERS[a], " = ",Groups[a]), line = -22-a, adj=0.1 ,cex.main = 0.6)

    shortNames <- LETTERS[1:length(Groups)]
    names(shortNames) <- Groups
    legend_text <- paste(shortNames, ": ", Groups, collapse = "\n")
    legend("topleft", legend = strsplit(legend_text, "\n")[[1]],
           xpd = TRUE, title = "",cex=0.4, bty="n")

    ## UP ##
    Venn <- vennDiagram(Compare[,Groups],include="up", circle.col = brewer.pal(length(Groups),"Set2"),
                        cex=1, names=LETTERS[1:length(Groups)],main="UP")
    for (a in 1:length(Groups))
      title(paste0(LETTERS[a], " = ",Groups[a]), line = -22-a, adj=0.1 ,cex.main = 0.6)
    legend("topleft", legend = strsplit(legend_text, "\n")[[1]],
           xpd = TRUE, title = "",cex=0.4, bty="n")

    ## DOWN ##
    Venn <- vennDiagram(Compare[,Groups],include="down", circle.col = brewer.pal(length(Groups),"Set2"),
                        cex=1, names=LETTERS[1:length(Groups)],main="DOWN")
    for (a in 1:length(Groups))
      title(paste0(LETTERS[a], " = ",Groups[a]), line = -22-a, adj=0.1 ,cex.main = 0.6)
    legend("topleft", legend = strsplit(legend_text, "\n")[[1]],
           xpd = TRUE, title = "",cex=0.4, bty="n")
      
      ## Create UpsetR diagramm
      
      shortNames <- LETTERS[1:length(Groups)]
      names(shortNames) <- Groups
      listInput <- lapply(Compare[, Groups], function(x) Compare$Gene_ID[which(abs(x) == 1)])
      names(listInput) <- shortNames
      
      colorPalette <- brewer.pal(length(Groups), "Set2")
      tmp <- lapply(seq_along(shortNames), function(i) {
        list(query = intersects, params = list(shortNames[i]), color = colorPalette[i], active = T)
      })

      upsetPlot <- upset(fromList(listInput), point.size = 1.5,sets = rev(shortNames),keep.order=TRUE,
            sets.x.label = "Contrasts",line.size = 0.5,text.scale = c(0.85, 1, 1, 1, 1, 0.75),
            order.by = "freq",queries = tmp) ##

      legend_text <- paste(shortNames, ": ", Groups, collapse = "\n")
      
      
      ggplot_upset <- ggplotify::as.ggplot(upsetPlot)
      ggplot_upset <- ggplot_upset + 
        annotate("text", x = Inf, y = Inf, label = legend_text, hjust = 1.1, vjust = 1.1, size = 3)
      
      ggsave(paste0(path,"/",Title,"/","UpSetR.pdf"),
             plot = ggplot_upset, width = 8, height = 6)
      dev.off()
    }
        
    pdf_combine(c(paste0(path,"/",Title,"/","UpSetR2.pdf")
              ,paste0(path,"/",Title,"/","UpSetR.pdf")), 
              output = paste0(path,"/",Title,"/","Venn_Diagram.pdf"))
        
    # Drop file "_Venn_Diagram.pdf" and "_UpSetR.pdf"
    file.remove(paste0(path,"/",Title,"/","UpSetR2.pdf"))
    file.remove(paste0(path,"/",Title,"/","UpSetR.pdf"))
}
