Load_Data_Files <- function(Data_Directory, Project_Name, Filter=NULL, Sep=Sep)
{
  ## Load Annotation file
  Annotation_File <- paste0(Data_Directory,"/",Project_Name,"_Annotation.csv")
  if (file.exists(Annotation_File)){
    Annotation <- read.csv2(Annotation_File, header = TRUE, sep=Sep, quote = "")
    Annotation <- Annotation[,c("Gene_ID","Gene_Name")]
  }
  else
  { Annotation = NULL }
  
  ## Load GO file
  Enrichment_File <- paste0(Data_Directory,"/",Project_Name,"_Enrichment.csv")
  if (file.exists(Enrichment_File)){
    Reference_Enrichment <- read.csv2(Enrichment_File, header = TRUE, sep=Sep, quote = "")
    Reference_Enrichment <- Reference_Enrichment[,c("Gene_ID","Annotation")]
    Reference_Enrichment <- unique(Reference_Enrichment)
  }
  else
  { Reference_Enrichment = NULL }
  
  ## Load Target table
  Target <- read.table(paste0(Data_Directory,"/",Project_Name,"_TARGET.csv"),header=TRUE,sep=Sep,row.names=1,colClasses = "factor")
  index=numeric()
  for (i in 1:ncol(Target))
  {
    if (nlevels(Target[,i])==1)
      index<-c(index,i) 
  }
  Target[,index]=NULL
  
  ##  Every first letter of each factor is a capital letter 
  for (i in 1:ncol(Target))
  {
    
    s <- strsplit(colnames(Target)[i], " ")[[1]]
    colnames(Target)[i] <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
                                 sep = "", collapse = " ")
  }
  
  ## Load Raw counts table
  Raw_Counts <- read.table(paste0(Data_Directory,"/",Project_Name,"_COUNTS.csv"),header=TRUE,sep=Sep,row.names=1)
  
  ## Control : labels of samples
  suppressMessages({
    if(!identical(colnames(Raw_Counts),rownames(Target))) {
      cat("Samples are not organized or named in the same manner in the expression file and the target file. The expression file is reorganized according to the target table \n")
      Raw_Counts<-Raw_Counts[,rownames(Target)]
    }
  })

  
  ## Filter : once the two files are organized in the same way
  if(!is.null(Filter))
  {
    Project_Name <- Filter[[length(Filter)]]
    Filter <- Filter[-length(Filter)]  
    
    # Capitalize the first letter of each factor.
    Filter <- sapply(Filter, function(x) {
      x[1] <- paste0(toupper(substring(x[1], 1, 1)), substring(x[1], 2))
      return(x)
    }, simplify = FALSE)
    
    
    unique_factors <- unique(sapply(Filter, `[`, 1))
    new_Filter <- list()
    
    #  Loop over each factor in the Filter list
    for (factor in unique_factors) {
      new_list_values <- list()
      factor = paste(toupper(substring(factor, 1, 1)),
                     substring(factor, 2),sep = "", collapse = " ")
      unique_values = names(table(Target[[factor]]))
      # Loop over each modality value
      for (modality in unique_values){
        #  Index of modality if present in Filter
        index<-which(sapply(Filter, function(x) x[2]) == modality)
        
        # If modality is present
        if (length(index)!=0){
          # If its value is FALSE, retrieve it to create the new Filter list
          if(unlist(Filter[index])[3]==FALSE){
            new_list_values <- c(new_list_values,modality)
          }
        }
        # Or if the modality is not present at all
        
        
        if(all(sapply(Filter, function(x) x[[2]] != modality))){
          new_list_values <- c(new_list_values,modality)
        }
        
        
      }
      new_list_values <- unlist(new_list_values)
      
      filtered_list <- Filter[sapply(Filter, function(x) x[1] == factor)]
      filter_bool <- sapply(filtered_list, function(x) x[[3]])
      all_false <- all(filter_bool == FALSE) 
      
      
      for (i in 1:length(new_list_values)){
        if(all_false){
          tmp <- c(factor,new_list_values[i],TRUE)
        }
        else{
          tmp <- c(factor,new_list_values[i],FALSE)
        }
        new_Filter <- c(new_Filter,list(tmp))
      }
    }
    
    Filter <- c(Filter,new_Filter,Project_Name)
    
    # Drop duplicates term in Filter
    tmp <- sapply(Filter, function(x) paste(x, collapse = "-----"))
    unique_term <- unique(tmp)
    Filter <- lapply(unique_term, function(x) strsplit(x, "-----")[[1]])
    
    
    # print(Filter)
    # CODE ML
    nb.rules=length(Filter)-1
    for (i in 1:nb.rules)
    {
      x=Filter[[i]]
      x[1]=paste(toupper(substring(x[1], 1, 1)),
                 substring(x[1], 2),sep = "", collapse = " ")
      # print(x)
      if(x[1]=="Sample")
      {
        index<-grep(x[2],rownames(Target))
        Raw_Counts=Raw_Counts[,-index]
        Target<-Target[-index,]
      }
      else
      {
        index.col<-grep(x[1],colnames(Target))
        # index<-grep(x[2],Target[,index.col])
        index <- which(Target[, index.col] == x[2])
        
        if(x[3]==FALSE)
        {
          #   Raw_Counts=Raw_Counts[,index]
          #   Target<-Target[index,-index.col]
          # }
          # else
          # {
          Raw_Counts=Raw_Counts[,-index]
          Target<-Target[-index,]
          Target[,index.col]=droplevels(Target[,index.col])
        }
      }
      Project_Name=Filter[[nb.rules+1]]
    }
    for (i in 1:ncol(Target))
      Target[,i]=droplevels(Target[,i])
  }
  
  
  # Remove factor with 1 modality
  uniques <- apply(Target, 2, function(x) length(unique(x)))
  cols_to_remove <- which(uniques == 1)
  indice <- which(cols_to_remove == "Replicate")
  if(length(indice)>0){
    cols_to_remove <- cols_to_remove[-indice]
  }
  if (length(cols_to_remove) > 0) {
    Target <- Target[,-cols_to_remove]
  }
  
  
  ## Number of biological factors
  NbBioFactors <- ncol(Target)-1
  Design_Ok=TRUE
  ## Control: Complete design
  if(sum(table(Target[,1:NbBioFactors]))==0)
  {
    cat("DiCoExpress is not suitable for this experimental design which is not complete \n")
    print(table(Target[,1:NbBioFactors]))
    Design_Ok=FALSE
    Data_Files=NULL
  }
  
  ## Control: only two biological factors and one Replicate factor
  if((colnames(Target)[ncol(Target)]!="Replicate"))
  {
    cat("Target is badly defined, the last column must be Replicate \n")
    Design_Ok=FALSE
    Data_Files=NULL
  }
  if( NbBioFactors>2)
  {
    cat("DiCoExpress is not suitable for this experimental design.\n
It must contain only 2 biological factors and the factor Replicate. \n
However you can perform the quality control step \n")
    Design_Ok=TRUE
  }
  ## Create Project Name directory
  if (!I( Project_Name %in% dir(paste0(Results_Directory)))) dir.create(paste0(Results_Directory,"/", Project_Name), showWarnings=FALSE)
  
  if(Design_Ok)
  {
  #   if(!is.null(Filter))
  #   {
  #     reducedTarget <-data.frame(Target,name=apply(Target,1,function(x) paste(x,collapse="_")))
  #     fileout=paste0(Results_Directory,"/", Project_Name,"/",Project_Name,"_reduced_TARGET.txt")
  #     write.table(reducedTarget,fileout,sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
  #   }
  #   rownames(Target)=apply(Target,1,function(x) paste(x,collapse="_"))
  #   colnames(Raw_Counts)=rownames(Target)
  }   
  
  Data_Files <- list("Target" = Target, "Raw_Counts" = Raw_Counts, "Project_Name"= Project_Name, "Annotation"= Annotation, "Reference_Enrichment"=Reference_Enrichment)
  
  invisible(Data_Files)
}