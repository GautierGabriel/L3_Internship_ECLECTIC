####################################################################
###                     Parameter Info                             ###
####################################################################
Save_Parameters=function()
{
Output = file(paste0(Results_Directory,"/",Project_Name,"/Parameter_Information.txt"), open="wt")
sink(Output)
sink(Output, type="message")
cat(date(),"\n")
cat("Project_Name:",Project_Name,"\n")
if(!is.null(Filter))
cat("Filter: ", Filter,"\n")
else
    cat("Filter: NULL \n")
cat("Filter_Strategy: filterByExpr of edger","\n")
cat("Min.count: ",min.count,"\n")
cat("Normalization_Method: ",Normalization_Method,"\n")
cat("Replicate: ", Replicate,"\n")
cat("Interaction: ", Interaction,"\n")
cat("NbGenes_Profiles: ",NbGenes_Profiles,"\n") 
cat("NbGenes_Clustering: ",NbGenes_Clustering,"\n") 
cat("Alpha_DiffAnalysis: ",Alpha_DiffAnalysis,"\n")
cat("Alpha_Enrichment: ",Alpha_Enrichment,"\n")
cat("\n")
cat("Groups (Contrats_Comparison): ",Groups_Contrast_Comparison,"\n")
cat("Groups (Coexpression): ",Groups_Coexpression,"\n")


sink(type="message")
sink()
####################################################################
###                     Session Info                             ###
####################################################################
writeLines(capture.output(sessionInfo()),
           paste0(Results_Directory,"/",Project_Name,"/SessionInfo.txt"))
}
