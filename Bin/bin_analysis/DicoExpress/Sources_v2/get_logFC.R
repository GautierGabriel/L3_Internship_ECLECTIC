#Author : Gabriel Gautier
#date : 19.06.2024
# create a file with logFC and gene ID for GSEA analysis, with the filter chosen

################# GET LOGFC FOR GSEA ANALYSIS ##################################

get_logFC <- function(filtre,contrast_files,contrast_comparison_path) {
  
  library(dplyr)
  
  annotations <- paste(filtre, collapse = "")
  
  # Séparer les lettres individuelles
  annotations <- unlist(strsplit(annotations, ""))
  
  # Retourner les lettres uniques
  annotations <- unique(annotations)
  
  if (!all(annotations %in% c("A", "B", "C", "D","E"))) {
    stop("Les annotations doivent être parmi 'A', 'B', 'C', 'D','E' ")
  }
  
  # Lire le fichier d'intersection pour obtenir les gènes correspondants à l'annotation
  results_df <-  read.csv(contrast_comparison_path,
                         header = TRUE,
                         sep = "\t") %>% 
    filter(DE_Group %in% filtre) %>% 
    dplyr::select(Gene_ID)
  
  # Lire les valeurs de logFC pour les contrastes correspondants à chaque annotation
  for (annot in annotations) {
    logFC_data <- read.csv(paste0("./Results/EXP3/DiffAnalysis/",contrast_files[[annot]],"/DEG.BH.txt"), header = T, sep = "\t") %>% 
      dplyr::select(Gene_ID,logFC)
    colnames(logFC_data)[2] <- paste0("logFC_",annot)
    
    results_df <- results_df %>%
      left_join(logFC_data, by = "Gene_ID")
  }
  results_df <- results_df %>% 
    mutate(logFC = rowMeans(abs(dplyr::select(., starts_with("logFC_"))), na.rm = TRUE)) %>% 
    dplyr::select(Gene_ID,logFC)
}
