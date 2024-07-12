#Author : Gabriel Gautier
#date : 10.06.2024
# transform Salmon counts into a dicoexpress readable COUNTS file


################## SALMON OUTPUT TOWARD CSV FORMATED FOR DICOEXPRESS ###########

clean_row_names <- function(Data_Directory, Project_Name, input_file) {
  # Importer le fichier avec les données
  data <- read.csv(paste0(Data_Directory,"/",input_file), header = TRUE, sep = "\t")
  
  # Nettoyer les noms de lignes
  data$Name <- gsub("^(mRNA:|ncRNA:)", "", data$Name) 
  
  # Afficher les premières lignes du nouveau fichier
  cat("Voici les premières lignes du nouveau fichier :\n")
  print(head(data))
  
  # Demander à l'utilisateur de confirmer avant de continuer
  confirmation <- readline("Voulez-vous continuer et sauvegarder le fichier ? (tapez 'yes' pour confirmer): ")
  
  # Vérifier si l'utilisateur a répondu "yes" pour continuer
  if (tolower(confirmation) == "yes") {
    # Sauvegarder les données dans un fichier CSV sans inclure les numéros de ligne et sans guillemets
    write.csv(data, file = paste0(Data_Directory,"/",Project_Name,"_COUNTS.csv"), row.names = FALSE, quote = FALSE)
    cat("Le fichier a été sauvegardé avec succès.\n")
  }
  else {
    cat("Opération annulée. Le fichier n'a pas été sauvegardé.\n")
  }
}
