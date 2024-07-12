#Author : Gabriel Gautier
#date : 10.06.2024
# transform ANR compilation into a dicoexpress readable TARGET file
################# EXCEL TO TARGET FILE ##########################################

create_target_file <- function(Data_Directory, Project_Name, input_file, sample_list_file) {
  # Charger les bibliothèques
  library(readxl)
  library(dplyr)
  
  # Fonction pour extraire le label
  extract_label <- function(id, sample_list) {
    label <- sample_list[grep(paste0("RD_", id, "_"), sample_list)]
    if (length(label) == 0) {
      return(NA)
    }
    return(label)
  }
  
  # Lire le fichier de labels une seule fois
  sample_list <- readLines(paste0(Data_Directory, "/", sample_list_file))

  
  # Charger et transformer les données

  data <- read_excel(paste0(Data_Directory, "/", input_file)) %>%
    select(ID, Population, Condition, Name) %>%
    mutate(labels = sapply(as.character(ID), extract_label, sample_list = sample_list)) %>%
    filter(!is.na(labels)) %>%
    mutate(Replicate = gsub(".*_(\\d+)$", "R\\1", Name)) %>%
    select(-Name) %>% 
    select(-ID)
  
  data <- data[, c("labels", setdiff(names(data), "labels"))]

  
  # Afficher les premières lignes du nouveau fichier
  cat("Voici les premières lignes du nouveau fichier :\n")
  print(head(data))
  
  # Boucle pour obtenir une réponse valide de l'utilisateur
  repeat {
    confirmation <- readline("Voulez-vous continuer et sauvegarder le fichier ? (tapez 'yes' pour confirmer) : ")
    if (tolower(confirmation) %in% c("yes", "no")) {
      break
    } else {
      cat("Réponse invalide. Veuillez répondre par 'yes' ou 'no'.\n")
    }
  }
  
  
  # Si l'utilisateur a répondu "yes" pour continuer
  if (tolower(confirmation) == "yes") {
    # Sauvegarder les données dans un fichier CSV
    output_file <- paste0(Data_Directory, "/", Project_Name, "_TARGET.csv")
    write.csv(data, file = output_file, row.names = FALSE, quote = FALSE)
    cat("Le fichier a été sauvegardé avec succès : ", output_file, "\n")
  } else {
    cat("Opération annulée. Le fichier n'a pas été sauvegardé.\n")
  }
}


