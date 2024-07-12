################# Infected notDEG Vienn diagramm ######################################
library(dplyr)
library(UpSetR)
library(RColorBrewer)
library(ggplot2)
library(ggplotify)

#new directory

dir.create("./Results/EXP3/Contrast_Comparison/new_infected_genepop")

Gene_pop <- read.delim("./Results/EXP3/Contrast_Comparison/infected_comparison/Union_Summary_Table.txt", header=T)

# Charger les fichiers
compare_table <- read.csv("./Results/EXP3/DiffAnalysis/Compare_table.txt", header = T,sep = "\t")[c(1,15:20)] %>% 
  filter(Gene_ID %in% Gene_pop$Gene_ID)

colnames(compare_table) <- c("Gene_ID","A-B","A-C","A-D","B-C","B-D","C-D")
  
# Create a data frame to store the results
all_genes <- compare_table$Gene_ID
results <- data.frame(gene = all_genes, group = "", stringsAsFactors = FALSE)

################# CONTRAST COMPARISON ON INFECTED TREES ########################

# Function to assign groups
assign_group <- function(row, gene_id) {
  groups <- c()
  for (col_name in names(row)[-1]) {
    if (row[col_name] == 1) {
      group <- gsub("-.*", "", col_name)  # Extract part before "-"
      groups <- c(groups, group)
    } else if (row[col_name] == -1) {
      group <- gsub(".*-", "", col_name)  # Extract part after "-"
      groups <- c(groups, paste0("", group))
    }
  }
  
  if (length(groups) == 0) {
    # Concatenate all unique letters without duplicates
    all_letters <- unique(unlist(strsplit(gsub("-+", "", names(row)[-1]), "")))
    groups <- paste(sort(all_letters), collapse = "")
  }
  
  # Remove duplicates
  groups <- paste(unique(unlist(strsplit(groups, ""))), collapse = "")
  
  results[results$gene == gene_id, "group"] <<- groups
}

# Apply the function to each row of the table
for (i in 1:nrow(compare_table)) {
  assign_group(compare_table[i, ], compare_table[i, "Gene_ID"])
}

colnames(results) <- c("Gene_ID","DE_Group")

# Fonction pour trier les lettres d'une chaîne par ordre alphabétique
sort_letters <- function(group) {
  paste(sort(strsplit(group, NULL)[[1]]), collapse = "")
}

# Appliquer la fonction pour trier les lettres des groupes dans results
results$DE_Group <- sapply(results$DE_Group, sort_letters)

########################## FILTERING DEG GENES #################################

DEG <- read.delim("./Results/EXP3/DiffAnalysis/[M. sylvestris_Ctrl-M. sylvestris_Inf]/Id_DEG.txt", header=F) %>% 
  mutate(DE_Group ="D") 

DEG_2 <- read.delim("./Results/EXP3/DiffAnalysis/[M. sieversii_Ctrl-M. sieversii_Inf]/Id_DEG.txt", header=F) %>% 
  mutate(DE_Group ="C") 
DEG_3 <-read.delim("./Results/EXP3/DiffAnalysis/[M. domestica dessert_Ctrl-M. domestica dessert_Inf]/Id_DEG.txt", header=F) %>% 
  mutate(DE_Group ="B") 
DEG_4 <-read.delim("./Results/EXP3/DiffAnalysis/[M. domestica cider_Ctrl-M. domestica cider_Inf]/Id_DEG.txt", header=F) %>% 
  mutate(DE_Group ="A") 

DEG <- rbind(DEG,DEG_2,DEG_3,DEG_4)


filter_DEG_by_population <- function(results_df, DEG_df, population_letter) {
  # Filtrer les gènes DEG pour la population spécifique
  DEG_filtered <- DEG_df %>%
    filter(DE_Group == population_letter)
  
  # Filtrer les groupes contenant la lettre spécifique
  results_with_letter <- results_df %>%
    filter(grepl(population_letter, DE_Group))
  
  # Filtrer les gènes de ces groupes pour exclure ceux présents dans DEG_filtered
  results_filtered <- results_with_letter %>%
    filter(!Gene_ID %in% DEG_filtered$V1)
  
  # Retourner les résultats filtrés et les résultats sans la lettre spécifique
  results_no_letter <- results_df %>%
    filter(!grepl(population_letter, DE_Group))
  
  return(rbind(results_no_letter, results_filtered))
}


# Initialiser un dataframe pour stocker les résultats filtrés
results_filtered_all <- results

# Lettres représentant les différentes populations
Populations <- c("A", "B", "C", "D")

# Itérer sur chaque population et appliquer le filtre
for (pop in Populations) {
  results_filtered_all <- filter_DEG_by_population(results_filtered_all, DEG, pop)
}

######################## PLOT UPSETR ###########################################

# Obtenir la liste des groupes uniques
Groups <- unique(unlist(strsplit(paste(results_filtered_all$DE_Group, collapse=""), "")))

# Préparer une liste des gènes pour chaque groupe
Compare <- sapply(Groups, function(g) {
  results_filtered_all$Gene_ID[grepl(g, results_filtered_all$DE_Group)]
}, simplify = FALSE)
names(Compare) <- Groups

# Utiliser les noms de groupes existants directement
listInput <- Compare

# Palette de couleurs
colorPalette <- brewer.pal(length(Groups), "Set2")

# Créer les requêtes pour les intersections
tmp <- lapply(seq_along(Groups), function(i) {
  list(query = intersects, params = list(Groups[i]), color = colorPalette[i], active = TRUE)
})

# Créer le diagramme UpSetR
upsetPlot <- upset(
  fromList(listInput),
  point.size = 1.5,
  sets = rev(Groups),
  keep.order = TRUE,
  sets.x.label = "Contrasts",
  line.size = 0.5,
  text.scale = c(0.85, 1, 1, 1, 1, 0.75),
  order.by = "freq",
  queries = tmp
)

# Créer le texte de la légende
legend_text <- paste(Groups, collapse = "\n")

# Convertir le diagramme en objet ggplot
ggplot_upset <- ggplotify::as.ggplot(upsetPlot)

# Ajouter la légende au diagramme
ggplot_upset <- ggplot_upset + 
  annotate("text", x = Inf, y = Inf, label = legend_text, hjust = 1.1, vjust = 1.1, size = 3)

# Enregistrer le fichier Presults
ggsave("./Results/EXP3/Contrast_Comparison/new_infected_genepop/UpSetR.pdf", plot = ggplot_upset, width = 8, height = 6)

################### PLOT VIENN DIAGRAMM ########################################