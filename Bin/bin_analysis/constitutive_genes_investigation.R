######################### CONSTITUTIVE GENE ANALYSIS ###########################
library(tidyverse)
library(UpSetR)
library(RColorBrewer)
DEG <- read.delim("./Results/EXP3/DiffAnalysis/[M. sieversii_Ctrl-M. sieversii_Inf]/Id_DEG.txt", header=F)

genes_constitutifs <- read.delim("./Results/EXP3/Contrast_Comparison/new_infected_genepop/Union_Summary_Table.txt",sep = ";", header=T) %>% 
  filter(DE_Group %in% c("CCC","CC","C")) %>% 
  select(Gene_ID,DE_Group) %>% 
  filter(!Gene_ID %in% DEG$V1)

A <- read.delim("./Results/EXP3/DiffAnalysis/[M. domestica cider_Ctrl-M. domestica cider_Inf]/Id_DEG.txt", header=FALSE) %>% 
  mutate(DE_Group = "A")
D <- read.delim("./Results/EXP3/DiffAnalysis/[M. sylvestris_Ctrl-M. sylvestris_Inf]/Id_DEG.txt", header=FALSE) %>% 
  mutate(DE_Group = "D")
B <- read.delim("./Results/EXP3/DiffAnalysis/[M. domestica dessert_Ctrl-M. domestica dessert_Inf]/Id_DEG.txt", header=FALSE) %>% 
  mutate(DE_Group = "B")

colnames(A)[1] <- "Gene_ID"
colnames(D)[1] <- "Gene_ID"
colnames(B)[1] <- "Gene_ID"

# Créer une liste nommée pour UpSetR
Table <- rbind(A,B,D) %>% 
  filter(Gene_ID %in% genes_constitutifs$Gene_ID)

C <- genes_constitutifs %>% 
  filter(!Gene_ID %in% Table$Gene_ID) %>%
  mutate(DE_Group = "Not DEG in any population")
  
Groups <- unique(unlist(strsplit(paste(Table$DE_Group, collapse=""), "")))

# Préparer une liste des gènes pour chaque groupe
Compare <- sapply(Groups, function(g) {
  Table$Gene_ID[grepl(g, Table$DE_Group)]
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
ggsave("./Results/EXP3/Contrast_Comparison/new_infected_genepop/UpSetR_not_DEG.pdf", plot = ggplot_upset, width = 8, height = 6)
