#Author : Gabriel Gautier
#date : 19.06.2024
#input : sample list,table of all traits
# output : table of phenotypic traits for WGCNA analysis
#################### CREATE PHENOTYPIC FILE ####################################
library(readxl)
library(dplyr)


# Lire le fichier de labels
sample_list <- readLines("./Data/sample_list_2.txt")


# Charger et transformer les données
extract_label <- function(id, sample_list) {
  label <- sample_list[grep(paste0("RD_", id, "_"), sample_list)]
  if (length(label) == 0) {
    return(NA)
  }
  return(label)
}


data <- read_excel("./Data/Expe3_ANR_Compilation_081123 (1).xlsx") %>%
  mutate(labels = sapply(as.character(ID), extract_label, sample_list = sample_list)) %>%
  filter(!is.na(labels)) %>%
  mutate(Replicate = gsub(".*_(\\d+)$", "R\\1", Name))

data <- data[, c("labels", setdiff(names(data), "labels"))] %>% 
  select(labels,`Stade (start)`, `Nb leaf (start)`, `Chl (start)`, `Chl_sd (start)`, `Flav (start)`,`Flav_sd (start)` , `NBI (start)`, `NBI_sd (start)`, 
         `aphid d2 (mother)`,`aphid d2 (larvea)`,
         `Stade (d4)`, `aphid d4 (mother)`, `aphid d4 (larvea)`, 
         `Stade (d8)`,`aphid d8 (all)`, 
         `Chl (d13)`, `Chl_sd (d13)`, `Flav (d13)`, `Flav_sd (d13)`, `NBI (d13)`,`NBI_sd (d13)`,
         `1_winged_Adult`, `2_apterous_Adult`, `3_nymph_(N4)`, `4_larvae_(L4)`,`5_larvae_(L1_L3)_nymph_(N1_N3)`, `6_molt`,Population,Condition,)

data <- data %>% 
  filter(!labels == "RD_106_B")

write.csv2(data, file= "./Data/Phenotypic_traits.csv",row.names = FALSE, quote = FALSE)
