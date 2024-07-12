# Author : adapted by Gabriel Gautier from Ronan Dadole
# Date : 17/10/2023 adapted 18.06.2024
# GO enrichment analysis EXP3

#################### GSEA and EnrichGO #########################################

################################################################################
#install.packages('GOplot')
#install.packages("./Data/org.MdomesticaGDDH13.eg.db", repos=NULL, type="sources")
library(GOplot)
require(org.MdomesticaGDDH13.eg.db)
library(GOSemSim)
library(enrichplot)
library(ggpubr)
library(stringr)
library(ggupset)
library(tidyverse)
library(patchwork)
library(clusterProfiler)

universe <- read.csv("./Data/universe.txt",sep = ",", header = F)[[1]]
semdata <-godata(  OrgDb = org.MdomesticaGDDH13.eg.db,
                   keytype = "GID",
                   ont = "BP"
)
################################################################################
# Gene  contrast strict
################################################################################

env_strict <- read.csv("./Results/EXP3/Contrast_Comparison/Union_20_to_23/Union_Summary_Table.txt",
                       header = TRUE, sep = "\t") %>%
  filter(DE_Group %in% c("B"))

env_strict <- env_strict[[1]]


# env_strict <- unique(c(env_strict, env_2))

# env_strict <- read.csv("./Results/EXP3/DiffAnalysis/[M. domestica dessert_Ctrl-M. domestica dessert_Inf]/Id_DEG.txt", header = F,sep = "\t" )[[1]]

pdf(file = "./Results/EXP3/Go/All_Des.pdf")
ego_env_strict <- enrichGO(gene          = env_strict,
                           universe      = universe,
                           OrgDb         = org.MdomesticaGDDH13.eg.db,
                           ont           = "BP",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE,
                           keyType = "GID")

results <- ego_env_strict@result
p1 <- dotplot(ego_env_strict, showCategory = 15)
p1
dev.off()

################################################################################
# Gene  Set Enrichment Analysis
###############################################################################
source("./Sources_v2/get_logFC.R")


# Correspondance entre les annotations et les contrastes
contrast_files <- list(
  "A" = "[M. domestica cider_Ctrl-M. domestica cider_Inf]",
  "B" = "[M. domestica dessert_Ctrl-M. domestica dessert_Inf]",
  "C" = "[M. sieversii_Ctrl-M. sieversii_Inf]",
  "D" = "[M. sylvestris_Ctrl-M. sylvestris_Inf]"
)

contrast_comparison_path <- "./Results/EXP3/Contrast_Comparison/Union_20_to_23/Union_Summary_Table.txt"

filtre = c("A","B","AB")

data <- get_logFC(filtre,contrast_files,contrast_comparison_path)
# 
# data <- read.csv("./Results/EXP3/DiffAnalysis/[M. domestica dessert_Ctrl-M. domestica dessert_Inf]/DEG.BH.txt", header = T, sep = "\t") %>% 
#   select(Gene_ID,logFC) %>% 
#   arrange(desc(logFC))

geneList = data[,2]
names(geneList) = as.character(data[,1])
geneList = sort(geneList, decreasing = TRUE)

ego_env_strict <- gseGO(
  geneList,
  ont = "BP",
  OrgDb = org.MdomesticaGDDH13.eg.db,
  keyType = "GID",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)

results <- ego_env_strict@result
p1 <- dotplot(ego_env_strict, showCategory = 15)

p <- p3 + p1 + p2 + plot_annotation(tag_levels = 'A')
ggsave("./Results/EXP3/Go/panel3.pdf", p,    width  = 20,
       height = 7)
