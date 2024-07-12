# Author : adapted by Gabriel Gautier from Ronan Dadole
# Date : 2023 adapted 20.06.2024
# Get cluster DEG EXP3

#################### EnrichGO and Cluster/DEG inters ###########################
library(GOplot)
require(org.MdomesticaGDDH13.eg.db)
library(GOSemSim)
library(enrichplot)
library(ggpubr)
library(stringr)
library(ggupset)
library(clusterProfiler)
library(tidyverse)
library(patchwork)

wd_modules = "./Results/EXP3/WGCNA/modules_filt/"
module_file = list.files(wd_modules, pattern = regex("module"))

# DEG <- read.delim("./Results/EXP3/DiffAnalysis/[M. sieversii_Ctrl-M. sieversii_Inf]/Id_DEG.txt", header=F)

Gene_pop <- read.delim("./Results/EXP3/Contrast_Comparison/Union_20_to_23/Union_Summary_Table.txt",sep = "\t", header=T) %>% 
  filter(DE_Group %in% c("C")) %>% 
  select(Gene_ID,DE_Group,Regulation) #%>% 
  # filter(!Gene_ID %in% DEG$V1)




name <- c()
perc_pop <- c()
nb_tot <- c()
nb_pop <- c()

for (i in seq(1, length(module_file))){
  file_name = module_file[i]
  module_color = substr(file_name,8,nchar(file_name)-4)
  name_variable = paste0("Genes_",module_color)
  assign(name_variable, read.csv(file.path(wd_modules,file_name), header=FALSE)[[1]])
  assign(paste0("it_pop_",module_color),intersect(get(name_variable),Gene_pop[[1]]))
  nb_tot <- append(nb_tot,length(get(name_variable)))
  nb_pop <- append(nb_pop,length(get(paste0("it_pop_",module_color))))
  name <- append(name,str_replace(module_color, "_filt$", ""))
  perc_pop <- append(perc_pop,(nb_pop/nb_tot*100))
}

info <- data.frame(name,nb_tot,nb_pop)
names(info) <- c("Module_Name","Total_Genes_Module","constitutive_resistance_genes")

info <- info %>% 
  mutate( Not_constitutive_resistance_genes = nb_tot - nb_pop) %>% 
  mutate( Percentage_constitutive_resistance_genes = (nb_pop/nb_tot*100)) %>% 
  mutate( Percentage_Not_constitutive_resistance_genes = 100-Percentage_constitutive_resistance_genes)
  

################################################################################
# library(ggplot2)
# library(ggpubr)
# 
# cluster_comparison <- function(colors) {
#   # Charger les bibliothèques nécessaires
#   library(UpSetR)
#   library(RColorBrewer)
#   library(ggplot2)
#   library(ggplotify)
#   
#   # creer le dataframe 
#   
#   selected_genes  <- read.csv("./Results/EXP3/WGCNA/gene_modules_filt.txt", header=T, sep = "\t")
#   colnames(selected_genes) <- c("Gene_ID","color")
#   
#   selected_genes <- selected_genes %>% 
#     left_join(Gene_pop) %>% 
#     filter(!is.na(DE_Group))
#   
#   if(colors != "all"){
#     selected_genes <- selected_genes %>% 
#       filter(color %in% colors)
#   }
#   
#   df <- selected_genes %>%
#     select(-color)
#   
#   # Obtenir la liste des groupes uniques
#   Groups <- unique(unlist(strsplit(paste(df$DE_Group, collapse=""), "")))
#   
#   # Préparer une liste des gènes pour chaque groupe
#   Compare <- sapply(Groups, function(g) {
#     df$Gene_ID[grepl(g, df$DE_Group)]
#   }, simplify = FALSE)
#   names(Compare) <- Groups
#   
#   # Utiliser les noms de groupes existants directement
#   listInput <- Compare
#   
#   # Palette de couleurs
#   colorPalette <- brewer.pal(length(Groups), "Set2")
#   
#   # Créer les requêtes pour les intersections
#   tmp <- lapply(seq_along(Groups), function(i) {
#     list(query = intersects, params = list(Groups[i]), color = colorPalette[i], active = TRUE)
#   })
#   
#   # Créer le diagramme UpSetR
#   upsetPlot <- upset(
#     fromList(listInput),
#     point.size = 1.5,
#     sets = rev(Groups),
#     keep.order = TRUE,
#     sets.x.label = "Contrasts",
#     line.size = 0.5,
#     text.scale = c(0.85, 1, 1, 1, 1, 0.75),
#     order.by = "freq",
#     queries = tmp
#   )
#   
#   # Créer le texte de la légende
#   legend_text <- paste(Groups, collapse = "\n")
#   
#   # Convertir le diagramme en objet ggplot
#   ggplot_upset <- ggplotify::as.ggplot(upsetPlot)
#   
#   # Ajouter la légende au diagramme
#   ggplot_upset <- ggplot_upset + 
#     annotate("text", x = Inf, y = Inf, label = legend_text, hjust = 1.1, vjust = 1.1, size = 3)
#   
#   # Enregistrer le fichier PDF
#   ggsave("./Results/EXP3/Go_WGCNA/UpSetR.pdf", plot = ggplot_upset, width = 8, height = 6)
#   
# }
# 
# # cluster_comparison(colors = "all")
# 
# ggplot(info, aes(x = Module_Name, y = Percentage_constitutive_resistance_genes, fill = Module_Name)) + 
#   geom_col(col = "black") + 
#   ylim(0, 100) + 
#   scale_fill_identity() + 
#   theme_bw() +
#   geom_text(aes(x = Module_Name, y = Percentage_constitutive_resistance_genes, label = round(Percentage_constitutive_resistance_genes, 1)), nudge_y = 3) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
#   labs(y = "Percentage constitutive resistance genes", x = "Gene Module")
# 
# ggsave(                                                                       
#   file   = paste0("./Results/EXP3/Go_WGCNA/Percentage_constitutive_resistance_genes_Module.pdf"),
#   plot   = last_plot(),
#   width  = 7,
#   height = 7
# )
# 
# data_long <- info %>% 
#   pivot_longer(cols = c(constitutive_resistance_genes, Not_constitutive_resistance_genes), names_to = "Genes", values_to = "nb")
# 
# data_long$Genes <- factor(data_long$Genes, levels = c("constitutive_resistance_genes", "Not_constitutive_resistance_genes"))
# 
# ggplot(data_long, aes(x = Module_Name, y = nb, fill = Genes)) + 
#   geom_col(position = "stack", col = "black") + 
#   theme_bw() + 
#   scale_fill_manual(values = c("red", "grey")) + 
#   theme_pubclean() +
#   labs(y = "Number of constitutive resistance genes", x = "Gene Module") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
#   theme(legend.title = element_blank(), legend.text = element_text(size = 12))
# 
# ggsave(                                                                       
#   file   = paste0("./Results/EXP3/Go_WGCNA/Number_DEG_module.pdf"),
#   plot   = last_plot(),
#   width  = 7,
#   height = 7
# )
# 
# data_long_p <- info %>% 
#   pivot_longer(cols = c(Percentage_constitutive_resistance_genes, Percentage_Not_constitutive_resistance_genes), names_to = "Genes", values_to = "nb")
# 
# data_long_p$Genes <- factor(data_long_p$Genes, levels = c("Percentage_constitutive_resistance_genes", "Percentage_Not_constitutive_resistance_genes"))
# 
# ggplot(data_long_p, aes(x = Module_Name, y = nb, fill = Genes)) + 
#   geom_col(position = "stack", col = "black") + 
#   theme_bw() + 
#   scale_fill_manual(values = c("Percentage_constitutive_resistance_genes" = "white", "Percentage_Not_constitutive_resistance_genes" = "grey")) + 
#   theme_pubclean() + 
#   labs(y = "Percentage of constitutive resistance genes", x = "Gene Module") +
#   theme(text = element_text(size = 13),
#         axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
#         axis.ticks.length.x = unit(-0.125, "cm"),
#         legend.title = element_blank(), legend.text = element_text(size = 12))
# 
# ggsave(                                                                       
#   file   = paste0("./Results/EXP3/Go_WGCNA/Percentage_constitutive_resistance_genes_Not_DEG_Module.pdf"),
#   plot   = last_plot(),
#   width  = 7,
#   height = 7
# )
# 

############################# GO enrichment ####################################

#run first the sript to create universe.txt or use your filtered gene list (after quality control filtering and normalization)

universe <- read.csv("./Data/universe.txt",sep = ",", header = F)[[1]]
semdata <-godata(  OrgDb = org.MdomesticaGDDH13.eg.db,
                   keytype = "GID",
                   ont = "BP"
)
################################################################################
# Gene  contrast strict
################################################################################

Cluster_GO <- function(color){
  
  eGO_global_strict <<- enrichGO(gene          = get(paste0("it_pop_",color,"_filt")),
                             universe      = universe,
                             OrgDb         = org.MdomesticaGDDH13.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "bonferroni",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE,
                             keyType = "GID")
  
  # results <- eGO_global_strict@result
  # ego <- pairwise_termsim(eGO_global_strict,method = "Wang", semData = semdata)
  # pdf(file = paste0("./Results/EXP3/Go_WGCNA/Global_",color,".pdf"),width = 11.7, height = 8.3)
  # p1 <- dotplot(eGO_global_strict, showCategory = 15)
  # p2 <- emapplot(ego,showCategory = 50,cex.params = list(category_label = 0.7))
  # print(p1)
  # print(p2)
  # dev.off()
  # write_delim(x = results, file = paste0("./Results/EXP3/Go_WGCNA/GO_",color,".tsv"), 
  #             quote = "none", delim = "\t")
  # 
  # results$Module = color
  # results <- filter(results, results$p.adjust < 0.05)
  # 
  # if(!exists("GO_global") || is.null(GO_global) || nrow(GO_global) == 0){
  #   GO_global <<- results
  # } else {
  #   GO_global <<- rbind(GO_global, results)
  # }
  # return(GO_global)
}


for (i in 1:nrow(info)) {
  if (info$Percentage_constitutive_resistance_genes[i] >= 14) {
    Cluster_GO(info$Module_Name[i])
  }
}

# Convert the character gene ratio for exemple "44/162" into a real ratio by separating both numbers and divide them
make_ratio_numeric = function(ratio) {
  ratio = as.numeric(unlist(str_split(ratio, "/"))[1])/as.numeric(unlist(str_split(ratio, "/"))[2])
  ratio
}

GO_global$`Gene Ratio` =  sapply(GO_global$GeneRatio, make_ratio_numeric)

order_GO <- GO_global$Description
# GO_global$Description <- factor(GO_global$Description, levels = make.unique(order_GO))

# GO_global_Limited <- GO_global %>% 
#   filter(Module %in% c("blue","yellow","red","turquoise","brown"))

ggplot(data = GO_global, aes(y = Description,
                             x = "", 
                          color = `p.adjust`,
                          size = `Gene Ratio`)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  
  
  xlab("") + 
  ggtitle("GO analysis") + 
  facet_wrap(~ Module, ncol = 10)
ggsave(                                                                       
  file   = paste0("./Results/EXP3/Go_WGCNA/Enrichment_C.pdf"),
  plot   = last_plot(),
  width  = 20,
  height = 15
)
write.csv(info,"./Results/EXP3/Go_WGCNA/info.csv", quote = F, row.names = F)

################## BUBBLE PLOT #################################################

library(GOplot)
library(clusterProfiler)
library(org.MdomesticaGDDH13.eg.db)
library(patchwork)



prep_gene <- function(str_obj){
    return(gsub("/", ", ", str_obj)[[1]])
  }


DEG <- read.csv2(file = "./Results/Exp3/DiffAnalysis/[M. sieversii_Ctrl-M. sieversii_Inf]/LRT_BH.txt", sep = "\t") %>% 
  filter(Gene_ID %in% it_pop_midnightblue_filt)
gene <- DEG[[1]]
DEG <- DEG[,c(1,3,6)]
names(DEG) <- c("ID","logFC","adj.P.Val")

ego <- enrichGO(gene          = gene,
                universe      = universe,
                OrgDb         = org.MdomesticaGDDH13.eg.db,
                ont           = "BP",
                pAdjustMethod = "bonferroni",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE, keyType = "GID")

ego2 <- clusterProfiler::simplify(ego, semData = semdata)
result = ego2@result
result= result[,c(1,2,6,8)]
names(result) <- c("ID","Term","adj_pval","Genes")
result$Category= "BP"
result$Genes <- unlist(lapply(result$Genes, prep_gene))
result <- result[which(result$adj_pval < 0.05) ,]
circ3 <- circle_dat(result, DEG) %>% 
  mutate(module_names = "midnightblue")

combined_data <- rbind(circ1, circ2, circ3,circ4)

p2 <- GOBubble_module(combined_data)

# p3.3 <- GOBubble(circ,labels = FALSE, table.col=FALSE, table.legend	= FALSE, ID = FALSE)
