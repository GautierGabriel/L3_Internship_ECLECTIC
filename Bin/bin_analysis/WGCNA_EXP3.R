# Author : adapted by Gabriel Gautier from Ronan Dadole
# Date : 30/11/2023 adapted 17.06.2024
# GO enrichment analysis EXP3

###############################  WGCNA #########################################


################################################################################
library(tidyverse)
library(magrittr)
library(WGCNA)
library(edgeR)
library(matrixStats)



input_mat = read.csv("./Results/EXP3/DiffAnalysis/NormCounts_log2.txt",row.names = 1, header= TRUE, sep = "\t") %>% 
  select('Gene_Name',infected_trees)

infected_trees <- read.csv("./Data/infected.txt", header= TRUE) %>% 
  t(.) %>% 
  as.vector(.)



input_mat <- input_mat[,-1]
input_mat <- 2^input_mat - 1
input_mat = as.data.frame(t(input_mat))

traitData <- read.csv2("./Data/Phenotypic_traits_infested.csv",row.names = 1, header = T, sep =  ";")
traitData <- traitData[,-c(1:15,17,19,21,28,29)]
names(traitData)
names(traitData) <- c("Chlf", "Flavf", "NBIf","Winged_Adult","Apterous_Adult","N4","L4","L2_L3_N2_N3","Molt")




allowWGCNAThreads() 
################################################################################

powers = c(1:25)
sft = pickSoftThreshold(
  input_mat,           
  powerVector = powers,
  verbose = 5,
  networkType = "unsigned"
  
  
)

par(mfrow = c(1,2));
cex1 = 0.9;
sft$fitIndices$Power <- as.numeric(sft$fitIndices$Power)
graphics::plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
graphics::plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

################################################################################

picked_power = 8

################################################################################

temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "unsigned",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 35000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.3,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
pdf(file = "./Results/EXP3/WGCNA/Cluster_Dendrogram.pdf")
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
dev.off()

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
################################################################################

length(unique(module_df$colors))


write_delim(module_df,
            file = "./Results/EXP3/WGCNA/gene_modules_filt.txt",
            delim = "\t")

################################################################################

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

pdf(file = "./Results/EXP3/WGCNA/Sample-trait_Relationships.pdf")
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90, size=8)) +  # Changez la taille à votre convenance
  labs(title = "Sample-trait Relationships", y = "Modules", fill="corr")
dev.off()

# # tidy & plot data
# mME = MEs0[,c(4,5,12,13,16,23,24,25,26,29,33,37,38,44)] %>%
#   pivot_longer(-treatment) %>%
#   mutate(
#     name = gsub("ME", "", name),
#     name = factor(name, levels = module_order)
#   )
# 
# pop <- read.csv("./Data/Phenotypic_traits_reduced.csv", sep = ";")
# merged <- mME %>%
#   inner_join(pop, by = c("treatment" = "labels")) %>%
#   mutate(treatment = factor(treatment, levels = unique(treatment[order(Population)])))
# 
# pdf(file = "./Results/EXP3/WGCNA/Sample-trait_Relationships_population.pdf")
# merged %>%
#   ggplot(aes(x = treatment, y = name, fill = value)) +
#   geom_tile() +
#   theme_bw() +
#   scale_fill_gradient2(
#     low = "blue",
#     high = "red",
#     mid = "white",
#     midpoint = 0,
#     limit = c(-1, 1)
#   ) +
#   theme(axis.text.x = element_text(angle = 90, size = 8)) +
#   labs(title = "Sample-trait Relationships_reduced_population", y = "Modules", fill = "corr") +
#   # Add population names below treatments
#   facet_wrap(~Population, scales = "free_x", strip.position = "bottom") +
#   theme(
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     strip.text.x = element_text(angle = 90, hjust = 1)
#   )
# dev.off()
# 
# 
# write.csv2(MEs0,file = "./Results/EXP3/WGCNA/Sample-trait.csv")

################################################################################

nGenes = ncol(input_mat)
nSamples = nrow(input_mat)

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(input_mat, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traitData, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(8,4)

#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

pdf(file = "./Results/EXP3/WGCNA/Module-trait_Relationship.pdf")
#Displaying the correlation values in a heatmap plot
p_trait_cor <- labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = row.names(moduleTraitCor),
               ySymbols = row.names(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(10),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.7,
               cex.lab = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               yColorWidth = 0.01)
p_trait_cor
dev.off()

write.table(moduleTraitCor,"./Results/EXP3/WGCNA/Trait_correlation_filt.txt", sep= "\t", quote = FALSE)

################################################################################

modules = unique(module_df$colors)
dir.create(file.path("./Results/EXP3/WGCNA/modules_filt"), showWarnings = FALSE)

for (nb in seq(1,length(modules))){
  mod = modules[nb]
  genes = filter(module_df, colors == mod)[[1]]
  write.table(genes,file=paste0("./Results/EXP3/WGCNA/modules_filt/module_",mod,"_filt.txt"),quote=FALSE,row.names = FALSE,col.names = FALSE)
}

################################################################################

chooseTopHubInEachModule(input_mat,
                         module_df$colors,
                         power = picked_power,
                         type = "unsigned")

################################################################################
################################################################################
# Cluster profiler analysis of gene modules
library(GOplot)
require(org.MdomesticaGDDH13.eg.db)
library(clusterProfiler)
library(GOSemSim)
library(enrichplot)
library(ggpubr)
library(stringr)
library(ggupset)

universe <- names(input_mat)
semdata <-godata(  OrgDb = org.MdomesticaGDDH13.eg.db,
                   keytype = "GID",
                   ont = "BP"
)

wd_modules = "./Results/EXP3/WGCNA/modules_filt/"
module_file = list.files(wd_modules, pattern = regex("module"))



enrich <- function(gene, universe, module,path){
  ego <- enrichGO(gene          = gene,
                  universe      = universe,
                  OrgDb         = org.MdomesticaGDDH13.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  minGSSize = 50,
                  readable      = TRUE,
                  keyType = "GID")
  
  ego2 <- simplify(ego, semData = semdata)
  ego3 <- pairwise_termsim(ego,method = "Wang", semData = semdata)
  ego4 <- pairwise_termsim(ego2,method = "Wang", semData = semdata)
  pdf(file.path(paste0(path,module,"_Overrepresentation.pdf")),width = 14,height = 11)
  x1 <- dotplot(ego, showCategory=20) 
  x2 <- dotplot(ego2, showCategory=20) 
  x3 <- emapplot(ego3,showCategory = 50,cex.params = list(category_label = 0.7))
  x4 <- emapplot(ego4,showCategory = 50,cex.params = list(category_label = 0.7))
  print(x1)
  print(x2)
  print(x3)
  print(x4)
  dev.off()
  write.table(ego@result,file.path(paste0(path,module,"_Overrepresentation.tsv")),quote = FALSE, sep="\t")
}



for (i in seq(1, length(module_file))){
  file_name = module_file[i]
  genes = read.csv(file.path(wd_modules,file_name), header=FALSE)[[1]]
  module_color = substr(file_name,8,nchar(file_name)-4)
  try(enrich(genes,universe,module_color,wd_modules))
}