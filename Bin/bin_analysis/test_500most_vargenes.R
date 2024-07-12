library(gplots)
library(factoextra)
library(DESeq2)
samples <- read.csv("./Data/EXP3_TARGET.csv", header = T, sep = ",")
counts <- read.csv("./Data/EXP3_COUNTS.csv", header = T, sep = ",", row.names = 1)
counts[] <- lapply(counts, function(x) {
  if (is.numeric(x)) {
    return(as.integer(x))
  } else {
    return(x)
  }
})


dds <- DESeqDataSetFromMatrix(counts,samples, design = ~Condition)
##stabilisation de la variance car RNAseq est biaisé vers les gènes de haut comptes ( cf exponentielle de la pcr) 
dds_counts <- assay(dds)


var_genes <- apply(dds_counts,1, FUN = var)
var_genes <- sort(var_genes,decreasing = TRUE)

high_var_genes <- var_genes[1:500]
Gene_names <- names(high_var_genes)

dds_vst_top_500 <- dds[Gene_names,]

palette <- hcl.colors(100, palette = "RdBu")

heatmap.2(assay(dds_vst_top_500), 
          trace = "none", 
          col = palette)
heatmap.2(t(scale(t(assay(dds_vst_top_500)))), trace = "none", col = palette)

# adding a colour bar to heatmap2

disease_color <- rep("red", dim(samples)[1])
disease_color[samples$Condition == "Ctrl"] <- "blue"

dds_scaled <- t(scale(t(assay(dds_vst_top_500))))
heatmap.2(dds_scaled,
          trace = "none",
          col = palette,
          cexRow = 0.15,
          cexCol = 0.25,
          ColSideColors = disease_color)


Ref = merge(x = read.csv2("./Data/gene_id.txt", header = T), y = Reference_Enrichment, by = "Gene_ID")