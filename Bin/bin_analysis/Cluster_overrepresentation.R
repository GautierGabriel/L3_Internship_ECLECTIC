library(ggplot2)   
library(tidyverse)


tot <- read.table("./Results/EXP3/Go_WGCNA/GO_saddlebrown.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "") %>% 
  mutate( "module" = "saddlebrown")

tot$Description <- factor(tot$Description, levels = rev(tot$Description))
tot<- head(tot, 20)

blue <- read.table("./Results/EXP3/Go_WGCNA/GO_turquoise.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
  mutate( "module" = "turquoise")
blue<- head(blue, 20)

tot <- rbind(tot,blue)

make_ratio_numeric = function(ratio) {
  ratio = as.numeric(unlist(str_split(ratio, "/"))[1])/as.numeric(unlist(str_split(ratio, "/"))[2])
  ratio
}

tot$`Gene Ratio` =  sapply(tot$GeneRatio, make_ratio_numeric)

ggplot(tot, aes(y = Description,
                             x = "", 
                             color = `p.adjust`,
                             size = `Gene Ratio`)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO analysis") + 
  facet_wrap(~ module, ncol = 10)

ggsave(                                                                       
  file   = "./Results/EXP3/Go_WGCNA/Presence_Sieversii.pdf",
  plot   = last_plot(),
  width  = 20,
  height = 15
)
