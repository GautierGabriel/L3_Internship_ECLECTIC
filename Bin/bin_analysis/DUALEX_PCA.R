#################### DUALEX PCA ################################################
library(ade4)
library(tidyverse)
traitData <- read.csv2("./Data/Phenotypic_traits.csv",row.names = 1, header = T, sep =  ";")[,c(22:29)] %>% 
  filter(Condition == "Inf")
  
traitData <- traitData[,-c(8)] %>% 
  mutate(somme = rowSums(across("winged adults":"molt")))

names(traitData)
names(traitData) <- c("winged adults", "apterous adults", "nymphs 4","larvae 4", "larvae and nymphs 1,3", "molt", "population")

traitData <- traitData %>% 
  mutate(somme = rowSums(across("winged adults":"molt")))

# Création du boxplot
ggplot(traitData, aes(x = population, y = somme)) +
  geom_boxplot() +
  labs(title = "Boxplot of Sum of Aphids by Population",
       x = "Population",
       y = "Sum of Aphids") +
  theme_minimal()

# ANOVA pour comparer plusieurs populations
modele <- lm(somme ~ population, data = traitData)
plot <- autoplot(modele)

shapiro.test(residuals(modele))

anova_result <- aov(somme ~ population, data = traitData)
summary(anova_result)
     
# Tests post-hoc pour identifier les différences spécifiques
dunnTest(SumAphids ~ Population, data = data_extended, method = "bonferroni")

# PCA
OutputPCA <- dudi.pca(traitData[,c(1:6,8)],center=T,scale=T,scannf=F,nf=5)
pve <- 100*OutputPCA$eig/sum(OutputPCA$eig)

fviz_pca_ind(OutputPCA,geom="point",
             habillage=traitData$population)
ggsave(                                                                       
  file   = paste0("./Results/EXP3/WGCNA/aphid_population.pdf"),
  plot   = last_plot(),
  width  = 7,
  height = 7
)

fviz_pca_var(OutputPCA, fill.var = "contrib")
ggsave(                                                                       
  file   = paste0("./Results/EXP3/WGCNA/aphid_contrib.pdf"),
  plot   = last_plot(),
  width  = 7,
  height = 7
)

#quantification des effets 
bca(OutputPCA,as.factor(traitData$population),scannf=F,nf=30)$ratio
# bca(OutputPCA,as.factor(traitData$condition),scannf=F,nf=30)$ratio
