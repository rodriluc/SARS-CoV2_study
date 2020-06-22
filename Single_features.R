setwd('~/Documents/Helsinki_COVID19/Olink_20200515')

######################################
#             Extract PP             #
######################################

Acute <- read.csv('PB_Covid-19_Plate1_filterA.csv', sep = ';')
head(Acute)
Recovered <- read.csv('PB_Covid-19_Plate1_filterR.csv', sep = ';')

######################################
#         Plot NPX vs. Days          #
######################################

library(grid)
library(directlabels)
library(ggplot2)
library(gtools)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggforce)
require(gridExtra)

y <- Acute[,-1]
x <- Recovered[,-1]

#Grid
setwd('~/Documents/Helsinki_COVID19')
cellA <- read.csv('meta_clinical.csv', sep = ';')
#y <- cellA[,-1]

#cellR <- read.csv('abundance_sample_lvl3_R.csv', sep = ';')
#x <- cellR[,-1]

#limits <- c(-0.008, 0.03)

pp = 'P/F (mmHg, lowest)'
p <- ggplot(cellA, aes(x=Days, y=P.F..mmHg..lowest., group=Subject_ID, colour=Subject_ID)) +  
  ggtitle(pp) +
  geom_line(aes(group = Subject_ID)) + labs(x= 'Days', y = 'P/F (mmHg, lowest)') +
  #scale_y_continuous(limits=limits) + 
  theme_bw()
p
q <- ggplot(x, aes(x=as.character(Days), y=Other.gdT)) +  
  ggtitle(pp) +
  geom_violin(trim=FALSE, colour = "3", alpha = .5, size=1) +
  geom_jitter(shape=16, position=position_jitter(0.15)) + labs(x= 'Days', y = 'Cell abundance') +
  scale_y_continuous(limits=limits) + theme_bw() 

grid.arrange(p, q, ncol=2)

#All feature plots
setwd('~/Documents/Helsinki_COVID19/Olink_20200515')
  
#Acute
plot1 <- y %>%
  gather(-Subject_ID, -ICU, -Days, key = "var", value = "value") %>% 
  ggplot(aes(x = Days, y = value)) +
  geom_line() +
  scale_x_continuous(breaks=seq(-1,16,3)) +
  facet_wrap(~ var, scales = "free")

#Recovered
plot2 <- x %>%
  gather(-Subject_ID, -Days, key = "var", value = "value") %>% 
  ggplot(aes(x = Days, y = value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  facet_wrap(~ var, scales = "free")

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))


######################################
#                PCA                 #
######################################
library(ggfortify)
library(devtools)

df <- Acute[,-c(1:4)]
head(df)

pca_res <- prcomp(df, scale. = TRUE)

autoplot(pca_res, data = Acute, colour = 'clinical_grade') + scale_color_gradient(low = "green",
                                                                        high = "red")
pairs(pca_res$loadings)

library(factoextra)
res.var <- get_pca_var(pca_res)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

write.csv(res.var$contrib, file = 'Contrib_PCs.csv')

# Dissimilarity matrix
Acute <- read.csv('PB_Covid-19_Plate1_filterA.csv', sep = '\t')
df <- Acute[,-c(1:2)]
df <- df[,-2]
head(df)

df <- df[,-1]
rownames(df) <- df[,1]

df <- na.omit(df)
df <- scale(df)
d <- dist(df, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)
