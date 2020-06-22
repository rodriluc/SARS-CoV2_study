library(ggplots)

setwd('~/Documents/Helsinki_COVID19/serology/')

##################################
#        Input and subset        #
##################################

#Load Recovered
d_R = read.csv('sero_mixedcorr.csv', sep=';', row.names = 1)
head(d_R)
#load Grid <=3
d_3 = read.csv('grid_4_2.csv', sep=';', row.names = 1) 
head(d_3)
#load Grid 4-7
d_7 = read.csv('8_2.csv', sep=';', row.names = 1) 
head(d_7)
#load Grid 8-14
d_14 = read.csv('grid_14_2.csv', sep=';', row.names = 1) 
head(d_14)

##################################
#    Spearman Correlation        #
##################################

# Spearman correlation matrix was undertaken and rounded to 2 decimal places
cormat_R <- round(cor(d_R, method = c("spearman")),2)
cormat_3 <- round(cor(d_3, method = c("spearman")),2)
cormat_7 <- round(cor(d_7, method = c("spearman")),2)
cormat_14 <- round(cor(d_14, method = c("spearman")),2)

##################################
#            Re-order            #
##################################
# Re-order the correlation matrix by using correlation between variables as distance
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

#cormat0 <- reorder_cormat(cormat0)
melted_cormatR <- reorder_cormat(cormat_R)

##################################
#      Prep for plotting         #
##################################

library(reshape2)
melted_cormatR <- melt(cormat_R)
head(melted_cormatR)

melted_cormat3 <- melt(cormat3)
head(melted_cormat3)

melted_cormat7 <- melt(cormat_7) 
head(melted_cormat7)

melted_cormat14 <- melt(cormat14)
head(melted_cormat14)

##################################
#              Plot              #
##################################
library(corrplot)
library("Hmisc")

# Plots Corr. matrix of AT 16 re-ordered
library(ggplot2)
library(viridis)

mc4 <- ggplot(data = melted_cormatR, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation")+
  #scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 3, hjust = 1))+
  theme(axis.text.y = element_text(size = 3))+

  coord_fixed() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 
mc4

# Set vector of levels you want
mylevels <- melted_cormat7$Var1
# Re-order factors
melted_cormat3$Var1 <- factor(melted_cormat3$Var1,levels=unique(mylevels))
melted_cormat3$Var2 <- factor(melted_cormat3$Var2,levels=unique(mylevels))

# Plots Corr. matrix of comat3 re-ordered to cormat7
mc0 <- ggplot(data = melted_cormat3, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#b2182b", high = "#2166ac", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1))+
  coord_fixed() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 
mc0
