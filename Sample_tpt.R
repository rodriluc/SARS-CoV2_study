setwd('~/Documents/Helsinki_COVID19')

######################################
#                Meta                #
######################################

meta <- read.csv('meta_Helsinki.csv', sep = ';')
head(meta)

library(tidyverse)
library(ggplots)
library(dplyr)
library(tidyr)
library(forcats)

meta %>% ggplot(aes(x = Sampling.date, y = Subject_ID, group = Subject_ID, col=ICU)) + 
  geom_line() + 
  geom_point(shape = 15) +
  geom_point(aes(x = Sampling.date.CyTOF, y = Subject_ID), shape=4, size=3, col="black") +
  geom_point(aes(x = Sampling.date.Olink, y = Subject_ID), shape=1, size=3,col="dark green") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) 
