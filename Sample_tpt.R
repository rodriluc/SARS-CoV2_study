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
  #scale_x_date(date_minor_breaks = "1 day") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(limits = c('2020-04-08','2020-04-09','2020-04-10','2020-04-11','2020-04-12','2020-04-13','2020-04-14','2020-04-15','2020-04-16','2020-04-17','2020-04-18','2020-04-19','2020-04-20','2020-04-21','2020-04-22','2020-04-23','2020-04-24','2020-04-25','2020-04-26','2020-04-27','2020-04-28','2020-04-29','2020-04-30','2020-05-01','2020-05-02','2020-05-03','2020-05-04','2020-05-05','2020-05-06','2020-05-07','2020-05-08'))
