library(MOFA)
library(MOFAdata)
library(MultiAssayExperiment)
library(reticulate)
library(mofapy)
library(rhdf5)

#Using a specific python binary
py_install("mofapy", envname = "r-reticulate", method="auto") 
use_python('/usr/local/bin/python3', required = TRUE) 
mofapy <- import('mofapy')


setwd('~/Documents/Helsinki_COVID19/MOFA_Helsinki') 
#Load Olink data
d_olink = read.csv('Olink_all.csv', sep='\t', row.names = 1, header = TRUE)
d_olink = as.data.frame(t(d_olink))
head(d_olink)
#Load Grid data (Z-score transformation)
d_grid = read.csv('grid_all.csv', sep=';', row.names = 1, header = TRUE)
d_grid = as.data.frame(t(d_grid))
head(d_grid)
#Load Metadata
d_meta = read.csv('meta_all_Rdays.csv', sep=';', row.names = 1)
d_meta = as.data.frame(d_meta)
head(d_meta)

#MOFA object
COV_data = list(d_olink, d_grid)
names(COV_data) <- c('Plasma Protein', 'Cell Population')
mae_COV <- MultiAssayExperiment(experiments = COV_data, colData = d_meta)
MOFAobject <- MOFA::createMOFAobject(mae_COV) 
MOFAobject 
MOFA::plotDataOverview(MOFAobject)


#Fit the MOFA model
DataOptions <- MOFA::getDefaultDataOptions()
DataOptions

ModelOptions <- MOFA::getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 10 
ModelOptions

TrainOptions <- MOFA::getDefaultTrainOptions()
#TrainOptions$DropFactorThreshold <- 0.02
TrainOptions

MOFAobject <- MOFA::prepareMOFA(
  MOFAobject, 
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)

#Run MOFA
MOFAobject_COV <- runMOFA(MOFAobject)
MOFAobject_COV
# Calculate the variance explained (R2) per factor in each view 
r2 <- calculateVarianceExplained(MOFAobject_COV)
r2$R2Total
plotVarianceExplained(MOFAobject_COV)

#plot the top loadings for a given factor and view
plotTopWeights(
  MOFAobject_COV, 
  view = "Cell Population", 
  factor = 2
)

plotTopWeights(
  MOFAobject_COV, 
  view = "Plasma Protein", 
  factor = 2
)

#correlation plot between factors -> should be uncorrelated
plotFactorCor(MOFAobject_COV, method = 'pearson')

#scatterplot between two factors, similar to PCA
fs <- plotFactorScatters(MOFAobject_COV, factors = 1:5, color_by = 'Overall.clinical.grade') #Day..

library(wesanderson)
wes_palette("Zissou1")
pal <- wes_palette("Zissou1", 3, type = "continuous")

#Key gradient
colfunc <- colorRampPalette(c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)

for(i in 1:fs$nrow) {
  for(j in 1:fs$ncol){
    fs[i,j] <- fs[i,j] + 
      geom_point(size=2.0) +
      scale_colour_gradient2(low='#3B9AB2', mid='#EBCC2A', high="#F21A00")
    
  }
}
fs  

MOFAweights <- getWeights(
  MOFAobject_COV, 
  views = "all", 
  factors = "all", 
  as.data.frame = TRUE    # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
)
MOFAweights
write.csv(MOFAweights, file='MOFAweights.csv')

factor2 <- sort(getFactors(MOFAobject_COV,"LF2")[,1])
order_samples <- names(factor2)
df <- data.frame(
  row.names = order_samples,
  factor = factor2
)

df1 <- getFactors(MOFAobject_COV, 'all')
write.csv(df1, file='AllFactors.csv')

#Correlation
library(ggplot2)
library(ggpubr)
ID <-getCovariates(MOFAobject_COV,'Subject_ID') 
cdr <-getCovariates(MOFAobject_COV,'Day..') 
ICU <-getCovariates(MOFAobject_COV,'ICU') 

factor2 <- getFactors(MOFAobject_COV,
                      factors=2)
COV_LF2 <- data.frame(factor = as.numeric(factor2), cdr = cdr)
ggplot(COV_LF2, aes_string(x = "factor", y = "cdr", color=ICU)) + 
  geom_point()+
  xlab("Latent Factor 2") +
  ylab("Days") +
  geom_smooth(se = FALSE, method = lm,aes(group=ID))+
  theme_bw() + coord_flip()
