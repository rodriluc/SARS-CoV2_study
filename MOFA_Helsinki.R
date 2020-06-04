library(MOFA)
library(MOFAdata)
library(MultiAssayExperiment)
library(reticulate)
library(mofapy)
library(rhdf5)

# Using a specific python binary
py_install("mofapy", envname = "r-reticulate", method="auto") 
use_python('/usr/local/bin/python3', required = TRUE) 
mofapy <- import('mofapy')


setwd('~/Documents/Helsinki_COVID19/MOFA_Helsinki') 
#Load Olink data
d_olink = read.csv('Olink_all.csv', sep='\t', row.names = 1, header = TRUE)
d_olink = as.data.frame(t(d_olink))
head(d_olink)
#load Grid data (Z-score transformation)
d_grid = read.csv('grid_all.csv', sep=';', row.names = 1, header = TRUE)
d_grid = as.data.frame(t(d_grid))
head(d_grid)
#Metadata
d_meta = read.csv('meta_all_Rdays.csv', sep=';', row.names = 1)
d_meta = as.data.frame(d_meta)
head(d_meta)
d_meta

#MOFA object
ME_data = list(d_olink, d_grid)
names(ME_data) <- c('Plasma Protein', 'Cell Population')
mae_ME <- MultiAssayExperiment(experiments = ME_data, colData = d_meta)
MOFAobject <- MOFA::createMOFAobject(mae_ME) 
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

#option to regress covariates HERE
#MOFAobject_H Last model
#Run MOFA
MOFAobject_fix <- runMOFA(MOFAobject)
MOFAobject_fix
# Calculate the variance explained (R2) per factor in each view 
r2 <- calculateVarianceExplained(MOFAobject_fix)
r2$R2Total
plotVarianceExplained(MOFAobject_fix)

#plot a heatmap of the loadings from multiple factors in a given view
plotWeightsHeatmap(
  MOFAobject_fix, 
  view = "Plasma Protein", 
  factors = 1:5,
  show_colnames = TRUE
)
#plot the top loadings for a given factor and view
plotTopWeights(
  MOFAobject_fix, 
  view = "Cell Population", 
  factor = 2
)

#plot all loadings for a given factor and view
plotWeights(
  MOFAobject_fix, 
  view = "Plasma Protein", 
  factor = 2
) 

#correlation plot between factors -> should be uncorrelated
plotFactorCor(MOFAobject_fix, method = 'pearson')

#scatterplot between two factors, similar to PCA
fs <- plotFactorScatters(MOFAobject_fix, factors = 1:5, color_by = 'Overall.clinical.grade') #Day..

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

set.seed(1234)
clusters <- clusterSamples(
  MOFAobject_fix, 
  k = 3,        # Number of clusters for the k-means function
  factors = 5   # factors to use for the clustering
)

plotFactorScatter(
  MOFAobject_fix, 
  factors = 4:5, 
  color_by = clusters
)

plotFactorBeeswarm(
  MOFAobject_fix,
  factors = 1,
  color_by = "Day.."
)

MOFAweights <- getWeights(
  MOFAobjectSS2, 
  views = "all", 
  factors = "all", 
  as.data.frame = TRUE    # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
)
MOFAweights
write.csv(MOFAweights, file='MOFAweights.csv')

factor2 <- sort(getFactors(MOFAobject_fix,"LF2")[,1])
order_samples <- names(factor2)
df <- data.frame(
  row.names = order_samples,
  factor = factor2
)

df1 <- getFactors(MOFAobject_fix, 'all')
write.csv(df1, file='AllFactors.csv')

plotDataHeatmap(
  object = MOFAobject_fix, 
  view = "Cell Population", 
  factor = "LF2", 
  features = 20, 
  transpose = TRUE, 
  show_colnames=FALSE, show_rownames=TRUE, # pheatmap options
  cluster_cols = FALSE, annotation_col=df # pheatmap options
)

plotWeights(
  object = MOFAobject_fix,
  view = "Plasma Protein", 
  factor = 3, 
  nfeatures = 0,
  abs = FALSE,
  scale = FALSE
)

#Correlation
library(ggplot2)
library(ggpubr)
ID <-getCovariates(MOFAobject_fix,'Subject_ID') 
cdr <-getCovariates(MOFAobject_fix,'Day..') 
ICU <-getCovariates(MOFAobject_fix,'ICU') 

factor2 <- getFactors(MOFAobject_fix,
                      factors=2)
foo <- data.frame(factor = as.numeric(factor2), cdr = cdr)
ggplot(foo, aes_string(x = "factor", y = "cdr", color=ICU)) + 
  geom_point()+
  xlab("Latent Factor 2") +
  ylab("Days") +
  geom_smooth(se = FALSE, method = lm,aes(group=ID))+
  theme_bw() + coord_flip()
