library(lme4)
library(ggplot2)
library(stargazer)
library(blme)
library(sjPlot)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

setwd('~/Documents/Mixed_Effects_Model/AT_rerun/')

##################################
#       Test trial for MEM       #
##################################

# Read in file, build dataframe (columns should be covariates and cell populations, 
# row should be ID and subsequent info.)
me_memG = read.csv('grid_ME12_new.csv', sep = ';')
head(me_memG)
me_memG_df = data.frame(me_memG)

## Construct Mixed-Effect model
# Linear MEM
me_lmer.model = lmer(IgD.pos..Memory.B ~ sex + age + symptom_score + Group + active_treatment +
                        (1|id/Group), data=me_memG)
# Partial Bayesian MEM
me_blmer.model = blmer(IgD.pos..Memory.B ~ sex + age + symptom_score + Group + active_treatment +
                        (1|id/Group), data=me_memG)

# Rescale and center continuous parameters
numcols <- grep("^c\\.",names(me_memG_df))
dfs <- me_memG_df
dfs[,numcols] <- scale(dfs[,numcols])
me_lmer.model <- update(me_lmer.model,data=dfs)

# Evaluates whther a fitted mixed model is singular, if singular TRUE then opt for solution 
# One of these solutions being a partial bayesion model especially with complex models
isSingular(me_lmer.model, tol=1e-05)

# Print model output
summary(me_lmer.model)

# Create table of MEM model either with tab_model() or stargazer()
library(sjPlot)
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table') 
tab_model(me_lmer.model,me_mem4.model,me_mem20.model, file = 'tableME_grid_trial.html')

stargazer(me_lmer.model,
          type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "", out = 'tableME_grid_trial.html')

## Statistical significance (ex. SS)

#1) construct null model first
me_lmer.null = lmer(IgD.pos..Memory.B ~ sex + age + (1|Group) +
                      (1|id), data=me_mem1, REML=FALSE)
#2) Re-do full model
me_lmer.model = lmer(IgD.pos..Memory.B ~ sex + age + symptom_score + (1|Group) + 
                       (1|id), data=me_mem1, REML=FALSE)
#3) Perform likelihood ratio test
anova(me_lmer.null, me_lmer.model)

# Looks at the coefficients of the model by subject and by item
coef(me_lmer.model)

##################################
#   Input and Modelling - GRID   #
##################################

# Read in file, build dataframe (columns should be covariates and cell populations, 
# row should be ID and subsequent info.)
me_memG = read.csv('grid_ME12_new.csv', sep = ';')
head(me_memG)
me_memG_df = data.frame(me_memG)

# List of cell population names in order to loop through, omitting first 6 columns which 
# includes ID and covariates
G_list = colnames(me_memG_df[-(1:6)]) 
head(G_list)

# Mixed-effect modelling (Partial-bayesian)
varlist=G_list 
blups.models_G <- lapply(varlist, function(x) {
  mod2 = try(blmer(substitute(i ~ sex + age + symptom_score + Group + active_treatment + 
                   (1|id/Group), list(i = as.name(x))), 
                   data = me_memG_df, na.action=na.exclude))
  if(isTRUE(class(mod2)=='try-error')) {return(NULL)} else{return(mod2)}
})

blups.models_nextG = as.list(blups.models_G)
# Remove NULL models that failed to converge
blups.models_nextG[sapply(blups.models_nextG, is.null)] <- NULL 

##################################
# Extract info from model - GRID #
##################################

library(predictmeans)

# Mixed-effect model expression
varlist1 <- lapply(blups.models_nextG, function(f) summary(f)$call[2])
varlist1

# Correlation Coefficient
estimate_AT <- lapply(blups.models_nextG, function(f) summary(f)$coefficients[6,1])
estimate_KOS <- lapply(blups.models_nextG, function(f) summary(f)$coefficients[5,1])
estimate_SS <- lapply(blups.models_nextG, function(f) summary(f)$coefficients[4,1])

# R^2 and adjusted
R2 <- lapply(blups.models_nextG, function(f) summary(f)$r.squared)
R2_adj <- lapply(blups.models_nextG, function(f) summary(f)$adj.r.squared)

# Median
med <- lapply(blups.models_nextG, function(f) summary(f)$residuals)
med_na <- lapply(med, function(f) na.exclude(f))
med_calc <- lapply(med_na, function(f) median(f))

## p-values for covariates (fixed-effect)
# SS
test_pSS <- lapply(blups.models_nextG, function(f) parameters::p_value_wald(f)[4,])
df.pvalueSS <- as.data.frame(test_pSS)
df.pvalueSS <-t(df.pvalueSS)
df.pvalueSS = df.pvalueSS[seq(0, nrow(df.pvalueSS), 2), ]
# KOS
test_pKOS <- lapply(blups.models_nextG, function(f) parameters::p_value_wald(f)[5,])
df.pvalueKOS <- as.data.frame(test_pKOS)
df.pvalueKOS <-t(df.pvalueKOS)
df.pvalueKOS = df.pvalueKOS[seq(0, nrow(df.pvalueKOS), 2), ]
# AT
test_pAT <- lapply(blups.models_nextG, function(f) parameters::p_value_wald(f)[6,])
df.pvalueAT <- as.data.frame(test_pAT)
df.pvalueAT <-t(df.pvalueAT)
df.pvalueAT = df.pvalueAT[seq(0, nrow(df.pvalueAT), 2), ]

# Prepare dataframe with extracted info. for downstream use
test_data = list(as.character(varlist1), med_calc, estimate_AT, estimate_KOS, estimate_SS, as.numeric(t(df.pvalueSS)), as.numeric(t(df.pvalueKOS)), as.numeric(t(df.pvalueAT)))
names(test_data) <- c('cells', 'Median', 'Estimate_AT', 'Estimate_KOS', 'Estimate_SS', 'pvalueSS', 'pvalueKOS', 'pvalueAT')

test_final <- as.data.frame(do.call(rbind, test_data))
t(test_final)
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table')
write.csv2(t(test_final), file='grid_MEMtable_newSS.csv')

##################################
#    Input and Modelling - PP    #
##################################

setwd('~/Documents/Mixed_Effects_Model/AT_rerun') 

# Read in file, build dataframe (columns should be covariates and plasma proteins, 
# row should be ID and subsequent info.)
me_mem1 = read.csv('Olink_ME12_new.csv', sep = ';')

setwd('~/Documents/Helsinki_COVID19/serology/MEM')
me_mem1 <- read.csv('MEM_RBD.csv', sep = ';')
#write.csv(me_mem1, file = 'abundance_TRIAL_lvl3_meta.csv')
head(me_mem1)
me_mem1_df = data.frame(me_mem1)

# List of plasma protein names in order to loop through, omitting first 6 columns which 
# includes ID and covariates
PP_list = colnames(me_mem1_df[-(1:4)])
head(PP_list)

# Mixed-effect modelling (Partial-bayesian)
varlist=PP_list 
blups.models_PP <- lapply(varlist, function(x) {
  mod2 = try(blmer(substitute(i ~ old.RBD + Days + 
                   (1|Subject_ID), list(i = as.name(x))), 
                   data = me_mem1_df, na.action=na.exclude))
  if(isTRUE(class(mod2)=='try-error')) {return(NULL)} else{return(mod2)}
})

blups.models_nextPP = as.list(blups.models_PP)
# Remove NULL models that failed to converge
blups.models_nextPP[sapply(blups.models_nextPP, is.null)] <- NULL

##################################
#  Extract info from model - PP  #
##################################
me_blmer.model = blmer(pDC ~ old.RBD + Days +
                         (1|Subject_ID), data=me_mem1_df)
summary(me_blmer.model)$coefficients[3,1]
parameters::p_value_wald(me_blmer.model)[3,]

# Mixed-effect model expression
varlist1 <- lapply(blups.models_nextPP, function(f) summary(f)$call[2])
varlist1

# Correlation Coefficient
#estimate_AT <- lapply(blups.models_nextPP, function(f) summary(f)$coefficients[6,1])
estimate_RBD <- lapply(blups.models_nextPP, function(f) summary(f)$coefficients[2,1])
estimate_DAYS <- lapply(blups.models_nextPP, function(f) summary(f)$coefficients[3,1])


# Median
med <- lapply(blups.models_nextPP, function(f) summary(f)$residuals)
med_na <- lapply(med, function(f) na.exclude(f))
med_calc <- lapply(med_na, function(f) median(f))

## p-values for covariates (fixed-effect)
# ICU
test_pSS <- lapply(blups.models_nextPP, function(f) parameters::p_value_wald(f)[2,])
df.pvalueSS <- as.data.frame(test_pSS)
df.pvalueSS <-t(df.pvalueSS)
df.pvalueRBD = df.pvalueSS[seq(0, nrow(df.pvalueSS), 2), ]
# Days
test_pKOS <- lapply(blups.models_nextPP, function(f) parameters::p_value_wald(f)[3,])
df.pvalueKOS <- as.data.frame(test_pKOS)
df.pvalueKOS <-t(df.pvalueKOS)
df.pvalueDays = df.pvalueKOS[seq(0, nrow(df.pvalueKOS), 2), ]


# Prepare dataframe with extracted info. for downstream use
test_data = list(as.character(varlist1), med_calc, estimate_RBD, estimate_DAYS, as.numeric(t(df.pvalueRBD)), as.numeric(t(df.pvalueDays)))
names(test_data) <- c('cells', 'Median', 'Estimate_RBD', 'Estimate_Days', 'pvalueRBD', 'pvalueDays')

test_final <- as.data.frame(do.call(rbind, test_data))
t(test_final)
#setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table')
write.csv2(t(test_final), file='RBD_MEMtable.csv')

##################################
#          Input - mRNA          #
##################################

setwd('~/Documents/Mixed_Effects_Model/AT_rerun/')

# Read in file, build dataframe (columns should be covariates and genes, 
# row should be ID and subsequent info.)
me_memRNA = read.csv('deseq_ME12.csv', sep = ';', header = FALSE, row.names=1)
head(me_memRNA[1:7])
me_memRNA_df = as.data.frame(t(me_memRNA))
head(me_memRNA_df)

# List of genes names in order to loop through, omitting first 6 columns which 
# includes ID and covariates
mRNA_list = colnames(me_memRNA_df[-(1:6)])
head(mRNA_list)

## If wish to convert Ensembl to HUGO should do so at this moment

##################################
#     Convert Ensembl to HUGO    #
##################################

library(biomaRt)

mart <- biomaRt::useDataset("hsapiens_gene_ensembl", biomaRt::useMart("ensembl"))
genes <- mRNA_list
df<-me_memRNA_df[,-1]
head(df)
G_list <- biomaRt::getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
head(G_list)
write.table(as.data.frame(G_list),file="gene_list.csv", quote=F,sep=";",row.names=F)

df_0 <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("gene", "ensembl_gene_id")
colnames(df_0) <- x
df_0 <- merge(df,G_list,by.df_0="gene",by.df_0="ensembl_gene_id") 
head(df_0)
write.table(as.data.frame(df_0),file="mRNA_HUGO.csv", quote=F,sep=";",row.names=F)

# Read in HUGO table now, and continue as before
me_memRNA = read.csv('mRNA_HUGO.csv', sep = ';', header = FALSE, row.names=1)
me_memRNA_df = as.data.frame(t(me_memRNA))
head(me_memRNA_df)

# List of genes names 
mRNA_list = colnames(me_memRNA_df[-(1:6)])
head(mRNA_list)

##################################
#        Modelling - mRNA        #
##################################

# Mixed-effect modelling (Partial-bayesian)
varlist=mRNA_list 
blups.models_1 <- lapply(varlist, function(x) {
  mod2 = try(blmer(substitute(i ~ sex + age + symptom_score + Group + active_treatment + 
                   (1|id/Group), list(i = as.name(x))), 
                   data = me_memRNA_df1, na.action=na.exclude))
  if(isTRUE(class(mod2)=='try-error')) {return(NULL)} else{return(mod2)}
})

blups.models_next1 = as.list(blups.models_1)
# Remove NULL models that failed to converge
blups.models_next1[sapply(blups.models_next1, is.null)] <- NULL

##################################
# Extract info from model - mRNA #
##################################

# Mixed-effect model expression
varlist1 <- lapply(blups.models_next1, function(f) summary(f)$call[2])
head(varlist1)

# Correlation Coefficient
estimate_AT <- lapply(blups.models_next1, function(f) summary(f)$coefficients[6,1])
estimate_KOS <- lapply(blups.models_next1, function(f) summary(f)$coefficients[5,1])
estimate_SS <- lapply(blups.models_next1, function(f) summary(f)$coefficients[4,1])

# Median
med <- lapply(blups.models_next1, function(f) summary(f)$residuals)
med_na <- lapply(med, function(f) na.exclude(f))
med_calc <- lapply(med_na, function(f) median(f))

## p-values for covariates (fixed-effect)
# SS
test_pSS <- lapply(blups.models_next1, function(f) parameters::p_value_wald(f)[4,])
df.pvalueSS <- as.data.frame(test_pSS)
df.pvalueSS <-t(df.pvalueSS)
df.pvalueSS = df.pvalueSS[seq(0, nrow(df.pvalueSS), 2), ]
# KOS
test_pKOS <- lapply(blups.models_next1, function(f) parameters::p_value_wald(f)[5,])
df.pvalueKOS <- as.data.frame(test_pKOS)
df.pvalueKOS <-t(df.pvalueKOS)
df.pvalueKOS = df.pvalueKOS[seq(0, nrow(df.pvalueKOS), 2), ]
# AT
test_pAT <- lapply(blups.models_next1, function(f) parameters::p_value_wald(f)[6,])
df.pvalueAT <- as.data.frame(test_pAT)
df.pvalueAT <-t(df.pvalueAT)
df.pvalueAT = df.pvalueAT[seq(0, nrow(df.pvalueAT), 2), ]

# Prepare dataframe with extracted info. for downstream use
test_data = list(as.character(varlist1), med_calc, estimate_AT, estimate_KOS, estimate_SS, as.numeric(t(df.pvalueSS)), as.numeric(t(df.pvalueKOS)), as.numeric(t(df.pvalueAT)))
names(test_data) <- c('genes', 'Median', 'Estimate_AT', 'Estimate_KOS', 'Estimate_SS', 'pvalueSS', 'pvalueKOS', 'pvalueAT')

test_final <- as.data.frame(do.call(rbind, test_data))
setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table')
write.csv2(t(test_final), file='gene2_MEMtable_newSS.csv')
