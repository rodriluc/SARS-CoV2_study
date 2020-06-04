library(lme4)
library(ggplot2)
library(stargazer)
library(blme)
library(sjPlot)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

setwd('~/Documents/Helsinki_COVID19/serology/MEM')

##################################
#       Test trial for MEM       #
##################################

# Read in file, build dataframe (columns should be covariates and cell populations, 
# row should be ID and subsequent info.)
me_memG <- read.csv('MEM_RBD.csv', sep = ';')
head(me_memG)
me_memG_df = data.frame(me_memG)

## Construct Mixed-Effect model
# Linear MEM
me_lmer.model = lmer(IgD.pos..Memory.B ~ old.RBD + Days + 
                   (1|Subject_ID), data=me_mem1, REML=FALSE), data=me_memG)
# Partial Bayesian MEM
me_blmer.model = blmer(IgD.pos..Memory.B ~ old.RBD + Days + 
                   (1|Subject_ID), data=me_memG)

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
tab_model(me_lmer.model,me_mem4.model,me_mem20.model, file = 'tableMEM_trial.html')

stargazer(me_lmer.model,
          type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "", out = 'table_trial.html')

## Statistical significance (ex. SS)

#1) construct null model first
me_lmer.null = lmer(IgD.pos..Memory.B ~ old.RBD + Days + 
                   (1|Subject_ID), data=me_mem1, REML=FALSE)
#2) Re-do full model
me_lmer.model = lmer(IgD.pos..Memory.B ~ Days + 
                   (1|Subject_ID), data=me_mem1, REML=FALSE)
#3) Perform likelihood ratio test
anova(me_lmer.null, me_lmer.model)

# Looks at the coefficients of the model by subject and by item
coef(me_lmer.model)

##################################
#    Input and Modelling         #
##################################

# Read in file, build dataframe (columns should be covariates and plasma proteins/cell abundance, 
# row should be ID and subsequent info.)

me_mem1 <- read.csv('MEM_RBD.csv', sep = ';')
#write.csv(me_mem1, file = 'abundance_TRIAL_lvl3_meta.csv')
head(me_mem1)
me_mem1_df = data.frame(me_mem1)

# List of plasma protein/cell names in order to loop through, omitting first 6 columns which 
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
#  Extract info from model       #
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
# RBD
test_pRBD <- lapply(blups.models_nextPP, function(f) parameters::p_value_wald(f)[2,])
df.pvalueRBD <- as.data.frame(test_pRBD)
df.pvalueRBD <-t(df.pvalueRBD)
df.pvalueRBD = df.pvalueRBD[seq(0, nrow(df.pvalueRBD), 2), ]
# Days
test_pDays <- lapply(blups.models_nextPP, function(f) parameters::p_value_wald(f)[3,])
df.pvalueDays <- as.data.frame(test_pDays)
df.pvalueDays <-t(df.pvalueDays)
df.pvalueDays = df.pvalueDays[seq(0, nrow(df.pvalueDays), 2), ]


# Prepare dataframe with extracted info. for downstream use
test_data = list(as.character(varlist1), med_calc, estimate_RBD, estimate_DAYS, as.numeric(t(df.pvalueRBD)), as.numeric(t(df.pvalueDays)))
names(test_data) <- c('cells', 'Median', 'Estimate_RBD', 'Estimate_Days', 'pvalueRBD', 'pvalueDays')

test_final <- as.data.frame(do.call(rbind, test_data))
t(test_final)
#setwd('~/Documents/Mixed_Effects_Model/AT_rerun/MEM_table')
write.csv2(t(test_final), file='RBD_MEMtable.csv')
