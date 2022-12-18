# This code fits a fractional polynomial mixed model for the Dialyzer data set,
# as presented in Garcia and Trinca ...

# We generaly use lme to fit the mixed models bu occasionally will use lme4 
# as well

set.seed(28102022)
library(nlme)

#-------------------------------------------------------------------------------
# Functions required
# Fit mixed FP models (mfp algorithm)
source("MFP_function_With_RandomSlops.R", encoding = 'UTF-8')  

# We will comparing a few tests including permutation tests 
# F test for RE (Demidenko)
source("FTest_RandomTerms_Function.R", encoding = 'UTF-8')  

# Permutation tests for RE 
# The routine for permutation tests are from Lee and Brown from 
# http://www-personal.umich.edu/~tombraun/LeeBraun/Permutation%20Test%20for%20Random%20Slope.R
# It uses lmer instead of lme
# We have modified the original in order to make it a bit more general
source("PermuTest_function.R", encoding = 'UTF-8')  

# Graph diagnostic tools
# Function for mixed model diagnostic tools 
# Downloaded from http://www.ime.usp.br/~jmsinger/lmmdiagnostics.zip
# Remove option LINPACK = F in line 82 
# Remove option contrasts.arg around line 110

source("residdiag_nlmeM.r")
#-------------------------------------------------------------------------------

# Plot for the response profiles
#---------------------------------------------------------------------------------
X11()
plot(Dialyzer, aspect = "fill", xlab="Transmembrane pressure (dmHg)", 
     ylab="Ultarfiltratio rate (ml/hr)", key = F, displayLevel="QB", 
     outer= ~ QB)

#-------------------------------------------------------------------------------
# Step 1: FP using a minimal variance-covariance structure for the random part
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Our code is not general enough, so please, **** name your grouped data frame 
# as "dadosag" **** before fitting the models and calling function MFP 
#-------------------------------------------------------------------------------

dadosag = Dialyzer
dadosag$QB = as.numeric(dadosag$QB)-1

# Dialyzer has "Subject" as ordered. We prefer to work with un-ordered Subject
# so that it is easier to keep track of the several matrices are being formed
# Here we generate a new column for Subject and group the data with 
# the order.group=FALSE

dadosag = cbind(Subject=factor(rep(1:20, rep(7, 20))), dadosag)
dadosag = groupedData(rate~pressure|Subject, order.groups= FALSE, data=dadosag)

# ------------------------------------------------------------------------------
# Pinheiro et al.(2000) - Conventional Polynomial Mixed Model
# ------------------------------------------------------------------------------

mod.CP = lme(rate ~ pressure + I(pressure^2) + I(pressure^3) + I(pressure^4) + 
             QB + pressure:QB, 
             random = pdSymm(~ 1 + pressure + I(pressure^2)), 
             weights = varPower(form = ~ pressure), method="ML", data=dadosag)
round(summary(mod.CP)$tTable, dig=2)
#-------------------------------------------------------------------------------

lmmLin = lme(rate ~ QB + pressure, method = "ML", random = ~1|Subject, 
             data=dadosag)
lmmNull = update(lmmLin, .~. -pressure)
RE1 = lmeStruct(reStruct(~1|Subject), corStruct=NULL, varStruct=varIdent())
RE2 = lmeStruct(reStruct(~1|Subject), corStruct=NULL, varStruct=varIdent())
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$pressure, RE1, RE2)

# Model_selected "FP2 Model"  "1" "1" ==> pressure +I(pressure*log(pressure))

# Fit the best FP model found above
lmmFP2 = lme(rate ~ QB + pressure + I(pressure*log(pressure)), method="ML",
             random= ~1|Subject, data=dadosag)
summary(lmmFP2)

#---------------------------------------------------------------------------------
# Including interactions (algorithm MFPI)
#-------------------------------------------------------------------------------

lmmFP2Main = lmmFP2

# Full model 
lmmFP2Full = update(lmmFP2Main, .~. + QB:pressure + QB:I(pressure*log(pressure)))
summary(lmmFP2Full)

# Global test for interactions
anova(lmmFP2Main, lmmFP2Full)

# Both interactions will be kept for the search of the best random structure

# ------------------------------------------------------------------------------
# Step 2: Selection of random effects structure
# ------------------------------------------------------------

# -------------------------------------------
# Fixed Effects Model F test
# -------------------------------------------

modfix = lm(rate ~ QB + pressure + I(pressure*log(pressure)) + QB:pressure + 
            QB:I(pressure*log(pressure)), data=dadosag)

# Subject effect (Random intercept)
modb0 = lm(rate ~ Subject + pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
             QB:I(pressure*log(pressure)), data=dadosag)

# Random slop for pressure
modb1 = lm(rate ~ Subject + pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
             QB:I(pressure*log(pressure)) + Subject:pressure, data=dadosag)

# Random intercept and both slops pressure and I(pressure*log(pressure))
modb2 = lm(rate ~ Subject + pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
             QB:I(pressure*log(pressure)) + Subject:pressure + 
             Subject:I(pressure*log(pressure)), data=dadosag)

Fanova = rbind(unlist(anova(modfix, modb2)[2,5:6]),
               unlist(anova(modfix, modb0)[2,5:6]),
               unlist(anova(modb0, modb1)[2,5:6]),
               unlist(anova(modb1, modb2)[2,5:6]))

colnames(Fanova) = c("F_FixedEff_Model", "Pr(>F)")
Fanova[,1] = round(Fanova[,1], dig=3)
Fanova = cbind(rbind(anova(modfix, modb2)[2,c(3,1)],
                             anova(modfix, modb0)[2,c(3,1)],
                             anova(modb0, modb1)[2,c(3,1)],
                             anova(modb1, modb2)[2,c(3,1)]), Fanova)
rownames(Fanova) = c("GlobalRandom", "Random_b0", "Random_b1", "Random_b2")                            
Fanova               
# ------------------------------------------------------------------------------
# Exact F test: Demidenko
# ------------------------------------------------------------------------------
library(MASS)
modfix = lm(rate ~ QB + pressure + I(pressure*log(pressure)) + QB:pressure + 
              QB:I(pressure*log(pressure)), data=dadosag)
modrandom = lm(rate ~ -1 + Subject + Subject:pressure + Subject:I(pressure*log(pressure)), data=dadosag)
FG = Ftest_Demidenko(modfix, modrandom, dadosag, dadosag$rate)

modrandom = lm(rate ~ -1 + Subject, data=dadosag)
Fb0 = Ftest_Demidenko(modfix, modrandom, dadosag, dadosag$rate)

modfix = lm(rate ~ Subject + QB + pressure + I(pressure*log(pressure)) + QB:pressure + 
              QB:I(pressure*log(pressure)), data=dadosag)
modrandom = lm(rate ~ -1 + Subject:pressure, data=dadosag)
Fb1 = Ftest_Demidenko(modfix, modrandom, dadosag, dadosag$rate)

modfix = lm(rate ~ Subject + Subject:pressure + QB + pressure + 
              I(pressure*log(pressure)) + QB:pressure + 
              QB:I(pressure*log(pressure)), data=dadosag)
modrandom = lm(rate ~ Subject:I(pressure*log(pressure)), data=dadosag)
Fb2 = Ftest_Demidenko(modfix, modrandom, dadosag, dadosag$rate)

FDemi = rbind(FG, Fb0, Fb1, Fb2)

cbind(Fanova, FDemi)

# The data is balanced in terms of Subject sizes, the tests considering the 
# fixed effects model and Demidenko Exact tests are  equivaletnt.

#-------------------------------------------------------------------------------
# LRT + adjusting p-values by the mixed chisquares (0.5 x 0.5 for 
# df=4 and df=5 (global) and df=0 and df=1, 1 and 2, 2 an 3 (sequential))
#-------------------------------------------------------------------------------

modf <- gls(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
             QB:I(pressure*log(pressure)), data=dadosag)

# Random intercept
modb0 = lme(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
             QB:I(pressure*log(pressure)), random=~1|Subject, data=dadosag)

# Random slop for pressure
modb1 = lme(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
             QB:I(pressure*log(pressure)), random=~1+pressure|Subject, data=Dialyzer)

# Random intercept and both slops pressure and I(pressure*log(pressure))
controlLME=list(opt =c("nlminb", "optim"), 
                optimMethod = "BFGS",
                maxIter=1000, msMaxIter=1000, 
                niterEM=500, msMaxEval=1000)
modb2 = lme(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
             QB:I(pressure*log(pressure)),  
             random=~1+pressure+I(pressure*log(pressure))|Subject, 
             control=do.call(lmeControl, controlLME), data=dadosag)

LRT = rbind(GlobalRandom = c(DF=5, unlist(anova(modf, modb2)[2,8:9])),
               RandomInt = c(1, unlist(anova(modf, modb0)[2,8:9])),
               RandomX = c(2, unlist(anova(modb0, modb1)[2,8:9])),
               RandomX2 =c(3, unlist(anova(modb1, modb2)[2,8:9])))

# Adjusting p-values for mixtures
library(emdbook)
my_chibar = function(x)
{
  p = x[1]; df =  x[2]
  pchibarsq(p, df = df, mix = 0.5, lower.tail=FALSE, log.p = FALSE)
}

pv_adj=apply(LRT[,c(2,1)], 1, my_chibar)

#-------------------------------------------------------------------------------
# Results of usual tests for random effects
#-------------------------------------------------------------------------------
cbind(cbind(LRT, pv_adj, Fanova), FDemi)

#-------------------------------------------------------------------------------
# Permutation tests
# Lee & Braun (2012)
# BLUP based permutation
# LR statistic based permutation
#-------------------------------------------------------------------------------
# The routine uses lme4 instead of nlme

library("lme4")
library(tidyverse)

set.seed(56575)
# Define the number of subjects and the number of observations per subject contained in dat.
numsub = nlevels(dadosag$Subject)

Nperms = 1000

# ------------------------------------------------------------------------------
# inclusion of random Intercept (b0)
# ------------------------------------------------------------------------------

# Test for the intercept, only BLUP based permutation
# Fitting the alternative and null models.

fit0 = lm(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
            QB:I(pressure*log(pressure)), data=dadosag)
fit1 = lmer(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
              QB:I(pressure*log(pressure)) + (1|Subject),
            data=dadosag) 

randpart0 = "1"
randpart1 = "(1|Subject)" 

controlMER = list(optimizer ='nloptwrap', optCtrl=list(method='NLOPT_LN_COBYLA', 
                                                       maxiter=500))

# Only BLUP version is performed
pValuesIntercept <- PermutationTest(Nperms, dadosag$rate, positionranef=1, 
                                    numsub, fit0, fit1, randpart0, randpart1, 
                                    controlMER, dadosag)

# ------------------------------------------------------------------------------
# inclusion of random slop for pressure (b1)
# ------------------------------------------------------------------------------

# Fitting the alternative and null models.
fit1 = fit1
fit2 = lmer(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
            QB:I(pressure*log(pressure)) + (1+pressure|Subject),
            data=dadosag) 

randpart0 = "(1|Subject)" 
randpart1 = "(1 + pressure|Subject)"

controlMER = list(optimizer ='nloptwrap', optCtrl=list(method='NLOPT_LN_COBYLA', 
                                                       maxiter=1000))

pValuesSlop1 <- PermutationTest(Nperms, dadosag$rate, positionranef=2, 
                                numsub, fit1, fit2, randpart0, randpart1, 
                                controlMER, dadosag)

# ------------------------------------------------------------------------------
# inclusion of random slop for I(pressure*log(pressure)) (b2)
# ------------------------------------------------------------------------------

# Fitting the alternative and null models.
fit2 = fit2
fit3 = lmer(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
            QB:I(pressure*log(pressure)) + (pressure + I(pressure*log(pressure))|Subject),
            control = do.call(lmerControl, controlMER),
            #   lmerControl(
            # optimizer ='optimx', optCtrl=list(method='L-BFGS-B')),
            data=dadosag) 

randpart0 = "(1+pressure|Subject)"
randpart1 = "(1 + pressure + I(pressure*log(pressure))|Subject)"
pValuesSlop2 = PermutationTest(Nperms, dadosag$rate, positionranef=3, 
                               numsub, fit2, fit3, randpart0, randpart1, 
                               controlMER, dadosag)

# P-values and NA's for the test statistics in the permutations
pValuesPerm <- rbind(b0=pValuesIntercept, b1=pValuesSlop1, b2=pValuesSlop2)
pValuesPerm

#-------------------------------------------------------------------------------
# STEP 2 Conclusions: all random terms contributes to the modelling
#-------------------------------------------------------------------------------

lmmFP2Step2 = lme(rate ~ pressure + I(pressure*log(pressure)) + QB + QB:pressure + 
                    QB:I(pressure*log(pressure)),  
                    random=~pressure + I(pressure*log(pressure))|Subject,
                    control = lmeControl(opt ='optim', optimMethod='BFGS'),
                    data=dadosag) 
summary(lmmFP2Step2)
round(summary(lmmFP2Step2)$tTable, dig=3)
VarCorr(lmmFP2Step2) 

# -------------------------------------------------------------------------
# STEP 3 - Variance-covariance structure for the error
# -------------------------------------------------------------------------

# -------------------------------------------------------------
# Graph diagnostic tools
# -------------------------------------------------------------

plot(lmmFP2Step2)

residStep2 = residdiag.nlme(fit=lmmFP2Step2,limit=2,plotid=1:6)

# As already explored in Pinheiro and Bates there some indication of variance 
# heterogeneity

#----------------------------------------------------------------------------------
# Selecting a variance function (REML)
#----------------------------------------------------------------------------------

lmmFP2_h1 = update(lmmFP2Step2, weights = varPower(form = ~ pressure),
                   control = lmeControl(opt ='nlminb')
                   )
summary(lmmFP2_h1)
anova(lmmFP2Step2, lmmFP2_h1)

residStep3 = residdiag.nlme(fit=lmmFP2_h1,limit=2,plotid=1:6)

# The graphs do not highlight problems with the assuptions
# Variance heterogeneity is well accommodated by the varPower function

lmmFP2Step3 = lmmFP2_h1

#----------------------------------------------------------------------------------
# STEP 4: Once a reasonable structure for the variance-covariance is found,
# we go back to the mean model
#----------------------------------------------------------------------------------

lmmFP2_h1ML = update(lmmFP2Step3, .~. -QB:pressure - QB:I(pressure*log(pressure)), 
                     method = "ML")

lmmLin = update(lmmFP2_h1ML, fixed=rate ~ QB + pressure)
controlLME=list(opt =c("nlminb", "optim"), 
           optimMethod = "BFGS",
           maxIter=10000, msMaxIter=10000, 
           niterEM=500, msMaxEval=10000)
lmmNull = update(lmmLin, .~. -pressure, 
                 control = do.call(lmeControl, controlLME))

RE1 = lmeStruct(reStruct(~x1|Subject),
                corStruct=NULL, varStruct=varPower(form = ~ pressure))
RE2 = lmeStruct(reStruct(~x1 + x2|Subject),
                corStruct=NULL, varStruct=varPower(form = ~ pressure))
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$pressure, RE1, RE2)

# $Model_selected: "FP2 Model"         "1"         "1" 

# Same FP2 is maintained!

#---------------------------------------------------------------------------------
# Including interactions (algorithm MFPI)
#-------------------------------------------------------------------------------

lmmFP2Main = lmmFP2_h1ML

# Full model 
lmmFP2Full = update(lmmFP2Main, .~. + QB:pressure + QB:I(pressure*log(pressure)))
summary(lmmFP2Full)

# Global test for interactions
anova(lmmFP2Main, lmmFP2Full)

# Hierarchical anova
anova(lmmFP2Full)

# The second fixed interaction term is not relevant

lmmFP2Step4 <- update(lmmFP2Full, .~. - QB:I(pressure*log(pressure)), method="REML")
round(summary(lmmFP2Step4)$tTable, dig=2)
VarCorr(lmmFP2Step4)
round(summary(lmmFP2Step4)$AIC, dig=2)
round(summary(lmmFP2Step4)$BIC, dig=2)
round(summary(lmmFP2Step4)$logLik, dig=2)

# -------------------------------------------------------------------------
# Final diagnostic graphs
# -------------------------------------------------------------------------

residSTep4 = residdiag.nlme(fit=lmmFP2Step4,limit=3,plotid=1:6)

