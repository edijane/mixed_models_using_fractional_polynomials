# This code fits a fractional polynomial mixed model for the Boston House data set,
# as presented in Garcia and Trinca ...             

# Package to fit the models
library(nlme)

# Package to get the corrected data set
library(mlbench)

# Function for the MFP algorithm for mixed models
source("MFP_function_No_RandomSlop.R")
#-------------------------------------------------------------------------------
# Function for mixed model diagnostic tools 
# Downloaded from http://www.ime.usp.br/~jmsinger/lmmdiagnostics.zip
# Remove option LINPACK = F in line 82 
source("residdiag_nlme.R") 
#-------------------------------------------------------------------------------
set.seed(25102022)
data(BostonHousing2)
dados = BostonHousing2
count = table(factor(dados$town, levels=unique(dados$town)))
dados$townid = rep(1:92, count)
dados$Y = log(dados$cmed)   # response variable

# re-scaling some variables 
dados$taxt = dados$tax/100
omega_crim =  0.2+((1-0.2)*(dados$crim - min(dados$crim))/(max(dados$crim) - min(dados$crim)))
dados$omegacrim = omega_crim
omega_zn   =  0.2+((1-0.2)*(dados$zn - min(dados$zn))/(max(dados$zn) - min(dados$zn)))
dados$omegazn = omega_zn

# -------------------------------------------------------------
# Fitting a preliminary LMM
# -------------------------------------------------------------
#-------------------------------------------------------------------------------
# Our code is not general enough, so please, name your grouped data frame 
# as "dadosag" before fitting the models and calling function MFP 
#-------------------------------------------------------------------------------

dadosag = groupedData(Y ~ lstat|townid, data = dados)

lmm0REML = lme(Y ~ lstat+rm+crim+nox+rad+dis+taxt+ptratio+age+zn+indus+chas, 
               method = "REML", random = ~1|townid, data=dadosag)
summary(lmm0REML)
VarCorr(lmm0REML)   # variance component estimates (and SD)
AIC(lmm0REML)
BIC(lmm0REML)
logLik(lmm0REML)

# -------------------------------------------------------------
###### DIAGNOSTIC tools -  lmm0REML
# load residdiag.nlme function (Singer et al., 2017)

X11()
resid0 = residdiag.nlme(fit=lmm0REML,limit=3,plotid=1:16)
#

# -------------------------------------------------------------
# Step 1: "Best" model for the fixed effects 
# -------------------------------------------------------------
# Algorithm MFP
# -------------------------------------------------------------
# 1) Use pre-transformed covariates omegacrim, taxt and omegazn 
# 2) Estimation method: ML

lmm0 = lme(Y ~ lstat+rm+omegacrim+nox+rad+dis+taxt+ptratio+age+omegazn+indus+chas, 
               method = "ML", random = ~1|townid, data=dadosag)

# Covariate visiting order
A = data.frame(round(summary(lmm0)$tTable, dig=3))
A = A[order(A$p.value),]
A

library(xtable)
xtable(A, digits = 3)   # Latex table

# Value Std.Error  DF t-value p-value
# 
# lstat       -0.023     0.002 407 -12.606   0.000
# rm           0.131     0.014 407   9.147   0.000
# omegacrim   -0.770     0.117 407  -6.567   0.000
# nox         -0.817     0.171 407  -4.789   0.000
# rad          0.019     0.005  86   3.942   0.000
# dis         -0.041     0.011 407  -3.683   0.000
# taxt        -0.078     0.025  86  -3.064   0.003
# ptratio     -0.029     0.011  86  -2.684   0.009
# age         -0.001     0.000 407  -2.412   0.016
# omegazn      0.130     0.098  86   1.330   0.187
# indus        0.003     0.005  86   0.598   0.552
# chas        -0.009     0.030 407  -0.303   0.762

# -------------------------------------------------------------
# 1st CYCLE
# -------------------------------------------------------------

# load function MFP

### 1) LSTAT (1st CYCLE)
lmmLin = lme(Y ~ lstat + rm + omegacrim + nox + rad + dis + taxt + ptratio + 
               age + omegazn + indus + chas,
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -lstat)

MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$lstat)

# Model_selected "FP1 Model" "0.5" =>  lstat^0.5

### 2) RM (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + rm + omegacrim + nox + rad + dis + taxt + ptratio + 
               age + omegazn + indus + chas,
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -rm)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$rm)

# Model_selected "FP2 Model" "0.5" "0.5"  =>   rm^0.5 + (rm^0.5)*log(rm)

### 3) OMEGACRIM (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + omegacrim + 
               nox + rad + dis + taxt + ptratio + age + omegazn + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -omegacrim)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$omegacrim)

# Model_selected: "FP1 Model" "-0.5"  => omegacrim^-0,5

### 4) NOX (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + dis + taxt + ptratio + age + omegazn + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -nox)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$nox)

# Model_selected "Linear Model" => nox

### 5) RAD (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + dis + taxt + ptratio + age + omegazn + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -rad)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$rad)

# Model_selected "Linear Model" => rad

### 6) DIS (1o CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + dis + taxt + ptratio + age + omegazn + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -dis)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$dis)

# Model_selected "Model FP1" "-0.5" => dis^-0.5

### 7) TAXT (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + age + omegazn + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -taxt)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$taxt)

# Model_selected "Linear Model" => taxt

### 8) PTRATIO (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + age + omegazn + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -ptratio)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$ptratio)

# Model_selected "Linear Model" => ptratio

### 9) AGE (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + age + omegazn + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -age)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$age)

# Model_selected: "Null Model" => age removed!

### 10) OMEGAZN (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + omegazn + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -omegazn)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$omegazn)

# Model_selected: "Null Model" => omegazn removed!

### 11) INDUS (1st CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + indus + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -indus)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$indus)

# Model_selected: "Null Model" => indus removed!

### 12) CHAS (qualitative variable) - 1st CYCLE
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + chas, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -chas)

# LRT
logLikdif.chas = -2*logLik(lmmNull) - (-2*logLik(lmmLin))
1 - pchisq(logLikdif.chas, 1)
anova(lmmLin, lmmNull)

# Model_selected: chas removed!!

# ------------------------------------------------------------------- 
# Final model at end of 1st cycle includes
# ------------------------------------------------------------------- 
# lstat^0.5, rm^0.5, rm^0.5*log(rm), omegacrim^-0.5, nox, rad, dis^-0.5, taxt e ptratio
# -------------------------------------------------------------
# 2nd CYCLE
# -------------------------------------------------------------
### 1) LSTAT (2nd CYCLE)
lmmLin = lme(Y ~ lstat + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
             nox + rad + I(dis^-0.5) + taxt + ptratio, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -lstat)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$lstat)

# Modelo_selecionado, lstat^0.5

### 2) RM (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + rm + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio, 
              method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -rm)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$rm)

# Model_selected "FP2 Model" "0.5" "0.5" => rm^0.5+rm^0.5*log(rm)

### 3) OMEGACRIM (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + omegacrim +
               nox + rad + I(dis^-0.5) + taxt + ptratio, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -omegacrim)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$omegacrim)

# Model_selected: "FP1 Model" "-0.5"  => omegacrim^-0.5 

### 4) NOX (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -nox)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$nox)

# Model_selected: Linear => nox

### 5) RAD (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -rad)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$rad)

# Model_selected: Linear => rad

### 6) DIS (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + dis + taxt + ptratio, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -dis)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$dis)

# Model_selected "Model FP1" "-0.5" => dis^-0.5

#-------------------------------------------------------------------------------
# THE PROCESS CAN FINISH HERE SINCE THERE WAS NO CHANGE IN THE PREVIOUS POWERS
# We leave the code below for completeness
#-------------------------------------------------------------------------------
### 7) TAXT (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -taxt)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$taxt)

# Model_selected: Linear => taxt

### 8) PTRATIO (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio, 
               method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -ptratio)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$ptratio)

# Model_selected: Linear => ptratio

### 9) AGE (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + age, 
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -age)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$age)

# Model_selected, Null model

### 11) omegazn (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + omegazn,
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -omegazn)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$omegazn)

# Model_selected, Null model

### 11) INDUS (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + indus,
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -indus)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$indus)

# Model_selected, Null model

### 13) CHAS (qualitative) (2nd CYCLE)
lmmLin = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
               nox + rad + I(dis^-0.5) + taxt + ptratio + chas,
             method = "ML", random = ~1|townid, data=dadosag)
lmmNull = update(lmmLin, .~. -chas)

# LRT
logLikdif.chas = -2*logLik(lmmNull) - (-2*logLik(lmmLin))
1 - pchisq(logLikdif.chas, 1)
anova(lmmLin, lmmNull)

# Model_selected, Null model

# ------------------------------------------------------------------- 
# Model at the end of 2nd CYCLE = final model at the end of 1st CYCLE, 
# end of MFP algorithm 
# ------------------------------------------------------------------- 

lmmFP2Step1 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio,
              method = "REML", random = ~1|townid, data=dadosag)

# Final model in Step 1
round(summary(lmmFP2Step1)$tTable, dig=3)
VarCorr(lmmFP2Step1) 
c(AIC(lmmFP2Step1), BIC(lmmFP2Step1), logLik(lmmFP2Step1))
# ------------------------------------------------------------------------------
# STEP 2: Random effect structure for the  MMPF
# ------------------------------------------------------------------------------

# The minimal random structure is the inclusion of random intercept that we think 
# should be in the model due to the structure of the data = grouped observations 
# within town. We are not really interested in hypothesis test for the variance 
# parameter, however, the procedures below give indication this parameter 
# contributes to the modeling.

mmpf1 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio,random = ~1|townid, data=dadosag)

mmpf1.1 = gls(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio, data=dadosag)

# LRT
VarTest = anova(mmpf1, mmpf1.1)
pchisq(VarTest$L.Ratio[2], df = 1, lower.tail=FALSE)

#        Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#mmpf1       1 12 -430.6105 -380.1316 227.3053                        
#mmpf1.1     2 11 -297.0647 -250.7924 159.5324 1 vs 2 135.5458  <.0001

# p-value correction using chi-square mixture (0 and 1 dfs, weights=0.5)
library(MASS)
library(emdbook)
pchibarsq(VarTest$L.Ratio[2], df = 1, mix = 0.5, lower.tail=FALSE, log.p = FALSE)
# 1.423777e-31

# Test F
modfix = lm(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
            nox + rad + I(dis^-0.5) + taxt + ptratio, data=dadosag)

modrandom = lm(Y ~ -1 + townid, data=dadosag)
Fb0 = Ftest_Demidenko(modfix, modrandom, dadosag, dadosag$Y)
# Fb0
# dfn         dfd          F_      pvalue 
# 88          408       5.269   2.88538e-31

#-------------------------------------------------------------------------------
# Permutation tests
# Lee & Braun (2012)
# BLUP based permutation
# LR statistic based permutation
#-------------------------------------------------------------------------------
# The routine uses lme4 instead of nlme

library("lme4")

set.seed(6547654)
numsub = nlevels(dadosag$townid)
Nperms = 1000
# ------------------------------------------------------------------------------
# inclusion of random Intercept (b0)
# ------------------------------------------------------------------------------

controlMER = list(optimizer ='nloptwrap', optCtrl=list(method='NLOPT_LN_COBYLA', 
                                                       maxiter=500))

# Test for the intercept, only BLUP based permutation
# Fitting the alternative and null models.
fit0 = lm(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
            nox + rad + I(dis^-0.5) + taxt + ptratio, data=dadosag)
fit1 = lmer(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio + (1|townid), data=dadosag) 

randpart0 = "1"
randpart1 = "(1|townid)" 
pValuesIntercept = PermutationTest(Nperms, dadosag$Y, numsub, positionranef=1, fit0, fit1, 
                                    randpart0, randpart1, controlMER, dadosag)

# Final model in STEP 2 = Final model in STEP 1
lmmFP2Step2 = lmmFP2Step1 


# ----------------------------------------------------------
# STEP 3: Variance structure for the random error
# ----------------------------------------------------------

# Use of Modifyed Lesaffre-Verbeke Index (MLVI) to explore variance 
# heterogeneity among towns
set.seed(25102022)
mmpf1 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) +  
              nox + rad + I(dis^-0.5) + taxt + ptratio, method = "REML", 
            random = ~1|townid, data=dadosag)

residmmpf1 = residdiag.nlme(fit=mmpf1,limit=2,plotid=1)
#
LS = residmmpf1$"lesaffreverbeke.measure"    
round(LS[order(LS[,2], decreasing = F), ],dig=2)

# MLVI >= 2.0 suggests towns with var-cov matrix not adequate
#               MLVI
# [81,]      79 2.08
# [82,]      83 2.45
# [83,]      89 3.07
# [84,]      77 4.35
# [85,]      68 4.50
# [86,]      76 4.61
# [87,]      33 4.92
# [88,]      81 5.71
# [89,]      80 6.23
# [90,]      69 6.69
# [91,]      67 6.92
# [92,]      78 7.80


# The idea is forming clusters of towns with the same vcov matrix,
# refitting the model, recalculate the index and relocate towns that flag out
# Lets try forming 2 clusters: grupo 1: MLVI <=2 grupo 2: MLVI > 2

dados$grupo = NULL
dados$grupo = 1
dados$grupo[dados$townid==78] = 2
dados$grupo[dados$townid==67] = 2
dados$grupo[dados$townid==69] = 2
dados$grupo[dados$townid==80] = 2
dados$grupo[dados$townid==81] = 2

dados$grupo[dados$townid==33] = 2
dados$grupo[dados$townid==76] = 2
dados$grupo[dados$townid==68] = 2
dados$grupo[dados$townid==77] = 2
dados$grupo[dados$townid==89] = 2
dados$grupo[dados$townid==83] = 2
dados$grupo[dados$townid==79] = 2

dados$group = factor(dados$grupo)
dadosag = groupedData(Y ~ lstat|group, data = dados)

mmpf2 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio, method = "REML", 
            random = ~1|townid, correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
            weights = varIdent(form = ~ 1 | group), 
            control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), data=dadosag)
summary(mmpf2)
# Not interested in test, only comparing AIC, BIC
# The clustering seems to improve
anova(mmpf1, mmpf2)

round(summary(mmpf2)$tTable, dig=3)
VarCorr(mmpf2)  

residmmpf2 = residdiag.nlme(fit=mmpf2,limit=2,plotid=1)
#
LS = residmmpf2$"lesaffreverbeke.measure"    
round(LS[order(LS[,2], decreasing = F), ],dig=2)

# Other towns flag out
# [87,]      66 2.23
# [88,]      87 2.38
# [89,]      77 2.72
# [90,]      76 2.92
# [91,]      82 4.62
# [92,]      90 8.39
dados$grupo = NULL
dados$grupo = 1
dados$grupo[dados$townid==78] = 2
dados$grupo[dados$townid==67] = 2
dados$grupo[dados$townid==69] = 2
dados$grupo[dados$townid==80] = 2
dados$grupo[dados$townid==81] = 2

dados$grupo[dados$townid==33] = 2
dados$grupo[dados$townid==76] = 2
dados$grupo[dados$townid==68] = 2
dados$grupo[dados$townid==77] = 2
dados$grupo[dados$townid==89] = 2
dados$grupo[dados$townid==83] = 2
dados$grupo[dados$townid==79] = 2

# grupo 3
dados$grupo[dados$townid==90] = 3
dados$grupo[dados$townid==82] = 3
dados$grupo[dados$townid==76] = 3
dados$grupo[dados$townid==77] = 3
dados$grupo[dados$townid==87] = 3
dados$grupo[dados$townid==66] = 3

dados$group = factor(dados$grupo)
dadosag = groupedData(Y ~ lstat|group, data = dados)

mmpf3 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio, method = "REML", 
            random = ~1|townid, correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
            weights = varIdent(form = ~ 1 | group), 
            control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), data=dadosag)
summary(mmpf3)
anova(mmpf1, mmpf2, mmpf3)
round(summary(mmpf3)$tTable, dig=3)
VarCorr(mmpf3)  

residmmpf3 = residdiag.nlme(fit=mmpf3,limit=2,plotid=1)
#
LS = residmmpf3$"lesaffreverbeke.measure"   
round(LS[order(LS[,2], decreasing = F), ],dig=2)

# [91,]      86 2.20
# [92,]      91 2.48

dados$grupo = NULL
dados$grupo = 1
dados$grupo[dados$townid==78] = 2
dados$grupo[dados$townid==67] = 2
dados$grupo[dados$townid==69] = 2
dados$grupo[dados$townid==80] = 2
dados$grupo[dados$townid==81] = 2

dados$grupo[dados$townid==33] = 2
dados$grupo[dados$townid==76] = 2
dados$grupo[dados$townid==68] = 2
dados$grupo[dados$townid==77] = 2
dados$grupo[dados$townid==89] = 2
dados$grupo[dados$townid==83] = 2
dados$grupo[dados$townid==79] = 2

# grupo 3
dados$grupo[dados$townid==90] = 3
dados$grupo[dados$townid==82] = 3
dados$grupo[dados$townid==76] = 3
dados$grupo[dados$townid==77] = 3
dados$grupo[dados$townid==87] = 3
dados$grupo[dados$townid==66] = 3

# grupo 4
dados$grupo[dados$townid==91] = 4
dados$grupo[dados$townid==86] = 4


dados$group = factor(dados$grupo)
dadosag = groupedData(Y ~ lstat|group, data = dados)

mmpf4 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio, method = "REML", 
            random = ~1|townid, correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
            weights = varIdent(form = ~ 1 | group), 
            control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), data=dadosag)
summary(mmpf4)
anova(mmpf1, mmpf2, mmpf3, mmpf4)
round(summary(mmpf4)$tTable, dig=3)
VarCorr(mmpf4)  

x11()
residmmpf4 = residdiag.nlme(fit=mmpf4,limit=2,plotid=1:8)
#
LS = residmmpf4$"lesaffreverbeke.measure"    
round(LS[order(LS[,2], decreasing = F), ],dig=2)

# [91,]      29 2.10
# [92,]      85 2.21

dados$grupo = NULL
dados$grupo = 1
dados$grupo[dados$townid==78] = 2
dados$grupo[dados$townid==67] = 2
dados$grupo[dados$townid==69] = 2
dados$grupo[dados$townid==80] = 2
dados$grupo[dados$townid==81] = 2

dados$grupo[dados$townid==33] = 2
dados$grupo[dados$townid==76] = 2
dados$grupo[dados$townid==68] = 2
dados$grupo[dados$townid==77] = 2
dados$grupo[dados$townid==89] = 2
dados$grupo[dados$townid==83] = 2
dados$grupo[dados$townid==79] = 2

# grupo 3
dados$grupo[dados$townid==90] = 3
dados$grupo[dados$townid==82] = 3
dados$grupo[dados$townid==76] = 3
dados$grupo[dados$townid==77] = 3
dados$grupo[dados$townid==87] = 3
dados$grupo[dados$townid==66] = 3

# grupo 4
dados$grupo[dados$townid==91] = 4
dados$grupo[dados$townid==86] = 4

# grupo 5

dados$grupo[dados$townid==29] =5
dados$grupo[dados$townid==85] =5

dados$group = factor(dados$grupo)
dadosag = groupedData(Y ~ lstat|group, data = dados)

mmpf5 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio, method = "REML", 
            random = ~1|townid, correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
            weights = varIdent(form = ~ 1 | group), 
            control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), data=dadosag)
summary(mmpf5)
anova(mmpf1, mmpf2, mmpf3, mmpf4, mmpf5)
round(summary(mmpf5)$tTable, dig=3)
VarCorr(mmpf5)  

residmmpf5 = residdiag.nlme(fit=mmpf5,limit=2,plotid=1)
#
LS = residmmpf5$"lesaffreverbeke.measure"    
round(LS[order(LS[,2], decreasing = F), ],dig=2)

dados$grupo = NULL
dados$grupo = 1
dados$grupo[dados$townid==78] = 2
dados$grupo[dados$townid==67] = 2
dados$grupo[dados$townid==69] = 2
dados$grupo[dados$townid==80] = 2
dados$grupo[dados$townid==81] = 2

dados$grupo[dados$townid==33] = 2
dados$grupo[dados$townid==76] = 2
dados$grupo[dados$townid==68] = 2
dados$grupo[dados$townid==77] = 2
dados$grupo[dados$townid==89] = 2
dados$grupo[dados$townid==83] = 2
dados$grupo[dados$townid==79] = 2

# grupo 3
dados$grupo[dados$townid==90] = 3
dados$grupo[dados$townid==82] = 3
dados$grupo[dados$townid==76] = 3
dados$grupo[dados$townid==77] = 3
dados$grupo[dados$townid==87] = 3
dados$grupo[dados$townid==66] = 3

# grupo 4
dados$grupo[dados$townid==91] = 4
dados$grupo[dados$townid==86] = 4

# grupo 5

dados$grupo[dados$townid==29] =5
dados$grupo[dados$townid==85] =5

# grupo 6

dados$grupo[dados$townid==4] = 6
dados$grupo[dados$townid==84] = 6
dados$grupo[dados$townid==31] = 6


dados$group = factor(dados$grupo)
dadosag = groupedData(Y ~ lstat|group, data = dados)

mmpf6 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio, method = "REML", 
            random = ~1|townid, correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
            weights = varIdent(form = ~ 1 | group), 
            control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), data=dadosag)

summary(mmpf6)
residmmpf6 = residdiag.nlme(fit=mmpf6,limit=2,plotid=1)
#
LS = residmmpf6$"lesaffreverbeke.measure"   
round(LS[order(LS[,2], decreasing = F), ],dig=2)

dados$grupo = NULL
dados$grupo = 1
dados$grupo[dados$townid==78] = 2
dados$grupo[dados$townid==67] = 2
dados$grupo[dados$townid==69] = 2
dados$grupo[dados$townid==80] = 2
dados$grupo[dados$townid==81] = 2

dados$grupo[dados$townid==33] = 2
dados$grupo[dados$townid==76] = 2
dados$grupo[dados$townid==68] = 2
dados$grupo[dados$townid==77] = 2
dados$grupo[dados$townid==89] = 2
dados$grupo[dados$townid==83] = 2
dados$grupo[dados$townid==79] = 2

# grupo 3
dados$grupo[dados$townid==90] = 3
dados$grupo[dados$townid==82] = 3
dados$grupo[dados$townid==76] = 3
dados$grupo[dados$townid==77] = 3
dados$grupo[dados$townid==87] = 3
dados$grupo[dados$townid==66] = 3

# grupo 4
dados$grupo[dados$townid==91] = 4
dados$grupo[dados$townid==86] = 4

# grupo 5

dados$grupo[dados$townid==29] =5
dados$grupo[dados$townid==85] =5

# grupo 6

dados$grupo[dados$townid==4] = 6
dados$grupo[dados$townid==84] = 6
dados$grupo[dados$townid==31] = 6

# adding towns
dados$grupo[dados$townid==15] = 6
dados$grupo[dados$townid==79] = 6

dados$group = factor(dados$grupo)
dadosag = groupedData(Y ~ lstat|group, data = dados)

mmpf7 = lme(Y ~ I(lstat^0.5) + I(rm^0.5) + I((rm^0.5)*log(rm)) + I(omegacrim^-0.5) + 
              nox + rad + I(dis^-0.5) + taxt + ptratio, method = "REML", 
            random = ~1|townid, correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
            weights = varIdent(form = ~ 1 | group), 
            control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), data=dadosag)

summary(mmpf7)

residmmpf7 = residdiag.nlme(fit=mmpf7,limit=2,plotid=1)
#
LS = residmmpf7$"lesaffreverbeke.measure"   
round(LS[order(LS[,2], decreasing = F), ],dig=2)

# Comparing improvements of using clusters 
anova(mmpf1, mmpf2, mmpf3, mmpf4, mmpf5, mmpf6, mmpf7)

# Big improvement, stop clustering, we found a reasonable variance-covariance 
# for the error.

lmmFP2Step3 = mmpf7
round(summary(lmmFP2Step3)$tTable, dig=3)
VarCorr(lmmFP2Step3)  

# --------------------------------------------------------------------------
# STEP 4: Model reduction in the fixed part 
# --------------------------------------------------------------------------

mmpf2.ML = update(lmmFP2Step3, method = "ML")
summary(mmpf2.ML)$tTable

# -------------------------------------------
# Obtaining F pvalue for the ordering, since some have 2df (FP2)

LSTAT <- anova(mmpf2.ML, Terms = 2)$"p-value"
RM <- anova(mmpf2.ML, Terms = 3:4)$"p-value"
OMEGACRIM <- anova(mmpf2.ML, Terms = 5)$"p-value"
NOX <- anova(mmpf2.ML, Terms = 6)$"p-value"
RAD <- anova(mmpf2.ML, Terms = 7)$"p-value"
DIS <- anova(mmpf2.ML, Terms = 8)$"p-value"
TAXT <- anova(mmpf2.ML, Terms = 9)$"p-value"
PTRATIO <- anova(mmpf2.ML, Terms = 10)$"p-value"

pvalues = sort(c(LSTAT=LSTAT, RM=RM, OMEGACRIM=OMEGACRIM, NOX=NOX, RAD=RAD, DIS=DIS, TAXT=TAXT, PTRATIO=PTRATIO))
as.data.frame(pvalues)

#           pvalues
# RM        1.314839e-57
# LSTAT     7.295392e-42
# OMEGACRIM 5.167342e-13
# PTRATIO   4.928659e-07
# RAD       6.592494e-07
# TAXT      3.406866e-04
# NOX       4.548672e-04
# DIS       4.282047e-03

# Visiting order
# RM, LSTAT, CRIM, RAD, PTRATIO, TAXT, NOX, DIS
# I(rm^0.5) + I((rm^0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-0.5) + rad + 
# ptratio + taxt + nox + I(dis^-0.5)

# ----------------------------------------------------------------
# 1st Cycle
# ----------------------------------------------------------------
# 1) RM (1st Cycle)
lmmLin = lme(Y ~ rm + I(lstat^0.5) + I(omegacrim^-0.5) + rad +
               ptratio + taxt + nox + I(dis^-0.5), 
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

lmmNull = update(lmmLin, .~. -rm)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$rm)

# Model_selected "FP2 Model" "-0.5"      "-0.5" => I(rm^-0.5) + I((rm^-0.5)*log(rm))

### 2) LSTAT (1st Cycle)
lmmLin = lme(Y ~ I(rm^-0.5) + I((rm^-0.5)*log(rm)) + lstat + I(omegacrim^-0.5) + rad + 
               ptratio + taxt + nox + I(dis^-0.5),  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)
lmmNull = update(lmmLin, .~. -lstat)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$lstat)

# Model_selected "FP1 Model" "0.5"  => lstat^0.5

### 3) OMEGACRIM (1st Cycle)
lmmLin = lme(Y ~  I(rm^-0.5) + I((rm^-0.5)*log(rm))  + I(lstat^0.5) + omegacrim + rad + 
               ptratio + taxt + nox + I(dis^-0.5),    
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)
lmmNull = update(lmmLin, .~. -omegacrim)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$omegacrim)

# Model_selected "FP1 Model" "-1" => I(omegacrim^-1)

### 4) RAD (1st Cycle)
lmmLin = lme(Y ~ I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + rad + 
               ptratio + taxt + nox + I(dis^-0.5),   
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

lmmNull = update(lmmLin, .~. -rad)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$rad)

# Model_selected Linear Model

### 5) PTRATIO (1st Cycle)
lmmLin = lme(Y ~ I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + 
               rad + ptratio + taxt + nox + I(dis^-0.5),  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)
lmmNull = update(lmmLin, .~. -ptratio)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$ptratio)

# Model_selected "Linear Model"

### 6) TAXT (1st Cycle)
lmmLin = lme(Y ~ I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + 
               rad + ptratio + taxt + nox + I(dis^-0.5),  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

lmmNull = update(lmmLin, .~. -taxt)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$taxt)

# Model_selected "FP1 Model" "-2" 

### 7) NOX (1st Cycle)
lmmLin = lme(Y ~ I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + 
               rad + ptratio + I(taxt^-2) + nox + I(dis^-0.5),  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)
lmmNull = update(lmmLin, .~. -nox)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$nox)

# Model_selected "Linear Model"

### 8) DIS (1st Cycle)
lmmLin = lme(Y ~ I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + 
               rad + ptratio + I(taxt^-2) + nox + dis,  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

lmmNull = update(lmmLin, .~. -dis)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$dis)

# Model_selected "Linear Model"

# --------------------------------------
# Variables selected 1st Cycle:
# --------------------------------------
# RM         |  (-0.5; -0.5)
# LSTAT      |    0.5
# OMEGACRIM  |    -1 
# RAD        |     1  
# PTRATIO    |     1
# TAXT       |    -2
# NOX        |     1
# DIS        |     1

# ----------------------------------------------------------------
# 2o CICLO
# ----------------------------------------------------------------
### 1) RM (2o ciclo)
lmmLin = lme(Y ~ rm + I(lstat^0.5) + I(omegacrim^-1) + 
               rad + ptratio + I(taxt^-2) + nox + dis,  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

lmmNull = update(lmmLin, .~. -rm)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$rm)

# Model_selected "FP2 Model" "-0.5"      "-0.5"     

### 2) LSTAT (2o ciclo)
lmmLin = lme(Y ~ I(rm^-0.5) + I((rm^-0.5)*log(rm)) + lstat + I(omegacrim^-1) + 
               rad + ptratio + I(taxt^-2) + nox + dis,  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)
lmmNull = update(lmmLin, .~. -lstat)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$lstat)

# Model_selected "FP1 Model" "0.5"     

### 3) OMEGACRIM (2o ciclo)
lmmLin = lme(Y ~  I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + omegacrim + 
               rad + ptratio + I(taxt^-2) + nox + dis,  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)


lmmNull = update(lmmLin, .~. -omegacrim)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$omegacrim)

# Model_selected "FP1 Model" "-1"    

### 4) RAD
lmmLin = lme(Y ~  I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + 
             rad + ptratio + I(taxt^-2) + nox + dis,    
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

lmmNull = update(lmmLin, .~. -rad)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$rad)

# Model_selected "Linear Model"

### 5) PTRATIO
lmmLin = lme(Y ~  I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + 
               rad + ptratio + I(taxt^-2) + nox + dis,  
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

lmmNull = update(lmmLin, .~. -ptratio)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$ptratio)

# Model_selected "Linear Model"

### 6) TAXT
lmmLin = lme(Y ~  I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + 
               rad + ptratio + taxt + nox + dis, 
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

lmmNull = update(lmmLin, .~. -taxt)
MPF(fit1=lmmLin, fit2=lmmNull, covarx=dadosag$taxt)

# Model_selected "FP1 Model" "-2"

######################### Finish here since there was no changes from cycle 1
# and the next variables entered as linear

lmmLin = lme(Y ~  I(rm^-0.5) + I((rm^-0.5)*log(rm)) + I(lstat^0.5) + I(omegacrim^-1) + 
               rad + ptratio + I(taxt^-2) + nox + dis, 
             method = "ML", random = ~1|townid, 
             correlation = corCompSymm(value = 0.3, form = ~ 1 | townid),
             weights = varIdent(form = ~ 1 | group), 
             control=lmeControl(opt = "optim", optimMethod = "BFGS", msMaxIter = 100), 
             data=dadosag)

# --------------------------------------
# Variable and powers in the 2nd cycle 
# --------------------------------------
# RM         |  (-0.5;-0.5)
# LSTAT      |     0.5
# OMEGACRIM  |     -1 
# RAD        |      1
# PTRATIO    |      1
# TAXT       |     -2
# NOX        |      1
# DIS        |      1

  
# ----------------------------------------------------------
# Step 4: FINAL MODEL (REML)
# ----------------------------------------------------------

mmpfStep4 = update(lmmLin, method="REML")
summary(mmpfStep4)
round(summary(mmpfStep4)$tTable, dig=3)
VarCorr(mmpfStep4)  


# Comparing the improvements from step to step
rbind(c(AIC(mmpfStep1), AIC(mmpfStep2), AIC(mmpfStep3), AIC(mmpfStep4)),
c(BIC(mmpfStep1), BIC(mmpfStep2), BIC(mmpfStep3), BIC(mmpfStep4)),
c(logLik(mmpfStep1),logLik(mmpfStep2),logLik(mmpfStep3),logLik(mmpfStep4)))

# Diagnostic graphs
residStep4 = residdiag.nlme(fit=mmpfStep4,limit=2,plotid=1)
#
LS = residStep4$"lesaffreverbeke.measure"    
round(LS[order(LS[,2], decreasing = F), ],dig=2)
# Largest MLV: #36 (2.01) 

# -------------------------------------------------------------
# Other diagnostic graphs

X11()
residStep4 = residdiag.nlme(fit=mmpfStep4,limit=2,plotid=1)
#

X11()
residStep4 = residdiag.nlme(fit=mmpfStep4,limit=2,plotid=3)
#

X11()
residStep4 = residdiag.nlme(fit=mmpfStep4,limit=2,plotid=6)
#
