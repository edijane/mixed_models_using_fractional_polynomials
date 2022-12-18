# ADAPTED from:
# http://www-personal.umich.edu/~tombraun/LeeBraun/Permutation%20Test%20for%20Random%20Slope.R
PermutationTest <- function(Nperms, Y, numsub, positionranef, fit0, fit1, 
                            randpart0, randpart1, controlMER, df)
{
  options(warn = 2)
  n = nrow(df)

  # --------------------------------------------------------------------------------
  # Adding random slopes for pressure
  # --------------------------------------------------------------------------------
  
  # Create objects that will contain the observed test statistics and 
  # permutation distributions.
  lrtest = c(rep(-9999,(Nperms+1)))
  t2     = c(rep(-9999,(Nperms+1)))
  
  # Recording the observed test statistics.
  t2[1]     = sum(matrix(unlist(ranef(fit1)),numsub,positionranef)[,positionranef]^2) 
  lrtest[1] = 2*(logLik(fit1, REML=TRUE)-logLik(fit0,REML=TRUE))
  
  # Calculating the marginal residuals
  MX      = model.matrix(fit1)
  eblue  = matrix(fixef(fit1), nc=1)
  errors = matrix(Y, nc=1) - MX%*%eblue
  
  # Generating the design matrix Z to get V = ZGt(Z)+R
  Z = getME(fit1, "Z")
  
  # Random effect covariance matrix under the null hypothesis.
  # Needs adaptation in case of more than 3 random effects
  if(positionranef==1) {
    G0 = 0
    R0 = diag(sigma(fit0)^2, n) 
  }
  if(positionranef==2) {
    G0 = diag(c(unlist(VarCorr(fit0))[1], 0))
    R0 = diag(rep(attr(VarCorr(fit0),"sc")^2, n))   
  }
  if(positionranef==3){
    G0 = cbind(rbind(matrix(unlist(VarCorr(fit0)[1]), nc=2), 
                                       rep(0,2)), rep(0,3))
    R0 = diag(rep(attr(VarCorr(fit0),"sc")^2, n))   
  }
  G0 = bdiag(replicate(numsub,G0,simplify=FALSE)) 
  
  V0 = Z%*%G0%*%t(Z) + R0
  invV0 = solve(V0,diag(1,dim(V0)[1]))
  U = chol(V0)
  A = solve(t(U),diag(dim(U)[1]))
  
  # Weighting the residuals.
  weightederrors = A%*%errors
  permweightederrors = matrix(weightederrors,n,Nperms+1)
  
  # Permuting the weighted residuals.
  for(i in 2:(Nperms+1)) {
    permweightederrors[,i] = sample(weightederrors, replace=FALSE)
  }
  
  # Unweighting the permuted residuals.
  weighted = (t(U)%*%permweightederrors)[,2:(Nperms+1)]
  
  # Calculating the restricted likelihood ratio test statistic for each permutation.  
  # try() used to replace estimates from warning message fits with NA.
  fitmodels = function(res)
  {
    fixpart0 = 'res ~ 1'
    blup = try(sum(matrix(unlist(ranef(suppressMessages(lmer(formula(paste(fixpart0, 
                   randpart1, sep = " + ")), data=df, control = do.call(lmerControl, controlMER))))),
      numsub,positionranef)[,positionranef]^2), silent=TRUE)
    if(positionranef>1) 
    {
     dev = try(2*(logLik(suppressMessages(lmer(formula(paste(fixpart0, randpart1, sep = " + ")),
                   REML=TRUE, data=df, control = do.call(lmerControl, controlMER)))) -
                  logLik(suppressMessages(lmer(formula(paste(fixpart0, randpart0, sep = " + ")), 
                   REML=TRUE, data=df, control = do.call(lmerControl, controlMER))))), silent=TRUE)
    } else 
   {
    dev = try(2*(logLik(suppressMessages(lmer(formula(paste(fixpart0, randpart1, 
                  sep = " + ")), REML=TRUE, data=df, control = do.call(lmerControl, 
                  controlMER)))) - logLik(lm(formula(paste(fixpart0)), data=df))), silent=TRUE)
    }
    
    lrt = ifelse(is.numeric(dev)==TRUE, dev, NA)
    blupt = ifelse(is.numeric(blup)==TRUE, blup, NA)
    return(cbind(lrt, blupt))
  }  
  
  resperm = matrix(apply(weighted[,2:Nperms],2,fitmodels), byrow=T, nc=2)
  lrtest2 = resperm[,1]
  
  lrtest = c(lrtest[1],lrtest2)
  lrtest       = ifelse(lrtest<0, 0, lrtest)
  permdistlpt1 = lrtest[-1]
  rejectLRp    = mean(ifelse(permdistlpt1>=lrtest[1],1,0),na.rm=TRUE)
  NA_CountsLRT=sum(is.na(lrtest))
  
  t22 = resperm[,2]
  t2     = c(t2[1],t22)
  permdistbpt2 = t2[-1]
  rejectBLUPp    = mean(ifelse(permdistbpt2>=t2[1], 1, 0),na.rm=TRUE)
  NA_CountsBLUP=sum(is.na(t2))
  
  return(c(pvPermLRT=rejectLRp, NA_CountsLRT= NA_CountsLRT, 
           pvPermBLUP=rejectBLUPp, NA_CountsBLUP=NA_CountsBLUP))
}
