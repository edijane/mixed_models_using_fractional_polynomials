MPF = function(fit1, fit2, covarx) 
  {
   dadosag$covarx = covarx
   devnull = as.numeric(-2*logLik(fit2))
   devlin = as.numeric(-2*logLik(fit1))
   PF2 = matrix(0,nr=36,nc=5)
   colnames(PF2) = c("p1", "p2", "dev.diff", "-2logLik", "AIC")
   pot = c(-2,-1,-0.5,0,0.5,1,2,3)
   i = 0
   for (k in 1:length(pot)) 
     {
      p1 = pot[k] 
      for (l in k:length(pot)) 
        {
         p2 = pot[l]
         i = i + 1
         if (p1==0) {dadosag$x1=log(covarx)} else {dadosag$x1=covarx^p1}
         if (p2==0) {dadosag$x2=log(covarx)} else {dadosag$x2=covarx^p2}
         if (p1==p2) {dadosag$x2=dadosag$x1*log(covarx)}
         mfp2 = try(update(fit2, .~. + x1 + x2), 
                    silent = FALSE)
         if(inherits(mfp2, "try-error")) 
           {PF2[i,] = c(p1,p2,NA,NA,NA)} else {
                  devdiff = -2*logLik(fit1) - (-2*logLik(mfp2))
                  PF2[i,] = c(p1,p2,devdiff,-2*logLik(mfp2),AIC(mfp2))
           }
      }
   }
   
   matrizPF2 = PF2[order(PF2[,4]), ]
   devpf2 = matrizPF2[1,4]
   logLikdif_null  = as.numeric(devnull - devpf2)
   valorp_null = as.numeric(1 - pchisq(logLikdif_null, 4))
   my_list = list(MatrizPF2=matrizPF2, Power_FP2=matrizPF2[1,1:2], 
   Deviance_FP2=devpf2, Deviance_Null=devnull, 
   Devdiff_Null_FP2=logLikdif_null, 
   PValue_Null_FP2=valorp_null, Model_selected="Null Model")
   
   if (valorp_null >= 0.05) 
     {
     return(my_list)
     } else {
             logLikdif_lin = as.numeric(devlin - devpf2)
             valorp_lin = as.numeric(1 - pchisq(logLikdif_lin, 3))
             my_list = list(MatrizPF2=matrizPF2, Power_FP2=matrizPF2[1,1:2], 
                            Power_Lin = 1, Deviance_FP2=devpf2, Deviance_Null=devnull, 
                            Deviance_Lin=devlin, DevDiff_Null_FP2=logLikdif_null, 
                            DevDiff_Lin_FP2=logLikdif_lin, PValue_Null_FP2=valorp_null, 
                            Pvalue_Lin_FP2=valorp_lin,
                            Model_selected = "Linear Model")
             if (valorp_lin >= 0.05) 
              {
               return(my_list)
              } else{
                     PF1  = matrix(0,nr=length(pot),nc=4)
                     colnames(PF1) = c("p1","dev.diff","-2logLik","AIC")
                     i = 0
                     for (k in 1:length(pot)) 
                      {
                       p1 = pot[k] 
                       i = i + 1
                       if (p1==0) 
                       {dadosag$x1 = log(covarx)} else {dadosag$x1 = covarx^p1}
                       mfp1 = try(update(fit2, .~. + x1), silent = FALSE)
                       if(inherits(mfp1, "try-error")) 
                        {
                         PF1[i,] = c(p1,NA,NA,NA)
                        } else {
                                devdiff = -2*logLik(fit1) - (-2*logLik(mfp1))
                                PF1[i,] = c(p1, devdiff, -2*logLik(mfp1), 
                                         AIC(mfp1))
                                }
                       }
                     matrizPF1 = PF1[order(PF1[,3]), ]
                     
                     devpf1 = matrizPF1[1,3]
                     logLikdif_pf1 = as.numeric(devpf1 - devpf2)
                     valorp_pf1 = as.numeric(1 - pchisq(logLikdif_pf1, 2))
                     my_list = list(MatrizPF2=matrizPF2, MatrizPF1=matrizPF1,
                              Power_FP2=matrizPF2[1,1:2], Power_Lin = 1, 
                              Power_FP1=matrizPF1[1,1], 
                              Deviance_FP2=devpf2, Deviance_Null=devnull, 
                              Deviance_Lin=devlin,
                              Deviance_FP1=devpf1, 
                              DevDiff_Null_FP2=logLikdif_null, 
                              DevDiff_Lin_FP2=logLikdif_lin, 
                              DevDiff_FP1_FP2=logLikdif_pf1, 
                              PValue_Null_FP2=valorp_null, PValue_Lin_FP2=valorp_lin, 
                              PValue_FP1_FP2=valorp_pf1, 
                              Model_selected=c("FP1 Model", matrizPF1[1,1]))
                     my_list2 = list(MatrizPF2=matrizPF2, MatrizPF1=matrizPF1,
                                     Power_FP2=matrizPF2[1,1:2], Power_Lin = 1, 
                                     Power_FP1=matrizPF1[1,1], 
                              Deviance_FP2=devpf2, Deviance_Null=devnull, 
                              Deviance_Lin=devlin,
                              Deviance_FP1=devpf1,
                              DevDiff_Null_FP2=logLikdif_null, DevDiff_Lin_FP2=logLikdif_lin, 
                              DevDiff_FP1_FP2=logLikdif_pf1, 
                              PValue_Null_FP2=valorp_null, PValue_Lin_FP2=valorp_lin, 
                              PValue_FP1_FP2=valorp_pf1, 
                              Model_selected=c("FP2 Model", 
                                               matrizPF2[1,1:2]))
                     if (valorp_pf1 >= 0.05) {return(my_list)} else{return(my_list2)}
              }
     }
            
}
