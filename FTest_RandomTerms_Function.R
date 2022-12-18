Ftest_Demidenko = function(modfix, modrandom, dataframe, Y)
{
  # Model matrix fixed part
  X = model.matrix(modfix, data=dataframe)
  
  # Z matrix random part
  Z = model.matrix(modrandom, data=dataframe)
  
  # Full model matrix
  W = cbind(X, Z)
  Y  = as.matrix(Y)
  X_inv = ginv(t(X)%*%X)
  pX = X%*%X_inv%*%t(X)
  W_invg = ginv(t(W)%*%W, tol = sqrt(.Machine$double.eps))
  pW = W%*%W_invg%*%t(W)
  n = nrow(dataframe)
  p = rankMatrix(X)[[1]]
  W_rank = rankMatrix(W)[[1]]
  In = diag(1, nrow = n)
  dfn = (W_rank - p)
  dfd = (n - W_rank)
  num_FD = (t(Y)%*%(pW - pX)%*%Y)/dfn
  den_FD = (t(Y)%*%(In - pW)%*%Y)/dfd
  FD = num_FD/den_FD
  round(FD, dig=2)
  pvalue = pf(FD, dfn, dfd, lower.tail = FALSE)
  c(dfn=dfn, dfd=dfd, F_=round(FD, dig=3), pvalue=pvalue )
}