statComp <- function(result, eps = 1e-04)
{
  nonZeroLASSO <- abs(result$betalasso) >= eps
  nonZeroMCP <- abs(result$betamcp) >= eps
  sum((!nonZeroLASSO[, 2]) & (non.zero.index)) / sum(non.zero.index)# type I ER
  sum((nonZeroLASSO[, 2]) & (!non.zero.index)) / sum(!non.zero.index) # type II ER
}
partEst = function(X, y)
{  
  estimation = list()
  N = nrow(X)
  P = ncol(X)
  gammaInd = seq(1:P)
  result = weighted.sparsenet(X = X, y = y, delta = rep(1, N), 
                              max.iter = 1000, kappa0 = 1 / 3, n_kappa = 1, 
                              method = "include", eps = 1e-03)
  CV.result = cv.sparsenet(result, X = X, y = y, delta = rep(1, N), 
                           max.iter = 1000, fold = 5, method = "include", 
                           kappa0 = 1 / 3)
  CV.result$betamcp != 0
  Xest.cv.mcp = as.matrix(X[, CV.result$betamcp != 0])
  Xest.cv.lasso = as.matrix(X[, CV.result$betalasso != 0])

  if (length(Xest.cv.mcp) > 0) 
    lmOLS_cv_mcp = lm(y ~ Xest.cv.mcp - 1) else
      lmOLS_cv_mcp = lm(y ~ 1)
  if (length(Xest.cv.lasso) > 0) 
    lmOLS_cv_lasso = lm(y ~ Xest.cv.lasso - 1) else
      lmOLS_cv_lasso = lm(y ~ 1)

  estimation$cv_mcp = as.numeric(lmOLS_cv_mcp$coef[1:n.beta0])
  estimation$cv_lasso = as.numeric(lmOLS_cv_lasso$coef[1:n.beta0])

  return(estimation)
}