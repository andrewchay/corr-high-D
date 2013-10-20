type1 <- function(a, b)
{
  if (length(a) != length(b)) stop("Dimensionalities are different.")
  sum((a) & (!b)) / sum(!b)
}
type2 <- function(a, b)
{
  if (length(a) != length(b)) stop("Dimensionalities are different.")
  sum((!a) & (b)) / sum(b)
}
FDR <- function(a, b)
{
  if (length(a) != length(b)) stop("Dimensionalities are different.")
  if (sum(a) > 0) return(sum((a) & (!b)) / sum(a)) else return(0)
}
FIR <- function(a, b)
{
  if (length(a) != length(b)) stop("Dimensionalities are different.")
  if (sum(!a) > 0) return(sum((!a) & (b)) / sum(!a)) else return(0)
}
Dist <- function(a, b, sp = 1)
{
  (sum(abs((a - b)) ^ sp)) ^ (1 / sp)
}
DistX <- function(a, b, X, sp = 1)
{
  (sum(abs(crossprod(t(X), (a - b))) ^ sp)) ^ (1 / sp)
}
statComp <- function(result, eps = 1e-04)
{
  nonZeroLASSO <- abs(result$betalasso) >= eps
  nonZeroMCP <- abs(result$betamcp) >= eps
  
  Ltype1 <- apply(nonZeroLASSO, 2, type1, non.zero.index) # type I ER
  Ltype2 <- apply(nonZeroLASSO, 2, type2, non.zero.index) # type II ER
  LFDR <- apply(nonZeroLASSO, 2, FDR, non.zero.index) # FDR
  LFIR <- apply(nonZeroLASSO, 2, FIR, non.zero.index) # FIR
  Mtype1 <- apply(nonZeroMCP, 2, type1, non.zero.index) # type I ER
  Mtype2 <- apply(nonZeroMCP, 2, type2, non.zero.index) # type II ER
  MFDR <- apply(nonZeroMCP, 2, FDR, non.zero.index) # FDR
  MFIR <- apply(nonZeroMCP, 2, FIR, non.zero.index) # FIR
  LL.5 <- apply(result$betalasso, 2, Dist, beta, sp = .5) # L1 norm (beta)
  LL1 <- apply(result$betalasso, 2, Dist, beta, sp = 1) # L1 norm (beta)
  LL1.5 <- apply(result$betalasso, 2, Dist, beta, sp = 1.5) # L1 norm (beta)
  LL2 <- apply(result$betalasso, 2, Dist, beta, sp = 2) # L2 norm (beta)
  LXL1 <- apply(result$betalasso, 2, DistX, beta, X, sp = 1) # L1 norm (X * beta)
  LXL2 <- apply(result$betalasso, 2, DistX, beta, X, sp = 2) # L2 norm (X * beta)
  ML.5 <- apply(result$betamcp, 2, Dist, beta, sp = .5) # L1 norm (beta)
  ML1 <- apply(result$betamcp, 2, Dist, beta, sp = 1) # L1 norm (beta)
  ML1.5 <- apply(result$betamcp, 2, Dist, beta, sp = 1.5) # L1 norm (beta)
  ML2 <- apply(result$betamcp, 2, Dist, beta, sp = 2) # L2 norm (beta)
  MXL1 <- apply(result$betamcp, 2, DistX, beta, X, sp = 1) # L1 norm (X * beta)
  MXL2 <- apply(result$betamcp, 2, DistX, beta, X, sp = 2) # L2 norm (X * beta)
  M <- cbind(Ltype1, Ltype2, LFDR, LFIR, Mtype1, Mtype2, MFDR, MFIR, LL.5, LL1, 
             LL1.5, LL2, LXL1, LXL2, ML.5, ML1, ML1.5, ML2, MXL1, MXL2, 
             lambda = result$lambda)
  return(M)
}
CDgrid = function(X, y, eps = 1e-04)
{  
  estimation = list()
  N = nrow(X)
  P = ncol(X)
  result = weighted.sparsenet(X = X, y = y, delta = rep(1, N), 
                              max.iter = 1000, kappa0 = 1 / 3, n_kappa = 1, 
                              method = "include", eps = 1e-03)
  CV.result = cv.sparsenet(result, X = X, y = y, delta = rep(1, N), 
                           max.iter = 1000, fold = 5, method = "include", 
                           kappa0 = 1 / 3, power = 2)
  CV.result$betamcp != 0
  mcp.nonzero = (CV.result$betamcp != 0)
  lasso.nonzero = (CV.result$betalasso != 0)
  Xest.cv.mcp = as.matrix(X[, mcp.nonzero])
  Xest.cv.lasso = as.matrix(X[, lasso.nonzero])

  if (length(Xest.cv.mcp) > 0) 
    lmOLS_cv_mcp = lm(y ~ Xest.cv.mcp - 1) else
      lmOLS_cv_mcp = lm(y ~ 1)
  if (length(Xest.cv.lasso) > 0) 
    lmOLS_cv_lasso = lm(y ~ Xest.cv.lasso - 1) else
      lmOLS_cv_lasso = lm(y ~ 1)

  tmp = rep(0, P)
  tmp[mcp.nonzero] = lmOLS_cv_mcp$coef
  estimation$beta_mcp = tmp
  tmp = rep(0, P)
  tmp[lasso.nonzero] = lmOLS_cv_lasso$coef
  estimation$beta_lasso = tmp
  estimation$lambda = result$lambda
  estimation$kappa = result$kappa
  estimation$lasso.lam = CV.result$lasso.lambda
  estimation$mcp.lam = CV.result$mcp.lambda
  estimation$Ltype1 = type1(abs(CV.result$betalasso) >= eps, non.zero.index)
  estimation$Ltype2 = type2(abs(CV.result$betalasso) >= eps, non.zero.index)
  estimation$LFDR = FDR(abs(CV.result$betalasso) >= eps, non.zero.index)
  estimation$LFIR = FIR(abs(CV.result$betalasso) >= eps, non.zero.index)
  estimation$Mtype1 = type1(abs(CV.result$betamcp) >= eps, non.zero.index)
  estimation$Mtype2 = type2(abs(CV.result$betamcp) >= eps, non.zero.index)
  estimation$MFDR = FDR(abs(CV.result$betamcp) >= eps, non.zero.index)
  estimation$MFIR = FIR(abs(CV.result$betamcp) >= eps, non.zero.index)
  estimation$LL1 = Dist(CV.result$betalasso, beta, 1)
  estimation$LL2 = Dist(CV.result$betalasso, beta, 2)
  estimation$LXL1 = DistX(CV.result$betalasso, beta, X, 1)
  estimation$LXL2 = DistX(CV.result$betalasso, beta, X, 2)
  estimation$ML1 = Dist(CV.result$betamcp, beta, 1)
  estimation$ML2 = Dist(CV.result$betamcp, beta, 2)
  estimation$MXL1 = DistX(CV.result$betamcp, beta, X, 1)
  estimation$MXL2 = DistX(CV.result$betamcp, beta, X, 2)
  estimation$stats = statComp(result)

  return(estimation)
}