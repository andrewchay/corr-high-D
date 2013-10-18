partEst = function(Z, y, R1 = NULL, X = NULL, part1 = NULL, R2 = NULL)
{  
  estimation = list()
  N = nrow(Z)
  P = ncol(Z)
  gammaInd = seq(1:P)
  if (is.null(part1)) part1 = sample(seq(1, N), size = N / 2)
  part2 = seq(1, N)[-part1]
  if (is.null(dim(X))) xdim = rep(1, length(part1)) else 
    xdim = rep(ncol(X), length(part1))
  if (xdim[1] == 1) X1 = X[part1] else X1 = X[part1,]
  if (xdim[1] == 1) X2 = X[part2] else X2 = X[part2,]
  if (is.null(R1)) 
    R1 = diag(rep(1, N / 2)) - 
    X1 %*% solve(crossprod(X1, X1)) %*% t(X1)
  if (is.null(R2)) 
    R2 = diag(rep(1, N / 2)) - 
    X2 %*% solve(crossprod(X2, X2)) %*% t(X2)
  yadj = R1 %*% y[part1]
  Zadj = R1 %*% Z[part1,]
  result = weighted.sparsenet(X = Zadj, y = yadj, delta = rep(1, N / 2), 
                              max.iter = 1000, kappa0 = 1 / 3, n_kappa = 1, 
                              method = "include", eps = 1e-03)
  CV.result = cv.sparsenet(result, X = Zadj, y = yadj, delta = rep(1, N / 2), 
                           max.iter = 1000, fold = 5, method = "include", 
                           kappa0 = 1 / 3)
  BIC.result = bic.sparsenet(result, Zadj, yadj)

  Zest.cv.mcp = as.matrix(Z[, CV.result$betamcp != 0])
  Zest.cv.lasso = as.matrix(Z[, CV.result$betalasso != 0])
  Zest.bic.mcp = as.matrix(Z[, BIC.result$betamcp != 0])
  Zest.bic.lasso = as.matrix(Z[, BIC.result$betalasso != 0])

  if (length(Zest.cv.mcp) > 0) 
    lmOLS_cv_mcp = lm(y[part2] ~ X2 + Zest.cv.mcp[part2,] - 1) else
      lmOLS_cv_mcp = lm(y[part2] ~ X2 - 1)
  if (length(Zest.bic.mcp) > 0) 
    lmOLS_bic_mcp = lm(y[part2] ~ X2 + Zest.bic.mcp[part2,] - 1) else
      lmOLS_bic_mcp = lm(y[part2] ~ X2 - 1)
  if (length(Zest.cv.lasso) > 0) 
    lmOLS_cv_lasso = lm(y[part2] ~ X2 + Zest.cv.lasso[part2,] - 1) else
      lmOLS_cv_lasso = lm(y[part2] ~ X2 - 1)
  if (length(Zest.bic.lasso) > 0) 
    lmOLS_bic_lasso = lm(y[part2] ~ X2 + Zest.bic.lasso[part2,] - 1) else
      lmOLS_bic_lasso = lm(y[part2] ~ X2 - 1)
  
  estimation$sd_cv_mcp = as.numeric(sqrt(diag(vcov(lmOLS_cv_mcp))[1:n.beta0]))
  estimation$sd_bic_mcp = as.numeric(sqrt(diag(vcov(lmOLS_bic_mcp))[1:n.beta0]))
  estimation$sd_cv_lasso = as.numeric(sqrt(diag(vcov(lmOLS_cv_lasso))[1:n.beta0]))
  estimation$sd_bic_lasso = as.numeric(sqrt(diag(vcov(lmOLS_bic_lasso))[1:n.beta0]))
  estimation$sigma_cv_mcp = sqrt(sum(lmOLS_cv_mcp$resid ^ 2) / 
                                   lmOLS_cv_mcp$df.resid)
  estimation$sigma_bic_mcp = sqrt(sum(lmOLS_bic_mcp$resid ^ 2) / 
                                    lmOLS_bic_mcp$df.resid)
  estimation$sigma_cv_lasso = sqrt(sum(lmOLS_cv_lasso$resid ^ 2) / 
                                     lmOLS_cv_lasso$df.resid)
  estimation$sigma_bic_lasso = sqrt(sum(lmOLS_bic_lasso$resid ^ 2) / 
                                      lmOLS_bic_lasso$df.resid)
  yadj = R2 %*% y[part2]
  Zadj = R2 %*% Z[part2,]
  result = weighted.sparsenet(X = Zadj, y = yadj, delta = rep(1, N / 2), 
                              max.iter = 1000, kappa0 = 1 / 3, n_kappa = 1, 
                              method = "include", eps = 1e-03)
  CV.result = cv.sparsenet(result, X = Zadj, y = yadj, delta = rep(1, N / 2), 
                           max.iter = 1000, fold = 5, method = "include", 
                           kappa0 = 1 / 3)
  BIC.result = bic.sparsenet(result, Zadj, yadj)
  
  Zest.cv.mcp = as.matrix(Z[, CV.result$betamcp != 0])
  Zest.cv.lasso = as.matrix(Z[, CV.result$betalasso != 0])
  Zest.bic.mcp = as.matrix(Z[, BIC.result$betamcp != 0])
  Zest.bic.lasso = as.matrix(Z[, BIC.result$betalasso != 0])
  
  if (length(Zest.cv.mcp) > 0) 
    lmOLS_cv_mcp = lm(y[part1] ~ X1 + Zest.cv.mcp[part1,] - 1) else
      lmOLS_cv_mcp = lm(y[part1] ~ X1 - 1)
  if (length(Zest.bic.mcp) > 0) 
    lmOLS_bic_mcp = lm(y[part1] ~ X1 + Zest.bic.mcp[part1,] - 1) else
      lmOLS_bic_mcp = lm(y[part1] ~ X1 - 1)
  if (length(Zest.cv.lasso) > 0) 
    lmOLS_cv_lasso = lm(y[part1] ~ X1 + Zest.cv.lasso[part1,] - 1) else
      lmOLS_cv_lasso = lm(y[part1] ~ X1 - 1)
  if (length(Zest.bic.lasso) > 0) 
    lmOLS_bic_lasso = lm(y[part1] ~ X1 + Zest.bic.lasso[part1,] - 1) else
      lmOLS_bic_lasso = lm(y[part1] ~ X1 - 1)
  
  estimation$sd_cv_mcp = (estimation$sd_cv_mcp + 
                            as.numeric(sqrt(diag(vcov(lmOLS_cv_mcp))[1:n.beta0]))) / 2
  estimation$sd_bic_mcp = (estimation$sd_bic_mcp + 
                             as.numeric(sqrt(diag(vcov(lmOLS_bic_mcp))[1:n.beta0]))) / 2
  estimation$sd_cv_lasso = (estimation$sd_cv_lasso + 
                              as.numeric(sqrt(diag(vcov(lmOLS_cv_lasso))[1:n.beta0]))) / 2
  estimation$sd_bic_lasso = (estimation$sd_bic_lasso +
                               as.numeric(sqrt(diag(vcov(lmOLS_bic_lasso))[1:n.beta0]))) / 2
  estimation$sigma_cv_mcp = (estimation$sigma_cv_mcp + 
                               sqrt(sum(lmOLS_cv_mcp$resid ^ 2) / 
                                   lmOLS_cv_mcp$df.resid)) / 2
  estimation$sigma_bic_mcp = (estimation$sigma_bic_mcp + 
                                sqrt(sum(lmOLS_bic_mcp$resid ^ 2) / 
                                    lmOLS_bic_mcp$df.resid)) / 2
  estimation$sigma_cv_lasso = (estimation$sigma_cv_lasso + 
                                 sqrt(sum(lmOLS_cv_lasso$resid ^ 2) / 
                                     lmOLS_cv_lasso$df.resid)) / 2
  estimation$sigma_bic_lasso = (estimation$sigma_bic_lasso + 
                                  sqrt(sum(lmOLS_bic_lasso$resid ^ 2) / 
                                      lmOLS_bic_lasso$df.resid)) / 2
  
  
  R = diag(rep(1, N)) - 
    X %*% solve(crossprod(X, X)) %*% t(X)
  Zadj = R %*% Z
  yadj = R %*% y
  result = weighted.sparsenet(X = Zadj, y = yadj, delta = rep(1, N), 
                              max.iter = 1000, kappa0 = 1 / 3, n_kappa = 1, 
                              method = "include", eps = 1e-03)
  CV.result = cv.sparsenet(result, X = Zadj, y = yadj, delta = rep(1, N), 
                           max.iter = 1000, fold = 5, method = "include", 
                           kappa0 = 1 / 3)
  BIC.result = bic.sparsenet(result, Zadj, yadj)
  Zest.cv.mcp = as.matrix(Z[, CV.result$betamcp != 0])
  Zest.cv.lasso = as.matrix(Z[, CV.result$betalasso != 0])
  Zest.bic.mcp = as.matrix(Z[, BIC.result$betamcp != 0])
  Zest.bic.lasso = as.matrix(Z[, BIC.result$betalasso != 0])
  if (length(Zest.cv.mcp) > 0) 
    lmOLS_cv_mcp = lm(y ~ X + Zest.cv.mcp - 1) else
      lmOLS_cv_mcp = lm(y ~ X - 1)
  if (length(Zest.bic.mcp) > 0) 
    lmOLS_bic_mcp = lm(y ~ X + Zest.bic.mcp - 1) else
      lmOLS_bic_mcp = lm(y ~ X - 1)
  if (length(Zest.cv.lasso) > 0) 
    lmOLS_cv_lasso = lm(y ~ X + Zest.cv.lasso - 1) else
      lmOLS_cv_lasso = lm(y ~ X - 1)
  if (length(Zest.bic.lasso) > 0) 
    lmOLS_bic_lasso = lm(y ~ X + Zest.bic.lasso - 1) else
      lmOLS_bic_lasso = lm(y ~ X - 1)
  
  estimation$beta_cv_mcp = as.numeric(CV.result$betamcp)
  estimation$beta_bic_mcp = as.numeric(BIC.result$betamcp)
  estimation$beta_cv_lasso = as.numeric(CV.result$betalasso)
  estimation$beta_bic_lasso = as.numeric(BIC.result$betalasso)
    
  estimation$gammaOLS_cv_mcp = as.numeric(lmOLS_cv_mcp$coef[1:n.beta0])
  estimation$gammaOLS_bic_mcp = as.numeric(lmOLS_bic_mcp$coef[1:n.beta0])
  estimation$gammaOLS_cv_lasso = as.numeric(lmOLS_cv_lasso$coef[1:n.beta0])
  estimation$gammaOLS_bic_lasso = as.numeric(lmOLS_bic_lasso$coef[1:n.beta0])
  tmp = rep(0, P)
  tmp[gammaInd[CV.result$betamcp != 0]] = lmOLS_cv_mcp$coef[-(1:n.beta0)]
  estimation$ZOLS_cv_mcp = as.numeric(tmp)
  tmp = rep(0, P)
  tmp[gammaInd[BIC.result$betamcp != 0]] = lmOLS_bic_mcp$coef[-(1:n.beta0)]
  estimation$ZOLS_bic_mcp = as.numeric(tmp)
  tmp = rep(0, P)
  tmp[gammaInd[CV.result$betalasso != 0]] = lmOLS_cv_lasso$coef[-(1:n.beta0)]
  estimation$ZOLS_cv_lasso = as.numeric(tmp)
  tmp = rep(0, P)
  tmp[gammaInd[BIC.result$betalasso != 0]] = lmOLS_bic_lasso$coef[-(1:n.beta0)]
  estimation$ZOLS_bic_lasso = as.numeric(tmp)
    
  if (length(Zest.cv.mcp) > 0)
    resp = y - 
    crossprod(t(Zest.cv.mcp), 
              CV.result$betamcp[CV.result$betamcp != 0]) else
                resp = y
  estimation$gammaILS_cv_mcp = as.numeric(lm(resp ~ X - 1)$coef)
  
  if (length(Zest.cv.lasso) > 0)
    resp = y - 
    crossprod(t(Zest.cv.lasso), 
              CV.result$betalasso[CV.result$betalasso != 0]) else
                resp = y
  estimation$gammaILS_cv_lasso = as.numeric(lm(resp ~ X - 1)$coef)
  
  if (length(Zest.bic.mcp) > 0)
    resp = y - 
    crossprod(t(Zest.bic.mcp), 
              BIC.result$betamcp[BIC.result$betamcp != 0]) else
                resp = y  
  estimation$gammaILS_bic_mcp = as.numeric(lm(resp ~ X - 1)$coef)
  
  if (length(Zest.bic.lasso) > 0)
    resp = y - 
    crossprod(t(Zest.bic.lasso), 
              BIC.result$betalasso[BIC.result$betalasso != 0]) else
                resp = y
  estimation$gammaILS_bic_lasso = as.numeric(lm(resp ~ X - 1)$coef)
  

  return(estimation)
}