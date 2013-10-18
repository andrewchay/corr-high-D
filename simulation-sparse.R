# This simulation compares all three methods. dlasso stands for double LASSO and post estimation method. scaled is the scaled LASSO by Zhang and partEst is our method.
setwd("D:/Lecture/Homework/291/corr-high-D")
library("CCSparsenet")
require("MASS")
library("ggplot2")

####simulate with 500 datesets ####
Oracle = c()
cv_mcp = c()
cv_lasso = c()

sd_Oracle = c()
sd_cv_mcp = c()
sd_cv_lasso = c()

sigma_cv_mcp = c()
sigma_cv_lasso = c()
sigma_Oracle = c()

for (i in 1:num.simu)
{
  e = rnorm(N, 0, sd(mu) / sqrt(sn))
  y = mu + e
  res0 = lm(y ~ X + Z[, 1:num.nzeros] - 1)
  res1 = CDgrid(Z, y, X = X, part1 = part1)
  
  Oracle = cbind(Oracle, res0$coef[1:n.beta0])
  cv_mcp = cbind(parEst_cv_mcp, res1$gammaOLS_cv_mcp)
  cv_lasso = cbind(parEst_O_cv_lasso, res1$gammaOLS_cv_lasso)
  sd_Oracle = cbind(sd_Oracle, sqrt(diag(vcov(res0))[1:n.beta0]))
  parEst_sd_cv_mcp = cbind(parEst_sd_cv_mcp, res1$sd_cv_mcp)
  parEst_sd_cv_lasso = cbind(parEst_sd_cv_lasso, res1$sd_cv_lasso)
  parEst_sd_bic_mcp = cbind(parEst_sd_bic_mcp, res1$sd_bic_mcp)
  parEst_sd_bic_lasso = cbind(parEst_sd_bic_lasso, res1$sd_bic_lasso)

  
  sigma_cv_mcp = cbind(sigma_cv_mcp, res1$sigma_cv_mcp)
  sigma_bic_mcp = cbind(sigma_bic_mcp, res1$sigma_bic_mcp)
  sigma_cv_lasso = cbind(sigma_cv_lasso, res1$sigma_cv_lasso)
  sigma_bic_lasso = cbind(sigma_bic_lasso, res1$sigma_bic_lasso)
  sigma_Oracle = cbind(sigma_Oracle, sqrt(sum(res0$resid ^ 2) / res0$df.resid))
  print(paste(i, " out of ", num.simu))
}
Oracle = Oracle - beta0
parEst_O_cv_mcp = parEst_O_cv_mcp - beta0
parEst_I_cv_mcp = parEst_I_cv_mcp - beta0
parEst_O_cv_lasso = parEst_O_cv_lasso - beta0
parEst_I_cv_lasso = parEst_I_cv_lasso - beta0



pmatrix1 = matrix(rep(0, 5 * 5), nrow = 5)
pmatrix1[1, 1] = round(mean(abs(parEst_O_cv_mcp)), 2)
pmatrix1[2, 1] = round(mean(abs(parEst_O_cv_lasso)), 2)
pmatrix1[3, 1] = round(mean(abs(parEst_I_cv_mcp)), 2)
pmatrix1[4, 1] = round(mean(abs(parEst_I_cv_lasso)), 2)
pmatrix1[5, 1] = round(mean(abs(Oracle)), 2)

pmatrix1[1, 2] = round(mean(parEst_sd_cv_mcp), 2)
pmatrix1[2, 2] = round(mean(parEst_sd_cv_lasso), 2)
pmatrix1[3, 2] = round(mean(parEst_sd_cv_mcp), 2)
pmatrix1[4, 2] = round(mean(parEst_sd_cv_lasso), 2)
pmatrix1[5, 2] = round(mean(sd_Oracle), 2)

pmatrix1[1, 3] = round(mean(parEst_sd_cv_mcp ^ 2 + parEst_O_cv_mcp ^ 2), 2)
pmatrix1[2, 3] = round(mean(parEst_sd_cv_lasso ^ 2 + parEst_O_cv_lasso ^ 2), 2)
pmatrix1[3, 3] = round(mean(parEst_sd_cv_mcp ^ 2 + parEst_I_cv_mcp ^ 2), 2)
pmatrix1[4, 3] = round(mean(parEst_sd_cv_lasso ^ 2 + parEst_I_cv_lasso ^ 2), 2)
pmatrix1[5, 3] = round(mean(sd_Oracle ^ 2 + Oracle ^ 2), 2)

zval = qnorm(0.975, 0, 1)
coverage = c(mean((parEst_O_cv_mcp + parEst_sd_cv_mcp * zval >= 0) & 
                     (parEst_O_cv_mcp - parEst_sd_cv_mcp * zval <= 0)),
             mean((parEst_O_cv_lasso + parEst_sd_cv_lasso * zval >= 0) & 
                     (parEst_O_cv_lasso - parEst_sd_cv_lasso * zval <= 0)),
             mean((Oracle + sd_Oracle * zval >= 0) & 
                    (Oracle - sd_Oracle * zval <= 0)))

pmatrix1[, 4] = t(coverage)
pmatrix1[1, 5] = round(sd(as.numeric(parEst_O_cv_mcp)), 2)
pmatrix1[2, 5] = round(sd(as.numeric(parEst_O_cv_lasso)), 2)
pmatrix1[5, 5] = round(sd(as.numeric(Oracle)), 2)
colnames(pmatrix1) = c("Bias", "SD", "MSE", 
                       "Coverage", "Simulation SD")
rownames(pmatrix1) = c("OLS CV MCP", "OLS CV LASSO",
                       "ILS CV MCP", "ILS CV LASSO",
                       "Oracle")



