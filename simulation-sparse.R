# This simulation compares all three methods. dlasso stands for double LASSO and post estimation method. scaled is the scaled LASSO by Zhang and partEst is our method.
library(CCSparsenet, lib.loc = "/mnt/nfs/netapp2/grad/hchai/Rpackages")
require("MASS")
library("ggplot2")
source("dlasso.R")
source("scaled.R")
source("partEst-modified.R")
source("scaled_XX.R")
source("dlasso_X.R")
#set.seed = 100
N = 200
P = 500
num.simu = 500
rho = 0.9
matrix_vec = function(P){
  tmp = matrix(rep(0, P * P), nrow = P)
  for (i in 1:P) 
    for (j in 1:P)
      if (i == j) tmp[i, j] = 1 else 
        tmp[i, j] = (rho) ^ abs(i - j)
        #tmp[i, j] = rho
  return(tmp)
}
lambdaUniv = sqrt(2 * log(P) / N)
beta0 = 1
n.beta0 = length(beta0)
g0 = c(1, 0.9, 0.8, 0.8, 0.9, 1)
g = c(g0, rep(0, P - length(g0)))
num.nzeros = length(g0)
num.zeros = P - num.nzeros
Z2 = mvrnorm(N, rep(0, num.zeros), Sigma = matrix(matrix_vec(num.zeros), ncol = num.zeros))
X = mvrnorm(N, rep(0, n.beta0 + num.nzeros), 
            Sigma = matrix(matrix_vec(n.beta0 + num.nzeros), ncol = n.beta0 + 
                             num.nzeros))
Z = cbind(X[, (n.beta0 + 1):(n.beta0 + num.nzeros)], Z2)
X = X[, 1:n.beta0]
mu = as.vector(crossprod(t(X), beta0) + crossprod(t(Z), g))
sn = 1.5
part1 = sample(seq(1, N), size = N / 2)
part2 = seq(1, N)[-part1]
R11 = list()
R11$R1 = diag(rep(1, N / 2)) - 
  X[part1] %*% solve(crossprod(X[part1], X[part1])) %*% t(X[part1])
R11$R2 = diag(rep(1, N / 2)) - 
  X[part2] %*% solve(crossprod(X[part2], X[part2])) %*% t(X[part2])
D = sqrt(diag(solve(crossprod(X, X))))
z = scaled_X(X = cbind(X, Z))
zx = dlasso_X(Z, X)



####simulate with 500 datesets ####
Oracle = c()
parEst_O_cv_mcp = c()
parEst_O_bic_mcp = c()
parEst_O_cv_lasso = c()
parEst_O_bic_lasso = c()
parEst_I_cv_mcp = c()
parEst_I_bic_mcp = c()
parEst_I_cv_lasso = c()
parEst_I_bic_lasso = c()
dlasso_cv_beta = c()
dlasso_bic_beta = c()
scaled_beta = c()

sd_Oracle = c()
parEst_sd_cv_mcp = c()
parEst_sd_bic_mcp = c()
parEst_sd_cv_lasso = c()
parEst_sd_bic_lasso = c()
dlasso_cv_sd = c()
dlasso_bic_sd = c()
scaled_sd = c()

beta_cv_mcp = c()
beta_bic_mcp = c()
beta_cv_lasso = c()
beta_bic_lasso = c()
betaO_cv_mcp = c()
betaO_bic_mcp = c()
betaO_cv_lasso = c()
betaO_bic_lasso = c()

sigma_cv_mcp = c()
sigma_bic_mcp = c()
sigma_cv_lasso = c()
sigma_bic_lasso = c()
sigma_scaled = c()
sigma_dlasso_cv = c()
sigma_dlasso_bic = c()
sigma_Oracle = c()

for (i in 1:num.simu)
{
  e = rnorm(N, 0, sd(mu) / sqrt(sn))
  y = mu + e
  res0 = lm(y ~ X + Z[, 1:num.nzeros] - 1)
  res1 = partEst(Z, y, X = X, part1 = part1)
  res2 = dlasso(Z, y, X, zx, part1 = part1, n_lambda = 1000, sigma0 = 1, 
                eps = 1e-3, max.iter = 1000)
  res3 = scaled(cbind(X, Z), y, z, n_lambda = 1000, sigma0 = sd(mu) / sqrt(sn), 
                eps = 1e-3, max.iter = 1000, contrast = c(1, rep(0, P)))
  beta_cv_mcp = cbind(beta_cv_mcp, res1$beta_cv_mcp)
  beta_bic_mcp = cbind(beta_bic_mcp, res1$beta_bic_mcp)
  beta_cv_lasso = cbind(beta_cv_lasso, res1$beta_cv_lasso)
  beta_bic_lasso = cbind(beta_bic_lasso, res1$beta_bic_lasso)
  betaO_cv_mcp = cbind(betaO_cv_mcp, res1$ZOLS_cv_mcp)
  betaO_bic_mcp = cbind(betaO_bic_mcp, res1$ZOLS_bic_mcp)
  betaO_cv_lasso = cbind(betaO_cv_lasso, res1$ZOLS_cv_lasso)
  betaO_bic_lasso = cbind(betaO_bic_lasso, res1$ZOLS_bic_lasso)
  
  Oracle = cbind(Oracle, res0$coef[1:n.beta0])
  parEst_O_cv_mcp = cbind(parEst_O_cv_mcp, res1$gammaOLS_cv_mcp)
  parEst_I_cv_mcp = cbind(parEst_I_cv_mcp, res1$gammaILS_cv_mcp)
  parEst_O_cv_lasso = cbind(parEst_O_cv_lasso, res1$gammaOLS_cv_lasso)
  parEst_I_cv_lasso = cbind(parEst_I_cv_lasso, res1$gammaILS_cv_lasso)
  parEst_O_bic_mcp = cbind(parEst_O_bic_mcp, res1$gammaOLS_bic_mcp)
  parEst_I_bic_mcp = cbind(parEst_I_bic_mcp, res1$gammaILS_bic_mcp)
  parEst_O_bic_lasso = cbind(parEst_O_bic_lasso, res1$gammaOLS_bic_lasso)
  parEst_I_bic_lasso = cbind(parEst_I_bic_lasso, res1$gammaILS_bic_lasso)
  sd_Oracle = cbind(sd_Oracle, sqrt(diag(vcov(res0))[1:n.beta0]))
  parEst_sd_cv_mcp = cbind(parEst_sd_cv_mcp, res1$sd_cv_mcp)
  parEst_sd_cv_lasso = cbind(parEst_sd_cv_lasso, res1$sd_cv_lasso)
  parEst_sd_bic_mcp = cbind(parEst_sd_bic_mcp, res1$sd_bic_mcp)
  parEst_sd_bic_lasso = cbind(parEst_sd_bic_lasso, res1$sd_bic_lasso)
  dlasso_cv_beta = cbind(dlasso_cv_beta, res2$beta_cv)
  dlasso_bic_beta = cbind(dlasso_bic_beta, res2$beta_bic)
  dlasso_cv_sd = cbind(dlasso_cv_sd, res2$sd_cv)
  dlasso_bic_sd = cbind(dlasso_bic_sd, res2$sd_bic)
  scaled_beta = cbind(scaled_beta, res3$beta)
  scaled_sd = cbind(scaled_sd, sqrt(res3$var))
  
  sigma_cv_mcp = cbind(sigma_cv_mcp, res1$sigma_cv_mcp)
  sigma_bic_mcp = cbind(sigma_bic_mcp, res1$sigma_bic_mcp)
  sigma_cv_lasso = cbind(sigma_cv_lasso, res1$sigma_cv_lasso)
  sigma_bic_lasso = cbind(sigma_bic_lasso, res1$sigma_bic_lasso)
  sigma_scaled = cbind(sigma_scaled, res3$sigma)
  sigma_dlasso_cv = cbind(sigma_dlasso_cv, res2$sigma_cv)
  sigma_dlasso_bic = cbind(sigma_dlasso_bic, res2$sigma_bic)
  sigma_Oracle = cbind(sigma_Oracle, sqrt(sum(res0$resid ^ 2) / res0$df.resid))
  print(paste(i, " out of ", num.simu))
}
Oracle = Oracle - beta0
parEst_O_cv_mcp = parEst_O_cv_mcp - beta0
parEst_I_cv_mcp = parEst_I_cv_mcp - beta0
parEst_O_cv_lasso = parEst_O_cv_lasso - beta0
parEst_I_cv_lasso = parEst_I_cv_lasso - beta0
parEst_O_bic_mcp = parEst_O_bic_mcp - beta0
parEst_I_bic_mcp = parEst_I_bic_mcp - beta0
parEst_O_bic_lasso = parEst_O_bic_lasso - beta0
parEst_I_bic_lasso = parEst_I_bic_lasso - beta0
dlasso_cv_beta = dlasso_cv_beta - beta0
dlasso_bic_beta = dlasso_bic_beta - beta0
scaled_beta = scaled_beta - beta0

pmatrix1 = matrix(rep(0, 5 * 12), nrow = 12)
pmatrix1[1, 1] = round(mean(abs(parEst_O_cv_mcp)), 2)
pmatrix1[2, 1] = round(mean(abs(parEst_O_bic_mcp)), 2)
pmatrix1[3, 1] = round(mean(abs(parEst_O_cv_lasso)), 2)
pmatrix1[4, 1] = round(mean(abs(parEst_O_bic_lasso)), 2)
pmatrix1[5, 1] = round(mean(abs(parEst_I_cv_mcp)), 2)
pmatrix1[6, 1] = round(mean(abs(parEst_I_bic_mcp)), 2)
pmatrix1[7, 1] = round(mean(abs(parEst_I_cv_lasso)), 2)
pmatrix1[8, 1] = round(mean(abs(parEst_I_bic_lasso)), 2)
pmatrix1[9, 1] = round(mean(abs(dlasso_cv_beta)), 2)
pmatrix1[10, 1] = round(mean(abs(dlasso_bic_beta)), 2)
pmatrix1[11, 1] = round(mean(abs(scaled_beta)), 2)
pmatrix1[12, 1] = round(mean(abs(Oracle)), 2)
pmatrix1[1, 2] = round(mean(parEst_sd_cv_mcp), 2)
pmatrix1[2, 2] = round(mean(parEst_sd_bic_mcp), 2)
pmatrix1[3, 2] = round(mean(parEst_sd_cv_lasso), 2)
pmatrix1[4, 2] = round(mean(parEst_sd_bic_lasso), 2)
pmatrix1[5, 2] = round(mean(parEst_sd_cv_mcp), 2)
pmatrix1[6, 2] = round(mean(parEst_sd_bic_mcp), 2)
pmatrix1[7, 2] = round(mean(parEst_sd_cv_lasso), 2)
pmatrix1[8, 2] = round(mean(parEst_sd_bic_lasso), 2)
pmatrix1[9, 2] = round(mean(dlasso_cv_sd), 2)
pmatrix1[10, 2] = round(mean(dlasso_bic_sd), 2)
pmatrix1[11, 2] = round(mean(scaled_sd), 2)
pmatrix1[12, 2] = round(mean(sd_Oracle), 2)
pmatrix1[1, 3] = round(mean(parEst_sd_cv_mcp ^ 2 + parEst_O_cv_mcp ^ 2), 2)
pmatrix1[2, 3] = round(mean(parEst_sd_bic_mcp ^ 2 + parEst_O_bic_mcp ^ 2), 2)
pmatrix1[3, 3] = round(mean(parEst_sd_cv_lasso ^ 2 + parEst_O_cv_lasso ^ 2), 2)
pmatrix1[4, 3] = round(mean(parEst_sd_bic_lasso ^ 2 + parEst_O_bic_lasso ^ 2), 2)
pmatrix1[5, 3] = round(mean(parEst_sd_cv_mcp ^ 2 + parEst_I_cv_mcp ^ 2), 2)
pmatrix1[6, 3] = round(mean(parEst_sd_bic_mcp ^ 2 + parEst_I_bic_mcp ^ 2), 2)
pmatrix1[7, 3] = round(mean(parEst_sd_cv_lasso ^ 2 + parEst_I_cv_lasso ^ 2), 2)
pmatrix1[8, 3] = round(mean(parEst_sd_bic_lasso ^ 2 + parEst_I_bic_lasso ^ 2), 2)
pmatrix1[9, 3] = round(mean(dlasso_cv_sd ^ 2 + dlasso_cv_beta ^ 2), 2)
pmatrix1[10, 3] = round(mean(dlasso_bic_sd ^ 2 + dlasso_bic_beta ^ 2), 2)
pmatrix1[11, 3] = round(mean(scaled_sd ^ 2 + scaled_beta ^ 2), 2)
pmatrix1[12, 3] = round(mean(sd_Oracle ^ 2 + Oracle ^ 2), 2)

zval = qnorm(0.975, 0, 1)
coverage = c(mean((parEst_O_cv_mcp + parEst_sd_cv_mcp * zval >= 0) & 
                     (parEst_O_cv_mcp - parEst_sd_cv_mcp * zval <= 0)),
             mean((parEst_O_bic_mcp + parEst_sd_bic_mcp * zval >= 0) & 
                     (parEst_O_bic_mcp - parEst_sd_bic_mcp * zval <= 0)),
             mean((parEst_O_cv_lasso + parEst_sd_cv_lasso * zval >= 0) & 
                     (parEst_O_cv_lasso - parEst_sd_cv_lasso * zval <= 0)),
             mean((parEst_O_bic_lasso + parEst_sd_bic_lasso * zval >= 0) & 
                     (parEst_O_bic_lasso - parEst_sd_bic_lasso * zval <= 0)),
             mean((parEst_I_cv_mcp + parEst_sd_cv_mcp * zval >= 0) & 
                     (parEst_I_cv_mcp - parEst_sd_cv_mcp * zval <= 0)),
             mean((parEst_I_bic_mcp + parEst_sd_bic_mcp * zval >= 0) & 
                     (parEst_I_bic_mcp - parEst_sd_bic_mcp * zval <= 0)),
             mean((parEst_I_cv_lasso + parEst_sd_cv_lasso * zval >= 0) & 
                     (parEst_I_cv_lasso - parEst_sd_cv_lasso * zval <= 0)),
             mean((parEst_I_bic_lasso + parEst_sd_bic_lasso * zval >= 0) & 
                     (parEst_I_bic_lasso - parEst_sd_bic_lasso * zval <= 0)),
             mean((dlasso_cv_beta + dlasso_cv_sd * zval >= 0) & 
                     (dlasso_cv_beta - dlasso_cv_sd * zval <= 0)),
             mean((dlasso_bic_beta + dlasso_bic_sd * zval >= 0) & 
                     (dlasso_bic_beta - dlasso_bic_sd * zval <= 0)),
             mean((scaled_beta + scaled_sd * zval >= 0) & 
                     (scaled_beta - scaled_sd * zval <= 0)),
             mean((Oracle + sd_Oracle * zval >= 0) & 
                    (Oracle - sd_Oracle * zval <= 0)))

pmatrix1[, 4] = t(coverage)
pmatrix1[1, 5] = round(sd(as.numeric(parEst_O_cv_mcp)), 2)
pmatrix1[2, 5] = round(sd(as.numeric(parEst_O_bic_mcp)), 2)
pmatrix1[3, 5] = round(sd(as.numeric(parEst_O_cv_lasso)), 2)
pmatrix1[4, 5] = round(sd(as.numeric(parEst_O_bic_lasso)), 2)
pmatrix1[5, 5] = round(sd(as.numeric(parEst_I_cv_mcp)), 2)
pmatrix1[6, 5] = round(sd(as.numeric(parEst_I_bic_mcp)), 2)
pmatrix1[7, 5] = round(sd(as.numeric(parEst_I_cv_lasso)), 2)
pmatrix1[8, 5] = round(sd(as.numeric(parEst_I_bic_lasso)), 2)
pmatrix1[9, 5] = round(sd(as.numeric(dlasso_cv_beta)), 2)
pmatrix1[10, 5] = round(sd(as.numeric(dlasso_bic_beta)), 2)
pmatrix1[11, 5] = round(sd(as.numeric(scaled_beta)), 2)
pmatrix1[12, 5] = round(sd(as.numeric(Oracle)), 2)
colnames(pmatrix1) = c("Bias", "SD", "MSE", 
                       "Coverage", "Simulation SD")
rownames(pmatrix1) = c("OLS CV MCP", "OLS BIC MCP", "OLS CV LASSO",
                       "OLS BIC LASSO", 
                       "ILS CV MCP", "ILS BIC MCP", 
                       "ILS CV LASSO", "ILS BIC LASSO", 
                       "DLASSO CV",
                       "DLASSO BIC", "SLASSO", "Oracle")

Xbeta_cv_mcp = X %*% solve(crossprod(X, X)) %*% t(X) %*% 
  crossprod(t(Z), beta_cv_mcp - g)
Xbeta_bic_mcp = X %*% solve(crossprod(X, X)) %*% t(X) %*%
  crossprod(t(Z), beta_bic_mcp - g)
Xbeta_cv_lasso = X %*% solve(crossprod(X, X)) %*% t(X) %*%
  crossprod(t(Z), beta_cv_lasso - g)
Xbeta_bic_lasso = X %*% solve(crossprod(X, X)) %*% t(X) %*%
  crossprod(t(Z), beta_bic_lasso - g)
XbetaO_cv_mcp = X %*% solve(crossprod(X, X)) %*% t(X) %*%
  crossprod(t(Z), betaO_cv_mcp - g)
XbetaO_bic_mcp = X %*% solve(crossprod(X, X)) %*% t(X) %*%
  crossprod(t(Z), betaO_bic_mcp - g)
XbetaO_cv_lasso = X %*% solve(crossprod(X, X)) %*% t(X) %*%
  crossprod(t(Z), betaO_cv_lasso - g)
XbetaO_bic_lasso = X %*% solve(crossprod(X, X)) %*% t(X) %*%
  crossprod(t(Z), betaO_bic_lasso - g)



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


xlim = 2
ylim = 2.5
tmp <- Xbeta_cv_mcp[-(1:num.nzeros), ]
Init <- data.frame(Error = c(as.vector(tmp[tmp != 0]), 
                             as.vector(Xbeta_cv_mcp[1:num.nzeros, ])))
tmp <- XbetaO_cv_mcp[-(1:num.nzeros), ]
OrdinaryD <- data.frame(Error = c(as.vector(tmp[tmp != 0]), 
                                 as.vector(XbetaO_cv_mcp[1:num.nzeros, ])))
Init$tag <- 'Initial'
OrdinaryD$tag <- 'Ordinary'
ErrorDiff <- rbind(OrdinaryD, Init)
p1 <- ggplot(ErrorDiff, aes(Error, fill = tag)) + 
  geom_density(alpha = 0.2) + labs(title = "CV-MCP") + 
  coord_cartesian(xlim=c(-xlim, xlim), ylim = c(0, ylim))

tmp <- Xbeta_bic_mcp[-(1:num.nzeros), ]
Init <- data.frame(Error = c(as.vector(tmp[tmp != 0]), 
                             as.vector(Xbeta_bic_mcp[1:num.nzeros, ])))
tmp <- XbetaO_bic_mcp[-(1:num.nzeros), ]
OrdinaryD <- data.frame(Error = c(as.vector(tmp[tmp != 0 ]),
                                 as.vector(XbetaO_bic_mcp[1:num.nzeros,])))
Init$tag <- 'Initial'
OrdinaryD$tag <- 'Ordinary'
ErrorDiff <- rbind(OrdinaryD, Init)
p2 <- ggplot(ErrorDiff, aes(Error, fill = tag)) + 
  geom_density(alpha = 0.2) + labs(title = "BIC-MCP") + 
  coord_cartesian(xlim=c(-xlim, xlim), ylim = c(0, ylim))

tmp <- Xbeta_cv_lasso[-(1:num.nzeros), ]
Init <- data.frame(Error = c(as.vector(tmp[tmp != 0]),
                             as.vector(Xbeta_cv_lasso[1:num.nzeros,])))
tmp <- XbetaO_cv_lasso[-(1:num.nzeros), ]
OrdinaryD <- data.frame(Error = c(as.vector(tmp[tmp != 0]),
                                 as.vector(XbetaO_cv_lasso[1:num.nzeros,])))
Init$tag <- 'Initial'
OrdinaryD$tag <- 'Ordinary'
ErrorDiff <- rbind(OrdinaryD, Init)
p3 <- ggplot(ErrorDiff, aes(Error, fill = tag)) + 
  geom_density(alpha = 0.2) + labs(title = "CV-LASSO") + 
  coord_cartesian(xlim=c(-xlim, xlim), ylim = c(0, ylim))

tmp <- Xbeta_bic_lasso[-(1:num.nzeros), ]
Init <- data.frame(Error = c(as.vector(tmp[tmp != 0]),
                             as.vector(Xbeta_bic_lasso[1:num.nzeros,])))
tmp <- XbetaO_bic_lasso[-(1:num.nzeros), ]
OrdinaryD <- data.frame(Error = c(as.vector(tmp[tmp != 0]),
                                 as.vector(XbetaO_bic_lasso[1:num.nzeros,])))
Init$tag <- 'Initial'
OrdinaryD$tag <- 'Ordinary'
ErrorDiff <- rbind(OrdinaryD, Init)
p4 <- ggplot(ErrorDiff, aes(Error, fill = tag)) + 
  geom_density(alpha = 0.2) + labs(title = "BIC-LASSO") + 
  coord_cartesian(xlim=c(-xlim, xlim), ylim = c(0, ylim))

multiplot(p1, p2, p3, p4, cols=2)
sigma_cv_mcp = sigma_cv_mcp / sd(mu) * sqrt(sn)
sigma_bic_mcp = sigma_bic_mcp / sd(mu) * sqrt(sn)
sigma_cv_lasso = sigma_cv_lasso / sd(mu) * sqrt(sn)
sigma_bic_lasso = sigma_bic_lasso / sd(mu) * sqrt(sn)
sigma_dlasso_cv = sigma_dlasso_cv / sd(mu) * sqrt(sn)
sigma_dlasso_bic = sigma_dlasso_bic / sd(mu) * sqrt(sn)
sigma_scaled = sigma_scaled / sd(mu) * sqrt(sn)
sigma_Oracle = sigma_Oracle / sd(mu) * sqrt(sn)
P_cv_mcp = data.frame(ratio = sigma_cv_mcp)
P_bic_mcp = data.frame(ratio = sigma_bic_mcp)
P_cv_lasso = data.frame(ratio = sigma_cv_lasso)
P_bic_lasso = data.frame(ratio = sigma_bic_lasso)
D_cv = data.frame(ratio = sigma_dlasso_cv)
D_bic = data.frame(ratio = sigma_dlasso_bic)
Scaled = data.frame(ratio = sigma_scaled)
OracleD = data.frame(ratio = sigma_Oracle)
P_cv_mcp$tag = "CV-MCP"
P_bic_mcp$tag = "BIC-MCP"
P_cv_lasso$tag = "CV-LASSO"
P_bic_lasso$tag = "BIC-LASSO"
D_cv$tag = "DLASSO-CV"
D_bic$tag = "DLASSO-BIC"
Scaled$tag = "Scaled"
OracleD$tag = "Oracle"
ErrorRatio = rbind(P_cv_mcp, P_bic_mcp, P_cv_lasso, P_bic_lasso, OracleD)
ggplot(ErrorRatio, aes(ratio, fill = tag)) + 
  geom_density(alpha = 0.4) + labs(title = "SigmaHat / Sigma") + 
  scale_fill_manual( values = c("red", "black", "yellow", "blue"))
ErrorRatio = rbind(D_cv, D_bic, Scaled, OracleD)
ggplot(ErrorRatio, aes(ratio, fill = tag)) + 
  geom_density(alpha = 0.4) + labs(title = "SigmaHat / Sigma") + 
  scale_fill_manual( values = c("red", "black", "yellow", "blue"))
MISS = c(cv_mcp_miss = mean(colSums(rbind(beta_cv_mcp[1:num.nzeros, ] == 0, 
              beta_cv_mcp[(num.nzeros + 1):P, ] != 0))) / P,
bic_mcp_miss = mean(colSums(rbind(beta_bic_mcp[1:num.nzeros, ] == 0, 
              beta_bic_mcp[(num.nzeros + 1):P, ] != 0))) / P,
cv_lasso_miss = mean(colSums(rbind(beta_cv_lasso[1:num.nzeros, ] == 0, 
              beta_cv_lasso[(num.nzeros + 1):P, ] != 0))) / P,
bic_lasso_miss = mean(colSums(rbind(beta_bic_lasso[1:num.nzeros, ] == 0, 
              beta_bic_lasso[(num.nzeros + 1):P, ] != 0))) / P)
save.image("simulation1.Rdata")
