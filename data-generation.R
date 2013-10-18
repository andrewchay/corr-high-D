# library(CCSparsenet, lib.loc = "/mnt/nfs/netapp2/grad/hchai/Rpackages")
require("MASS")
set.seed = 123

N = 200
P = 120
num.simu = 80
rho = 0.5
matrix_vec = function(P){
  tmp = matrix(rep(0, P * P), nrow = P)
  for (i in 1:P) 
    for (j in 1:P)
      if (i == j) tmp[i, j] = 1 else 
        tmp[i, j] = (rho) ^ abs(i - j)
  #tmp[i, j] = rho
  return(tmp)
}
beta0 = c(1, -1)
n.beta0 = length(beta0)
num.zeros = P - n.beta0
beta = c(beta0, rep(0, num.zeros))
non.zero.index = (beta != 0)
X = mvrnorm(N * P, rep(0, P), Sigma = matrix(matrix_vec(P), ncol = P))
mu = as.vector(crossprod(t(X), beta))
center = "Old"
sn = 1.5

y = mu + rnorm(N, 0, sd(mu) / sqrt(sn))