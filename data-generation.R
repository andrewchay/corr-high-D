# library(CCSparsenet, lib.loc = "/mnt/nfs/netapp2/grad/hchai/Rpackages")
require("MASS")
set.seed = 100

N = 200
P = 500
num.simu = 20
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
beta0 = c(1, -1)
n.beta0 = length(beta0)
g0 = c(1, 0.9, 0.8, 0.8, 0.9, 1)
g = c(g0, rep(0, P - length(g0)))
nonzero.ind = g != 0
num.nzeros = length(g0)
num.zeros = P - num.nzeros
Z2 = mvrnorm(N, rep(0, num.zeros), Sigma = matrix(matrix_vec(num.zeros), ncol = num.zeros))
X = mvrnorm(N, rep(0, n.beta0 + num.nzeros), 
            Sigma = matrix(matrix_vec(n.beta0 + num.nzeros), ncol = n.beta0 + 
                             num.nzeros))
Z = cbind(X[, (n.beta0 + 1):(n.beta0 + num.nzeros)], Z2)
X = X[, 1:n.beta0]
mu = exp(as.vector(crossprod(t(X), beta0) + crossprod(t(Z), g)))
part1 = sample(seq(1, N), size = N / 2)
part2 = seq(1, N)[-part1]
center = "New"

