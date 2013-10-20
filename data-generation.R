# library(CCSparsenet, lib.loc = "/mnt/nfs/netapp2/grad/hchai/Rpackages")
require("MASS")
library("CCSparsenet")
set.seed = 123
setwd("D:/Lecture/Homework/291/corr-high-D")
matrix_vec = function(P, rho){
  tmp = matrix(rep(0, P * P), nrow = P)
  for (i in 1:P) 
    for (j in 1:P)
      if (i == j) tmp[i, j] = 1 else 
        tmp[i, j] = (rho) ^ abs(i - j)
  #tmp[i, j] = rho
  return(tmp)
}
matrix_move <- function(N, P, q, ratio)
# ratio controls the number of columns to be shifted
# q is the bandwidth
{
  pr <- round(P * ratio)
  ind = sort(sample(P)[1:pr])
  tmp <- matrix(rnorm(N * length(ind)), nrow = N)
  M <- matrix(0, N, P)
  M[, ind] <- tmp
  tmp <- M
  M <- cbind(matrix(rep(0, N * q), nrow = N), M, matrix(rep(0, N * q), nrow = N))
  for (i in c(seq(-q, -1), seq(1, q)))
    tmp <- tmp + M[, (q + i + 1):(q + i + P)]
  return(tmp)
}

cov0 <- function(N, P)
{
  matrix(rnorm(N * P, 0, 1), nrow = N)
}
cov1 <- function(N, P, rho = 0.5) # generate the covariates from multivariate normal
{
  X = mvrnorm(N, rep(0, P), Sigma = matrix(matrix_vec(P, rho), ncol = P))
}
cov2 <- function(N, P, ratio = 0.1, a = 2, q = 3) # generate the covariates from latent variables
{
  matrix(rnorm(N * P, 0, 1), nrow = N) + a * matrix_move(N, P, q, ratio)
}

method = "cov2"
N = 100
P = 500
num.simu = 400
beta0 = c(3, 1, -1, -3)
n.beta0 = length(beta0)
num.zeros = P - n.beta0
beta = c(beta0, rep(0, num.zeros))
non.zero.index = (beta != 0)
if (method == "cov0") X = cov0(N, P)
if (method == "cov1") X = cov1(N, P, rho = 0.8)
if (method == "cov2") X = cov2(N, P, a = 2, q = 3)
mu = as.vector(crossprod(t(X), beta))
center = "Old"
sn = 1.5
