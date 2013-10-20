# This simulation compares all three methods. dlasso stands for double LASSO and post estimation method. scaled is the scaled LASSO by Zhang and partEst is our method.
library("ggplot2")
library("reshape")
setwd("D:/Lecture/Homework/291/corr-high-D")
source("data-generation.R")
source("CDA.R")
####simulate with 500 datesets ####
Oracle = c()
cv_mcp = c()
cv_lasso = c()
stat.names = c("Ltype1", "Ltype2", "LFDR", "LFIR", "LL1", "LL2", "LXL1", "LXL2",
               "Mtype1", "Mtype2", "MFDR", "MFIR", "ML1", "ML2", "MXL1", "MXL2")

stats = c()
lasso.lam = c()
mcp.lam = c()
Ltype1 = c()
Ltype2 = c()
LFDR = c()
LFIR = c()
Mtype1 = c()
Mtype2 = c()
MFDR = c()
MFIR = c()
LL1 = c()
LL2 = c()
LXL1 = c()
LXL2 = c()
ML1 = c()
ML2 = c()
MXL1 = c()
MXL2 = c()

for (i in 1:num.simu)
{
  e = rnorm(N, 0, sd(mu) / sqrt(sn))
  y = mu + e
  res0 = lm(y ~ X[, non.zero.index] - 1)
  res1 = CDgrid(X, y)
  
  Oracle = cbind(Oracle, res0$coef[1:n.beta0])
  cv_mcp = cbind(cv_mcp, res1$beta_mcp)
  cv_lasso = cbind(cv_lasso, res1$beta_lasso)
  
  lasso.lam = c(lasso.lam, res1$lasso.lam)
  mcp.lam = c(mcp.lam, res1$mcp.lam)
  Ltype1 = c(Ltype1, res1$Ltype1)
  Ltype2 = c(Ltype2, res1$Ltype2)
  LFDR = c(LFDR, res1$LFDR)
  LFIR = c(LFIR, res1$LFIR)
  Mtype1 = c(Mtype1, res1$Mtype1)
  Mtype2 = c(Mtype2, res1$Mtype2)
  MFDR = c(MFDR, res1$MFDR)
  MFIR = c(MFIR, res1$MFIR)
  LL1 = c(LL1, res1$LL1)
  LL2 = c(LL2, res1$LL2)
  LXL1 = c(LXL1, res1$LXL1)
  LXL2 = c(LXL2, res1$LXL2)
  ML1 = c(ML1, res1$ML1)
  ML2 = c(ML2, res1$ML2)
  MXL1 = c(MXL1, res1$MXL1)
  MXL2 = c(MXL2, res1$MXL2)
  stats = rbind(stats, res1$stats)

  print(paste(i, " out of ", num.simu))
}
Oracle = Oracle - beta0
cv_mcp = cv_mcp - beta
cv_lasso = cv_lasso - beta
stats = data.frame(stats)


pmatrix1 = matrix(rep(0, 5 * 5), nrow = 5)
pmatrix1[1, 1] = round(mean(abs(cv_mcp)), 2)
pmatrix1[2, 1] = round(mean(abs(cv_lasso)), 2)
pmatrix1[5, 1] = round(mean(abs(Oracle)), 2)

pmatrix1[1, 5] = round(sd(as.numeric(cv_mcp)), 2)
pmatrix1[2, 5] = round(sd(as.numeric(cv_lasso)), 2)
pmatrix1[5, 5] = round(sd(as.numeric(Oracle)), 2)
colnames(pmatrix1) = c("Bias", "SD", "MSE", 
                       "Coverage", "Simulation SD")
rownames(pmatrix1) = c("OLS CV MCP", "OLS CV LASSO",
                       "ILS CV MCP", "ILS CV LASSO",
                       "Oracle")

# plot(stats$Mtype1 ~ stats$lambda, pch = 1, col = 4, pty = 1, cex = .5)
# points(stats$Ltype1 ~ stats$lambda, pch = ".", col = 2, pty = 1, cex = 3)
# plot(stats$Mtype2 ~ stats$lambda, pch = ".", col = 4, pty = 1, cex = 3)
# points(stats$Ltype2 ~ stats$lambda, pch = ".", col = 2, pty = 1, cex = 3)
# plot(stats$MFDR ~ stats$lambda, pch = ".", col = 4, pty = 1, cex = 3)
# points(stats$LFDR ~ stats$lambda, pch = ".", col = 2, pty = 1, cex = 3)
# plot(stats$MFIR ~ stats$lambda, pch = ".", col = 4, pty = 1, cex = 3)
# points(stats$LFIR ~ stats$lambda, pch = ".", col = 2, pty = 1, cex = 3)
# 
# plot(stats$ML1 ~ stats$lambda, pch = 1, col = 4, pty = 1, cex = .5)
# points(stats$LL1 ~ stats$lambda, pch = ".", col = 2, pty = 1, cex = 3)
# plot(stats$ML2 ~ stats$lambda, pch = ".", col = 4, pty = 1, cex = 3)
# points(stats$LL2 ~ stats$lambda, pch = ".", col = 2, pty = 1, cex = 3)
# plot(stats$MXL1 ~ stats$lambda, pch = ".", col = 4, pty = 1, cex = 3)
# points(stats$LXL1 ~ stats$lambda, pch = ".", col = 2, pty = 1, cex = 3)
# plot(stats$MXL2 ~ stats$lambda, pch = ".", col = 4, pty = 1, cex = 3)
# points(stats$LXL2 ~ stats$lambda, pch = ".", col = 2, pty = 1, cex = 3)

tmp.data = melt(stats[, c(1:4, 21)], id = c("lambda"))
ggplot(data = tmp.data, aes(x = lambda, y = value, group = variable)) + 
  geom_smooth(method = "loess", size = 1.5, aes(colour = variable), size = .5) + 
  geom_vline(xintercept = mean(lasso.lam)) + ggtitle("LASSO Selection") + 
  xlab("lambda") + ylab("rate")
tmp.data = melt(stats[, c(5:8, 21)], id = c("lambda"))
ggplot(data = tmp.data, aes(x = lambda, y = value, group = variable)) + 
  geom_smooth(method = "loess", size = 1.5, aes(colour = variable), size = .5) + 
  geom_vline(xintercept = mean(mcp.lam)) + ggtitle("MCP Selection") + 
  xlab("lambda") + ylab("rate")
tmp.data = melt(stats[, c(9:14, 21)], id = c("lambda"))
ggplot(data = tmp.data, aes(x = lambda, y = value, group = variable)) + 
  geom_smooth(method = "loess", size = 1.5, aes(colour = variable), size = .5) + 
  geom_vline(xintercept = mean(lasso.lam)) + ggtitle("LASSO Distance") + 
  xlab("lambda") + ylab("rate")
tmp.data = melt(stats[, c(10:12, 14, 21)], id = c("lambda"))
ggplot(data = tmp.data, aes(x = lambda, y = value, group = variable)) + 
  geom_smooth(method = "loess", size = 1.5, aes(colour = variable), size = .5) + 
  geom_vline(xintercept = mean(lasso.lam)) + ggtitle("MCP Distance") + 
  xlab("lambda") + ylab("rate")
tmp.data = melt(stats[, c(15:20, 21)], id = c("lambda"))
ggplot(data = tmp.data, aes(x = lambda, y = value, group = variable)) + 
  geom_smooth(method = "loess", size = 1.5, aes(colour = variable), size = .5) + 
  geom_vline(xintercept = mean(mcp.lam)) + ggtitle("LASSO Distance") + 
  xlab("lambda") + ylab("rate")
tmp.data = melt(stats[, c(16:18, 20, 21)], id = c("lambda"))
ggplot(data = tmp.data, aes(x = lambda, y = value, group = variable)) + 
  geom_smooth(method = "loess", size = 1.5, aes(colour = variable), size = .5) + 
  geom_vline(xintercept = mean(mcp.lam)) + ggtitle("MCP Distance") + 
  xlab("lambda") + ylab("rate")


comb.M <- matrix(0, 8, 3)
rownames(comb.M) <- stat.names[c(1:4, 9:12)]
colnames(comb.M) <- c("cov0", "cov1", "cov2")
for (i in stat.names[c(1:4, 9:12)])
{
  comb.M[i, method] = eval(parse(text = paste("mean(", i, ")", sep = "")))
}
hist(cor(X), main = paste("Histogram of ", method, sep = ""))

save.image(paste("./", method, ".Rdata", sep = ""))