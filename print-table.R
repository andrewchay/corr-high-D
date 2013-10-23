load("S:/research/291simu/cov0.Rdata")
library(xtable)
table.M <- matrix(nrow = 8, ncol = 3)
table.D <- matrix(nrow = 8, ncol = 3)
table.M[,1] <- paste(round(comb.M[, 1], 3), "(", round(comb.Msd[, 1], 3), ")", sep = "")
rownames(table.M) <- rownames(comb.M)
table.D[,1] <- paste(round(comb.D[, 1], 3), "(", round(comb.Dsd[, 1], 3), ")", sep = "")
rownames(table.D) <- rownames(comb.D)

load("S:/research/291simu/cov1.Rdata")

table.M[,2] <- paste(round(comb.M[, 2], 3), "(", round(comb.Msd[, 2], 3), ")", sep = "")
rownames(table.M) <- rownames(comb.M)
table.D[,2] <- paste(round(comb.D[, 2], 3), "(", round(comb.Dsd[, 2], 3), ")", sep = "")
rownames(table.D) <- rownames(comb.D)

load("S:/research/291simu/cov2.Rdata")

table.M[,3] <- paste(round(comb.M[, 3], 3), "(", round(comb.Msd[, 3], 3), ")", sep = "")
rownames(table.M) <- rownames(comb.M)
table.D[,3] <- paste(round(comb.D[, 3], 3), "(", round(comb.Dsd[, 3], 3), ")", sep = "")
rownames(table.D) <- rownames(comb.D)

colnames(table.M) <- c("cov0", "cov1", "cov2")
colnames(table.D) <- c("cov0", "cov1", "cov2")
xtable(table.M)
xtable(table.D)