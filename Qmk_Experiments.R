#############################
# Consistency Factor of Qm^k
#############################

# Author: Santiago Ortiz (saortizar@unal.edu.co)
# Date: 30/12/2025

# Load Package
library(MASS)
set.seed(123)

# An Example 1
p = 2
n = 1000
rho = 0.5
X = mvrnorm(n, rep(0, p), matrix(c(1, rho, rho, 1), p, p))
pairs(X, pch = 20)
S = cov(X[,1], X[,2])
Qnm = biCovQn(X[,1], X[,2], 1)

# Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load codes and auxiliary programms
source("Code_Qmk.R")

########################
# Simulation Experiment
########################

# Example 2
rhos = seq(-1, 1, by = 0.01)
Tn0 = Tn1 = Tn2 = Tn3 = Tn4 = rep(0, length(rhos))
for (i in 1:length(rhos)) {
  X = mvrnorm(n, rep(0, p), matrix(c(1, rhos[i], rhos[i], 1), p, p))
  S = cov(X[,1], X[,2])
  Qnm0 = biCovQn(X[,1], X[,2], 0)
  Qnm1 = biCovQn(X[,1], X[,2], 1)
  Qnm2 = biCovQn(X[,1], X[,2], 2)
  Qnm3 = biCovQn(X[,1], X[,2], 3)
  Qnm4 = biCovQn(X[,1], X[,2], 4)
  Tn0[i] = Qnm0/S
  Tn1[i] = Qnm1/S
  Tn2[i] = Qnm2/S
  Tn3[i] = Qnm3/S
  Tn4[i] = Qnm4/S
}
plot(rhos, Tn0, pch = 20, xlab = "Correlation", ylab = "Consistency Relation")
grid()
lines(rhos, Tn0, lwd = 2)
#abline(h = median(Tn0), lwd = 2, col = "red")
plot(rhos, Tn1, pch = 20, xlab = "Correlation", ylab = "Consistency Relation")
grid()
lines(rhos, Tn1, lwd = 2)
plot(rhos, Tn2, pch = 20, xlab = "Correlation", ylab = "Consistency Relation")
grid()
lines(rhos, Tn2, lwd = 2)
plot(rhos, Tn3, pch = 20, xlab = "Correlation", ylab = "Consistency Relation")
grid()
lines(rhos, Tn3, lwd = 2)
plot(rhos, Tn4, pch = 20, xlab = "Correlation", ylab = "Consistency Relation")
grid()
lines(rhos, Tn4, lwd = 2)

# Example 3
rhos = seq(-1, 1, by = 0.01)
S = Qnm.0 = Qnm.1 = Qnm.2 = Qnm.3 = Qnm.4 = rep(0, length(rhos))
for (i in 1:length(rhos)) {
  X = mvrnorm(n, rep(0, p), matrix(c(1, rhos[i], rhos[i], 1), p, p))
  S[i] = cov(X[,1], X[,2])
  Qnm.0[i] = biCovQn(X[,1], X[,2], 0)
  Qnm.1[i] = biCovQn(X[,1], X[,2], 1)
  Qnm.2[i] = biCovQn(X[,1], X[,2], 2)
  Qnm.3[i] = biCovQn(X[,1], X[,2], 3)
  Qnm.4[i] = biCovQn(X[,1], X[,2], 4)
}
plot(Qnm.0, S, pch = 20, xlim = c(-1,1), ylim = c(-1,1))
grid()
abline(reg = c(0,1), lwd = 2, col = "deeppink")
plot(Qnm.1, S, pch = 20, xlim = c(-1,1), ylim = c(-1,1))
grid()
abline(reg = c(0,1), lwd = 2, col = "deeppink")
plot(Qnm.2, S, pch = 20, xlim = c(-1,1), ylim = c(-1,1))
grid()
abline(reg = c(0,1), lwd = 2, col = "deeppink")
plot(Qnm.3, S, pch = 20, xlim = c(-1,1), ylim = c(-1,1))
grid()
abline(reg = c(0,1), lwd = 2, col = "deeppink")
plot(Qnm.4, S, pch = 20, xlim = c(-1,1), ylim = c(-1,1))
grid()
abline(reg = c(0,1), lwd = 2, col = "deeppink")

##########
# Figures
##########

par(mar = c(5,5,2,3))
#par(mar = c(5,4,2,2))
par(mfrow = c(1,2))
plot(rhos, Tn2, pch = 19, xlab = "Correlation Value", ylab = "Consistency Quotient", cex.axis = 1.5, cex.lab = 1.5)
grid()
plot(Qnm.2, S, pch = 19, xlim = c(-1,1), ylim = c(-1,1),
     xlab = expression(Q[m]^(n/4)), ylab = "Sample Covariance", cex.axis = 1.5, cex.lab = 1.5)
grid()
abline(reg = c(0,1), lwd = 3, col = "deeppink")

plot(rhos, Tn3, pch = 19, xlab = "Correlation Value", ylab = "Consistency Quotient", cex.axis = 1.5, cex.lab = 1.5)
grid()
plot(Qnm.3, S, pch = 19, xlim = c(-1,1), ylim = c(-1,1),
     xlab = expression(Q[m]^(n/2)), ylab = "Sample Covariance", cex.axis = 1.5, cex.lab = 1.5)
grid()
abline(reg = c(0,1), lwd = 3, col = "deeppink")

plot(rhos, Tn4, pch = 19, xlab = "Correlation Value", ylab = "Consistency Quotient", cex.axis = 1.5, cex.lab = 1.5)
grid()
plot(Qnm.4, S, pch = 19, xlim = c(-1,1), ylim = c(-1,1),
     xlab = expression(Q[m]^(3*n/4)), ylab = "Sample Covariance", cex.axis = 1.5, cex.lab = 1.5)
grid()
abline(reg = c(0,1), lwd = 3, col = "deeppink")

plot(rhos, Tn0, pch = 19, xlab = "Correlation Value", ylab = "Consistency Quotient", cex.axis = 1.5, cex.lab = 1.5)
grid()
plot(Qnm.1, S, pch = 19, xlim = c(-1,1), ylim = c(-1,1),
     xlab = expression(Q[m]^{(n/4 * "," * 3*n/4)}), ylab = "Sample Covariance", cex.axis = 1.5, cex.lab = 1.5)
grid()
abline(reg = c(0,1), lwd = 3, col = "deeppink")