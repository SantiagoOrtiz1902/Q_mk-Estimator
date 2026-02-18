############################
# Independence of Q_m^(n\2)
############################

# Author: Santiago Ortiz (saortizar@unal.edu.co)
# Date: 30/12/2025

########################
# Simulation Experiment
########################

# Required Packages
my_packages = c("boot", "MASS", "scales")
not_installed = my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if (length(not_installed)) install.packages(not_installed, dependencies = TRUE)
for (q in 1:length(my_packages)) {library(my_packages[q], character.only = TRUE)}

# Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load codes and auxiliary programms
source("Code_Qmk.R")

# Experiment
set.seed(2025)
n = 1000
mu = c(0,0)
S = diag(2)
m = 500
Qm = rep(0,m)
for (i in 1:m) {
  X.l = mvrnorm(n, mu, S)
  Qm[i] = biCovQn(X.l[,1], X.l[,2], type = 3)
}
media = function(X, idx){df = mean(X[idx]); return(df)}
boot.res = boot(Qm, media, 10000)$t
boot.Qm = round(as.vector(quantile(boot.res, c(0.025, 0.975))), 4)
par(mfrow = c(1,2))
plot(Qm, pch = 20, ylim = c(-0.045,0.045), xlab = "Iteration Index", ylab = "Qm^(n/2)",
     col = alpha("black", 1.0))
grid()
lines(Qm, type = "h", lwd = 2, col = alpha("black", 0.4))
abline(h = 0, col = "deeppink", lwd = 3)
hist(boot.res, freq = F, col = alpha("lightblue4", 0.4), xlab = "Bootstrapped Mean Qm^(n/2)",
     main = "Bootstrapped distribution of mean of Qm^(n/2)\nwith 95% confidence interval")
abline(v = boot.Qm, lty = 3, lwd = 5, col = "deeppink")
legend("topright", legend = boot.Qm, col = "deeppink", lty = 3, lwd = 5)

##############
# Other Graphs
##############

par(mar = c(5, 6, 3, 2))
plot(Qm, pch = 20, ylim = c(-0.045,0.045), xlab = "Iteration Index", ylab = expression(Q[m]^(n/2)), col = alpha("black", 1.0), cex.lab = 1.5, cex.axis = 1.5)
grid()
lines(Qm, type = "h", lwd = 2, col = alpha("black", 0.4))
abline(h = 0, col = "deeppink", lwd = 4)

par(mar = c(5, 6, 3, 2))
hist(boot.res, freq = F, col = alpha("lightblue4", 0.4), xlab = expression(paste("Bootstrapped Mean of ", Q[m]^(n/2))), main = "", cex.lab = 1.5, cex.axis = 1.5)
abline(v = boot.Qm, lty = 3, lwd = 5, col = "deeppink")
legend("topright", legend = boot.Qm, col = "deeppink", lty = 3, lwd = 5)

