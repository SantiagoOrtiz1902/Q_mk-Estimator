#######################
### Q_n^k Applications
#######################

# Author: Santiago Ortiz (saortizar@unal.edu.co)
# Date: 30/12/2025

# Load Packages
my_packages = c("robustbase", "MASS", "doParallel", "foreach", "rrcov", "parallel", "DetMCD", "R.matlab", "depthTools", "caret", "e1071")
not_installed = my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if (length(not_installed)) install.packages(not_installed, dependencies = TRUE)
for (q in 1:length(my_packages)) {library(my_packages[q], character.only = TRUE)}

# Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load outlier detection codes and auxiliary routines
source("Code_Qmk.R")
source("AtipTestv7_2.R")
source("stahel_donoho_modified.R")

kn.par = rbind(c(1,1,0,1,0,1), c(1,1,-1,-1/4,0,2), c(1,1,-1,0,0,2))

################################################
# Application in Outlier Detection (Glass Data)
################################################

datos = read.csv("glass.csv", header = F)
Ind = c(which(datos$V8 == 1), which(datos$V8 == 7))
X = datos[Ind,1:7]
Y = datos$V8[Ind]
#Y[Y == 1] = 0; Y[Y == 2] = 1#; Y[Y == 3] = 1

S = qnCov(X, type = 3)
m = apply(X, 2, median)
mbd.x = MBD(X, plotting = F)$ordering == nrow(X)
m2 = unlist(X[mbd.x,])
qn.mah = as.numeric(mahalanobis(X, m, S))
qn.mah.mbd = as.numeric(mahalanobis(X, m2, S))
mah = as.numeric(mahalanobis(X, apply(X, 2, mean), cov(X)))
plot(mah, type = "h", lwd = 2, col = Y)
abline(h = qchisq(0.975, ncol(X)), col = "green3", lty = 2, lwd = 3)
plot(qn.mah, type = "h", lwd = 2, col = Y)
abline(h = quantile(qn.mah, 0.9), col = "deeppink", lty = 2, lwd = 3)
abline(h = qchisq(0.975, ncol(X)), col = "green3", lty = 2, lwd = 3)
plot(qn.mah.mbd, type = "h", lwd = 2, col = Y)
abline(h = quantile(qn.mah.mbd, 0.9), col = "deeppink", lty = 2)
Ind.out = which(qn.mah.mbd > quantile(qn.mah.mbd, 0.9))
abline(h = qchisq(0.975, ncol(X)), col = "green3", lty = 2, lwd = 3)

## Mz Qm-Cov
est.qn.mz = mz.Cov(X, 1)
mah.qn.mz = mahalanobis(X, est.qn.mz$m, est.qn.mz$S)
plot(mah.qn.mz, type = "h", col = Y)
outl.qn.mz = cut.Sn(mah.qn.mz, ncol(X))

colores = rep("black", length(Y))
colores[Y != 1] = "deepskyblue1"
par(mar = c(5,5,2,8))
plot(mah.qn.mz, mah, pch = 19, col = colores,
     main = "DD-Plot of Glass Data",
     xlab = "Robust Mahalanobis Distances",
     ylab = "Mahalanobis Distances")
grid()
abline(h = qchisq(0.975, ncol(X)), col = "deeppink", lty = 3, lwd = 3)
abline(v = outl.qn.mz$cutoff, col = "deeppink", lty = 2, lwd = 3)
par(xpd = TRUE)
legend("topright", inset = c(-0.25, 0), box.lty = 1, lwd = 2,
       legend = c("k = n/4","k = n/2","k = 3n/4", "k1 = n/4\nk2 = 3n/4"),
       pch = c(20, 20),
       lty = c(0, 0, 3, 4),
       col = c("black", "deepskyblue1", "deeppink", "deeppink"),
       horiz = F, text.font = 1)
par(xpd = FALSE)

## DETMCD
outl.detmcd = DetMCD(X)$flag
outl.detmcd[outl.detmcd == 1] = 2
outl.detmcd[outl.detmcd == 0] = 1
outl.detmcd[outl.detmcd == 2] = 0

## MM-Cov
outl.mm = as.numeric(getFlag(rrcov::CovMMest(X)))
outl.mm[outl.mm == 1] = 2
outl.mm[outl.mm == 0] = 1
outl.mm[outl.mm == 2] = 0

## Sn-Cov
est.sn = mz.Cov(X, 2)
mah.sn = mahalanobis(X, est.sn$m, est.sn$S)
outl.sn = cut.Sn(mah.sn, ncol(X))

## KASP
outl.kasp = KurMain(scale(X), kn.par, mode.sim = 1)$lbl

## SDE-SSD
outl.SDESSD = as.numeric(stahel_donoho_estimator(scale(X), return_mah = T)$flag)
outl.SDESSD[outl.SDESSD == 1] = 2
outl.SDESSD[outl.SDESSD == 0] = 1
outl.SDESSD[outl.SDESSD == 2] = 0

Y2 = Y
Y2[Y2 == 1] = 2
Y2[Y2 == 7] = 1
Y2[Y2 == 2] = 0
confusionMatrix(as.factor(outl.qn.mz$Flag), as.factor(Y2))
confusionMatrix(as.factor(outl.detmcd), as.factor(Y2))
confusionMatrix(as.factor(outl.mm), as.factor(Y2))
confusionMatrix(as.factor(outl.sn$Flag), as.factor(Y2))
confusionMatrix(as.factor(outl.kasp), as.factor(Y2))
confusionMatrix(as.factor(outl.SDESSD), as.factor(Y2))

### FIGURES 

# REAL
colores = rep("black", nrow(X)); colores[Y != 1] = "deeppink"
mark = rep(20, nrow(X)); mark[Y != 1] = 19
pairs(X[,1:4], pch = mark, col = colores, horOdd = TRUE, diag.panel = panel.hist, oma = c(3,3,7,16),
      #main = "Scatterplot Matrix for First Four Variables: Glass-Identification Data",
      labels = c("Variable 1", "Variable 2", "Variable 3", "Variable 4"), cex.lab = 1.5, cex.axis = 2.0)
par(xpd = TRUE)
legend("topright", legend = c("Real\nOutliers","Data"), col = c("deeppink","black"), pch = 19, box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# Qm(n/2)
colores = rep("black", nrow(X)); colores[outl.qn.mz$Flag == 1] = "deepskyblue1"
mark = rep(20, nrow(X)); mark[outl.qn.mz$Flag == 1] = 19
pairs(X[,1:4], pch = mark, col = colores, horOdd = TRUE, diag.panel = panel.hist, oma = c(3,3,7,16),
      #main = "Scatterplot Matrix for First Four Variables: Glass-Identification Data",
      labels = c("Variable 1", "Variable 2", "Variable 3", "Variable 4"), cex.lab = 1.5, cex.axis = 2.0)
par(xpd = TRUE)
legend("topright", legend = c(expression(atop(Q[m]^(n/2), "Outliers")),"Data"), col = c("deepskyblue1","black"), pch = 19, box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# DetMCD
colores = rep("black", nrow(X)); colores[outl.detmcd == 1] = "green4"
mark = rep(20, nrow(X)); mark[outl.detmcd == 1] = 19
pairs(X[,1:4], pch = mark, col = colores, horOdd = TRUE, diag.panel = panel.hist, oma = c(3,3,7,16),
      #main = "Scatterplot Matrix for First Four Variables: Glass-Identification Data",
      labels = c("Variable 1", "Variable 2", "Variable 3", "Variable 4"), cex.lab = 1.5, cex.axis = 2.0)
par(xpd = TRUE)
legend("topright", legend = c("DetMCD\nOutliers","Data"), col = c("green4","black"), pch = 19, box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# Sn-Cov
colores = rep("black", nrow(X)); colores[outl.sn$Flag == 1] = "firebrick1"
mark = rep(20, nrow(X)); mark[outl.sn$Flag == 1] = 19
pairs(X[,1:4], pch = mark, col = colores, horOdd = TRUE, diag.panel = panel.hist, oma = c(3,3,7,16),
      #main = "Scatterplot Matrix for First Four Variables: Glass-Identification Data",
      labels = c("Variable 1", "Variable 2", "Variable 3", "Variable 4"), cex.lab = 1.5, cex.axis = 2.0)
par(xpd = TRUE)
legend("topright", legend = c("Sn-Cov\nOutliers","Data"), col = c("firebrick1","black"), pch = 19, box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# MM-Cov
colores = rep("black", nrow(X)); colores[outl.mm == 1] = "darkorange1"
mark = rep(20, nrow(X)); mark[outl.mm == 1] = 19
pairs(X[,1:4], pch = mark, col = colores, horOdd = TRUE, diag.panel = panel.hist, oma = c(3,3,7,16),
      #main = "Scatterplot Matrix for First Four Variables: Glass-Identification Data",
      labels = c("Variable 1", "Variable 2", "Variable 3", "Variable 4"), cex.lab = 1.5, cex.axis = 2.0)
par(xpd = TRUE)
legend("topright", legend = c("MM-Cov\nOutliers","Data"), col = c("darkorange1","black"), pch = 19, box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# KASP

colores = rep("black", nrow(X)); colores[outl.kasp == 1] = "purple"
mark = rep(20, nrow(X)); mark[outl.kasp == 1] = 19
pairs(X[,1:4], pch = mark, col = colores, horOdd = TRUE, diag.panel = panel.hist, oma = c(3,3,7,16),
      #main = "Scatterplot Matrix for First Four Variables: Glass-Identification Data",
      labels = c("Variable 1", "Variable 2", "Variable 3", "Variable 4"), cex.lab = 1.5, cex.axis = 2.0)
par(xpd = TRUE)
legend("topright", legend = c("KASP\nOutliers","Data"), col = c("purple","black"), pch = 19, box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# SDE-SSD

colores = rep("black", nrow(X)); colores[outl.SDESSD == 1] = "gold3"
mark = rep(20, nrow(X)); mark[outl.SDESSD == 1] = 19
pairs(X[,1:4], pch = mark, col = colores, horOdd = TRUE, diag.panel = panel.hist, oma = c(3,3,7,16),
      #main = "Scatterplot Matrix for First Four Variables: Glass-Identification Data",
      labels = c("Variable 1", "Variable 2", "Variable 3", "Variable 4"), cex.lab = 1.5, cex.axis = 2.0)
par(xpd = TRUE)
legend("topright", legend = c("SDE-SSD\nOutliers","Data"), col = c("gold3","black"), pch = 19, box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

################################
# Application in PCA (Bus Data)
################################

data("bus")
dat.bus = bus[,-c(6,9)]
n = nrow(dat.bus)
p = ncol(dat.bus)
bus.Qn = apply(dat.bus, 2, sd)
dat.bus.st = scale(dat.bus)
# dat.bus.st = sweep(dat.bus, 2, bus.Qn, "/", check.margin = FALSE)
apply(dat.bus.st, 2, skewness)

# MM Cov
MM.Cov = CovMMest(x = dat.bus)
## DetMCD
DetMCD.Cov = DetMCD(dat.bus)

pairs(dat.bus.st[,1:4], pch = mark, horOdd = TRUE, diag.panel = panel.hist,
      #main = "Scatterplot Matrix for First Four Variables: Bus Data",
      labels = c("Variable 1", "Variable 2", "Variable 3", "Variable 4"), cex.lab = 1.5, cex.axis = 2.0)
pairs(dat.bus.st[,5:8], pch = mark, horOdd = TRUE, diag.panel = panel.hist,
      #main = "Scatterplot Matrix for Second Four Variables: Bus Data",
      labels = c("Variable 5", "Variable 6", "Variable 7", "Variable 8"), cex.lab = 1.5, cex.axis = 2.0)

S.qm = mz.Cov(datos = dat.bus, type = 1)
S.sn = mz.Cov(datos = dat.bus, type = 2)
qm.Dis = SD.OD(dat.bus.st, S.qm$S, S.qm$m)
sn.Dis = SD.OD(dat.bus.st, S.sn$S, S.sn$m)
MM.Dis = SD.OD(dat.bus.st, MM.Cov$cov, MM.Cov$center)
DetMCD.Dis = SD.OD(dat.bus.st, DetMCD.Cov$cov, DetMCD.Cov$center)

#plot(sort(log(OD.S)), sort(log(OD.qm)), xlab="classical", ylab="robust")
#lines(sort(log(OD.S)), sort(log(OD.qm)))
#abline(reg = c(0,1), col = "red")

pca3 <- PcaClassic(dat.bus.st) # classical
hpca3 <- PcaHubert(dat.bus.st) # Hubert
pca.OD <- pca3@od
pca.SD <- pca3@sd
robpca.OD <- hpca3@od
robpca.SD <- hpca3@sd

### Outlier maps

# Standard PCA
cut.sd = sqrt(qchisq(0.975, 16))
cut.OD = quantile(pca.OD, 0.975)
outl1 = unique(c(which(pca.SD >= cut.sd), which(pca.OD >= cut.OD)))
col1 = rep("black", n); col1[outl1]= "deeppink"
mark1 = rep(20, n); mark1[outl1] = 19
par(mar = c(5,5,2,14))
plot(pca.SD, pca.OD, pch = mark1, col = col1,
     #main = "Outlier Map Standard PCA: Bus Data",
     ylab = "Orthogonal Distance", xlab = "Score Distance", cex.lab = 1.5, cex.axis = 1.5)
grid()
abline(v = cut.sd, col = "gray50", lwd = 4, lty = 2)
abline(h = cut.OD, col = "gray50", lwd = 4, lty = 3)
par(xpd = TRUE)
legend("topright", inset = c(-0.5, 0),
       legend = c("Standard PCA\nOutliers","Data", "Cut-off\nScore Dist.", "Cut-off\nOrthogonal Dist."),
       col = c("deeppink","black","gray50","gray50"), pch = c(19,19, NA, NA), lty = c(0, 0, 2, 3), lwd = c(NA, NA, 3, 3),
       box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# ROBPCA
cut.sd = sqrt(qchisq(0.975, 16))
cut.OD = quantile(robpca.OD, 0.975)
outl1 = unique(c(which(robpca.SD >= cut.sd), which(robpca.OD >= cut.OD)))
col2 = rep("black", n); col2[outl1]= "purple"
mark2 = rep(20, n); mark2[outl1] = 19
par(mar = c(5,5,2,14))
plot(robpca.SD, robpca.OD, pch = mark2, col = col2,
     #main = "Outlier Map ROBPCA: Bus Data",
     ylab = "Orthogonal Distance", xlab = "Score Distance", cex.lab = 1.5, cex.axis = 1.5)
grid()
abline(v = cut.sd, col = "gray50", lwd = 4, lty = 2)
abline(h = cut.OD, col = "gray50", lwd = 4, lty = 3)
par(xpd = TRUE)
legend("topright", inset = c(-0.5, 0),
       legend = c("ROBPCA\nOutliers","Data", "Cut-off\nScore Dist.", "Cut-off\nOrthogonal Dist."),
       col = c("purple","black","gray50","gray50"), pch = c(19,19, NA, NA), lty = c(0, 0, 2, 3), lwd = c(NA, NA, 3, 3),
       box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# Qm^(n/2)
cut.sd = sqrt(qchisq(0.975, 16))
cut.OD = quantile(qm.Dis$OD, 0.975)
outl1 = unique(c(which(qm.Dis$SD >= cut.sd), which(qm.Dis$OD >= cut.OD)))
col3 = rep("black", n); col3[outl1]= "deepskyblue1"
mark3 = rep(20, n); mark3[outl1] = 19
par(mar = c(5,5,2,14))
plot(qm.Dis$SD, qm.Dis$OD, pch = mark3, col = col3,
     #main = "Outlier Map Robust Qm(n/2) PCA: Bus Data",
     ylab = "Orthogonal Distance", xlab = "Score Distance", cex.lab = 1.5, cex.axis = 1.5)
grid()
abline(v = cut.sd, col = "gray50", lwd = 4, lty = 2)
abline(h = cut.OD, col = "gray50", lwd = 4, lty = 3)
par(xpd = TRUE)
legend("topright", inset = c(-0.5, 0),
       legend = c(expression(atop(Q[m]^(n/2), "Outliers")),"Data", "Cut-off\nScore Dist.", "Cut-off\nOrthogonal Dist."),
       col = c("deepskyblue1","black","gray50","gray50"), pch = c(19,19, NA, NA), lty = c(0, 0, 2, 3), lwd = c(NA, NA, 3, 3),
       box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# Sn-Cov
cut.sd = sqrt(qchisq(0.975, 16))
cut.OD = quantile(sn.Dis$OD, 0.975)
outl1 = unique(c(which(sn.Dis$SD >= cut.sd), which(sn.Dis$OD >= cut.OD)))
col4 = rep("black", n); col4[outl1]= "firebrick1"
mark4 = rep(20, n); mark4[outl1] = 19
par(mar = c(5,5,2,14))
plot(sn.Dis$SD, sn.Dis$OD, pch = mark4, col = col4,
     #main = "Outlier Map Robust Sn-Cov PCA: Bus Data",
     ylab = "Orthogonal Distance", xlab = "Score Distance", cex.lab = 1.5, cex.axis = 1.5)
grid()
abline(v = cut.sd, col = "gray50", lwd = 4, lty = 2)
abline(h = cut.OD, col = "gray50", lwd = 4, lty = 3)
par(xpd = TRUE)
legend("topright", inset = c(-0.5, 0),
       legend = c("Sn-Cov\nOutliers","Data", "Cut-off\nScore Dist.", "Cut-off\nOrthogonal Dist."),
       col = c("firebrick1","black","gray50","gray50"), pch = c(19,19, NA, NA), lty = c(0, 0, 2, 3), lwd = c(NA, NA, 3, 3),
       box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# MM Estimators
cut.sd = sqrt(qchisq(0.975, 16))
cut.OD = quantile(MM.Dis$OD, 0.975)
outl1 = unique(c(which(MM.Dis$SD >= cut.sd), which(MM.Dis$OD >= cut.OD)))
col5 = rep("black", n); col5[outl1]= "darkorange1"
mark5 = rep(20, n); mark5[outl1] = 19
par(mar = c(5,5,2,14))
plot(MM.Dis$SD, MM.Dis$OD, pch = mark5, col = col5,
     #main = "Outlier Map Robust MM-Cov PCA: Bus Data",
     ylab = "Orthogonal Distance", xlab = "Score Distance", cex.lab = 1.5, cex.axis = 1.5)
grid()
abline(v = cut.sd, col = "gray50", lwd = 4, lty = 2)
abline(h = cut.OD, col = "gray50", lwd = 4, lty = 3)
par(xpd = TRUE)
legend("topright", inset = c(-0.5, 0),
       legend = c("MM-Cov\nOutliers","Data", "Cut-off\nScore Dist.", "Cut-off\nOrthogonal Dist."),
       col = c("darkorange1","black","gray50","gray50"), pch = c(19,19, NA, NA), lty = c(0, 0, 2, 3), lwd = c(NA, NA, 3, 3),
       box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)

# DETMCD
cut.sd = sqrt(qchisq(0.975, 16))
cut.OD = quantile(DetMCD.Dis$OD, 0.975)
outl1 = unique(c(which(DetMCD.Dis$SD >= cut.sd), which(DetMCD.Dis$OD >= cut.OD)))
col6 = rep("black", n); col6[outl1]= "green4"
mark6 = rep(20, n); mark6[outl1] = 19
par(mar = c(5,5,2,14))
plot(DetMCD.Dis$SD, DetMCD.Dis$OD, pch = mark6, col = col6,
     #main = "Outlier Map Robust DetMCD PCA: Bus Data",
     ylab = "Orthogonal Distance", xlab = "Score Distance", cex.lab = 1.5, cex.axis = 1.5)
grid()
abline(v = cut.sd, col = "gray50", lwd = 4, lty = 2)
abline(h = cut.OD, col = "gray50", lwd = 4, lty = 3)
par(xpd = TRUE)
legend("topright", inset = c(-0.5, 0),
       legend = c("DetMCD\nOutliers","Data", "Cut-off\nScore Dist.", "Cut-off\nOrthogonal Dist."),
       col = c("green4","black","gray50","gray50"), pch = c(19,19, NA, NA), lty = c(0, 0, 2, 3), lwd = c(NA, NA, 3, 3),
       box.lty = 1, horiz = FALSE, text.font = 1, cex = 1.4)
par(xpd = FALSE)