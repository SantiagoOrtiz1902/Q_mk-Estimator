# Load required packages
my_packages = c("robustbase", "MASS", "doParallel", "foreach", "rrcov", "parallel", "DetMCD", "SpatialNP", "fastM")
not_installed = my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if (length(not_installed)) install.packages(not_installed, dependencies = TRUE)
for (q in 1:length(my_packages)) {
  library(my_packages[q], character.only = TRUE)
}

# Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load codes and auxiliary programms
source("Code_Qmk.R")
source("AtipTestv7_2.R")
source("stahel_donoho_modified.R")

# KASP parameters
kn.par = rbind(c(1,1,0,1,0,1), c(1,1,-1,-1/4,0,2), c(1,1,-1,0,0,2))

#################################
# EXPERIMENT 2: Breakdown Point #
#################################

set.seed(123)

m = 100
alpha = c(0.1,0.2, 0.3, 0.4)
delta = c(4, 6)
p = c(10, 20, 40)
lambda = 0.5

resultado = matrix(, 0, 9)

for (i1 in 1:length(p)) {
  n = 10*p[i1]
  for (i2 in 1:length(alpha)) {
    for (i3 in 1:length(delta)) {
      MFE = matrix(NA, m, 9)
      for (i4 in 1:m) {
        # Clean Data
          # Generate data
          Xu = mvrnorm(n, rep(0, p[i1]), diag(p[i1]))
          # Computation of methods
          Qmk1.k2 = qnCov(Xu, type = 1)
          Qmk.50 = qnCov(Xu, type = 3)
          detmcd = DetMCD(Xu)$cov
          Sn.cov = snCov(Xu)
          mm.cov = CovMMest(Xu)$cov
          ssrc = SCov(Xu)
          tyler.DNS = TYLERshape(Xu)$Sigma
          KASP = KurMain(Xu, kn.par, mode.sim = 1)$cov.rob
          SDE.SSD = stahel_donoho_estimator(Xu, return_mah = F)$S_SD
        # Contamination
          outs = sample(1:n, round((1-alpha[i2])*n), replace = F)
          Xo = Xu
          Xo[outs,] = mvrnorm(length(outs), delta[i3]*rep(1,p[i1]), lambda*diag(p[i1]))
          # Computation of methods under contamination
          Qmk1.k2.c = qnCov(Xo, type = 1)
          Qmk.50.c = qnCov(Xo, type = 3)
          detmcd.c = DetMCD(Xo)$cov
          Sn.cov.c = snCov(Xo)
          mm.cov.c = CovMMest(Xo)$cov
          ssrc.c = SCov(Xo)
          tyler.DNS.c = TYLERshape(Xo)$Sigma
          KASP.c = KurMain(Xo, kn.par, mode.sim = 1)$cov.rob
          SDE.SSD.c = stahel_donoho_estimator(Xo, return_mah = F)$S_SD
        # Metrics
        MFE[i4,] = c(norm(Qmk1.k2 - Qmk1.k2.c, type = "F"), norm(Qmk.50 - Qmk.50.c, type = "F"), norm(detmcd - detmcd.c, type = "F"),
                     norm(Sn.cov - Sn.cov.c, type = "F"), norm(mm.cov - mm.cov.c, type = "F"), norm(ssrc - ssrc.c, type = "F"),
                     norm(tyler.DNS - tyler.DNS.c, type = "F"), norm(KASP - KASP.c, type = "F"), norm(SDE.SSD - SDE.SSD.c, type = "F"))
        cat("p =", p[i1], ";", "alpha =", alpha[i2], ";", "delta =", delta[i3], ";", "m =", i4, "\n")
      }
      resultado = rbind(resultado, apply(MFE, 2, mean))
      #cat("p=", p[i1], ";", "alpha", alpha[i2], ";", "delta", delta[i3], "\n")
    }
  }
}

save(resultado, file = "MFE_cov.RData")

###################################
# EXPERIMENT 1: Sample Efficiency #
###################################

set.seed(123)

m = 500
alpha = c(0.1, 0.15, 0.2,0.25)
delta = c(0.5, 0.8, 1.2)
n = c(10, 50, 100)
lambda = 0.1
rho = 0.8
p = 2

resultado = matrix(, 0, 9)

for (i1 in 1:length(n)) {
  for (i2 in 1:length(alpha)) {
    nu = round((1-alpha[i2])*n[i1])
    no = n[i1] - nu
    for (i3 in 1:length(delta)) {
      MAE = matrix(NA, m, 9)
      for (i4 in 1:m) {
        # Generate data
        Xu = mvrnorm(nu, rep(0,p), matrix(c(1, rho, rho, 1), p))
        Xo = mvrnorm(no, delta[i3]*c(min(Xu[,1]),max(Xu[,2])), lambda*diag(p))
        X = rbind(Xu,Xo)
        # Computation of methods
        Qmk1.k2 = biCovQn(X[,1], X[,2], type = 1)
        Qmk.50 = biCovQn(X[,1], X[,2], type = 3)
        detmcd = DetMCD(X)$cov[1,2]
        Sn.cov = biCovSn(X[,1], X[,2])
        mm.cov = CovMMest(X)$cov[1,2]
        ssrc = SCov(X)[1,2]
        tyler.DNS = TYLERshape(X)$Sigma[1,2]
        KASP = KurMain(X, kn.par, mode.sim = 1)$cov.rob[1,2]
        SDE.SSD = stahel_donoho_estimator(X, return_mah = F)$S_SD[1,2]
        # Metrics
        MAE[i4,] = abs( c(Qmk1.k2, Qmk.50, detmcd, Sn.cov, mm.cov, ssrc, tyler.DNS, KASP, SDE.SSD) - rho )
      }
      resultado = rbind(resultado, apply(MAE, 2, mean))
      cat("n=", n[i1], ";", "alpha", alpha[i2], ";", "delta", delta[i3], "\n")
    }
  }
}
