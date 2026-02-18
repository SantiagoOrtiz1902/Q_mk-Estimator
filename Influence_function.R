###########################################
# Empirical Robust Properties of Q_m^(n\2)
###########################################

# Author: Santiago Ortiz (saortizar@unal.edu.co)
# Date: 30/12/2025

# Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load codes and auxiliary programms
source("Code_Qmk.R")

library(doParallel)
library(MASS)

biCovQn2 <- function(x, y, type){
  frame.x = outer(x, x, "-")
  frame.y = outer(y, y, "-")
  frame.pairs = frame.x[lower.tri(frame.x)] * frame.y[lower.tri(frame.y)]
  if (missing(type)) {type = 0.5}
  if (type == 2) {return(0.355*(quantile(frame.pairs , probs = 0.75) + quantile(frame.pairs , probs = 0.25)))} # k1 = n/4 and k2 = 3n/4
  return(quantile(frame.pairs , probs = type))
}

calcular_influencia_empirica_biCovQn <- function(estimador, datos_x, datos_y, ...) {
  h1 <- seq(-5, 5, by = 0.5)  # Tamaño del paso en X
  h2 <- seq(-5, 5, by = 0.5)  # Tamaño del paso en Y
  n <- length(datos_y)  # Longitud de los datos
  
  infl_emp <- c()  # Inicializar el vector de influencia empírica
  
  for (i in 1:length(h1)) {
    for (j in 1:length(h2)) {
      datos_x_perturbados <- datos_x
      datos_y_perturbados <- datos_y
      datos_x_perturbados[n] = datos_x_perturbados[n] + h1[i]  # Perturbar los datos de salida en la n-ésima posición
      datos_y_perturbados[n] = datos_y_perturbados[n] + h2[j]  # Perturbar los datos de salida en la i-ésima posición
      
      influencia_perturbada = estimador(datos_x_perturbados, datos_y_perturbados, ...)  # Calcular la influencia con los datos perturbados
      infl_emp = c(infl_emp, influencia_perturbada)  # Almacenar la influencia empírica para la i-ésima observación
      print(c(h1[i],h2[j]))
    }
  }
  m.infl_emp = matrix(infl_emp, length(h1), length(h2), byrow = T)
  return(m.infl_emp)  # Devolver el vector de influencia empírica
}

###############################
# Empirical Influence Function
###############################

set.seed(123)
mis_datos_x <- c(rnorm(999),0)
mis_datos_y <- c(rnorm(999),0)

# Q_m^k for k = n/4

  # Calcular la influencia empírica para el estimador biCovQn
  influencia_Cov <- calcular_influencia_empirica_biCovQn(cov, mis_datos_x, mis_datos_y)
  influencia_biCovQn1 <- calcular_influencia_empirica_biCovQn(biCovQn2, mis_datos_x, mis_datos_y, type = 0.25)
  
  par(mfrow = c(1,3))
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn1,
        phi = 35, theta = 45, col = alpha("deeppink",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn1,
        phi = 0, theta = 70, col = alpha("deeppink",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn1,
        phi = 0, theta = 160, col = alpha("deeppink",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)

# Q_m^k for k = n/2

  # Calcular la influencia empírica para el estimador biCovQn
  influencia_biCovQn2 <- calcular_influencia_empirica_biCovQn(biCovQn2, mis_datos_x, mis_datos_y, type = 0.5)
  
  par(mfrow = c(1,3))
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn2,
        phi = 35, theta = 45, col = alpha("turquoise",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn2,
        phi = 0, theta = 70, col = alpha("turquoise",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn2,
        phi = 0, theta = 160, col = alpha("turquoise",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)
  
# Q_m^k for k = 3n/4

  # Calcular la influencia empírica para el estimador biCovQn
  influencia_biCovQn3 <- calcular_influencia_empirica_biCovQn(biCovQn2, mis_datos_x, mis_datos_y, type = 0.75)
  
  par(mfrow = c(1,3))
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn3,
        phi = 35, theta = 45, col = alpha("gold",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn3,
        phi = 0, theta = 70, col = alpha("gold",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn3,
        phi = 0, theta = 160, col = alpha("gold",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 2, cex.axis = 2)

# Q_m^(k1,k2) for  k1 = n/4 and k2 = 3n/4
  
  # Calcular la influencia empírica para el estimador biCovQn
  influencia_biCovQn4 <- calcular_influencia_empirica_biCovQn(biCovQn2, mis_datos_x, mis_datos_y, type = 2)
  
  par(mfrow = c(1,3))
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn4,
        phi = 35, theta = 45, col = alpha("blue",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 1.5, cex.axis = 2)
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn4,
        phi = 0, theta = 70, col = alpha("blue",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 1.5, cex.axis = 2)
  par(mar = c(0,3,0,0))
  persp(seq(-5, 5, by = 0.5), seq(-5, 5, by = 0.5), influencia_biCovQn4,
        phi = 0, theta = 160, col = alpha("blue",0.5), shade = 0.5, expand = 0.7,
        ticktype = "detailed", xlab = "X", ylab = "Y", zlab = "", cex.lab = 1.5, cex.axis = 2)
  

####################################
# Empirical Gross-error Sensitivity
####################################

rhos = seq(-1, 1, by = 0.1)
Ge.S = rep(0, length(rhos))
p = 2
for (i in 1:length(rhos)) {
  X = rbind(mvrnorm(999, rep(0, p), matrix(c(1, rhos[i], rhos[i], 1), p, p)), c(0,0))
  Ge.S[i] = max(abs(calcular_influencia_empirica_biCovQn(biCovQn2, X[,1], X[,2])))
  print(i)
}
Ge.S25 = c(2.570402181, 2.521543263, 2.287276074, 2.150430288, 1.756508314, 1.489117217, 1.548881201,
           1.158070532, 1.012445470, 0.855794948, 0.785515803, 0.623553789, 0.432973652, 0.354508519,
           0.264754960, 0.153447653, 0.059972539, 0.002348876, 0.039068594, 0.125755792, 0.198284887)
Ge.S50 = c(0.896612776, 0.813137677, 0.656265083, 0.595004391, 0.397930533, 0.340730343, 0.277169700,
           0.143979807, 0.101854003, 0.046966424, 0.002368253, 0.034751167, 0.119862180, 0.152260278,
           0.235159587, 0.311598283, 0.486795553, 0.523000636, 0.660706029, 0.770154650, 0.828043610)
Ge.S75 = c(0.216326930, 0.098843644, 0.050960270, 0.001867563, 0.060921095, 0.157680578, 0.250437673,
           0.353706959, 0.480863521, 0.551289771, 0.691045154, 0.846359827, 1.188443412, 1.356328391,
           1.397852278, 1.657736441, 1.922309840, 2.060623087, 2.169174159, 2.569501611, 2.481273547)
Ge.S2575 = c(0.9699367, 0.8951087, 0.8004339, 0.7951412, 0.6150631, 0.5267708, 0.4338911,
             0.3382764, 0.1704625, 0.1284554, 0.0190859, 0.1720768, 0.2032217, 0.3051420,
             0.4154203, 0.4837311, 0.6730289, 0.7274655, 0.7901317, 0.9028783, 1.0212076)

# Fig 1
par(mar = c(5,5,2,9), xpd = F)
plot(rhos, Ge.S25, type = "l", lwd = 2, xlab = "Correlation Value",
     ylab = "Gross-error Senstivity", cex.lab = 1.5, cex.axis = 1.5)
grid()
par(mar = c(5,5,2,8), xpd = T)
points(rhos, Ge.S25, pch = 20)
lines(rhos, Ge.S50, lwd = 2, col = "purple")
points(rhos, Ge.S50, pch = 20, col = "purple")
lines(rhos, Ge.S75, lwd = 2, col = "turquoise3")
points(rhos, Ge.S75, pch = 20, col = "turquoise3")
legend("topright", inset = c(-0.3, 0), box.lty = 1,
       legend = c("k = n/4","k = n/2","k = 3n/4"),
       lty = 1, pch = 20, lwd = 2,
       col = c("black", "purple", "turquoise"),
       horiz = F, text.font = 1, cex = 1.4)

# Fig 2
par(mar = c(5,5,2,9.2), xpd = F)
plot(rhos, Ge.S25, type = "l", lwd = 2, xlab = "Correlation Value",
     ylab = "Gross-error Senstivity", cex.lab = 1.5, cex.axis = 1.5)
grid()
par(mar = c(5,5,2,8), xpd = T)
points(rhos, Ge.S25, pch = 20)
lines(rhos, Ge.S50, lwd = 2, col = "purple")
points(rhos, Ge.S50, pch = 20, col = "purple")
lines(rhos, Ge.S75, lwd = 2, col = "turquoise3")
points(rhos, Ge.S75, pch = 20, col = "turquoise3")
lines(rhos, Ge.S2575, lwd = 2, col = "deeppink")
points(rhos, Ge.S2575, pch = 20, col = "deeppink")
legend("topright", inset = c(-0.32, 0), box.lty = 1,
       legend = c("k = n/4","k = n/2","k = 3n/4", "k1 = n/4\nk2 = 3n/4"),
       lty = 1, pch = 20, lwd = 2,
       col = c("black", "purple", "turquoise", "deeppink"),
       horiz = F, text.font = 1, cex = 1.4)