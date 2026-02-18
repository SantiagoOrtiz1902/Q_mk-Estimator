###############################################################
# Q_m^k Estimator and other methods for simulation experiments
###############################################################

# Author: Santiago Ortiz (saortizar@unal.edu.co)
# Date: 30/12/2025

######################################
# Generate random covariance matrixes
######################################

generate_cor <- function(cond, p){
    if (p < 2){
        print("p must be larger that 2")
        stop("Error in the parameters")
    }else{
        maxits = 300
        tol = 1/(10^(5))
        ramdon = floor(runif(p-2,0,cond))
        lambda = sort(c(1,ramdon, cond))
        x = mvrnorm(n=p, mu= rep(0,p), Sigma = diag(1,p))
        Sigma = x %*% t(x);
        Q = eigen(Sigma)$vectors
        ratio = 0
        iter = 0
        
        while ((abs(ratio - cond) > tol) & (iter < maxits)){
            iter = iter + 1;
            Sigma = Q %*% diag(lambda) %*% t(Q);
            factor = diag(diag(Sigma)^(-1/2))
            Sigma = factor %*% Sigma %*% factor;
            Q = eigen(Sigma)$vectors;
            lambda = eigen(Sigma)$values
            ratio = lambda[1]/lambda[p];
            lambda[p] = lambda[1]/cond;
        }
        return(Sigma)
    }
}

###############
# Pn Estimator
############### 

Pn = function(y){
  n=length(y)
  if(n <= 2){
    warning("Your sample size is too small")
    return()
  }
  y.pairs = outer(y,y,"+")
  y.pairs = y.pairs[lower.tri(y.pairs)]/2
  const = 1/0.9539 # asymptotic correction factor
  scale.est=const*as.numeric(diff(quantile(y.pairs,c(1/4,3/4),type=1)))
  # Correction factors obtained through simulation over 1 million replications
  correction.factors =c(1.128,1.303,1.109,1.064,1.166,1.103,1.087,1.105,1.047,1.063,1.057,1.040,1.061,1.047,1.043,1.048,1.031,1.037,1.035,1.028,1.036,1.030,1.029,1.032,1.023,1.025,1.024,1.021,1.026,1.022,1.021,1.023,1.018,1.020,1.019,1.017,1.020,1.018,1.017,1.018,1.015,1.016,1.016,1.014,1.016,1.015,1.014,1.015)
  if(n <= 40){
    scale.est = scale.est*correction.factors[n-2]
    # n-2 as the first element of the correction.factors vector is for n=3
    }else if(n > 40) scale.est = scale.est*n/(n-0.7)
  return(scale.est)
}

biCovPn <- function(X, Y){
 alpha <- 1/Pn(X)
 beta <- 1/Pn(Y)
 Cov_PN = (Pn(alpha*X + beta*Y)^2 - Pn(alpha*X - beta*Y)^2)*(1/(beta*alpha))
 return(0.25*Cov_PN)
}

PnCov <- function(data){
  p <- ncol(data)
  covqn = matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      if (i == j){
        X <- data[,i]
        covqn[i, j] = Pn(X)^2 # <-- usamos el factor que ya tiene 
      }else{
        X <- data[,i]
        Y <- data[,j]
        covqn[i, j] <- biCovPn(X, Y) 
      }
    }
  }
  return(covqn)
}

###############
# Sn-Cov Estimator
############### 

biCovSn = function(X, Y) {
  p = length(X)
  # Create matrices of differences
  diff_X = matrix(rep(X, each = p), ncol = p, byrow = TRUE)
  diff_Y = matrix(rep(Y, each = p), ncol = p, byrow = TRUE)
  # Exclude diagonal elements
  diff_X = diff_X[row(diff_X) != col(diff_X)]
  diff_Y = diff_Y[row(diff_Y) != col(diff_Y)]
  # Compute the product and reshape
  product_mat = matrix(diff_X * diff_Y, ncol = p-1)
  # Compute medians
  medians = apply(product_mat, 2, median)
  # Compute the median of medians
  result = median(medians)
  return(result)
}

snCov = function(data) {
  p = ncol(data)
  covsn = matrix(0, p, p)
  for (i in 1:(p - 1)) {
    X = data[, i]
    covsn[i, i] = Sn(X)^2  # Use Sn scale estimator for diagonal elements
    for (j in (i + 1):p) {
      Y <- data[, j]
      covsn[i, j] = covsn[j, i] = biCovSn(X, Y)
    }
  }
  return((1.1926^2) * covsn)
}

##################
# Q_m^k Estimator
##################

biCovQn <- function(x, y, type){
  frame.x = outer(x, x, "-")
  frame.y = outer(y, y, "-")
  frame.pairs = frame.x[lower.tri(frame.x)] * frame.y[lower.tri(frame.y)]
  if (type == 0){
    return(as.numeric(1*(quantile(frame.pairs , probs = 0.75) + quantile(frame.pairs , probs = 0.25))))
  } else if (type == 1){
    return(as.numeric(0.355*(quantile(frame.pairs , probs = 0.75) + quantile(frame.pairs , probs = 0.25))))
  } else if (type == 2){
    return(as.numeric(1*(quantile(frame.pairs , probs = 0.25))))
  } else if (type == 3){
    return(as.numeric(1*(quantile(frame.pairs , probs = 0.5))))
  } else if (type == 4){
    return(as.numeric(1*(quantile(frame.pairs , probs = 0.75))))
  }
}

qnCov <- function(data, type) {
  p <- ncol(data)
  D = diag(apply(data, 2, Qn))
  covqn <- matrix(0, p, p) + D
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      covqn[i, j] = covqn[j, i] = biCovQn(data[, i], data[, j], type)
    }
  }
  return(covqn)
}

###########################################################
# Positive-definite Procedure for Scatter Matrix Estimator
###########################################################

mz.Cov = function(datos, type){
  if (type == 1){
    #1
    D = diag(apply(datos, 2, Qn))
    Y = solve(D) %*% as.matrix(t(datos))
    #2
    U = qnCov(t(Y), 3)
    diag(U) = 1
    #3
    EIG = eigen(U)
    E = EIG$vectors
    A = diag(EIG$values)
    #4
    A2 = D %*% E
    Z = t(E) %*% Y
    V.X = A2 %*% diag(apply(t(Z), 2, Qn)^2) %*% t(A2)
    t.X = as.vector(A2 %*% apply(t(Z), 2, median))
    return(list(S = V.X, m = t.X))
  } else {
    #1
    D = diag(apply(datos, 2, Sn))
    Y = solve(D) %*% as.matrix(t(datos))
    #2
    U = snCov(t(Y))
    diag(U) = 1
    #3
    EIG = eigen(U)
    E = EIG$vectors
    A = diag(EIG$values)
    #4
    A2 = D %*% E
    Z = t(E) %*% Y
    V.X = A2 %*% diag(apply(t(Z), 2, Sn)^2) %*% t(A2)
    t.X = as.vector(A2 %*% apply(t(Z), 2, median))
    return(list(S = V.X, m = t.X))
    }
}

#########################################
# Kullback-Liebler for Matrix Comaprison
#########################################

KLD = function(Sigma, S){
  # This function computes de Kullback-Liebler between Sigma and S
  p = nrow(S)
  M = S %*% solve(Sigma)
  result = sum(diag(M)) - log(det(M)) - p
  return(result)
}

KLD.2 = function(Sigma, S){
  # This function computes de Kullback-Liebler between Sigma and S
  result = norm(Sigma - S, type = "F")
  return(result)
}

#########################################
# Panel histogram for scatterplot matrix
#########################################

panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, border = "black")#, col = "navyblue", border = "white", ...)
}

#################################################
# Cut-off of Sn-Covariance for Outlier Detection
#################################################

cut.Sn = function(d, p){
  ind = rep(0, length(d))
  if (p < 15){
    if (p <= 5){b = 1} else {b = 2.5}
    cv = b*qchisq(0.95, p)
  } else {
    cv = (qchisq(0.95, p) / qchisq(0.5, p)) * median(d)
  }
  ind[d > cv] = 1
  return(list(cutoff = cv, Flag = ind))
}

###############
# Outlier Map
###############

SD.OD = function(datos, S, m){
  n = nrow(datos)
  pca.qm = eigen(S)
  var.expl = cumsum(pca.qm$values / sum(pca.qm$values))
  k = min(which(var.expl >= 0.95))
  P.p = pca.qm$vectors[,1:k]
  T.k = as.matrix(dat.bus - rep(1, n) %*% t(m)) %*% P.p
  S.Qm.2 = P.p %*% diag(pca.qm$values[1:k]) %*% t(P.p)
  OD.qm = rep(0, n)
  for (i in 1:n) {
    OD.qm[i] = sqrt(sum(as.vector(dat.bus[i,] - m - P.p %*% T.k[i,])^2))
  }
  SD.qm = rep(0, n)
  for (i in 1:n) {
    SD.qm[i] = sqrt(sum(T.k[i,]^2/pca.qm$values[1:k]))
  }
  return(list(SD = SD.qm, OD = OD.qm, loading = k))
}
