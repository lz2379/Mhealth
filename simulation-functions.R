##########Simulation functions######################

library(MASS)

#1-dim Kernel function, h is the bandwidth
K = function(x, h){
  return(1 / (h * (sqrt(2 * pi))) * exp( - (x / h)^2 / 2))
}


generate_sample = function(seed = 1,  
                           n, t, d = 1, 
                           beta10 = 0, beta11 = 2, 
                           eta1 =  - 0.2, eta2 = 0.2,  
                           theta1 = 0.8, theta2 = 0,  
                           gamma1 = 1,  gamma2 =  - 0.5,  
                           cc = sqrt(0.5), sigma = 0.5, 
                           print="FALSE"){
  set.seed(seed)
  
  #generate dataset
  eps_cov = (sigma^2) * exp(as.matrix(log(cc) * dist(1:(t + 1), upper = TRUE, diag = TRUE)))
  epsilon = as.matrix(mvrnorm(n = n, mu = rep(0, t + 1), Sigma = eps_cov))
  A = matrix(NA, nrow = n, ncol = t)
  X = array(NA, dim = c(n, t, d + 1)) # sample size, time index, covariates
  Y = matrix(NA, nrow = n, ncol = t + 1)
  #t=1
  A[, 1] = runif(n = n, min = 0, max = 1)
  X[, , 1] = 1 
  X[, 1, 2] = rnorm(n = n, mean = 0, sd = 0.5)
  for (i in 2:t){
    X[, i, 2] = rnorm(n = n, mean = eta1 * X[, i - 1, 2] + eta2 * A[, i - 1], sd = 0.5)
    A[, i] = rnorm(n = n, mean = gamma1 * X[, i, 2] + gamma2 * A[, i - 1], sd = 0.5)
    
    if (i>2) {
      Y[, i] = theta1 * X[, i - 1, 2] + theta2 * (A[, i - 2]) - A[, i - 1] * (A[, i - 1] - beta10 - beta11 * X[, i - 1, 2]) + epsilon[, i]
    }
  }
  Y[, t + 1] = theta1 * (X[, t, 2]) + theta2 * (A[, t - 1]) - A[, t] * (A[, t] - beta10 - beta11 * X[, t, 2]) + epsilon[, t + 1]
  return(list(A = A, X = X, Y = Y, n = n, t = t, d = d))
}


weightedK = function(S, A, Y, 
                     lag = 3, 
                     n, t, d = 1, 
                     starttime = 2, # the first time point when Y is available
                     bootstrap = 0, 
                     w = NA, #weights
                     print = "FALSE"){ #k is lag
  
  #if weights are not given, we assume equal weights for all lags
  if (sum(!is.na(w)) == 0){
    w = rep(1 / lag, lag)
  }
  
  Ahat = matrix(NA, nrow = n, ncol = t - lag + 1)
  Yhat = matrix(NA, nrow = n, ncol = t - lag + 1)
  AAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
  firstpart = matrix(0, nrow = d + 2, ncol = d + 2)
  secondpart = rep(0, d + 2)
  KKtt = array(NA, dim = c(t - lag + 1, n, n))
  
  wY = matrix(0, nrow = n, ncol = t + 1)
  for (j in 1:lag){
    wY[, (lag):(t + 1)] = wY[, (lag):(t + 1)] + (w[lag + 1 - j] / sum(w)) * Y[, (lag - (j - 1)):(t + 1 - (j - 1))]
  }
  for (tt in ((max(starttime, 1)):(t - lag + 1))){
    if (d == 1){
      ht = const * n^{ - 1 / expo} * sd(S[, tt, 2:(d + 1)])#bw.nrd(S[, tt, 2:(d + 1)])#
      KK = K(as.matrix(dist(S[, tt, 2:(d + 1)])), h = ht)
    } else {
      ht = const * n^{ - 1 / expo} * apply(S[, tt, 2:(d + 1)], 2, sd)
      KK = K(as.matrix(dist(S[, tt, 2:(d + 1)] %*% diag(1 / ht))), h = 1)
    }
    KKtt[tt, , ] = KK
    Ahat[, tt] = rowMeans(KK * matrix(rep(A[, tt], each = n), ncol = n)) / rowMeans(KK)
    Yhat[, tt] = rowMeans(KK * matrix(rep(wY[, tt + lag], each = n), ncol = n)) / rowMeans(KK)
    AAhat[, tt] = rowMeans(KK * matrix(rep(A[, tt]^2, each = n), ncol = n)) / rowMeans(KK)
    temp =  rbind(t(A[, tt]^2 - AAhat[, tt]), 
                  t(S[, tt, 1:(d + 1)]) %*% 
                    diag(A[, tt] - Ahat[, tt]))
    firstpart = firstpart + temp %*% t(temp)
    secondpart = secondpart + 
      rowSums(temp %*% diag(wY[, tt + lag] - Yhat[, tt]))
  }
  betahat = solve(as.matrix(firstpart)) %*%  secondpart
  if (print){
    cat("\n", "beta:",  betahat,  "\n")
  }
  
  L1hat_phi = matrix(0, nrow = n, ncol = d + 2)
  L2hat = matrix(0, nrow = n, ncol = d + 2)
  for (tt in ((max(starttime, 1)):(t - lag + 1))){
    L1hat_phi[, 1] = L1hat_phi[, 1] + (A[, tt]^2 - AAhat[, tt])^2 * betahat[1] + 
      (A[, tt]^2 - AAhat[, tt]) *  (A[, tt] - Ahat[, tt])  *  (S[, tt, 1:(d + 1)] %*% betahat[2:(d + 2)])
    L1hat_phi[, 2:(d + 2)] = L1hat_phi[, 2:(d + 2)] +  diag((A[, tt] - Ahat[, tt]) * (A[, tt]^2 - AAhat[, tt]) * betahat[1]) %*% S[, tt, 1:(d + 1)] + 
      diag(as.numeric((A[, tt] - Ahat[, tt])^2  *  (S[, tt, 1:(d + 1)] %*% betahat[2:(d + 2)]))) %*%  S[, tt, 1:(d + 1)]
    
    L2hat[, 1] = L2hat[, 1] + (A[, tt]^2 - AAhat[, tt])  * (wY[, tt + lag] - Yhat[, tt])
    L2hat[, 2:(d + 2)] = L2hat[, 2:(d + 2)] + diag((wY[, tt + lag] - Yhat[, tt]) * (A[, tt] - Ahat[, tt])) %*% S[, tt, 1:(d + 1)]
  }
  estimate_var = solve(firstpart / n) %*%  var(L1hat_phi - L2hat) %*% solve(firstpart / n) / n
  
  if (bootstrap>0){
    cat("BOOT: ")
    boot = matrix(NA, nrow = bootstrap, ncol = d + 2)
    for (i in 1:bootstrap){
      cat(i, ", ")
      v = sample(n,  n, replace = TRUE)
      bootS = S[v, , ]
      bootA = A[v, ]
      bootY = Y[v, ]
      bootAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootYhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootAAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootfirstpart = matrix(0, nrow = d + 2, ncol = d + 2)
      bootsecondpart = rep(0, d + 2)
      for (tt in ((max(starttime - lag + 1, 1)):(t - lag + 1))){
        if (d == 1){
          hboot = const * n^{ - 1 / expo} * sd(bootS[, tt, 2:(d + 1)])#bw.nrd(bootS[, tt, 2:(d + 1)])
          KK = K(as.matrix(dist(bootS[, tt, 2:(d + 1)])), h = hboot)
        } else {
          hboot = const * n^{ - 1 / expo} * apply(bootS[, tt, 2:(d + 1)], 2, sd)#bw.nrd(bootS[, tt, 2:(d + 1)])
          KK = K(as.matrix(dist(bootS[, tt, 2:(d + 1)]) %*%  diag(1 / hboot)), h = 1)
        }
        bootAhat[, tt] = rowMeans(KK * matrix(rep(bootA[, tt], each = n), ncol = n)) / rowMeans(KK)
        bootYhat[, tt] = rowMeans(KK * matrix(rep(bootY[, tt + lag], each = n), ncol = n)) / rowMeans(KK)
        bootAAhat[, tt] = rowMeans(KK * matrix(rep(bootA[, tt]^2, each = n), ncol = n)) / rowMeans(KK)
        temp =  rbind(t(bootA[, tt]^2 - bootAAhat[, tt]), 
                      t(bootS[, tt, 1:(d + 1)]) %*% 
                        diag(bootA[, tt] - bootAhat[, tt]))
        bootfirstpart = bootfirstpart + temp %*% t(temp)
        bootsecondpart = bootsecondpart + rowSums(temp %*% diag(bootY[, tt + lag] - bootYhat[, tt]))
      }
      boot[i, ] = solve(as.matrix(bootfirstpart)) %*%  bootsecondpart
    }
    result = list(beta = betahat, sd = apply(boot, 2, sd), estimate_sd = sqrt(diag(estimate_var)))
    cat("\n",  "sd:",  result$sd,  "\n")
  } else{ 
    result = list(beta = betahat, sd = NA, estimate_sd = sqrt(diag(estimate_var)))
  }
  return(result)
}

#Calculate the weighted advantage of the suggested dose
AdvK = function(betahat,  lag = 3, n = 5000, t = 50, 
                   d = 1, starttime = 2, w = c(1 / 3, 1 / 3, 1 / 3), 
                   beta10 = 0, beta11 = 2, 
                   eta1 =  - 0.2, eta2 = 0.2,  
                   theta1 = 0.8, theta2 = 0,  
                   gamma1 = 1,  gamma2 =  - 0.5,  
                   cc = sqrt(0.5), sigma = 0.5, 
                   truebeta){
  alphaK = betahat[1]
  betaK = betahat[2:length(betahat)]
  #generate testing datasets
  
  eps_cov = (sigma^2) * exp(as.matrix(log(cc) * dist(1:(t + 4), upper = TRUE, diag = TRUE)))
  epsilon = as.matrix(mvrnorm(n = n, mu = rep(0, t + 4), Sigma = eps_cov))
  testA = matrix(NA, nrow = n, ncol = t + 3)
  testX = array(NA, dim = c(n, t + 3, d + 1))
  testY = matrix(NA, nrow = n, ncol = t + 4)
  testA[, 1] = runif(n = n, min = 0, max = 1)
  testX[, , 1] = 1
  testX[, 1, 2] = rnorm(n = n, mean = 0, sd = 0.5)  
  
  suggestYw = array(NA, dim = c(n, t))
  originalYw = array(NA, dim = c(n, t))
  suggestA = array(NA, dim = c(n, t))
  
  suggestA[, 1] =  - testX[, 1, ] %*% betaK / (2 * alphaK)
  set.seed(100)
  for (i in (2:t)){
    randomerrorX = rnorm(n = n, mean = 0, sd = 0.5)
    testX[, i, 2] = eta1 * testX[, i - 1, 2] + eta2 * testA[, i - 1] + randomerrorX
    randomerrorA = rnorm(n = n, mean = 0, sd = 0.5)
    testA[, i] = gamma1 * testX[, i, 2] + gamma2 * testA[, i - 1] + randomerrorA
    
    suggestA[, i] =  - testX[, i, ] %*% betaK / (2 * alphaK)
    
    if (i>2){
      testY[, i] = theta1 * testX[, i - 1, 2] + theta2 * (testA[, i - 2])
      - testA[, i - 1] * (testA[, i - 1] - beta10 - beta11 * testX[, i - 1, 2]) + epsilon[, i]
    }
    suggestYw[, i] = truebeta[1] * suggestA[, i]^2 + truebeta[2] * suggestA[, i] + truebeta[3] * suggestA[, i] * testX[, i, 2]
    originalYw[, i] = truebeta[1] * testA[, i]^2   + truebeta[2] * testA[, i]    + truebeta[3] * testA[, i] * testX[, i, 2]
  }
  suggestMeanAdvantage = colMeans(suggestYw[, 1:t])
  originalMeanAdvantage = colMeans(originalYw[, 1:t])
  
  return(list(suggestA = suggestMeanAdvantage, 
              originalA = originalMeanAdvantage))
}






















weightedK_corrected = function(X, A, Y, 
                               lag = 3, 
                               n, t, d0 = 1, starttime = 2, 
                               bootstrap = 0, w = NA,
                               print = "FALSE"){ #k is lag
  if (sum(!is.na(w)) == 0){
    w = rep(1 / lag, lag)
  }
  d = d0 + 1
  Ahat = matrix(NA, nrow = n, ncol = t - lag + 1)
  Yhat = matrix(NA, nrow = n, ncol = t - lag + 1)
  AAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
  firstpart = matrix(0, nrow = d + 2, ncol = d + 2)
  secondpart = rep(0, d + 2)
  KKtt = array(NA, dim = c(t - lag + 1, n, n))
  
  wY = matrix(0, nrow = n, ncol = t + 1)
  for (j in 1:lag){
    wY[, (lag):(t + 1)] = wY[, (lag):(t + 1)] + (w[lag + 1 - j] / sum(w)) * Y[, (lag - (j - 1)):(t + 1 - (j - 1))]
  }
  S = array(NA, dim = c(n, t, d + 1))
  S[, , 1:d] = X
  S[, 2:t, d + 1] = A[, 1:t - 1]
  for (tt in ((max(starttime, 1)):(t - lag + 1))){
    ht = const * n^{ - 1 / expo} * apply(S[, tt, 2:(d + 1)], 2, sd)
    KK = K(as.matrix(dist(S[, tt, 2:(d + 1)] %*% diag(1 / ht))), h = 1)
    KKtt[tt, , ] = KK
    Ahat[, tt] = rowMeans(KK * matrix(rep(A[, tt], each = n), ncol = n)) / rowMeans(KK)
    Yhat[, tt] = rowMeans(KK * matrix(rep(wY[, tt + lag], each = n), ncol = n)) / rowMeans(KK)
    AAhat[, tt] = rowMeans(KK * matrix(rep(A[, tt]^2, each = n), ncol = n)) / rowMeans(KK)
    temp =  rbind(t(A[, tt]^2 - AAhat[, tt]), 
                  t(S[, tt, 1:(d + 1)]) %*% 
                    diag(A[, tt] - Ahat[, tt]))
    firstpart = firstpart + temp %*% t(temp)
    secondpart = secondpart + 
      rowSums(temp %*% diag(wY[, tt + lag] - Yhat[, tt]))
  }
  betahat = solve(as.matrix(firstpart)) %*%  secondpart
  if (print){
    cat("\n", "beta:",  betahat,  "\n")
  }
  
  L1hat_phi = matrix(0, nrow = n, ncol = d + 2)
  L2hat = matrix(0, nrow = n, ncol = d + 2)
  for (tt in ((max(starttime, 1)):(t - lag + 1))){
    L1hat_phi[, 1] = L1hat_phi[, 1] + (A[, tt]^2 - AAhat[, tt])^2 * betahat[1] + 
      (A[, tt]^2 - AAhat[, tt]) *  (A[, tt] - Ahat[, tt])  *  (S[, tt, 1:(d + 1)] %*% betahat[2:(d + 2)])
    L1hat_phi[, 2:(d + 2)] = L1hat_phi[, 2:(d + 2)] +  diag((A[, tt] - Ahat[, tt]) * (A[, tt]^2 - AAhat[, tt]) * betahat[1]) %*% S[, tt, 1:(d + 1)] + 
      diag(as.numeric((A[, tt] - Ahat[, tt])^2  *  (S[, tt, 1:(d + 1)] %*% betahat[2:(d + 2)]))) %*%  S[, tt, 1:(d + 1)]
    
    L2hat[, 1] = L2hat[, 1] + (A[, tt]^2 - AAhat[, tt])  * (wY[, tt + lag] - Yhat[, tt])
    L2hat[, 2:(d + 2)] = L2hat[, 2:(d + 2)] + diag((wY[, tt + lag] - Yhat[, tt]) * (A[, tt] - Ahat[, tt])) %*% S[, tt, 1:(d + 1)]
  }
  estimate_var = solve(firstpart / n) %*%  var(L1hat_phi - L2hat) %*% solve(firstpart / n) / n
  
  if (bootstrap>0){
    cat("BOOT: ")
    boot = matrix(NA, nrow = bootstrap, ncol = d + 2)
    for (i in 1:bootstrap){
      cat(i, ", ")
      v = sample(n,  n, replace = TRUE)
      bootS = S[v, , ]
      bootA = A[v, ]
      bootY = Y[v, ]
      bootAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootYhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootAAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootfirstpart = matrix(0, nrow = d + 2, ncol = d + 2)
      bootsecondpart = rep(0, d + 2)
      for (tt in ((max(starttime - lag + 1, 1)):(t - lag + 1))){
        hboot = const * n^{ - 1 / expo} * apply(bootS[, tt, 2:(d + 1)], 2, sd)#bw.nrd(bootS[, tt, 2:(d + 1)])
        KK = K(as.matrix(dist(bootS[, tt, 2:(d + 1)]) %*%  diag(1 / hboot)), h = 1)
        bootAhat[, tt] = rowMeans(KK * matrix(rep(bootA[, tt], each = n), ncol = n)) / rowMeans(KK)
        bootYhat[, tt] = rowMeans(KK * matrix(rep(bootY[, tt + lag], each = n), ncol = n)) / rowMeans(KK)
        bootAAhat[, tt] = rowMeans(KK * matrix(rep(bootA[, tt]^2, each = n), ncol = n)) / rowMeans(KK)
        temp =  rbind(t(bootA[, tt]^2 - bootAAhat[, tt]), 
                      t(bootS[, tt, 1:(d + 1)]) %*% 
                        diag(bootA[, tt] - bootAhat[, tt]))
        bootfirstpart = bootfirstpart + temp %*% t(temp)
        bootsecondpart = bootsecondpart + rowSums(temp %*% diag(bootY[, tt + lag] - bootYhat[, tt]))
      }
      boot[i, ] = solve(as.matrix(bootfirstpart)) %*%  bootsecondpart
    }
    result = list(beta = betahat, sd = apply(boot, 2, sd), estimate_sd = sqrt(diag(estimate_var)))
    cat("\n",  "sd:",  result$sd,  "\n")
  } else{ 
    result = list(beta = betahat, sd = NA, estimate_sd = sqrt(diag(estimate_var)))
  }
  return(result)
}



























#Generate samples simulation
generate_sample4 = function(seed = 1,  
                            n, t, d = 1, 
                            beta10 = 0, beta11 = 2, 
                            eta1 =  - 0.2, eta2 = 0.2,  
                            theta1 = 0.8, theta2 = 0,  
                            gamma1 = 1,  gamma2 =  - 0.5,  
                            cc = sqrt(0.5), sigma = 0.5){
  set.seed(seed)
  #correlation of the errors
  eps_cov = (sigma^2) * exp(as.matrix(log(cc) * dist(1:(t + 1), upper = TRUE, diag = TRUE)))
  epsilon = as.matrix(mvrnorm(n = n, mu = rep(0, t + 1), Sigma = eps_cov))
  #generate data
  A = matrix(NA, nrow = n, ncol = t + 1) #treatment
  X = array(NA, dim = c(n, t + 1, d + 1)) #covariates dim d
  Y = matrix(NA, nrow = n, ncol = t + 1) #outcome
  A[, 1] = runif(n = n, min = 0, max = 1)
  X[, , 1] = 1
  X[, 1, 2] = rnorm(n = n, mean = 0, sd = 0.5)
  for (i in 2:t){
    X[, i, 2] = rnorm(n = n, mean = eta1 * X[, i - 1, 2] + eta2 * A[, i - 1], sd = 0.5)
    A[, i] = rnorm(n = n, mean = gamma1 * X[, i, 2] + gamma2 * A[, i - 1], sd = 0.5)
    if (i>2) {
      Y[, i] = theta1 * X[, i - 1, 2] + theta2 * (A[, i - 2]) - A[, i - 1] * (A[, i - 1] - beta10 - beta11 * X[, i - 1, 2]) + epsilon[, i]
    }
  }
  Y[, t + 1] = theta1 * (X[, t, 2]) + theta2 * (A[, t - 1]) - A[, t] * (A[, t] - beta10 - beta11 * X[, t, 2]) + epsilon[, t + 1]
  X[, t + 1, 2] = rnorm(n = n, mean = eta1 * X[, t, 2] + eta2 * A[, t], sd = 0.5)
  A[, t + 1] = NA
  Xstacked = X[1, , ]
  Astacked = A[1, ]
  Ystacked = Y[1, ]
  for (i in 2:n){
    Xstacked = rbind(Xstacked, X[i, , ])
    Astacked = c(Astacked, A[i, ])
    Ystacked = c(Ystacked, Y[i, ])
  }
  id = rep(c(1:n), each = t + 1)
  Xstacked = as.data.frame(Xstacked)
  names(Xstacked) = c("Intercept", "Covariate")
  return(list(A = Astacked, X = Xstacked, Y = Ystacked, id = id, n = n, t = t + 1, d = d))
}
#not now,  it is a bit wierd this way
#using all decision points for kernel estimation
weightedK4 = function(X, A, Y, decision = NA, id, 
                    lag, n, t, d = 1, starttime = 2, 
                    bootstrap = 2, w = NA, 
                    seed = 100, constX = 0.305, expo = 3, 
                    decision_var = NA){
  if (sum(!is.na(decision_var)) == 0){
    decision_var = names(X)
  }
  decision_var_index = sapply(decision_var,  function(x) which(names(X) == x))
  n_decision_var = length(decision_var)
  if (sum(!is.na(decision)) == 0){
    decision = rep(TRUE, length(Y))
  }
  set.seed(seed)
  total = n * t - lag
  #first calculate weighted outcome
  if (sum(!is.na(w)) == 0){
    w = c(rep(0, lag - 1), 1)
  }
  wY = rep(0, n * t)
  for (j in 1:lag){
    if (w[lag + 1 - j]>0){
      wY[(lag + 1):(lag + total)] = wY[(lag + 1):(lag + total)] + (w[lag + 1 - j] / sum(w)) * as.matrix(Y[(lag + 1 - (j - 1)):(lag + total - (j - 1))])
    }
  }
  missingY = !complete.cases(wY[(lag + 1):(lag + total)])
  missing = (!complete.cases(X[1:total, 2:(d + 1)])) | missingY
  decision_nonmiss = (!missing)&decision[1:total]
  decision_nonmiss_index = c(1:total)[decision_nonmiss]
  decision_number = sum(decision_nonmiss)
  
  Ahat = rep(NA, decision_number)
  Yhat = rep(NA, decision_number)
  AAhat = matrix(NA, decision_number)
  #######################################
  ##calculate betahat
  if (d>1){
    ht = constX * decision_number^{ - 1 / expo} * sqrt(sum(apply(X[decision_nonmiss_index, 2:(d + 1)], 2, FUN = function(x) var(x, na.rm = TRUE))))
  }else{
    ht = constX * decision_number^{ - 1 / expo} * sd(X[decision_nonmiss_index, 2:(d + 1)], na.rm = TRUE)
  }
    
  Ahat = rowMeans(KK * matrix(rep(A[decision_nonmiss_index], each = decision_number), ncol = decision_number), na.rm = TRUE)  /  rowMeans(KK, na.rm = TRUE)
  Yhat = rowMeans(KK * matrix(rep(wY[decision_nonmiss_index + lag], each = decision_number), ncol = decision_number), na.rm = TRUE)  /  rowMeans(KK, na.rm = TRUE)
  AAhat = rowMeans(KK * matrix(rep(A[decision_nonmiss_index]^2, each = decision_number), ncol = decision_number), na.rm = TRUE) / rowMeans(KK, na.rm = TRUE)
  temp =  rbind(t(A[decision_nonmiss_index]^2 - AAhat), 
              t(X[decision_nonmiss_index, decision_var_index]) %*% 
                makediag(A[decision_nonmiss_index] - Ahat))
  firstpart = temp %*% t(temp)
  secondpart = rowSums(temp %*% makediag(wY[decision_nonmiss_index + lag] - Yhat))
  betahat = solve(as.matrix(firstpart),  tol  =  1e-20) %*%  secondpart
  cat("beta:",  betahat,  "\n")
  ########################################
  ##estimate sd
  L1hat_phi = matrix(0, nrow = n, ncol = n_decision_var + 1)
  L2hat = matrix(0, nrow = n, ncol = n_decision_var + 1)
  id_nonmiss = id[1:total][decision_nonmiss]
  for (i in unique(id_nonmiss)){
    tempi = t(temp)[id_nonmiss == i, ]
    if (sum(id_nonmiss == i)<2){
      tempi = t(tempi)
    }
    L1hat_phi[i, ] = t(t(tempi) %*%  tempi  %*%  betahat)
    L2hat[i, ] = t(t(tempi) %*%  ((wY[(lag + 1):(lag + total)][decision_nonmiss] - Yhat)[id_nonmiss == i]))
  }
  estimate_var = solve(firstpart / n) %*%  var(L1hat_phi - L2hat) %*% solve(firstpart / n) / n
  
  cat('est_sd:', sqrt(diag(estimate_var)), '\n')
  #####################################################
  ##bootstrap to estimate sd
  if (bootstrap>0){
    cat("BOOT: ")
    boot = matrix(NA, nrow = bootstrap, ncol = n_decision_var + 1)
    X_index = as.numeric(sapply(c(1:n), FUN = function(x) (x - 1) * t + 1:(t - lag)))
    Y_index = X_index + lag
    for (i in 1:bootstrap){
      cat(i, ", ")
      v = sample(n,  n, replace = TRUE)
      #find the index that has the days chosen by bootstrap
      index = as.numeric(sapply(v, FUN = function(x) (x - 1) * t + 1:t))
      bootX = X[index, ]
      bootA = A[index]
      bootY = wY[index]
      bootdecision = decision[index]
      bootid = id[index]
      
      
      missingY = !complete.cases(bootY[Y_index])
      missing = (!complete.cases(bootX[X_index, 2:(d + 1)])) | missingY
      decision_nonmiss = X_index[(!missing)&bootdecision[X_index]]
      
      decision_number = length(decision_nonmiss)
      
      Ahat = rep(NA, decision_number)
      Yhat = rep(NA, decision_number)
      AAhat = matrix(NA, decision_number)
      
      ht = constX * decision_number^{ - 1 / expo} * sqrt(sum(apply(bootX[decision_nonmiss, 2:(d + 1)], 2, FUN = function(x) var(x, na.rm = TRUE))))#bw.nrd(X[, tt, 2:(d + 1)])#
      KK = K(as.matrix(dist(bootX[decision_nonmiss, 2:(d + 1)])), h = ht)
      Ahat = rowMeans(KK * matrix(rep(bootA[decision_nonmiss], each = decision_number), ncol = decision_number), na.rm = TRUE)  /  rowMeans(KK, na.rm = TRUE)
      Yhat = rowMeans(KK * matrix(rep(bootY[decision_nonmiss + lag], each = decision_number), ncol = decision_number), na.rm = TRUE)  /  rowMeans(KK, na.rm = TRUE)
      AAhat = rowMeans(KK * matrix(rep(bootA[decision_nonmiss]^2, each = decision_number), ncol = decision_number), na.rm = TRUE) / rowMeans(KK, na.rm = TRUE)
      temp =  rbind(t(bootA[decision_nonmiss]^2 - AAhat), 
                  t(bootX[decision_nonmiss, decision_var_index]) %*% 
                    makediag(bootA[decision_nonmiss] - Ahat))
      firstpart = temp %*% t(temp)
      secondpart = rowSums(temp %*% makediag(bootY[decision_nonmiss + lag] - Yhat))
      tryCatch({boot[i, ] = solve(as.matrix(firstpart),  tol  =  1e-20) %*%  secondpart}, 
               error = function(err){
                 print(err)
               })
    }
    result = list(beta = betahat, sd = apply(boot, 2, sd, na.rm = TRUE), estimate_sd = sqrt(diag(estimate_var)), boot = boot)
    cat("\n",  "sd:",  result$sd,  "\n")
  } else{ 
    result = list(beta = betahat, sd = NA, estimate_sd = sqrt(diag(estimate_var)),  boot = NA)
  }
  return(result)
}


##not now a bit wierd
valueK4 = function(betahat,  X, A, Y, 
                 decision = NA,  id,   #id is day index,  decision T / F is whether this is a decision point
                 lag = 3, n, t, 
                 d = 1, starttime = 2, w = c(1, 0.8, 0.8^2), 
                 maxdose = Inf, 
                 constX = 0.305, constA = 3.05, 
                 expo = 3, 
                 decision_var){ 
  decision_var_index = sapply(decision_var,  function(x) which(names(X) == x))
  
  if (sum(!is.na(decision)) == 0){
    decision = rep(TRUE, length(Y))
  }
  total = n * t - lag
  
  wY = rep(0, n * t)
  for (j in 1:lag){
    if (w[lag + 1 - j]>0){
      wY[(lag + 1):(lag + total)] = wY[(lag + 1):(lag + total)] + (w[lag + 1 - j] / sum(w)) * as.matrix(Y[(lag + 1 - (j - 1)):(lag + total - (j - 1))])
    }
  }
  
  missingY = !complete.cases(wY[(lag + 1):(lag + total)])
  missing = (!complete.cases(X[1:total, 2:(d + 1)])) | missingY
  decision_nonmiss = (!missing)&decision[1:total]
  decision_nonmiss_index = c(1:total)[decision_nonmiss]
  decision_number = sum(decision_nonmiss)
  
  ht = constX * decision_number^{ - 1 / expo} * sqrt(sum(apply(X[decision_nonmiss, 2:(d + 1)], 2, FUN = function(x) var(x, na.rm = TRUE))))#bw.nrd(X[, tt, 2:(d + 1)])#
  hta = constA * decision_number^{ - 1 / expo} * sd(A[decision_nonmiss_index])#bw.nrd()
  KK = K(as.matrix(dist(X[decision_nonmiss_index, 2:(d + 1)])), h = ht)
  Ahat = as.matrix(X[decision_nonmiss_index, decision_var_index]) %*% betahat[ - 1] / ( - 2 * betahat[1])
  Ahat[Ahat<0] = 0
  Ahat[Ahat>maxdose] = maxdose
  KKa = K(as.matrix(pdist(as.matrix(Ahat), as.matrix(A[decision_nonmiss_index]))), h = hta)
  KKa_orig = K(as.matrix(dist(as.matrix(A[decision_nonmiss_index]))), h = hta)
  wYhat = rowSums(KK * KKa * 
                  matrix(rep(wY[decision_nonmiss_index + lag], each = decision_number), ncol = decision_number), 
                na.rm = TRUE) / 
    rowSums(KK * KKa, na.rm = TRUE) 
  wYhat_orig = rowSums(KK * KKa_orig * 
                       matrix(rep(wY[decision_nonmiss_index + lag], each = decision_number), ncol = decision_number), 
                     na.rm = TRUE) / 
    rowSums(KK * KKa_orig, na.rm = TRUE) 
  result = list(wY = wY[decision_nonmiss_index + lag],    #original wY
              wYhat = wYhat,   #estimated wY with suggested dose
              wYhat_orig = wYhat_orig,  #estimated wY with original dose
              dosehat = Ahat,  #suggested dose
              dose = A,   #original dose
              id = id[decision_nonmiss_index])
  return(result)
}



lagk = function(X, A, Y, lag, n, t, d = 1, starttime = 2, bootstrap = 0,
                print = "FALSE"){ #k is lag
  Ahat = matrix(NA, nrow = n, ncol = t - lag + 1)
  Yhat = matrix(NA, nrow = n, ncol = t - lag + 1)
  AAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
  
  firstpart = matrix(0, nrow = d + 2, ncol = d + 2)
  secondpart = rep(0, d + 2)
  KKtt = array(NA, dim = c(t - lag + 1, n, n))
  
  for (tt in ((max(starttime - lag + 1, 1)):(t - lag + 1))){
    ht = const * n^{ - 1 / expo} * sd(X[, tt, 2:(d + 1)])#bw.nrd(X[, tt, 2:(d + 1)])#
    KK = K(as.matrix(dist(X[, tt, 2:(d + 1)])), h = ht)
    KKtt[tt, , ] = KK
    Ahat[, tt] = rowMeans(KK * matrix(rep(A[, tt], each = n), ncol = n)) / rowMeans(KK)
    Yhat[, tt] = rowMeans(KK * matrix(rep(Y[, tt + lag], each = n), ncol = n)) / rowMeans(KK)
    AAhat[, tt] = rowMeans(KK * matrix(rep(A[, tt]^2, each = n), ncol = n)) / rowMeans(KK)
    temp =  rbind(t(A[, tt]^2 - AAhat[, tt]), 
                t(X[, tt, 1:(d + 1)]) %*% 
                  diag(A[, tt] - Ahat[, tt]))
    firstpart = firstpart + temp %*% t(temp)
    secondpart = secondpart + 
      rowSums(temp %*% diag(Y[, tt + lag] - Yhat[, tt]))
  }
  betahat = solve(as.matrix(firstpart)) %*%  secondpart
  
  L1hat_phi = matrix(0, nrow = n, ncol = d + 2)
  L2hat = matrix(0, nrow = n, ncol = d + 2)
  for (tt in ((max(starttime - lag + 1, 1)):(t - lag + 1))){
    L1hat_phi[, 1] = L1hat_phi[, 1] + (A[, tt]^2 - AAhat[, tt])^2 * betahat[1] + 
      (A[, tt]^2 - AAhat[, tt]) *  (A[, tt] - Ahat[, tt])  *  (X[, tt, 1:(d + 1)] %*% betahat[2:(d + 2)])
    L1hat_phi[, 2:(d + 2)] = L1hat_phi[, 2:(d + 2)] +  diag((A[, tt] - Ahat[, tt]) * (A[, tt]^2 - AAhat[, tt]) * betahat[1]) %*% X[, tt, 1:(d + 1)] + 
      diag(as.numeric((A[, tt] - Ahat[, tt])^2  *  (X[, tt, 1:(d + 1)] %*% betahat[2:(d + 2)]))) %*%  X[, tt, 1:(d + 1)]
    
    L2hat[, 1] = L2hat[, 1] + (A[, tt]^2 - AAhat[, tt])  * (Y[, tt + lag] - Yhat[, tt])
    L2hat[, 2:(d + 2)] = L2hat[, 2:(d + 2)] + diag((Y[, tt + lag] - Yhat[, tt]) * (A[, tt] - Ahat[, tt])) %*% X[, tt, 1:(d + 1)]
  }
  estimate_var = solve(firstpart / n) %*%  var(L1hat_phi - L2hat) %*% solve(firstpart / n) / n
  
  if (print){
    cat("\n", "beta:",  betahat,  "\n")
  }
  if (bootstrap>0){
    cat("BOOT: ")
    boot = matrix(NA, nrow = bootstrap, ncol = d + 2)
    for (i in 1:bootstrap){
      cat(i, ", ")
      v = sample(n,  n, replace = TRUE)
      bootX = X[v, , ]
      bootA = A[v, ]
      bootY = Y[v, ]
      bootAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootYhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootAAhat = matrix(NA, nrow = n, ncol = t - lag + 1)
      bootfirstpart = matrix(0, nrow = d + 2, ncol = d + 2)
      bootsecondpart = rep(0, d + 2)
      for (tt in ((max(starttime - lag + 1, 1)):(t - lag + 1))){
        hboot = const * n^{ - 1 / expo} * sd(bootX[, tt, 2:(d + 1)])#bw.nrd(bootX[, tt, 2:(d + 1)])
        KK = K(as.matrix(dist(bootX[, tt, 2:(d + 1)])), h = hboot)
        bootAhat[, tt] = rowMeans(KK * matrix(rep(bootA[, tt], each = n), ncol = n)) / rowMeans(KK)
        bootYhat[, tt] = rowMeans(KK * matrix(rep(bootY[, tt + lag], each = n), ncol = n)) / rowMeans(KK)
        bootAAhat[, tt] = rowMeans(KK * matrix(rep(bootA[, tt]^2, each = n), ncol = n)) / rowMeans(KK)
        temp =  rbind(t(bootA[, tt]^2 - bootAAhat[, tt]), 
                    t(bootX[, tt, 1:(d + 1)]) %*% 
                      diag(bootA[, tt] - bootAhat[, tt]))
        bootfirstpart = bootfirstpart + temp %*% t(temp)
        bootsecondpart = bootsecondpart + rowSums(temp %*% diag(bootY[, tt + lag] - bootYhat[, tt]))
      }
      boot[i, ] = solve(as.matrix(bootfirstpart)) %*%  bootsecondpart
    }
    result = list(beta = betahat, sd = apply(boot, 2, sd), estimate_sd = sqrt(diag(estimate_var)))
    cat("\n",  "sd:",  result$sd,  "\n")
  } else{ 
    result = list(beta = betahat, sd = NA, estimate_sd = sqrt(diag(estimate_var)))
  }
  return(result)
}






