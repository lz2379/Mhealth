#-----------------------------------
# Functions for running model with real data
#-----------------------------------

library(MASS)
library(pdist)

#kernel function
K = function(x, h){
  return(1 / (h * (sqrt(2 * pi))) * exp( - (x / h)^2 / 2))
}

#create diagonal matrix considering one dimensional vector
makediag = function(x){
  if (length(x) == 1) {
    return(as.matrix(x))
  } else{
    return(diag(x))
  }
}

#calculate the lag K weighted advantage
weightedK = function(S, A, Y, 
                     decision = NA, 
                     # decision is a boolean variable indicating whether this time point is taken as a decision point, 
                     # used if we are not treating all timepoints as decision point
                     id, #id is day index
                     lag = 3, 
                     n, t, d = 1, 
                     lagmax = 2, 
                     bootstrap = 0, w = c(1, 0.8, 0.64), 
                     seed = 100,  
                     constS = 0.305,  expo = 3, 
                     decision_var
                     ){ 
  
  decision_var_index = sapply(decision_var,  function(x) which(names(S) == x))
  n_decision_var = length(decision_var)
  
  if (sum(!is.na(decision)) == 0){
    decision = rep(TRUE, length(Y))
  }
  
  set.seed(seed)
  total = n * t - lag
  
  #first calculate weighted outcome
  wY = rep(0, n * t)
  for (j in 1:lag){
    if (w[lag + 1 - j] > 0){
      wY[(lag + 1):(lag + total)] = wY[(lag + 1):(lag + total)] + (w[lag + 1 - j] / sum(w)) * as.matrix(Y[(lag + 1 - (j - 1)):(lag + total - (j - 1))])
    }
  }
  
  #exclude non decision points and missing points
  missingY = !complete.cases(wY[(lag + 1):(lag + total)])
  missing = (!complete.cases(S[1:total, 2:(d + 1)])) | missingY
  decision_nonmiss = (!missing)&decision[1:total]
  decision_nonmiss_index = c(1:total)[decision_nonmiss]
  decision_number = sum(decision_nonmiss)
  
  Ahat = rep(NA, decision_number)
  Yhat = rep(NA, decision_number)
  AAhat = matrix(NA, decision_number)
  
  ##calculate betahat
  ht = constS * decision_number^{ - 1 / expo} * sqrt(sum(apply(S[decision_nonmiss_index, 2:(d + 1)], 2, FUN = function(x) var(x, na.rm = TRUE))))#bw.nrd(S[, tt, 2:(d + 1)])#
  KK = K(as.matrix(dist(S[decision_nonmiss_index, 2:(d + 1)])), h = ht)
  Ahat = rowMeans(KK * matrix(rep(A[decision_nonmiss_index], each = decision_number), ncol = decision_number), na.rm = TRUE)  /  rowMeans(KK, na.rm = TRUE)
  Yhat = rowMeans(KK * matrix(rep(wY[decision_nonmiss_index + lag], each = decision_number), ncol = decision_number), na.rm = TRUE)  /  rowMeans(KK, na.rm = TRUE)
  AAhat = rowMeans(KK * matrix(rep(A[decision_nonmiss_index]^2, each = decision_number), ncol = decision_number), na.rm = TRUE) / rowMeans(KK, na.rm = TRUE)
  temp =  rbind(t(A[decision_nonmiss_index]^2 - AAhat), 
              t(S[decision_nonmiss_index, decision_var_index]) %*% 
                makediag(A[decision_nonmiss_index] - Ahat))
  firstpart = temp %*% t(temp)
  secondpart = rowSums(temp %*% makediag(wY[decision_nonmiss_index + lag] - Yhat))
  betahat = solve(as.matrix(firstpart),  tol  =  1e-20) %*%  secondpart
  cat("beta:",  betahat,  "\n")
  
  ##estimate sd
  L1hat_phi = matrix(0, nrow = n, ncol = n_decision_var + 1)
  L2hat = matrix(0, nrow = n, ncol = n_decision_var + 1)
  id_nonmiss = id[1:total][decision_nonmiss]
  for (i in unique(id_nonmiss)){
    tempi = t(temp)[id_nonmiss == i, ]
    if (sum(id_nonmiss == i) < 2){
      tempi = t(tempi)
    }
    L1hat_phi[i, ] = t(t(tempi) %*%  tempi  %*%  betahat)
    L2hat[i, ] = t(t(tempi) %*%  ((wY[(lag + 1):(lag + total)][decision_nonmiss] - Yhat)[id_nonmiss == i]))
  }
  estimate_var = solve(firstpart / n) %*%  var(L1hat_phi - L2hat) %*% solve(firstpart / n) / n
  
  cat('est_sd:', sqrt(diag(estimate_var)), '\n')
  
  ##bootstrap to estimate sd
  if (bootstrap  >  0){
    cat("BOOT: ")
    boot = matrix(NA, nrow = bootstrap, ncol = n_decision_var + 1)
    S_index = as.numeric(sapply(c(1:n), FUN = function(x) (x - 1) * t + 1:(t - lag)))
    Y_index = S_index  +  lag
    for (i in 1:bootstrap){
      cat(i, ", ")
      v = sample(n,  n, replace = TRUE)
      #find the index that has the days chosen by bootstrap
      index = as.numeric(sapply(v, FUN = function(x) (x - 1) * t + 1:t))
      bootS = S[index, ]
      bootA = A[index]
      bootY = wY[index]
      bootdecision = decision[index]
      bootid = id[index]
      
      
      missingY = !complete.cases(bootY[Y_index])
      missing = (!complete.cases(bootS[S_index, 2:(d + 1)])) | missingY
      decision_nonmiss = S_index[(!missing)&bootdecision[S_index]]
      
      decision_number = length(decision_nonmiss)
      
      Ahat = rep(NA, decision_number)
      Yhat = rep(NA, decision_number)
      AAhat = matrix(NA, decision_number)
      
      ht = constS * decision_number^{ - 1 / expo} * sqrt(sum(apply(bootS[decision_nonmiss, 2:(d + 1)], 2, FUN = function(x) var(x, na.rm = TRUE))))#bw.nrd(S[, tt, 2:(d + 1)])#
      KK = K(as.matrix(dist(bootS[decision_nonmiss, 2:(d + 1)])), h = ht)
      Ahat = rowMeans(KK * matrix(rep(bootA[decision_nonmiss], each = decision_number), ncol = decision_number), na.rm = TRUE)  /  rowMeans(KK, na.rm = TRUE)
      Yhat = rowMeans(KK * matrix(rep(bootY[decision_nonmiss + lag], each = decision_number), ncol = decision_number), na.rm = TRUE)  /  rowMeans(KK, na.rm = TRUE)
      AAhat = rowMeans(KK * matrix(rep(bootA[decision_nonmiss]^2, each = decision_number), ncol = decision_number), na.rm = TRUE) / rowMeans(KK, na.rm = TRUE)
      temp =  rbind(t(bootA[decision_nonmiss]^2 - AAhat), 
                  t(bootS[decision_nonmiss, decision_var_index]) %*% 
                    makediag(bootA[decision_nonmiss] - Ahat))
      firstpart = temp %*% t(temp)
      secondpart = rowSums(temp %*% makediag(bootY[decision_nonmiss + lag] - Yhat))
      tryCatch({boot[i, ] = solve(as.matrix(firstpart),  tol  =  1e-20) %*%  secondpart}, 
               error = function(err){
                 print(err)
               })
    }
    result = list(beta = betahat, 
                  sd = apply(boot, 2, sd, na.rm = TRUE), 
                  estimate_sd = sqrt(diag(estimate_var)), 
                  boot = boot)
    cat("\n",  "sd:",  result$sd,  "\n")
  } else{ 
    result = list(beta = betahat, 
                  sd = NA, 
                  estimate_sd = sqrt(diag(estimate_var)),  
                  boot = NA)
  }
  return(result)
}

#print the estimated parameters
printresult = function(model){
  a = data.frame(beta = model$result$beta,  
                 estimate_sd = model$result$estimate_sd,  
                 boot_sd = model$result$sd)
  row.names(a) = c("dose^2", "dose", model$decision_var)
  print(a)
}

#Calculate the weighted advantage of the suggested dose
AdvK  = function(betahat,  S, A, Y, 
                 decision = NA,  id,   #id is day index,  decision T / F is whether this is a decision point
                 lag = 3, n, t, 
                 d = 1, 
                 lagmax = 2, 
                 w = c(1, 0.8, 0.8^2), 
                 maxdose = Inf, 
                 constS = 0.305, constA = 3.05, 
                 expo = 3, 
                 decision_var, 
                 times){ 
  decision_var_index = sapply(decision_var,  function(x) which(names(S) == x))
  
  if (sum(!is.na(decision)) == 0){
    decision = rep(TRUE, length(Y))
  }
  total = n * t - lag
  missingY = rep(FALSE,  total)
  missing = (!complete.cases(S[1:total, 2:(d + 1)]))
  decision_nonmiss = (!missing)&decision[1:total]
  decision_nonmiss_index = c(1:total)[decision_nonmiss]
  decision_number = sum(decision_nonmiss)
  
  betaS = as.matrix(S[decision_nonmiss_index, decision_var_index]) %*% betahat[ - 1]
  alpha = betahat[1]
  Ahat = betaS / ( - 2 * alpha)
  Ahat[Ahat < 0] = 0
  Ahat[Ahat > maxdose] = maxdose
  Aorig = A[decision_nonmiss_index]
  advantage_suggest = alpha * Ahat^2 + betaS * Ahat
  advantage_orig = alpha * Aorig^2 + betaS * Aorig
  result = list(advantage_suggest = advantage_suggest, 
              advantage_orig = advantage_orig, 
              dosehat = Ahat,  #suggested dose
              dose = Aorig,   #original dose
              id = id[decision_nonmiss_index], 
              times = times[decision_nonmiss_index], 
              nonmiss_index = decision_nonmiss_index)
  return(result)
}




#use all decisions
runmodel = function(personID, 
                   outcome = "icg", 
                   covariates, 
                   dose = "bolus", 
                   decision_var = NA,  #variables for deciding dosages
                   day = "day", 
                   time = "tt", 
                   LAG, 
                   interval = 30, 
                   weights = NA, 
                   bootstrap = 0, 
                   test_index = NA, 
                   constX = 0.305,  
                   constA = 3.05, 
                   expo = 3){
  cat("Person       :",  personID, "\n")
  cat("Model        :",  outcome,  "~",  covariates,  "| Dose = ", dose, "\n")
  cat("Decision Var :",  decision_var, "\n")
  cat("ID           :",  day,  "\n")
  cat("Time_interval:",  interval,  "minutes \n")
  cat("LAG          :",  LAG,  "\n")
  data = data_full[[personID]]
  day_index = which(names(data) == day)
  time_index = which(names(data) == time)
  outcome_index = which(names(data) == outcome)
  covariates_index = sapply(covariates,  function(x) which(names(data) == x))
  dose_index = which(names(data) == dose)
  
  Y = data[, outcome_index]
  A = data[, dose_index]
  id = data[, day_index]
  
  n = length(unique(data[, day_index])) #number of days
  TT = length(unique(data[, time_index])) #number of timepoints each day
  d = length(covariates)
  X = cbind(1,  data[, covariates_index])
  names(X) = c("Intercept", covariates)
  
  if (sum(!is.na(test_index)) == 0){ 
    #if test data is not specified,  just let the original dataset as testdata
    testY = Y
    testX = X
    testA = A
    testid = id
    testn = n
  } else{
    testn = length(test_index)
    test_index = as.numeric(sapply(test_index,  FUN = function(x) (x - 1) * TT + 1:TT))
    testY = Y[test_index]
    testX = X[test_index, ]
    testA = A[test_index]
    testid = id[test_index]
    Y = Y[ - test_index]
    X = X[ - test_index, ]
    A = A[ - test_index]
    id = id[ - test_index]
  }
  
  if (sum(!is.na(weights)) == 0){
    #result = lagk4(X, A, Y, decision, id = id,  lag = LAG, n = n, t = TT, d = d, lagmax = 1, bootstrap = bootstrap)
    weights = c(rep(0, LAG - 1), 1)
  }
  if (sum(!is.na(decision_var)) == 0){
    decision_var = covariates
  }
  result = weightedK(S = X, A, Y, 
                     decision = NA, 
                     id = id,  
                     lag = LAG, 
                     n = n, t = TT, d = d, 
                     lagmax = 1, 
                     w = weights, 
                     bootstrap = bootstrap, 
                     constS = constX, expo = expo, 
                     decision_var = c("Intercept", decision_var))
  betahat = result$beta
  maxdose = max(A)
  value = AdvK(betahat, 
               S = testX, A = testA, Y = testY, 
               decision = NA,  
               id = testid , 
               lag = LAG, 
               n = testn, t = TT, d = d, 
               lagmax = 1, 
               w = weights, 
               maxdose = maxdose, 
               constS = constX, constA = constA, expo = expo, 
              decision_var = c("Intercept", decision_var),  times = data[, time_index])
  
  return(list(result = result, 
              Y = Y, 
              X = X, 
              A = A, 
              covariates = covariates, 
              value = value, 
              decision_var = decision_var))
}