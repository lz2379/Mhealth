##########Simulation code######################

############################
############################
# Scenario 1 with theta2 = 0
############################

source("simulation-functions.R")

#constants for bandwidths
const = 0.305
expo = 3

d = 1    #dimension for X
it = 200 #iteration
t = 50   #total time points
LAGMAX = 3

#define true parameters
true_parameters = matrix(NA, nrow = LAGMAX + 1, ncol = d + 2)
true_parameters[1, ] = c(-1, 0, 2)
true_parameters[2, ] = c(-0.21, 0.16, -0.08)
true_parameters[3, ] = c(-0.0125, -0.08, -0.03)
true_parameters[4, ] = colMeans(true_parameters[1:3, ])

finalresult = list()
for (n in c(100, 200, 400)){ #sample size

  #initialize datasets for saving results
  finalresult[[n]] = list()
  for (lag in 1:LAGMAX){
    finalresult[[n]][[lag]] = list(beta = matrix(NA, nrow = it, ncol = d + 2), 
                                   est_sd = matrix(NA, nrow = it, ncol = d + 2))  
  }
  finalresult[[n]][[LAGMAX + 1]] = list(beta = matrix(NA, nrow = it, ncol = d + 2), 
                                   est_sd = matrix(NA, nrow = it, ncol = d + 2),
                                   advantage = matrix(NA, nrow = it, ncol = t))
  
  cat("Sample size:",  n,  "\n")
  cat("Generating samples:")
  for (i in 1:it){
    
    #generate dataset
    set.seed(20000 + i)
    cat(i, ", ")
    data = generate_sample(seed = i, n = n, t = t, d = d)
    A = data$A
    X = data$X
    Y = data$Y
    
    #estimate parameters for each lag
    for (lag in 1:LAGMAX){
      weights=rep(0, LAGMAX)
      weights[lag] = 1
      res = weightedK(S = X, 
                      A = A, 
                      Y = Y, 
                      lag = lag, w = weights, 
                      n = n, t = t, d = d, 
                      starttime = 2, 
                      bootstrap = 0)
      finalresult[[n]][[lag]]$beta[i, ] = res$beta
      finalresult[[n]][[lag]]$est_sd[i, ] = res$estimate_sd
    }
    
    #estimate parameters for weighted result
    weights=rep(1/LAGMAX, LAGMAX)
    res = weightedK(S = X, A, Y, 
                    lag = lag, w = weights, 
                    n = n, t = t, d = d, 
                    starttime = 2, 
                    bootstrap = 0)
    finalresult[[n]][[lag + 1]]$beta[i, ] = res$beta
    finalresult[[n]][[lag + 1]]$est_sd[i, ] = res$estimate_sd
  }
  cat("\n")
  
  #Print out summary of results:
  for (lag in 1:LAGMAX){
    beta = finalresult[[n]][[lag]]$beta
    est_sd = finalresult[[n]][[lag]]$est_sd
    
    #coverage of confidence intervals
    lower = beta - 1.96 * est_sd
    upper = beta + 1.96 * est_sd
    coverage = colMeans( (lower - matrix(rep(true_parameters[lag, ], each = it), nrow = it)) * 
                           (upper - matrix(rep(true_parameters[lag, ], each = it), nrow = it)) < 0 )   
    cat("--- Lag", lag, "---\n")
    cat("Means of Beta_hat:",  colMeans(beta), "\n")
    cat("Sd    of Beta_hat:",  apply(beta, 2, sd),  "\n")
    cat("SE               :", colMeans(est_sd),  "\n")
    cat("Coverage of 95%CI:", coverage, "\n" )

  }
  
  beta = finalresult[[n]][[LAGMAX + 1]]$beta
  est_sd = finalresult[[n]][[LAGMAX + 1]]$est_sd
  lower = beta - 1.96 * est_sd
  upper = beta + 1.96 * est_sd
  coverage = colMeans( (lower - matrix(rep(true_parameters[LAGMAX + 1, ], each = it), nrow = it)) * 
                         (upper - matrix(rep(true_parameters[LAGMAX + 1, ], each = it), nrow = it)) < 0 )   
  cat("--- Weighted ---\n")
  cat("Means of Beta_hat:",  colMeans(beta), "\n")
  cat("Sd    of Beta_hat:",  apply(beta, 2, sd),  "\n")
  cat("SE               :", colMeans(est_sd),  "\n")
  cat("Coverage of 95%CI:", coverage, "\n" )
  
  ####### Calculate weighted advantage of the suggested dose ########
  for (i in 1:it){
    v = AdvK(betahat = beta[i, ],
               truebeta = true_parameters[4, ])
    finalresult[[n]][[LAGMAX + 1]]$advantage[i, ] = v$suggestA
  }
  optimal_advantage = AdvK(betahat = true_parameters[LAGMAX + 1, ],
                             truebeta = true_parameters[4, ])
  cat("--- Advantage ---\n")
  cat("Mean of Advantage:", mean(finalresult[[n]][[LAGMAX + 1]]$advantage[, 2:t]), "\n")
  cat("Optimal Advantage:", mean(optimal_advantage$suggestA[2:t]), "\n")
}






############################
############################
# Scenario 2 with theta2 = -0.1 
# Biased Estimation S_t=X_t
############################


theta2 = -0.1 

#constants for bandwidths
const = 0.305
expo = 3

d = 1    #dimension for S
it = 200 #iteration
t = 50   #total time points
LAGMAX = 3

#define true parameters
true_parameters = matrix(NA, nrow = LAGMAX + 1, ncol = d + 2)
true_parameters[1, ] = c(-1, 0, 2)
true_parameters[2, ] = c(-0.21, 0.06, -0.08)
true_parameters[3, ] = c(-0.0125, -0.05, -0.03)
true_parameters[4, ] = colMeans(true_parameters[1:3, ])

finalresult = list()
for (n in c(100, 200, 400)){ #sample size
  
  #initialize datasets for saving results
  finalresult[[n]] = list()
  for (lag in 1:LAGMAX){
    finalresult[[n]][[lag]] = list(beta = matrix(NA, nrow = it, ncol = d + 2), 
                                   est_sd = matrix(NA, nrow = it, ncol = d + 2))  
  }
  finalresult[[n]][[LAGMAX + 1]] = list(beta = matrix(NA, nrow = it, ncol = d + 2), 
                                        est_sd = matrix(NA, nrow = it, ncol = d + 2),
                                        advantage = matrix(NA, nrow = it, ncol = t))
  
  cat("Sample size:",  n,  "\n")
  cat("Generating samples:")
  for (i in 1:it){
    
    #generate dataset
    set.seed(20000 + i)
    cat(i, ", ")
    data = generate_sample(seed = i, n = n, t = t, d = d, 
                            theta2 = theta2)
    A = data$A
    X = data$X
    Y = data$Y
    
    #estimate parameters for each lag
    for (lag in 1:LAGMAX){
      weights=rep(0, LAGMAX)
      weights[lag] = 1
      res = weightedK(S = X, 
                      A = A, 
                      Y = Y, 
                      lag = lag, w = weights, 
                      n = n, t = t, d = d, 
                      starttime = 2, 
                      bootstrap = 0)
      finalresult[[n]][[lag]]$beta[i, ] = res$beta
      finalresult[[n]][[lag]]$est_sd[i, ] = res$estimate_sd
    }
    
    #estimate parameters for weighted result
    weights=rep(1/LAGMAX, LAGMAX)
    res = weightedK(S = X, 
                    A = A, 
                    Y = Y, 
                    lag = lag, w = weights, 
                    n = n, t = t, d = d, 
                    starttime = 2, 
                    bootstrap = 0)
    finalresult[[n]][[lag + 1]]$beta[i, ] = res$beta
    finalresult[[n]][[lag + 1]]$est_sd[i, ] = res$estimate_sd
  }
  cat("\n")
  
  #Print out summary of results:
  for (lag in 1:LAGMAX){
    beta = finalresult[[n]][[lag]]$beta
    est_sd = finalresult[[n]][[lag]]$est_sd
    
    #coverage of confidence intervals
    lower = beta - 1.96 * est_sd
    upper = beta + 1.96 * est_sd
    coverage = colMeans( (lower - matrix(rep(true_parameters[lag, ], each = it), nrow = it)) * 
                           (upper - matrix(rep(true_parameters[lag, ], each = it), nrow = it)) < 0 )   
    cat("--- Lag", lag, "---\n")
    cat("Means of Beta_hat:",  colMeans(beta), "\n")
    cat("Sd    of Beta_hat:",  apply(beta, 2, sd),  "\n")
    cat("SE               :", colMeans(est_sd),  "\n")
    cat("Coverage of 95%CI:", coverage, "\n" )
    
  }
  
  beta = finalresult[[n]][[LAGMAX + 1]]$beta
  est_sd = finalresult[[n]][[LAGMAX + 1]]$est_sd
  lower = beta - 1.96 * est_sd
  upper = beta + 1.96 * est_sd
  coverage = colMeans( (lower - matrix(rep(true_parameters[LAGMAX + 1, ], each = it), nrow = it)) * 
                         (upper - matrix(rep(true_parameters[LAGMAX + 1, ], each = it), nrow = it)) < 0 )   
  cat("--- Weighted ---\n")
  cat("Means of Beta_hat:",  colMeans(beta), "\n")
  cat("Sd    of Beta_hat:",  apply(beta, 2, sd),  "\n")
  cat("SE               :", colMeans(est_sd),  "\n")
  cat("Coverage of 95%CI:", coverage, "\n" )
  
  ####### Calculate weighted advantage of the suggested dose ########
  for (i in 1:it){
    v = AdvK(betahat = beta[i, ], 
               truebeta = true_parameters[4, ])
    finalresult[[n]][[LAGMAX + 1]]$advantage[i, ] = v$suggestA
  }
  optimal_advantage = AdvK(betahat = true_parameters[LAGMAX + 1, ],
                             truebeta = true_parameters[4, ])
  cat("--- Advantage ---\n")
  cat("Mean of Advantage:", mean(finalresult[[n]][[LAGMAX + 1]]$advantage[, 2 : t]),"\n")
  cat("Optimal Advantage:", mean(optimal_advantage$suggestA[2 : t]),"\n")
}



#################################
#
#Bias-corrected version with S_t=(X_t, A_{t-1})
#################################

theta2 = -0.1 

#constants for bandwidths
const = 1.2
expo = 4

d = 1    #dimension for X
it = 200 #iteration
t = 50   #total time points
LAGMAX = 3

#define true parameters
true_parameters = matrix(NA, nrow = LAGMAX + 1, ncol = d + 3)
true_parameters[1, ] = c(-1, 0, 2, 0)
true_parameters[2, ] = c(-0.21, 0.06, -0.08, 0)
true_parameters[3, ] = c(-0.0125, -0.05, -0.03, 0)
true_parameters[4, ] = colMeans(true_parameters[1:3, ])

finalresult = list()
for (n in c(100, 200, 400)){ #sample size
  
  #initialize datasets for saving results
  finalresult[[n]] = list()
  for (lag in 1:LAGMAX){
    finalresult[[n]][[lag]] = list(beta = matrix(NA, nrow = it, ncol = d + 3), 
                                   est_sd = matrix(NA, nrow = it, ncol = d + 3))  
  }
  finalresult[[n]][[LAGMAX + 1]] = list(beta = matrix(NA, nrow = it, ncol = d + 3), 
                                        est_sd = matrix(NA, nrow = it, ncol = d + 3),
                                        advantage = matrix(NA, nrow = it, ncol = t))
  
  cat("Sample size:",  n,  "\n")
  cat("Generating samples:")
  for (i in 1:it){
    
    #generate dataset
    set.seed(20000 + i)
    cat(i, ", ")
    data = generate_sample(seed = i, n = n, t = t, d = d, 
                            theta2 = theta2)
    A = data$A
    X = data$X
    Y = data$Y
    
    #Create S_t which include both X_t and A_{t-1}
    S = array(NA, dim = c(n, t, d + 2))
    S[, , 1:(d + 1)] = X
    S[, 2:t, d + 2] = A[, 1:(t - 1)]
    
    #estimate parameters for each lag
    for (lag in 1:LAGMAX){
      weights=rep(0, LAGMAX)
      weights[lag] = 1
      res = weightedK(S = S, 
                      A = A, 
                      Y = Y, 
                      lag = lag, 
                      w = weights, 
                      n = n, t = t, d = d + 1, 
                      starttime = 2, 
                      bootstrap = 0)
      finalresult[[n]][[lag]]$beta[i, ] = res$beta
      finalresult[[n]][[lag]]$est_sd[i, ] = res$estimate_sd
    }
    
    #estimate parameters for weighted result
    weights=rep(1 / LAGMAX, LAGMAX)
    res = weightedK(S = S, 
                    A = A, 
                    Y = Y, 
                    lag = lag, w = weights, 
                    n = n, t = t, d = d + 1, 
                    starttime = 2, 
                    bootstrap = 0)
    finalresult[[n]][[lag + 1]]$beta[i, ] = res$beta
    finalresult[[n]][[lag + 1]]$est_sd[i, ] = res$estimate_sd
  }
  cat("\n")
  
  #Print out summary of results:
  for (lag in 1:LAGMAX){
    beta = finalresult[[n]][[lag]]$beta
    est_sd = finalresult[[n]][[lag]]$est_sd
    
    #coverage of confidence intervals
    lower = beta - 1.96 * est_sd
    upper = beta + 1.96 * est_sd
    coverage = colMeans( (lower - matrix(rep(true_parameters[lag, ], each = it), nrow = it)) * 
                           (upper - matrix(rep(true_parameters[lag, ], each = it), nrow = it)) < 0 )   
    cat("--- Lag", lag, "---\n")
    cat("Means of Beta_hat:",  colMeans(beta), "\n")
    cat("Sd    of Beta_hat:",  apply(beta, 2, sd),  "\n")
    cat("SE               :", colMeans(est_sd),  "\n")
    cat("Coverage of 95%CI:", coverage, "\n" )
    
  }
  
  beta = finalresult[[n]][[LAGMAX + 1]]$beta
  est_sd = finalresult[[n]][[LAGMAX + 1]]$est_sd
  lower = beta - 1.96 * est_sd
  upper = beta + 1.96 * est_sd
  coverage = colMeans( (lower - matrix(rep(true_parameters[LAGMAX + 1, ], each = it), nrow = it)) * 
                       (upper - matrix(rep(true_parameters[LAGMAX + 1, ], each = it), nrow = it)) < 0 )   
  cat("--- Weighted ---\n")
  cat("Means of Beta_hat:",  colMeans(beta), "\n")
  cat("Sd    of Beta_hat:",  apply(beta, 2, sd),  "\n")
  cat("SE               :", colMeans(est_sd),  "\n")
  cat("Coverage of 95%CI:", coverage, "\n" )

}