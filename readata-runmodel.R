#-----------------------------------------------
# Run model for the Type 1 Diabetes dataset
load("data_preprocessed-halfhour.Rdata")
source("realdata-runmodel-functions.R")

nperson = 6

#Extract the days for test dataset
testday = list()
for (i in 1:nperson){
  testday[[i]] = as.numeric(unique(data_test[[i]]$day))
}


# set up bandwidths and covariates
constX=0.305
constA=3.05
expo=8
outcome="icg"
covariates = c("meal", "glucose", "heart", "basal", "basal_4_8", 
               "meal_future_halfhour", "bolus_lag_halfhour")
decision_var = c("meal","meal_future_halfhour","basal_4_8","bolus_lag_halfhour")
d = length(decision_var)

#------------------------------------------------------
# run model for person 6
# for each lag, just change the weights
personID = 6
p6_model1 = runmodel(personID = personID,
                      outcome = outcome,
                      covariates = covariates,
                      LAG = 4,
                      interval = 30, 
                      bootstrap = 0,
                      weights = c(1/4, 1/4, 1/4, 1/4),
                      test_index = testday[[personID]],
                      decision_var = decision_var,
                      constX = constX,
                      constA = constA,
                      expo = expo)

#print the average weighted advantage of the suggested doses
mean(p6_model1$value$advantage_suggest)
#print the average weighted advantage of the original doses
mean(p6_model1$value$advantage_orig)

#plot the suggested dose
par(mfrow = c(1, 1))
plot(p6_model1$value$dosehat, type = "l", col = "red")
lines(p6_model1$value$dose, type = "l")
#plot the weighted advantage at each time point
par(mfrow = c(1, 1))
plot(p6_model1$value$advantage_orig, type = "l", ylab = "Weighted Outcome", xlab = "Timestamp", main = "Estimated Weighted Outcome With Suggested Doses in Test Data")
lines(p6_model1$value$advantage_suggest, type = "l", col = "red")
legend("bottomright",  col = c(1, 2), legend = c("Weighted Advantage by Original Dose", "Weighted Advantage by Suggested Dose"), lty = 1)


#---------------------------------------------------
#printing result for all persons into latex
LAGMAX = 4
for (i in 1:nperson){
  #cat("\n", dataset, "\n")
  beta_final = matrix(NA, ncol = d + 2,  nrow = LAGMAX + 1)
  sd_final = matrix(NA, ncol = d + 2, nrow = LAGMAX + 1)
  for (LAG in 1:LAGMAX){
    output_now = invisible(runmodel(personID = i, 
                                    outcome = outcome, 
                                    covariates = covariates, 
                                    LAG = LAG, 
                                    interval = 30, 
                                    bootstrap = 0, 
                                    test_index = testday[[i]], 
                                    decision_var = decision_var, 
                                    constX = constX, constA = constA, expo = expo))
    beta_final[LAG, ] = round(output_now$result$beta * 100, digits = 1)
    sd_final[LAG, ] = round(output_now$result$estimate_sd * 100, digits = 1)
  }
  output_now = invisible(runmodel(personID = i, 
                                  outcome = outcome, 
                                  covariates = covariates, 
                                  LAG = 4, 
                                  interval = 30, 
                                  bootstrap = 0, 
                                  weights = c(1/4, 1/4, 1/4, 1/4), 
                                  test_index = testday[[i]], 
                                  decision_var = decision_var, 
                                  constX = constX, constA = constA, expo = expo))
  beta_final[LAG + 1, ] = round(output_now$result$beta * 100, digits = 1)
  sd_final[LAG + 1, ] = round(output_now$result$estimate_sd * 100, digits = 1)
  for (j in 1:(d + 2)){
    if (j == 1){
      cat("$\\alpha_k$")
    }else{
      cat("$\\beta_{k, ", j-2, "}$", sep = "")
    }
    for (k in 1:(LAG + 1)){
      cat("&",  beta_final[k, j], "(", sd_final[k, j], ")", sep = "")
    }
    cat("\\\\")
    cat("\n")
  }
}

