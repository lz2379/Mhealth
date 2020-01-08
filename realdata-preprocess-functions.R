

#------------------------------------------------------------------------------
# 1. Functions for data preprocessing with Ohio Type 1 Diabetes dataset
#------------------------------------------------------------------------------


library(XML) #package for reading XML data
library(lubridate) #package for dealing with date and time

#Functions for measuring the effect of blood glucose control
ICG = function(x, a = 1.35, b = 2){
  return( - (x < 80) * ((abs(80 - x))^b) / 30 - (x > 140) * ((abs(x - 140))^a) / 30)
}

HyperIndex = function(x, a = 1.35){
  return( - (x > 140) * ((abs(x - 140))^a) / 30)
}

HypoIndex = function(x, b = 2){
  return( - (x < 80) * ((abs(80 - x))^b) / 30)
}
M100 = function(x){
  return( - 1000 * (abs(log10(x / 100)))^3)
}
InRange = function(x){
  return(as.numeric((x <= 140)&(x >= 70)))
}
#transform hour and minutes into minutes in a day
daytime = function(x){
  return(60 * hour(x) + minute(x))
}

#Put data into a list
preprocess0 = function(data0){
  data = list(glucose = data0[[1]], 
            basal = data0[[3]],  
            tempbasal = data0[[4]],  
            bolus = data0[[5]][, c(1, 4)], 
            meal = data0[[6]], 
            heart = data0[[13]], 
            gsr = data0[[14]], 
            skintemp = data0[[15]], 
            airtemp = data0[[16]], 
            steps = data0[[17]])
  return(data)
}

preprocess1 = function(data0, interval_min = 60,  startdate = NA){
  #start date is the first date of subject's participation in the experiment
  if (is.na(startdate)){
    startdate = date(dmy_hms(data0$glucose[1, 1]))
  }
  #tt is the time index in the day
  #if interval_min=60 than all the time between 0:00-1:00 will have tt=0, all the time between 1:00-2:00 have tt=1, etc. 
  data = list(glucose  = data.frame( day = as.numeric(date(dmy_hms(data0$glucose[, 1])) - startdate),  
                                    tt = floor(daytime(dmy_hms(data0$glucose[, 1])) / interval_min), 
                                    daytime = daytime(dmy_hms(data0$glucose[, 1])), 
                                    glucose = as.numeric(data0$glucose[, 2])) , 
              basal    = data.frame(day = as.numeric(date(dmy_hms(data0$basal[, 1])) - startdate),  
                                    daytime = daytime(dmy_hms(data0$basal[, 1])), 
                                    rate = as.numeric(data0$basal[, 2])), 
              bolus    = data.frame(day = as.numeric(date(dmy_hms(data0$bolus[, 1])) - startdate), 
                                    daytime = daytime(dmy_hms(data0$bolus[, 1])), 
                                    tt = floor(daytime(dmy_hms(data0$bolus[, 1])) / interval_min), 
                                    bolus = as.numeric(data0$bolus[, 2])),  
              meal     = data.frame(day = as.numeric(date(dmy_hms(data0$meal[, 1])) - startdate),  
                                    daytime = daytime(dmy_hms(data0$meal[, 1])), 
                                    tt = floor(daytime(dmy_hms(data0$meal[, 1])) / interval_min), 
                                    meal = as.numeric(data0$meal[, 3])), 
              heart    = data.frame(day = as.numeric(date(dmy_hms(data0$heart[, 1])) - startdate),  
                                    tt = floor(daytime(dmy_hms(data0$heart[, 1])) / interval_min), 
                                    heart = as.numeric(data0$heart[, 2])), 
              gsr      = data.frame(day = as.numeric(date(dmy_hms(data0$gsr[, 1])) - startdate), 
                                    tt = floor(daytime(dmy_hms(data0$gsr[, 1])) / interval_min), 
                                    gsr = as.numeric(data0$gsr[, 2])), 
              skintemp = data.frame(day = as.numeric(date(dmy_hms(data0$skintemp[, 1])) - startdate), 
                                    tt = floor(daytime(dmy_hms(data0$skintemp[, 1])) / interval_min), 
                                    skintemp = as.numeric(data0$skintemp[, 2])), 
              airtemp  = data.frame(day = as.numeric(date(dmy_hms(data0$airtemp[, 1])) - startdate), 
                                    tt = floor(daytime(dmy_hms(data0$airtemp[, 1])) / interval_min), 
                                    airtemp = as.numeric(data0$airtemp[, 2])), 
              steps    = data.frame(day = as.numeric(date(dmy_hms(data0$steps[, 1])) - startdate), 
                                    tt = floor(daytime(dmy_hms(data0$steps[, 1])) / interval_min), 
                                    steps = as.numeric(data0$steps[, 2]))
              )
  if (length(data0$tempbasal) > 0){
    data[[10]] = data.frame(start_day = as.numeric(date(dmy_hms(data0$tempbasal[, 1])) - startdate), 
                            start_daytime = daytime(dmy_hms(data0$tempbasal[, 1])), 
                            end_day = as.numeric(date(dmy_hms(data0$tempbasal[, 2])) - startdate), 
                            end_daytime = daytime(dmy_hms(data0$tempbasal[, 2])), 
                            rate = as.numeric(data0$tempbasal[, 3])  )
    names(data)[10] = "tempbasal"
  }
  return(data)
}

#turn basal rate change information into basal rate for each interval
preprocess2_basal = function(data,  interval_min = 60, startbasal = NA){
  data1 = data
  firstday = min(data$glucose$day)
  lastday = max(data$glucose$day)
  nday = lastday - firstday + 1
  
  #initialize
  ntt = 1440 / interval_min
  datett = data.frame(day = rep(seq(firstday, lastday), each = ntt), 
                      tt  = rep(seq(0, ntt - 1), nday))
  tbasal = datett
  tbasal$rate = NA
  ob = 1
  trate = startbasal
  data1$basal = rbind(data1$basal,  c(lastday + 1, 0, 0)) #set an end point
  
  for (i in 1:(ntt * nday)){
    if ((data1$basal$day[ob] - firstday) * ntt + floor(data1$basal$daytime[ob] / interval_min) > i - 1) { 
      #basal rate not changing at this time interval
      tbasal$rate[i] = trate
    } else{
      #basal rate changing at this time interval
      begin_interval = 0
      tbasal$rate[i] = 0
      while ((data1$basal$day[ob] - firstday) * ntt + floor(data1$basal$daytime[ob] / interval_min) <= i - 1) {
        tbasal$rate[i] = tbasal$rate[i] +  trate *  (data1$basal$daytime[ob] %% interval_min - begin_interval) / interval_min
        begin_interval = data1$basal$daytime[ob] %% interval_min
        trate = data1$basal$rate[ob]
        ob = ob + 1
      }
      tbasal$rate[i] = tbasal$rate[i] +  trate *  (interval_min - begin_interval) / interval_min
    }
  }
  
  #adjust with temp basal information
  if (length(which(names(data1) == "tempbasal")) > 0){
    for (i in 1:nrow(data1$tempbasal)){
      begin_int_whole  = (data1$tempbasal$start_day[i] - firstday) * ntt + 
                         floor(data1$tempbasal$start_daytime[i] / interval_min) + 1
      end_int_whole    = (data1$tempbasal$end_day[i] - firstday) * ntt + 
                         floor(data1$tempbasal$end_daytime[i] / interval_min) + 1
      begin_int_rest   = data1$tempbasal$start_daytime[i] %% interval_min
      end_int_rest     = data1$tempbasal$end_daytime[i] %% interval_min
      trate = data1$tempbasal$rate[i]
      while (begin_int_whole  <  end_int_whole){
        tbasal[begin_int_whole, ]$rate = tbasal[begin_int_whole, ]$rate * (begin_int_rest) / interval_min + 
          trate * (interval_min - begin_int_rest) / interval_min
        begin_int_whole = begin_int_whole + 1
        begin_int_rest = 0
      }
      tbasal[ begin_int_whole, ]$rate = tbasal[ begin_int_whole, ]$rate * (begin_int_rest + interval_min - end_int_rest) / interval_min  + 
                                        trate * (end_int_rest - begin_int_rest) / interval_min
    }
  }
  names(tbasal)[which(names(tbasal) == "rate")] = "basal"
  return(tbasal)
}


#
#aggregated variables into time intervals
#Each variable has its own function for aggregating (sum or mean)
preprocess2_aggregate = function(data,  
                                 variables = c("meal", "bolus",  "glucose", "heart", "gsr", "skintemp", "airtemp", "steps"), 
                                 FUN = list(sum,  sum, mean, mean, mean, mean, mean, sum), 
                                 glucose_eval = list(ICG, HyperIndex, HypoIndex, M100, InRange), 
                                 eval_names = c("icg", "hyper", "hypo", "m100", "in_range"), 
                                 missing_as_zero = c("meal", "bolus"), 
                                 interval_min = 60){
  firstday = min(data$glucose$day)
  lastday = max(data$glucose$day)
  nday = lastday - firstday + 1
  ntt = 1440 / interval_min
  datett = data.frame(day = rep(seq(firstday, lastday), each = ntt), 
                    tt  = rep(seq(0, ntt - 1), nday))
  #find the index of the variables
  var_index = as.numeric(sapply(variables, function(x) which(names(data) == x)))
  data1 = datett
  for (i in 1:length(variables)){
    temp = aggregate(as.formula(paste(variables[i], "~day + tt", sep = "")), 
                   data = data[[var_index[i]]], 
                   FUN = FUN[[i]], 
                   na.rm = TRUE)
    data1 = merge(data1, temp, by = c("day", "tt"), all = TRUE)
  }
  for (i in 1:length(missing_as_zero)){
    m_index =  which(names(data1) == missing_as_zero[i])
    data1[is.na(data1[, m_index]), m_index] = 0
  }
  #transformations of glucose
  glucoseindex = which(names(data) == "glucose")
  gl = data[[glucoseindex]]$glucose
  for (i in 1:length(glucose_eval)){
    data[[glucoseindex]] = cbind(data[[glucoseindex]],  glucose_eval[[i]](gl))
    names(data[[glucoseindex]])[dim(data[[glucoseindex]])[2]] = eval_names[i]
    temp = aggregate(as.formula(paste(eval_names[i], "~day + tt", sep = "")), 
                   data = data[[glucoseindex]], 
                   FUN = mean, 
                   na.rm = TRUE)
    data1 = merge(data1, temp, by = c("day", "tt"), all = TRUE)
  }
  return(data1)
}

#calculate the last observed glucose level within the time interval (vs average glucose level)
preprocess2_lastglucose = function(data, interval_min = 30){
  tglucose = data$glucose
  C = 100000
  tglucose$daytime_glucose = tglucose$daytime * C + tglucose$glucose
  tglucose2 = aggregate(daytime_glucose~tt + day, data = tglucose, FUN = max)
  tglucose2$lastglucose = tglucose2$daytime_glucose %% C
  tglucose2$daytime_glucose = NULL
  return(tglucose2)
}

preprocess2 = function(data, interval_min = 60, startbasal = NA){
  tbasal = preprocess2_basal(data, 
                           interval_min = interval_min, 
                           startbasal = startbasal)
  temp = preprocess2_aggregate(data, interval_min = interval_min)
  lastglucose = preprocess2_lastglucose(data, interval_min = interval_min)
  data1 = merge(tbasal, temp, by = c("day", "tt"), all = TRUE)
  data1 = merge(data1, lastglucose, by = c("day", "tt"), all = TRUE)
  data1 = data1[order(data1$day, data1$tt), ]
  return(data1)
}

#calculate past variable information
preprocess3_past = function(data,  variable = "basal",  interval_min = 60,  
                            lag_time_range = c(120, 240), 
                            new_var_name = "basal_2_4", 
                            FUN = mean){
  lag_time = floor(lag_time_range / interval_min)
  lagmin = lag_time[1]
  lagmax = lag_time[2] - 1
  data1 = data
  var_index = which(names(data) == variable)
  #data1$newvar = 0
  index = (1 + lagmax):(dim(data)[1])
  temp = matrix(NA, nrow = length(index), ncol = 0)
  for (i in lagmin:lagmax){
    temp = cbind(temp, data[index - i, var_index])
  }
  data1$newvar = NA
  data1$newvar[index] = apply(temp, 1, FUN = FUN)
  names(data1)[which(names(data1) == "newvar")] = new_var_name
  return(data1)
}

#calculate future variable information
preprocess3_future = function(data,  variable = "meal",  interval_min = 30,  
                          lag_time_range = c(0, 30), 
                          new_var_name = "meal_future_halfhour", 
                          FUN = sum){
  lag_time = floor(lag_time_range / interval_min)
  lagmin = lag_time[1] + 1
  lagmax = lag_time[2]
  data1 = data
  var_index = which(names(data) == variable)
  #data1$newvar = 0
  index = (1):(dim(data)[1] - lagmax)
  temp = matrix(NA, nrow = length(index), ncol = 0)
  for (i in lagmin:lagmax){
    temp = cbind(temp, data[index + i, var_index])
  }
  data1$newvar = NA
  data1$newvar[index] = apply(temp, 1, FUN = FUN)
  names(data1)[which(names(data1) == "newvar")] = new_var_name
  return(data1)
}





