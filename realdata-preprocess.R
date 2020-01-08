#-----------------------------------------------
# Data proproecssing
#-----------------------------------------------
source("realdata-preprocess-functions.R")

#---------------------------------
# 1. Extract data from xml

data_path = "raw_data/"

#training data
raw_train = list()
file_names = c("559-ws-training.xml", 
               "563-ws-training.xml", 
               "570-ws-training.xml", 
               "575-ws-training.xml", 
               "588-ws-training.xml", 
               "591-ws-training.xml")
nperson = length(file_names)
for (n in 1:nperson){
  data = xmlParse(paste(data_path, file_names[n],sep=""))
  rootnode = xmlRoot(data)
  variables = xmlChildren(rootnode)
  nvar = xmlSize(rootnode)
  raw_train[[n]] = list()
  varnames = rep("",nvar)
  for (i in 1:nvar){
    vari = variables[[i]]
    varnames[i] = xmlName(variables[[i]])
    data_var = t(as.matrix(xmlSApply(vari, xmlAttrs)))
    raw_train[[n]][[i]] = data_var
  }
}

#test data
raw_test = list()
file_names = c("559-ws-testing.xml", 
               "563-ws-testing.xml", 
               "570-ws-testing.xml", 
               "575-ws-testing.xml", 
               "588-ws-testing.xml", 
               "591-ws-testing.xml")
nperson = length(file_names)
for (n in 1:nperson){
  data = xmlParse(paste(data_path, file_names[n],sep=""))
  rootnode = xmlRoot(data)
  variables = xmlChildren(rootnode)
  nvar = xmlSize(rootnode)
  raw_test[[n]] = list()
  varnames = rep("",nvar)
  for (i in 1:nvar){
    vari = variables[[i]]
    varnames[i] = xmlName(variables[[i]])
    data_var = t(as.matrix(xmlSApply(vari, xmlAttrs)))
    raw_test[[n]][[i]] = data_var
  }
}

#-----------------------------------------------
# 2. Preprocess data

# Into 30 min intervals 
interval = 30 #separate by half an hour
data_train = list()
data_test  = list()
data_full  = list()
for (n in 1:nperson){
  #training data
  train0 = preprocess0(raw_train[[n]])
  train1 = preprocess1(train0, interval_min = interval)
  data_train[[n]] = preprocess2(train1, interval_min = interval)
  
  #testing data
  startdate = date(dmy_hms(train0$glucose[1, 1]))
  startbasal = data_train[[n]]$basal[length(data_train[[n]]$basal)]
  test0 = preprocess0(raw_test[[n]])
  test1 = preprocess1(test0, interval_min = interval, startdate = startdate)
  data_test[[n]] = preprocess2(test1, interval_min = interval, startbasal = startbasal)
  
  #combine training and testing  
  data_full[[n]] = rbind(data_train[[n]], data_test[[n]])
  
  #generate extra variables
  #average basal rate of 4-8 hours before the current time point
  data_full[[n]] = preprocess3_past(data_full[[n]], 
                                    variable = "basal",
                                    interval_min = interval,
                                    lag_time_range = c(240, 480),
                                    new_var_name = "basal_4_8")
  #total amount of bolus from last time interval 
  data_full[[n]] = preprocess3_past(data_full[[n]],
                                    variable = "bolus",
                                    interval_min = interval,
                                    lag_time_range = c(30, 60),
                                    new_var_name = "bolus_lag_halfhour",
                                    FUN = sum)
  #total amount of carb input in the next half an hour
  data_full[[n]] = preprocess3_future(data_full[[n]],
                                      variable = "meal",
                                      interval_min = interval,
                                      lag_time_range = c(0,30),
                                      new_var_name = "meal_future_halfhour",
                                      FUN = sum)
}
save.image("data_preprocessed-halfhour.Rdata")

#--------
#Another option is to preprocess into 5min intervals
interval = 5 #separate by half an hour
data_train = list()
data_test  = list()
data_full  = list()
for (n in 1:nperson){
  #training data
  train0 = preprocess0(raw_train[[n]])
  train1 = preprocess1(train0, interval_min = interval)
  data_train[[n]] = preprocess2(train1, interval_min = interval)
  
  #testing data
  startdate = date(dmy_hms(train0$glucose[1, 1]))
  startbasal = data_train[[n]]$basal[length(data_train[[n]]$basal)]
  test0 = preprocess0(raw_test[[n]])
  test1 = preprocess1(test0, interval_min = interval, startdate = startdate)
  data_test[[n]] = preprocess2(test1, interval_min = interval, startbasal = startbasal)
  
  #combine training and testing  
  data_full[[n]] = rbind(data_train[[n]], data_test[[n]])
}
save.image("data_preprocessed-5min.Rdata")


#--------------------------------------------
# Plot part of the 5min dataset
#--------------------------------------------

load("data_preprocessed-5min.Rdata")
i = 1
dataset = data_full[[i]]
startdate = 7
enddate = 13
nday = enddate - startdate + 1
ntt = 60 * 24 / 5
temp = subset(dataset, (day >= startdate) & (day <= enddate))
temp$timestamp = (temp$day - startdate) * ntt + temp$tt

#plot the glucose
par(mar=c(5, 4, 4, 6) + 0.1)
plot(temp$timestamp, temp$glucose, 
     ylim = c(0, 400),
     xlim = c(0, ntt * nday), type = "l", axes = FALSE, 
     xlab = "", ylab = "",
     main = paste("Mobile Health Data for Diabetes Patients", i, "in", nday, "days"))
labels = rep("", 2 * nday + 1)
labels[(1:nday) * 2] = paste("Day ", c(1 : nday) )
axis(1, at = (ntt / 2 * c(0:(nday * 2))), labels = labels, 
     col="black")
abline(h = 80, col = "green")
abline(h = 140, col = "green")
mtext("Timestamp with 5 Min Intervals", side=1, line=2.5)
#axis(2,ylim=c(0,400),col="black")
#mtext("Blood Glucose (mg/dL)",side=2,line=2.5)
box()
par(new=TRUE)

#plot the meals
temp2 = subset(temp, meal > 0)
plot(temp2$timestamp, temp2$meal,
     xlim = c(0, ntt * nday), pch = 20, 
     xlab = "", ylab = "",ylim = c(0, 200), 
     axes = FALSE, col = "red", col.axis="red")
arrows(temp2$timestamp, 0, temp2$timestamp, temp2$meal, length = 0,col = "red")
mtext("Carb Intake (g)", side = 4,col = "red",line = 2.5)
axis(4, ylim = c(0, 200), col = "red", col.ticks = "red")

#plot insulin dosage
par(new = TRUE)
temp3 = subset(temp, bolus > 0)
plot(temp3$timestamp, temp3$bolus, 
     xlim = c(0, ntt * nday), pch = 20,
     xlab = "", ylab = "", col = "blue", ylim = c(0, 20), axes = FALSE)
arrows(temp3$timestamp, 0, temp3$timestamp, temp3$bolus, length = 0, col = "blue")
axis(2, ylim = c(0, 20), col = "blue", col.ticks = "blue")
mtext("Insulin Dosage (Unit)", side = 2, col = "blue", line = 2.5)

legend("topleft", col = c("black","green"), 
       lty = c(1, 1), legend = c("Blood Glucose (mg/dL)", "Safe Blood Glucose Range"))


