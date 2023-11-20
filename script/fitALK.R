library(spatioTemporalALK)
rm(list=ls())

#Read raw data
dataAge = readRDS("data/catch_at_age_data_ex_rus.rds")
#dataAge = dataAge[order(dataAge$startdatetime,dataAge$age),]

#Define configurations
conf = defConf_alk(years = 2018:2020,
                   maxAge = 10,
                   minAge = 3,
                   spatioTemporal = 0,
                   spatial = 0,
                   betaLength = 1,
                   cutoff =80, cbound = 130, 
                   rwBeta0 = 1,
                   readability = 0)

#Set up data
data = setUpData_alk(dataAge,conf)

#Define parameters
par = defpar_alk(data,conf)

#Fit model
startTime = Sys.time()
run = fitALK(data,par,conf)
endTime = Sys.time()
endTime-startTime

plotALK(run,year = 2019)

#OSA residuals
conditional = 1:(length(run$data$age)-50)
ageRange = seq(min(run$data$age), max(run$data$age))
res <- oneStepPredict(run$obj, observation.name ="age",
                      data.term.indicator="keep",
                      discrete = TRUE,
                      discreteSupport = ageRange,
                      method = "cdf",
                      conditional = conditional)
plot(res$residual)




#Example of plot of ALK wihthout spatio-temporal effect
x11(width = 20, height = 15)
par(oma=c(3,3,2,0.5),mar=c(2,2,2,0.5),mfrow = c(5,6))
for(year in 1994:2020){
  plotALK(run,year = year)
}


#Example of plot of ALK with spatial effect
x11()
year = 2019
age = 4
par(mfrow = c(3,2))
for(length in c(30,40,45,50,55,60)){
  plotALK(run,year = year,length = length,age = age)
}







#Go a newton step
he <- function(par){ optimHess(par, run$obj$fn, run$obj$gr) }
newtonsteps = 2
for(i in seq_len(newtonsteps)) { 
  g <- as.numeric( run$obj$gr(run$opt$par) )
  h <- stats::optimHess(run$opt$par, run$obj$fn, run$obj$gr)
  run$opt$par <- run$opt$par- solve(h, g)
  run$opt$objective <- run$obj$fn(run$opt$par)
}
run$opt$he <- optimHess(run$opt$par, run$obj$fn, run$obj$gr)





