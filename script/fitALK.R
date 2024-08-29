library(spatioTemporalALK)
rm(list=ls())

#Read raw data
dataAge = readRDS("data/catch_at_age.rds")

#Define configurations
conf = defConf_alk(years = 2018:2020,
                   maxAge = 10,
                   minAge = 3,
                   spatioTemporal = 0,
                   spatial = 0,
                   betaLength = 1,
                   cutoff =40, cbound = 130, 
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

#Plot ALK with spatial effect set to zero
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






