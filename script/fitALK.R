library(spatioTemporalALK)
rm(list=ls())

#Read raw data
dataAge = readRDS("catch_at_age_data_ex_rus.rds")

#TODO: do not care about distinguishing between aged less than min_age


#Define configurations
conf = defConf_alk(years = 1994:2020,maxAge = 10,spatioTemporal = 2,cutoff = 100, rwBeta0 = 1)

#Set up data
data = setUpData_alk(dataAge,conf)

#Define parameters
par = defpar_alk(data,conf)

#Fit model
startTime = Sys.time()
run = fitALK(data,par,conf)
endTime = Sys.time()
endTime-startTime


run

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

