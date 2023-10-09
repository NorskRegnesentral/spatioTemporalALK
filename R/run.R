##' Set up data
##' @param data data
##' @param conf  configurations
##' @details define parameters
##' @return Prameters to be provided to TMB
##' @useDynLib spatioTemporalALK
##' @export
fitALK = function(data,par, conf){
  
  map = setMap_alk(conf,par)
  
  random = c("xST_alk","xS_alk")
  profile = numeric(0)
  if(conf$rwBeta0==1){
    random = c(random,"beta0_alk")
    if(conf$betaLength == 2){
      random = c(random,"betaLength_alk")
    }
  }else{
    profile = c("beta0_alk","betaLength_alk")
  }
  if(conf$spatioTemporal==0 & conf$spatial ==0 & conf$rwBeta0==0){
    profile = c("beta0_alk") #Need one parameter that is not profiled
  }
  
  if(length(profile)>0){
    obj = MakeADFun(data,par,random = random,profile = profile,DLL = "spatioTemporalALK", map = map)
  }else{
    obj = MakeADFun(data,par,random = random,DLL = "spatioTemporalALK", map = map)
  }
  opt = nlminb(obj$par,obj$fn,obj$gr, control = list(trace = 1, iter.max = 1000, eval.max =1000))
  
  rep = obj$report()
  sdrep = sdreport(obj)
  pl = as.list(sdrep, "Est")
  plsd = as.list(sdrep, "Std")
  
  run = list(obj = obj,opt = opt, data = data, rep = rep, sdrep = sdrep, pl = pl, plsd = plsd, conf = conf)

  class(run) <- "alk"
  return(run)
}
