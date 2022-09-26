##' Set up data
##' @param data data
##' @param conf  configurations
##' @details define parameters
##' @return Prameters to be provided to TMB
##' @useDynLib spatioTemporalALK
##' @export
fitALK = function(data,par, conf){
  
  map = setMap_alk(conf,par)
  
  if(conf$spatioTemporal!=0){
    if(conf$rwBeta0==0){
      obj = MakeADFun(data,par,random = c("xST_alk","xS_alk"),profile = c("beta0_alk","betaLength_alk"),DLL = "spatioTemporalALK", map = map)
    }else{
      obj = MakeADFun(data,par,random = c("xST_alk", "beta0_alk","xS_alk"),profile = c("betaLength_alk"),DLL = "spatioTemporalALK", map = map)
    }
  }else{#Quick-fix, needs minimum 1 parameter which is not in "random" or "profile"
    if(conf$rwBeta0==0){
      obj = MakeADFun(data,par,random = c("xS_alk"),profile = c("beta0_alk","betaLength_alk"),DLL = "spatioTemporalALK", map = map)
    }else{
      obj = MakeADFun(data,par,random = c("xS_alk","beta0_alk"),profile = c("betaLength_alk"),DLL = "spatioTemporalALK", map = map)
    }
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
