##' Set up conf
##' @param years years to include
##' @param maxAge max age
##' @param spatioTemporal include spatio-tempora: 0-no 1-yes 2- yes but without correlation structure in time
##' @param cutoff mesh spesific, min distance between points
##' @param cbound mesh spesific, size of boundary meshSimilar
##' @param meshSimilar If true we apply the same mesh as in the package spatitemporalIndices
##' @param zone Zone used when converting to UTM coordinates
##' @param readability If 1: Utilize age reading quality, if 0: Do not utilize age reading quality
##' @details defines the configurations
##' @return Configurations to set up the model
##' @export
defConf_alk = function(years= NULL,maxAge, spatioTemporal = 1,rwBeta0 = 1, 
                       cutoff = 100,cbound = 200,meshSimilar = FALSE,zone = NULL,
                       readability = 1){
  conf = list()
  conf$maxAge = maxAge
  if(!is.null(zone)){
    conf$zone = zone
  }else{
    conf$zone = 35 #TODO make more dynamic
  }
  conf$spatioTemporal = spatioTemporal
  conf$rwBeta0 = rwBeta0
  conf$cutoff =cutoff
  conf$cbound = cbound
  conf$meshSimilar = meshSimilar
  conf$years = years
  conf$readability = readability
  return(conf)
}



##' Set up map variable
##' @param conf Configurations
##' @param par Parameters
##' @details Sets up the map-variable
##' @return Returns the map-variable
##' @export
setMap_alk = function(conf,par){
  map = list()
  if(conf$spatioTemporal==2){#No use of ar1 structure in time
    map$transRho_alk = as.factor(NA) 
  }else if(conf$spatioTemporal==0){#Turn off spatio-temporal structures
    map$transRho_alk = as.factor(NA) 
    map$xST_alk = as.factor(rep(NA,length(par$xST)))
    map$logKappa_alk = as.factor(NA)
    map$logSigma_alk = as.factor(NA) 
  }
  
  if(conf$rwBeta0==0){
    map$log_sigma_beta0_alk = as.factor(NA)
  }
  return(map)
}
