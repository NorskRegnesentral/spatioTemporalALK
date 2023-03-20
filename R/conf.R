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
defConf_alk = function(years= NULL,minAge,maxAge, spatioTemporal = 0,spatial = 0,rwBeta0 = 1, 
                       cutoff = 100,cbound = 200,meshSimilar = FALSE,zone = NULL,
                       readability = 1, 
                       usePCpriorsALK = 0, pcPriorsALKRange = c(300,0.1), pcPriorsALKSD = c(1,0.1)){
  conf = list()
  conf$minAge = minAge
  conf$maxAge = maxAge
  if(!is.null(zone)){
    conf$zone = zone
  }else{
    conf$zone = 36 #TODO make more dynamic
  }
  conf$rwBeta0 = rwBeta0
  conf$cutoff =cutoff
  conf$cbound = cbound
  conf$meshSimilar = meshSimilar
  conf$years = years
  conf$readability = readability
  conf$usePCpriorsALK = usePCpriorsALK
  conf$pcPriorsALKRange = pcPriorsALKRange
  conf$pcPriorsALKSD = pcPriorsALKSD
  conf$spatioTemporal = spatioTemporal
  conf$spatial = spatial
  
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
  map$logKappa_alk = c(1,2)
  map$logSigma_alk = c(1,2)
  
  if(conf$spatial==0){
    map$xS_alk = as.factor(rep(NA,length(par$xS_alk)))
    map$logKappa_alk[1] = NA
    map$logSigma_alk[1] = NA 
  }
  
  if(conf$spatioTemporal==2){#No use of ar1 structure in time
    map$transRho_alk = as.factor(NA) 
  }else if(conf$spatioTemporal==0){#Turn off spatio-temporal structures
    map$transRho_alk = as.factor(NA) 
    map$xST_alk = as.factor(rep(NA,length(par$xST_alk)))
    map$logKappa_alk[2] = NA
    map$logSigma_alk[2] = NA 
  }

  map$logKappa_alk = as.factor(map$logKappa_alk)
  map$logSigma_alk = as.factor(map$logSigma_alk)
  
  if(conf$rwBeta0==0){
    map$log_sigma_beta0_alk = as.factor(NA)
  }
  return(map)
}
