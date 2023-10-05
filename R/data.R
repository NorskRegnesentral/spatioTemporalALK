##' Set up data
##' @param df raw data
##' @param conf  configurations
##' @details Prepare data
##' @return Data to be provided to TMB
##' @export
setUpData_alk = function(d, conf_alk,conf_l = NULL){
  #This function will change a lot when data is read more dynamically
  d$year = as.integer(format(d$startdatetime, format = "%Y"))
  
  d = d[which(d$year %in%conf_alk$years),] 
  
  if(conf_alk$readability==0){#Do not utilize age reading quality
    d$readability[d$readability==5 | d$readability==6] = 1
  }
  
  d$ageNotTruncated = d$age
  d$age[d$age>conf_alk$maxAge] = conf_alk$maxAge
  
  if(length(d$age<conf_alk$minAge)>0){
    d$age[d$age<conf_alk$minAge] = conf_alk$minAge-1
  }
  d$age = d$age - min(d$age)+1 #First age is set to 1
  
  
  ageRange = c(min(d$age), max(d$age))
  nAge = ageRange[2]-ageRange[1] + 1

  
  uniqueYears = unique(d$year)
  idx1 = rep(0,length(uniqueYears)); idx2 = idx1; #Bookkeeping, observations within years
  for(y in 1:length(uniqueYears)){
    idx1[y] = min(which(d$year==uniqueYears[y]))-1
    idx2[y] = max(which(d$year==uniqueYears[y]))-1
  }
  
  
  #Convert to UTM coordinates-------------------------------------------------------------------------
  loc = data.frame(d$longitude,d$latitude)
  names(loc) = c("X","Y")
  attr(loc, "projection") = "LL"
  attr(loc, "zone") = conf_alk$zone
  locUTM <- PBSmapping::convUL(loc)
  colnames(locUTM) = c("UTMX", "UTMY")
  
  if(conf_alk$meshSimilar){
    mesh = spatioTemporalIndices::createMesh(conf_l)$mesh
  }else{
    boundary.loc <- SpatialPoints(as.matrix(locUTM))
    boundary <- list(
      inla.nonconvex.hull(coordinates(boundary.loc), 20,resolution = 100),
      inla.nonconvex.hull(coordinates(boundary.loc), conf_alk$cbound))
    
    max.edge = c(conf_alk$cutoff,4*conf_alk$cutoff)
    mesh <- inla.mesh.2d(boundary=boundary,
                         max.edge=max.edge,
                         cutoff=conf_alk$cutoff)
    plot(mesh)
    print(paste("Mesh points:",mesh$n))
  }

  
  spde = inla.spde2.matern(mesh, alpha=2)
  spdeMatricesST = spde$param.inla[c("M0","M1","M2")]
  
  A = inla.spde.make.A(mesh,as.matrix(locUTM))
  
  A_list =list()
  for(i in 1:length(uniqueYears)){
    A_list[[i]] = inla.spde.make.A(mesh,loc = as.matrix(locUTM[which(d$year == uniqueYears[i]),]))
  }
  #----------------------------------------------------------------------------------------------------
  
  data = list(age = d$age, 
              ageNotTruncated = d$ageNotTruncated,
              length = d$length,
              spdeMatricesST_alk = spdeMatricesST,
              A_alk_list = A_list,
              readability = d$readability,
              idx1 = idx1,
              idx2 = idx2,
              ageRange = ageRange,
              rwBeta0_alk = conf_alk$rwBeta0,
              maxAge = conf_alk$maxAge,
              minAge = conf_alk$minAge,
              usePCpriorsALK = conf_alk$usePCpriorsALK,
              pcPriorsALKRange = conf_alk$pcPriorsALKRange,
              pcPriorsALKSD = conf_alk$pcPriorsALKSD,
              spatioTemporalALK = conf_alk$spatioTemporal,
              spatialALK = conf_alk$spatial)
  
  attributes(data)$uniqueYears = uniqueYears 
  attributes(data)$loc = loc 
  attributes(data)$locUTM = locUTM 
  attributes(data)$years = d$year 
  attributes(data)$mesh = mesh 
  
  return(data)
}
