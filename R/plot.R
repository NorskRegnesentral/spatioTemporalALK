##' Plot ALK 
##' @param run object of class stALK
##' @param year  Year of interest
##' @param length if != NULL, plot spatial ALK based on length and year 
##' @details Function to plot ALK. 
##' @return 
##' @export
plotALK = function(run, year, length = NULL, age = NULL,...){
  
  if(is.null(run$conf)){
    run$conf = run$conf_alk #Called conf_alk in spatioTemporalIndices
  }
  years = run$conf$years
  
  yearIndex = which(years == year)-1
  nAges = run$data$ageRange[2]- run$data$ageRange[1] + 1
  
  if(is.null(length)){
    lengthInt = seq(min(run$data$length), max(run$data$length), length.out = 200)
    
    linPredMatrix = matrix(0,length(lengthInt), nAges-1)
    for(a in 1:(nAges-1)){
      linPredMatrix[,a] = 
#        run$pl$beta0[yearIndex*(nAges-1) + a] + run$pl$betaLength[a]*lengthInt
      run$pl$beta0_alk[yearIndex+1, a] + run$pl$betaLength_alk[a]*lengthInt
    }
    
    ALK = matrix(0,length(lengthInt), nAges)
    for(a in 1:nAges){
      probLess = rep(0,dim(ALK)[1])
      if(a>1){
        for(b in 1:(a-1)){
          tmp = ALK[,b];
          probLess = probLess + tmp;
        }
      }
      if(a <(nAges)){
        tmp2 = linPredMatrix[,a];
        ALK[,a] = plogis(tmp2)*(1-probLess);
      }else{
        ALK[,nAges] = (1-probLess);
      }
    }
    
    plot(lengthInt,ALK[,1], type = "l", ylab = "Prob", xlab = "Length",ylim = c(0,1), main = paste0("ALK year ",year),...)
    for(l in 2:nAges){
      lines(lengthInt,ALK[,l],col = l,...)
    }
  }else{
    linPredMatrix = matrix(0,dim(run$pl$xST_alk)[1], nAges-1)
    scaleST = 1/((4*3.14159265)*exp(run$pl$logKappa_alk*2)); 
    
    for(a in 1:(nAges-1)){
      linPredMatrix[,a] = 
        run$pl$beta0_alk[yearIndex+1, + a] + run$pl$betaLength_alk[a]*length +   run$pl$xST_alk[,yearIndex+1,a]/sqrt(scaleST) *exp(run$pl$logSigma_alk)
#       run$pl$beta0[yearIndex*(nAges-1) + a] + run$pl$betaLength[a]*length +   run$pl$xST[,yearIndex+1,a]/sqrt(scaleST) *exp(run$pl$logSigma)
    }
    
    ALK = matrix(0,dim(run$pl$xST_alk)[1], nAges)
    for(a in 1:nAges){
      probLess = rep(0,dim(ALK)[1])
      if(a>1){
        for(b in 1:(a-1)){
          tmp = ALK[,b];
          probLess = probLess + tmp;
        }
      }
      if(a <(nAges)){
        tmp2 = linPredMatrix[,a];
        ALK[,a] = plogis(tmp2)*(1-probLess);
      }else{
        ALK[,nAges] = (1-probLess);
      }
    }
    
    
    #Spatial effect
  
    mesh = attributes(data)$mesh
    proj = inla.mesh.projector(mesh)
    xlim = c(min(mesh$loc[,1]),max(mesh$loc[,1]))
    ylim = c(min(mesh$loc[,2]),max(mesh$loc[,2]))
    latentFieldMAP = rowSums(ALK[,1:age])
    fields::image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","yellow", "red"))(12),
                       xlab = 'Easting', ylab = 'Northing',
                       main = paste0("Year ", year,": Prob of age <", age+1 , " given length ", length, "cm"),
                       cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1,
                       xlim = xlim,
                       ylim = ylim,
                       zlim = c(0,1),...)
    contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldMAP) ,add = T,labcex  = 1,cex = 1,levels = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95),...)
    
    newmap <- maps::map("world", c("Norway","Sweden","Finland","Russia"),fill = TRUE,plot = FALSE, col = "transparent")
    mapTmp = data.frame(newmap$x,newmap$y)
    mapTmp[which(is.na(mapTmp[,1])),] = 3.141592 #Need something different from NA
    names(mapTmp) = c("X","Y")
    attr(mapTmp, "projection") = "LL"
    attr(mapTmp, "zone") = 35
    ddpcr::quiet(mapTmp <- PBSmapping::convUL(mapTmp))
    colnames(mapTmp) = c("UTMX", "UTMY")
    mapTmp[which(is.na(newmap$x)),] = NA
    polygon(mapTmp,col = 'lightgrey')
  }
}

