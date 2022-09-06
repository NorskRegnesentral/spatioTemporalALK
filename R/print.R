##' Print partabld
##' @param  run the fitted object as returned from the fitALK function
##' @details Print ...
##' @return Partable..
##' @export
partable<-function(run){
  pl = run$pl
  plsd = run$plsd
  
  sigmaS = sqrt(c(exp(pl$logSigma), exp(pl$logSigma -1.96*plsd$logSigma), exp(pl$logSigma +1.96*plsd$logSigma) ))
  spatialS = round(c(sqrt(8)/exp(pl$logKappa), sqrt(8)/exp(pl$logKappa +1.96*plsd$logKappa) , sqrt(8)/exp(pl$logKappa -1.96*plsd$logKappa)),1 )
  rhoT = c(2*plogis(pl$transRho)-1, 2*plogis(pl$transRho- 1.96*plsd$transRho)-1, 2*plogis(pl$transRho+ 1.96*plsd$transRho)-1)

  parTab = round(matrix(c(sigmaS,spatialS,rhoT),ncol = 3, byrow = TRUE),3)
  
  parTab = as.data.frame(parTab)
  rownames(parTab) = c("Variance space", "Spatial range", "Rho time")
  colnames(parTab) = c("Mode", "0.025Q", "0.975Q")
  
  return(parTab)
}
