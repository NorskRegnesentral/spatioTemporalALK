##' Set up data
##' @param data data
##' @param conf  configurations
##' @details define parameters
##' @return Prameters to be provided to TMB
##' @export
defpar_alk = function(data, conf){
  
  nAge = data$ageRange[2]-data$ageRange[1] + 1
  
  xS_alk = array(0.0, dim = c(dim(data$A_alk_list[[1]])[2],nAge-1))
  xST_alk = array(0.0, dim = c(dim(data$A_alk_list[[1]])[2],length(conf$years),nAge-1))
  
  par = list(beta0_alk = array(0,dim=c(length(conf$years),nAge-1)),
             log_sigma_beta_alk = c(0,0),
             betaLength_alk = rep(0,nAge-1),
             logSigma_alk = c(2,2),
             logKappa_alk = c(-4,-4),
             transRho_alk = 0,
             xS_alk = xS_alk,
             xST_alk = xST_alk)
  
  if(conf$betaLength == 2){
    par$betaLength_alk = rep(0,(nAge-1)*length(conf$years))
  }
  
  return(par)
}

