##' Print alk object 
##' @method print alk 
##' @param  x the fitted object as returned from the fitALK function
##' @details prints the log-likelihood and the main convergence criteria
##' @export
print.alk<-function(x, ...){
  cat("ALK model: log likelihood is", logLik.alk(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

##' Log likelihood of stALK object 
##' @method logLik alk 
##' @param  object fitted object as returned from the fitALK function
##' @details log likelihood of a ALK model run 
##' @export
logLik.alk<-function(object, ...){
  ret<- -object$opt$objective 
  
  beta0 = 0
  if(object$conf$rwBeta0==0){
    beta0 = beta0+ length(object$conf$minAge:object$conf$maxAge )*length(object$conf$years)
  }else{
    beta0 = beta0+ length(object$conf$minAge:object$conf$maxAge) + length(object$conf$years) -1
  }
  
  attr(ret,"df")<-length(object$opt$par) + length(object$pl$betaLength_alk) + beta0
  class(ret)<-"logLik"
  ret
}
