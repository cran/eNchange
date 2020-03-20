#' A bootstrap method to calculate the threshold (stopping rule) in the BS or EBS segmentation.
#' @references Cho, Haeran, and Karolos Korkas. "High-dimensional GARCH process segmentation with an application to Value-at-Risk." arXiv preprint <arXiv:1706.01155> (2018).
#' @rdname boot_thresh-methods
#' @description A bootstrap method to calculate the threshold (stopping rule) in the BS or EBS segmentation described in Cho and Korkas (2018) and adapted for irregularly time series in Korkas (2020).
#' @param H The input irregular time series.
#' @param q The bootstrap distribution quantile. Default is 0.75.
#' @param r The number of bootrstap simulations. Default is 100.
#' @param p The support of the CUSUM statistic. Default is 1.
#' @param start.values Warm starts for the optimizers of the likelihood functions.
#' @param process Choose between acd or hawkes. Default is acd.
#' @param do.parallel Choose the number of cores for parallel computation. If 0 no parallelism is done. Default is 2.
#' @param dampen.factor The dampen factor in the denominator of the residual process. Default is "auto".
#' @param epsilon A parameter added to ensure the boundness of the residual process. Default is  1e-5.
#' @param LOG Take the log of the residual process. Default is TRUE.
#' @param acd_p The p order of the ACD model. Default is 0.
#' @param acd_q The q order of the ACD model. Default is 1.
#' @examples
#' pw.acd.obj <- new("simACD")
#' pw.acd.obj@cp.loc <- c(0.25,0.75)
#' pw.acd.obj@lambda_0 <- c(1,2,1)
#' pw.acd.obj@alpha <- rep(0.2,3)
#' pw.acd.obj@beta <- rep(0.7,3)
#' pw.acd.obj@N <- 3000
#' pw.acd.obj <- pc_acdsim(pw.acd.obj)
#' boot_thresh(pw.acd.obj@x,r=20)
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats rnorm rgeom runif optim quantile
#' @importFrom ACDm acdFit
#' @importFrom  hawkes simulateHawkes
#' @useDynLib eNchange, .registration = TRUE
#' @export
#' @return Returns the threshold \code{C}.
#' @docType methods
#' @aliases boot_thresh boot_thresh-class boot_thresh-methods
setGeneric(name="boot_thresh",
           def=function(H,q=.75,r=100,p=1,start.values=c(0.9,0.6),process="acd",do.parallel=2,dampen.factor="auto",epsilon= 1e-5,LOG=TRUE,acd_p=0,acd_q=1)
           {
             standardGeneric("boot_thresh")
           }
)
#' @rdname boot_thresh-methods
setMethod(f="boot_thresh", definition = function(H,q=.75,r=100,p=1,start.values=c(0.9,0.6),process="acd",do.parallel=2,dampen.factor="auto",epsilon= 1e-5,LOG=TRUE,acd_p=0,acd_q=1) {
  if (do.parallel>0) {
    cl=makeCluster(do.parallel)
    registerDoParallel(cl)
  }
  `%parDo%` <- ifelse(do.parallel > 0, `%dopar%`, `%do%`)
  if (process=="hawkes") {
    estpar = optim(par=start.values,hawkesLike,t_i=cumsum(H),lower=c(0.001,0.001),method="L-BFGS-B")$par
    if (estpar[2]>estpar[3]) estpar[2]=.9*estpar[3]
    n=ceiling(tail(cumsum(H),1))
    X=list()
    M=matrix(0,2,r)
    i=1
    while(i != (r+1)){
      d=diff(simulateHawkes(lambda0=estpar[1],alpha=estpar[2],beta=estpar[3],horizon=n)[[1]])
      if ((length(d) > .95*length(H)) && (length(d) < 1.05*length(H))) {
        X[[i]]=d
        i=i+1
      }
    }
    for (i in 1:r){
      x=X[[i]]
      z=Z_trans(H=x,start.values = start.values,dampen.factor = dampen.factor,epsilon = epsilon,LOG=LOG,process = "hawkes")
      M[1,i] = abs(finner_prod_maxp(cusum(z),p=p)[[1]]) * (log(length(z))^(-1))
      M[2,i] = length(z)
    }
    return (quantile(M[1,],probs=q))
  } else if (process=="acd") {
    if (!is.null(start.values) & length(start.values) != sum(1,acd_p,acd_q)) stop("starting values do not match the ACD order")
    ACD.obj = tryCatch({acdFit(H,startPara=start.values,order=c(acd_p,acd_q),output=FALSE)},
                       error=function(e) {message("...Optimizer not converging: using start.values");
                         acd.est(H,lambda_0 = start.values[1],alpha=switch (acd_p ==0,NULL,start.values[2:(acd_p+1)]),beta=switch(acd_q==0,NULL,start.values[(2+acd_p):(2+acd_p+acd_q)])
                         )})
    
    if (is.null(ACD.obj$convergence)) ACD.obj = acd.est(H,lambda_0 = start.values[1],alpha=switch (acd_p ==0,NULL,start.values[2:(acd_p+1)]),beta=switch(acd_q==0,NULL,start.values[(2+acd_p):(2+acd_p+acd_q)]))
    boots.resids=ACD.obj$residuals
    estpar=ACD.obj$mPara
    n=length(H) 
    
    M <- foreach(i=iter(1:r),.combine = c,.packages = c('Rcpp','ACDm','eNchange')) %parDo% {
      x=sim.ACD(N=n,lambda_0 = estpar['omega'],alpha = ifelse(is.na(estpar['alpha1']),0,estpar['alpha1']),beta = ifelse(is.na(estpar['beta1']),0,estpar['beta1']),resids=boots.resids)
      z=Z_trans(H=x[[1]],process = "acd",start.values = start.values,dampen.factor = dampen.factor,epsilon = epsilon,LOG=LOG,acd_p = acd_p,acd_q = acd_q)
      abs(finner_prod_maxp(cusum(z),p=p)[[1]])*(log(length(z))^(-1))
    }
  }
  if (do.parallel) stopCluster(cl)
  return(quantile(M,probs = q))
})