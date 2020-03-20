#' Transformation of an irregularly spaces time series.
#' @references Korkas Karolos. "Ensemble Binary Segmentation for irregularly spaced data with change-points" Preprint <arXiv:2003.03649>.
#' @rdname Z_trans-methods
#' @description Transformation of a irregularly spaces time series. For the tvACD model, we calculate
#' \eqn{U_t = g_0(x_t, \psi_t) = \frac{x_t}{{\psi}_t}}, where
#' \eqn{{\psi}_t = C_0 + \sum_{j=1}^p C_j x_{t-j} + \sum_{k=1}^q C_{p+k} \psi_{t-k}+\epsilon x_t}.
#' where the last term \eqn{\epsilon x_t} is added to ensure the boundness of \eqn{U_t}.
#' @param H The input irregular time series.
#' @param start.values Warm starts for the optimizers of the likelihood functions.
#' @param dampen.factor The dampen factor in the denominator of the residual process. Default is "auto".
#' @param epsilon A parameter added to ensure the boundness of the residual process. Default is 1e-6.
#' @param LOG Take the log of the residual process. Default is TRUE.
#' @param process Choose between acd or hawkes. Default is acd.
#' @param acd_p The p order of the ACD model. Default is 0.
#' @param acd_q The q order of the ACD model. Default is 1.
#' @examples
#' pw.acd.obj <- new("simACD")
#' pw.acd.obj@cp.loc <- c(0.25,0.75)
#' pw.acd.obj@lambda_0 <- c(1,2,1)
#' pw.acd.obj@alpha <- rep(0.2,3)
#' pw.acd.obj@beta <- rep(0.7,3)
#' pw.acd.obj@N <- 1000
#' pw.acd.obj <- pc_acdsim(pw.acd.obj)
#' ts.plot(Z_trans(pw.acd.obj@x))
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats rnorm rgeom runif optim
#' @importFrom ACDm acdFit
#' @useDynLib eNchange, .registration = TRUE
#' @export
#' @return Returns the transformed residual series.
#' @docType methods
#' @aliases Z_trans Z_trans-class Z_trans-methods
setGeneric(name="Z_trans",
           def=function(H,start.values=c(0.9,0.6),dampen.factor="auto",epsilon=1e-5,LOG=TRUE,process="acd",acd_p=0,acd_q=1)
           {
             standardGeneric("Z_trans")
           }
)
#' @rdname Z_trans-methods
setMethod(f="Z_trans", definition = function(H,start.values=c(0.9,0.6),dampen.factor="auto",epsilon=1e-5,LOG=TRUE,process="acd",acd_p=0,acd_q=1) {
  if (process == "hawkes") {
    estpar = tryCatch({optim(par=start.values,hawkesLike,t_i=cumsum(H),lower=c(0.001,0.001,0.001),method="L-BFGS-B")$par},
                      error=function(e) {message("...Optimizer not converging: using start.values");start.values})
    if (estpar[2]>estpar[3]) estpar[2]=.9*estpar[3] ##this needs to be fixed later
    if (dampen.factor=="auto") dampen.factor=dyn_dampen(estpar[2],estpar[3],type="div")
    estpar[2]=estpar[2]/dampen.factor
    estpar[3]=estpar[3]*dampen.factor
    Z=timeChangeTrans(cumsum(H),estpar)
    if (LOG){
      return(log(Z+epsilon))
    } else return(Z)
  } else if (process=="acd"){
    if (!is.null(start.values) & length(start.values) != sum(1,acd_p,acd_q)) stop("starting values do not match the acd order")
    estpar = tryCatch({acdFit(H,startPara=start.values,order=c(acd_p,acd_q),output=FALSE)},
                      error=function(e) {message("...Optimizer not converging: using start.values");
                        acd.est(H,lambda_0 = start.values[1],alpha=switch(acd_p==0,NULL,start.values[2:(acd_p+1)]),beta=switch(acd_q==0,NULL,start.values[(2+acd_p):(1+acd_p+acd_q)]))})
    #outside unit root test
    if (is.null(estpar)){
      estpar = acd.est(H,lambda_0 = start.values[1],alpha = switch (acd_p==0,NULL,start.values[2:(acd_p+1)]),beta = switch(acd_q==0,NULL,start.values[(2+acd_p):(1+acd_p+acd_q)]))
      estpar$residual=H/estpar$muHats
    } else if(min(sign(estpar$mPara)) == -1 || sum(estpar$mPara[-1]) > (1-.0001) || sum(estpar$order)==0){
      estpar = acd.est(H,lambda_0 = start.values[1],alpha = switch(acd_p==0,NULL,start.values[2:(acd_p+1)]),beta = switch(acd_q==0,NULL,start.values[(2+acd_p):(1+acd_p+acd_q)]))
      estpar$residuals=H/estpar$muHats
      message("...estimated parameters add up to more than 1 or a negative parameter found...")
    }
    order.of.estpar=sum(estpar$order+1)
    psi = estpar$muHats
    if(dampen.factor == "auto") {
      if (estpar$order[1] != 0){
        paraA = estpar$mPara[paste("alpha",1:estpar$order[1],sep="")]
      } else paraA = 0
      if (estpar$order[2] != 0) {
        paraB = estpar$mPara[paste("beta",1:estpar$order[2],sep = "")]
      } else paraB = 0
      dampen.factor=dyn_dampen(paraA,paraB,type="sum")
    }
    if (estpar$order[1]==0) estpar$mPara["alpha1"] = 0
    if (estpar$order[2]==0) estpar$mPara["beta1"] = 0
    dampPara = c(estpar$mPara[1],estpar$mPara[-1]/dampen.factor)
    Z=H
    Z[1] = H[1]/(dampPara['omega']+dampPara['alpha1'] * estpar$residual[1]+epsilon)
    for (i in 2:length(H)) Z[i] = (H[i]/dampPara['omega']+dampPara['alpha1'] * H[i-1]+dampPara['beta1']*psi[i-1]+epsilon)
    if (LOG) {
      return(log(Z+epsilon))
    } else return(Z)
  }
})