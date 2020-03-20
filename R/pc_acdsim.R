#' A method to simulate nonstationary ACD models.
#' @name pc_acdsim-class
#' @param object a simACD object
#' @description A S4 method that takes as an input a \code{simACD} object and outputs a simulated nonstationary ACD(1,1) model. The formulation of the of the 
#' piecewise constant ACD model is given in the \code{simACD} class.
#' @references
#' Korkas Karolos. "Ensemble Binary Segmentation for irregularly spaced data with change-points" Preprint.
#' @examples
#' pw.acd.obj <- new("simACD")
#' pw.acd.obj@cp.loc <- c(0.25,0.75)
#' pw.acd.obj@lambda_0 <- c(1,2,1)
#' pw.acd.obj@alpha <- rep(0.2,3)
#' pw.acd.obj@beta <- rep(0.7,3)
#' pw.acd.obj@N <- 3000
#' pw.acd.obj <- pc_acdsim(pw.acd.obj)
#' ts.plot(pw.acd.obj@x)
#' ts.plot(pw.acd.obj@psi)
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats rnorm rgeom runif rexp
#' @importFrom utils head tail 
#' @useDynLib eNchange, .registration = TRUE
#' @export
#' @return Returns an object of \code{simACD} class containing a simulated piecewise constant ACD time series.
#' @docType methods
#' @rdname pc_acdsim-methods
#' @aliases pc_acdsim pc_acdsim-class pc_acdsim-methods
setGeneric(name="pc_acdsim",
           def=function(object)
           {
             standardGeneric("pc_acdsim")
           }
)
#' @rdname pc_acdsim-methods
setMethod(f="pc_acdsim", signature= "simACD", definition = function(object) {
    if (is.null(object@cp.loc)){
      temp.H=sim.ACD(object@N,object@lambda_0,object@alpha,object@beta,object@BurnIn)
      object@x=temp.H[[1]]
      object@psi=temp.H[[2]]
      return(object)
    } else {
      num.of.cp=length(object@cp.loc)
      if (num.of.cp==1){
        temp.H1=sim.ACD(N=floor(object@N*object@cp.loc),lambda_0=object@lambda_0[1],alpha = object@alpha[1],beta = object@beta[1],BurnIn = object@BurnIn)
        temp.H2=sim.ACD(N=floor(object@N*(1-object@cp.loc)),lambda_0=object@lambda_0[2],alpha = object@alpha[2],beta = object@beta[2],BurnIn = NULL)
        object@x=c(temp.H1[[1]],temp.H2[[1]])
        object@psi=c(temp.H1[[2]],temp.H2[[2]])
        return(object)
      } else {
        res=c()
        res2=c()
        for (i in 1:num.of.cp){
          if (i==1){
            temp.H=sim.ACD(N=floor(object@N*object@cp.loc[i]),lambda_0=object@lambda_0[i],alpha=object@alpha[i],beta=object@beta[i],BurnIn = object@BurnIn)
          } else temp.H=sim.ACD(N=floor(object@N*(object@cp.loc[i]-object@cp.loc[i-1])),lambda_0=object@lambda_0[i],alpha=object@alpha[i],beta=object@beta[i],BurnIn = NULL)
          res=c(res,temp.H[[1]])
          res2=c(res2,temp.H[[2]])
        }
        temp.H=sim.ACD(N=floor(object@N*(1-object@cp.loc[num.of.cp])),lambda_0=object@lambda_0[num.of.cp+1],alpha=object@alpha[num.of.cp+1],beta=object@beta[num.of.cp+1],BurnIn = object@BurnIn)
        res=c(res,temp.H[[1]])
        res2=c(res2,temp.H[[2]])
        object@x = res
        object@psi = res2
        return(object)
    }
  }
}
)

