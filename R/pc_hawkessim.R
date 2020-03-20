#' A method to simulate nonstationary Hawkes models.
#' @name pc_hawkessim-class
#' @param object a simHawkes object
#' @description A S4 method that takes as an input a \code{simHawkes} object and outputs a simulated nonstationary Hawkes model. The formulation of the of the 
#' piecewise constant ACD model is given in the \code{simHawkes} class.
#' @references
#' Korkas Karolos. "Ensemble Binary Segmentation for irregularly spaced data with change-points" <arXiv:2003.03649>.
#' @examples
#' pw.hawk.obj <- new("simHawkes")
#' pw.hawk.obj@cp.loc <- c(0.5)
#' pw.hawk.obj@lambda_0 <- c(1,2)
#' pw.hawk.obj@alpha <- c(0.2,0.2)
#' pw.hawk.obj@beta <- c(0.7,0.7)
#' pw.hawk.obj@horizon <- 1000
#' pw.hawk.obj <- pc_hawkessim(pw.hawk.obj)
#' ts.plot(pw.hawk.obj@H)
#' ts.plot(pw.hawk.obj@cH)
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom hawkes simulateHawkes
#' @importFrom utils head tail 
#' @importFrom stats rnorm rgeom runif
#' @useDynLib eNchange, .registration = TRUE
#' @export
#' @return Returns an object of \code{simHawkes} class containing a simulated piecewise constant Hawkes series.
#' @docType methods
#' @rdname pc_hawkessim-methods
#' @aliases pc_hawkessim pc_hawkessim-class pc_hawkessim-methods
setGeneric(name="pc_hawkessim",
           def=function(object)
           {
             standardGeneric("pc_hawkessim")
           }
)
#' @rdname pc_hawkessim-methods
setMethod(f="pc_hawkessim", signature= "simHawkes", definition = function(object) {
  if (is.null(object@cp.loc)){
    temp.H = diff(simulateHawkes(object@lambda_0,object@alpha,object@beta,object@horizon)[[1]])
    object@H=temp.H
    object@cH=cumsum(temp.H)
    object@N = length(object@H)
    return(object)
  } else {
    num.of.cp = length(object@cp.loc)
    if (num.of.cp==1) {
      temp.H1 = diff(simulateHawkes(object@lambda_0[1],object@alpha[1],object@beta[1],horizon=floor(object@horizon*object@cp.loc))[[1]])
      temp.H2 = diff(simulateHawkes(object@lambda_0[2],object@alpha[2],object@beta[2],horizon=floor(object@horizon*(1-object@cp.loc)))[[1]])
      object@H=c(temp.H1,temp.H2)
      object@cH=cumsum(object@H)
      object@N = length(object@H)
      return(object)
    } else {
      res=c()
      res.length=c()
      for (i in 1:(num.of.cp)){
        temp.H = diff(simulateHawkes(object@lambda_0[i],object@alpha[i],object@beta[i],horizon=floor(object@horizon*object@cp.loc[i]))[[1]])
        temp.H.length = length(temp.H.length)
        res = c(res,temp.H)
        res.length=c(res.length,temp.H.length)
      }
      temp.H=diff(simulateHawkes(object@lambda_0[num.of.cp+1],object@alpha[num.of.cp+1],object@beta[num.of.cp+1],horizon=floor(object@horizon*(1-object@cp.loc[num.of.cp])))[[1]])
      temp.H.length=length(temp.H)
      object@H=res
      object@cH=cumsum(object@H)
      object@N = length(object@H)
      return(object)
    }
  }
}
)

