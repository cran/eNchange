#' An S4 class for a nonstationary ACD model.
#' @name simHawkes-class
#' @rdname simHawkes-class
#' @description A specification class to create an object of a simulated piecewise constant Hawkes model of order (1,1). 
#' We consider the following time-varying piecewise constant Hawkes process (which we term tvHawkes)
#' \eqn{\lambda({\upsilon}) = \lambda_0({\upsilon}) +\sum_{{\upsilon}_t < s} \alpha({\upsilon})e^{-\beta({\upsilon}) ({\upsilon}-{\upsilon}_t)}, \ \mbox{for} \ {\upsilon} = 1, \ldots,T}.
#' @slot H The durational time series.
#' @slot cH The psi time series.
#' @slot horizon The time horizon of a Hawkes process typically expressed in seconds. Effective sample size will differ depending on the size of the parameters.
#' @slot N Effective sample size which differs depending on the size of the parameters.
#' @slot cp.loc The vector with the location of the changepoints. Takes values from 0 to 1 or NULL if none. Default is NULL.
#' @slot lambda_0 The vector of the parameters lambda_0 in the Hawkes model as in the above formula.
#' @slot alpha The vector of the parameters alpha in the Hawkes model as in the above formula.
#' @slot beta The vector of the parameters beta in the Hawkes model as in the above formula.
#' @references
#' Korkas Karolos. "Ensemble Binary Segmentation for irregularly spaced data with change-points" Preprint.
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
#' @import Rcpp foreach doParallel parallel iterators methods
#' @importFrom  hawkes simulateHawkes
#' @importFrom stats rnorm rgeom runif
#' @useDynLib eNchange, .registration = TRUE
#' @export
#' @return Returns an object of \code{simHawkes} class.
setClass("simHawkes", 
         slots = c(
           H = "numeric",
           cH = "numeric",
           horizon = "numeric",
           N = "numeric",
           cp.loc = "numeric",
           lambda_0 = "numeric",
           alpha = "numeric",
           beta = "numeric"
         )
)