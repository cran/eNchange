#' An S4 class for a nonstationary ACD model.
#' @name simACD-class
#' @rdname simACD-class
#' @description A specification class to create an object of a simulated piecewise constant conditional duration model of order (1,1). 
#' \eqn{x_t / \psi_t = \varepsilon_t \; \sim \mathcal{G}(\theta_2)}
#' \eqn{\psi_t = \omega(t) + \sum_{j=1}^p \alpha_{j}(t)x_{t-j} + \sum_{k=1}^q \beta_{k}(t)\psi_{t-k}.}
#' where \eqn{\psi_{t} = \mathcal{E} [x_t | x_t,\ldots,x_1| \theta_1]} is the conditional mean duration of the \eqn{t}-th event with parameter vector \eqn{\theta_1} and \eqn{\mathcal{G}(.)}
#' is a general distribution over \eqn{(0,+\infty)} with mean equal to 1 and parameter vector 
#' \eqn{\theta_2}. In this work we assume that \eqn{\varepsilon_t \; \sim \exp(1)}.
#' @slot x The durational time series.
#' @slot psi The psi time series.
#' @slot N Sample sze of the time series.
#' @slot cp.loc The vector with the location of the changepoints. Takes values from 0 to 1 or NULL. Default is NULL.
#' @slot lambda_0 The vector of the parameters lambda_0 in the ACD series as in the above formula.
#' @slot alpha The vector of the parameters alpha in the ACD series as in the above formula.
#' @slot beta The vector of the parameters beta in the ACD series as in the above formula.
#' @slot BurnIn The size of the burn-in sample. Note that this only applies at the first simulated segment. Default is 500.
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
#' @import Rcpp foreach doParallel parallel iterators methods
#' @importFrom stats rnorm rgeom runif rexp
#' @useDynLib eNchange, .registration = TRUE
#' @export
#' @return Returns an object of \code{simACD} class.
setClass("simACD", 
         slots = c(
           x = "numeric",
           psi = "numeric",
           N = "numeric", 
           cp.loc = "numeric",
           lambda_0 = "numeric",
           alpha = "numeric",
           beta = "numeric",
           BurnIn = "numeric"
         ), 
         prototype = list(
           cp.loc=NULL,
           BurnIn=500
         )
)