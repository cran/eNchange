% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simHawkes.R
\docType{class}
\name{simHawkes-class}
\alias{simHawkes-class}
\title{An S4 class for a nonstationary ACD model.}
\value{
Returns an object of \code{simHawkes} class.
}
\description{
A specification class to create an object of a simulated piecewise constant Hawkes model of order (1,1).
We consider the following time-varying piecewise constant Hawkes process (which we term tvHawkes)
\eqn{\lambda({\upsilon}) = \lambda_0({\upsilon}) +\sum_{{\upsilon}_t < s} \alpha({\upsilon})e^{-\beta({\upsilon}) ({\upsilon}-{\upsilon}_t)}, \ \mbox{for} \ {\upsilon} = 1, \ldots,T}.
}
\section{Slots}{

\describe{
\item{\code{H}}{The durational time series.}

\item{\code{cH}}{The psi time series.}

\item{\code{horizon}}{The time horizon of a Hawkes process typically expressed in seconds. Effective sample size will differ depending on the size of the parameters.}

\item{\code{N}}{Effective sample size which differs depending on the size of the parameters.}

\item{\code{cp.loc}}{The vector with the location of the changepoints. Takes values from 0 to 1 or NULL if none. Default is NULL.}

\item{\code{lambda_0}}{The vector of the parameters lambda_0 in the Hawkes model as in the above formula.}

\item{\code{alpha}}{The vector of the parameters alpha in the Hawkes model as in the above formula.}

\item{\code{beta}}{The vector of the parameters beta in the Hawkes model as in the above formula.}
}}

\examples{
pw.hawk.obj <- new("simHawkes")
pw.hawk.obj@cp.loc <- c(0.5)
pw.hawk.obj@lambda_0 <- c(1,2)
pw.hawk.obj@alpha <- c(0.2,0.2)
pw.hawk.obj@beta <- c(0.7,0.7)
pw.hawk.obj@horizon <- 1000
pw.hawk.obj <- pc_hawkessim(pw.hawk.obj)
ts.plot(pw.hawk.obj@H)
ts.plot(pw.hawk.obj@cH)
}
\references{
Korkas, K.K., 2022. Ensemble binary segmentation for irregularly spaced data with change-points. Journal of the Korean Statistical Society, 51(1), pp.65-86.
}
