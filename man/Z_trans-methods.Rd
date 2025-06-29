% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Z_trans.R
\docType{methods}
\name{Z_trans}
\alias{Z_trans}
\alias{Z_trans-class}
\alias{Z_trans-methods}
\alias{Z_trans,ANY-method}
\title{Transformation of an irregularly spaces time series.}
\usage{
Z_trans(
  H,
  start.values = c(0.9, 0.6),
  dampen.factor = "auto",
  epsilon = 1e-05,
  LOG = TRUE,
  process = "acd",
  acd_p = 0,
  acd_q = 1
)

\S4method{Z_trans}{ANY}(
  H,
  start.values = c(0.9, 0.6),
  dampen.factor = "auto",
  epsilon = 1e-05,
  LOG = TRUE,
  process = "acd",
  acd_p = 0,
  acd_q = 1
)
}
\arguments{
\item{H}{The input irregular time series.}

\item{start.values}{Warm starts for the optimizers of the likelihood functions.}

\item{dampen.factor}{The dampen factor in the denominator of the residual process. Default is "auto".}

\item{epsilon}{A parameter added to ensure the boundness of the residual process. Default is 1e-6.}

\item{LOG}{Take the log of the residual process. Default is TRUE.}

\item{process}{Choose between acd or hawkes. Default is acd.}

\item{acd_p}{The p order of the ACD model. Default is 0.}

\item{acd_q}{The q order of the ACD model. Default is 1.}
}
\value{
Returns the transformed residual series.
}
\description{
Transformation of a irregularly spaces time series. For the tvACD model, we calculate
\eqn{U_t = g_0(x_t, \psi_t) = \frac{x_t}{{\psi}_t}}, where
\eqn{{\psi}_t = C_0 + \sum_{j=1}^p C_j x_{t-j} + \sum_{k=1}^q C_{p+k} \psi_{t-k}+\epsilon x_t}.
where the last term \eqn{\epsilon x_t} is added to ensure the boundness of \eqn{U_t}.
}
\examples{
pw.acd.obj <- new("simACD")
pw.acd.obj@cp.loc <- c(0.25,0.75)
pw.acd.obj@lambda_0 <- c(1,2,1)
pw.acd.obj@alpha <- rep(0.2,3)
pw.acd.obj@beta <- rep(0.7,3)
pw.acd.obj@N <- 1000
pw.acd.obj <- pc_acdsim(pw.acd.obj)
ts.plot(Z_trans(pw.acd.obj@x))
}
\references{
Korkas, K.K., 2022. Ensemble binary segmentation for irregularly spaced data with change-points. Journal of the Korean Statistical Society, 51(1), pp.65-86.
}
