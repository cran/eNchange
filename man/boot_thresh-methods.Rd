% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_thresh.R
\docType{methods}
\name{boot_thresh}
\alias{boot_thresh}
\alias{boot_thresh-class}
\alias{boot_thresh-methods}
\alias{boot_thresh,ANY-method}
\title{A bootstrap method to calculate the threshold (stopping rule) in the BS or EBS segmentation.}
\usage{
boot_thresh(
  H,
  q = 0.75,
  r = 100,
  p = 1,
  start.values = c(0.9, 0.6),
  process = "acd",
  do.parallel = 2,
  dampen.factor = "auto",
  epsilon = 1e-05,
  LOG = TRUE,
  acd_p = 0,
  acd_q = 1
)

\S4method{boot_thresh}{ANY}(
  H,
  q = 0.75,
  r = 100,
  p = 1,
  start.values = c(0.9, 0.6),
  process = "acd",
  do.parallel = 2,
  dampen.factor = "auto",
  epsilon = 1e-05,
  LOG = TRUE,
  acd_p = 0,
  acd_q = 1
)
}
\arguments{
\item{H}{The input irregular time series.}

\item{q}{The bootstrap distribution quantile. Default is 0.75.}

\item{r}{The number of bootrstap simulations. Default is 100.}

\item{p}{The support of the CUSUM statistic. Default is 1.}

\item{start.values}{Warm starts for the optimizers of the likelihood functions.}

\item{process}{Choose between acd or hawkes. Default is acd.}

\item{do.parallel}{Choose the number of cores for parallel computation. If 0 no parallelism is done. Default is 2.}

\item{dampen.factor}{The dampen factor in the denominator of the residual process. Default is "auto".}

\item{epsilon}{A parameter added to ensure the boundness of the residual process. Default is  1e-5.}

\item{LOG}{Take the log of the residual process. Default is TRUE.}

\item{acd_p}{The p order of the ACD model. Default is 0.}

\item{acd_q}{The q order of the ACD model. Default is 1.}
}
\value{
Returns the threshold \code{C}.
}
\description{
A bootstrap method to calculate the threshold (stopping rule) in the BS or EBS segmentation described in Cho and Korkas (2022) and adapted for irregularly time series in Korkas (2022).
}
\examples{
pw.acd.obj <- new("simACD")
pw.acd.obj@cp.loc <- c(0.25,0.75)
pw.acd.obj@lambda_0 <- c(1,2,1)
pw.acd.obj@alpha <- rep(0.2,3)
pw.acd.obj@beta <- rep(0.7,3)
pw.acd.obj@N <- 3000
pw.acd.obj <- pc_acdsim(pw.acd.obj)
boot_thresh(pw.acd.obj@x,r=20)
}
\references{
Cho, H. and Korkas, K.K., 2022. High-dimensional GARCH process segmentation with an application to Value-at-Risk. Econometrics and Statistics, 23, pp.187-203.
}
