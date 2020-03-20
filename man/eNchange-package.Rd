\name{eNchange-package}
\alias{eNchange-package}
\alias{eNchange}
\docType{package}
\title{
 Ensemble Methods for Multiple Change-Point Detection
}
\description{
Implements a segmentation algorithm for multiple change-point detection in univariate time series using the Ensemble Binary Segmentation of Korkas (2020) <arXiv:2003.03649>.
}
\details{
We propose a new technique for
consistent estimation of the number and locations of the
change-points in the structure of an irregularly spaced time series. The
core of the segmentation procedure is the Ensemble Binary Segmentation
method (EBS), a technique in which  a
large number of multiple change-point detection tasks using the Binary Segmentation (BS) method are applied on sub-samples of the
data of differing lengths, and then the results are combined to create
an overall answer. This methodology is applied to irregularly time series models such as the time-varying Autoregressive Conditional Duration model or the time-varying  Hawkes process.
}
\author{
Karolos K. Korkas <kkorkas@yahoo.co.uk>.

Maintainer: Karolos K. Korkas <kkorkas@yahoo.co.uk>
}
\references{
Korkas Karolos. "Ensemble Binary Segmentation for irregularly spaced data with change-points" <arXiv:2003.03649>.
}
\keyword{ eNchange }

\examples{
  \dontrun{
 pw.acd.obj <- new("simACD")
 pw.acd.obj@cp.loc <- seq(0.1,0.95,by=0.025)
 pw.acd.obj@lambda_0 <- rep(c(0.5,2),1+length(pw.acd.obj@cp.loc)/2)
 pw.acd.obj@alpha <- rep(0.2,1+length(pw.acd.obj@cp.loc))
 pw.acd.obj@beta <- rep(0.4,1+length(pw.acd.obj@cp.loc))
 pw.acd.obj@N <- 5000
 pw.acd.obj <- pc_acdsim(pw.acd.obj)
 ts.plot(pw.acd.obj@x,main="Ensemble BS");abline(v=EnBinSeg(pw.acd.obj@x)[[1]],col="red")
 #real change-points in grey
 abline(v=floor(pw.acd.obj@cp.loc*pw.acd.obj@N),col="grey",lty=2) 
 ts.plot(pw.acd.obj@x,main="Standard BS");abline(v=BinSeg(pw.acd.obj@x)[[1]],col="blue")
 #real change-points in grey
 abline(v=floor(pw.acd.obj@cp.loc*pw.acd.obj@N),col="grey",lty=2)

  }
}
