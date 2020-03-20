#' An S4 method to detect the change-points in an irregularly spaced time series using Binary Segmentation.
#' @references Korkas Karolos. "Ensemble Binary Segmentation for irregularly spaced data with change-points" Preprint <arXiv:2003.03649>.
#' @rdname BinSeg-methods
#' @description An S4 method to detect the change-points in an irregularly spaced time series using the Binary Segmentation methodology described in Korkas (2020).
#' @param H The input irregular time series.
#' @param thresh The threshold parameter which acts as a stopping rule to detect further change-points and has the form C log(sample). If "universal" then C is data-independent and preselected using the approach described in Korkas (2020). If "boot" it uses the data-dependent method \code{boot_thresh}. Default is "universal".
#' @param q The universal threshold simulation quantile or the bootstrap distribution quantile. Default is 0.99.
#' @param p The support of the CUSUM statistic. Default is 1.
#' @param start.values Warm starts for the optimizers of the likelihood functions.
#' @param dampen.factor The dampen factor in the denominator of the residual process. Default is "auto".
#' @param epsilon A parameter added to ensure the boundness of the residual process. Default is  1e-5.
#' @param LOG Take the log of the residual process. Default is TRUE.
#' @param process Choose between acd or hawkes. Default is acd.
#' @param acd_p The p order of the ACD model. Default is 0.
#' @param acd_q The q order of the ACD model. Default is 1.
#' @param z Transform the time series to use for post-processing. If NULL this is done automatically. Default is NULL.
#' @param do.parallel Choose the number of cores for parallel computation. If 0 no parallelism is done. Default is 2. (Only applies if thresh = "boot").
#' @examples
#' pw.acd.obj <- new("simACD")
#' pw.acd.obj@cp.loc <- seq(0.1,0.95,by=0.025)
#' pw.acd.obj@lambda_0 <- rep(c(0.5,2),1+length(pw.acd.obj@cp.loc)/2)
#' pw.acd.obj@alpha <- rep(0.2,1+length(pw.acd.obj@cp.loc))
#' pw.acd.obj@beta <- rep(0.4,1+length(pw.acd.obj@cp.loc))
#' pw.acd.obj@N <- 5000
#' pw.acd.obj <- pc_acdsim(pw.acd.obj)
#' ts.plot(pw.acd.obj@x,main="Standard BS");abline(v=BinSeg(pw.acd.obj@x)[[1]],col="blue")
#' #real change-points in grey
#' abline(v=floor(pw.acd.obj@cp.loc*pw.acd.obj@N),col="grey",lty=2)
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats rnorm rgeom runif
#' @useDynLib eNchange, .registration = TRUE
#' @export
#' @return Returns a list with the detected change-points and the transformed series.
#' @docType methods
#' @aliases BinSeg BinSeg-class BinSeg-methods
setGeneric(name="BinSeg",
           def=function(H, thresh="universal", q=0.99, p= 1, z=NULL,start.values=c(0.9,0.6),dampen.factor="auto",epsilon= 0.00001,LOG=TRUE,process="acd",acd_p=0,acd_q=1,do.parallel=2)
           {
             standardGeneric("BinSeg")
           }
)
#' @rdname BinSeg-methods
setMethod(f="BinSeg", definition = function(H, thresh="universal", q=0.99, p= 1, z=NULL,start.values=c(0.9,0.6),dampen.factor="auto",epsilon= 0.00001,LOG=TRUE,process="acd",acd_p=0,acd_q=1,do.parallel=2) {
    if (thresh == "universal"){
        thresh = pi_thresh(N=length(H),q=q,process=process)*log(length(H))
    } else if (thresh == "boot"){
        thresh = boot_thresh(H=H,q=q,r=100,p=p,start.values=start.values,process=process,do.parallel=do.parallel,dampen.factor=dampen.factor,epsilon= epsilon,LOG=LOG,acd_p=acd_p,acd_q=acd_q)
    }
    
    if (is.null(z)) z=Z_trans(H=H,start.values = start.values,dampen.factor = dampen.factor,epsilon = epsilon,LOG = LOG,process = process,acd_p = acd_p,acd_q = acd_q)
    cp.est = BinSegTree(z,thresh=thresh,p=p)$breakpoints
    out = list()
    out[[1]] = cp.est
    out[[2]] = z
    return(out)
})