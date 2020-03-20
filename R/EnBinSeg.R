#' An S4 method to detect the change-points in an irregularly spaced time series using Ensemble Binary Segmentation.
#' @references Korkas Karolos. "Ensemble Binary Segmentation for irregularly spaced data with change-points" Preprint <arXiv:2003.03649>.
#' @rdname EnBinSeg-methods
#' @description An S4 method to detect the change-points in an irregularly spaced time series using the Ensemble Binary Segmentation methodology described in Korkas (2020).
#' @param H The input irregular time series.
#' @param thresh The threshold parameter which acts as a stopping rule to detect further change-points and has the form C log(sample). If "universal" then C is data-independent and preselected using the approach described in Korkas (2020). If "boot" it uses the data-dependent method \code{boot_thresh}. Default is "universal".
#' @param q The universal threshold simulation quantile or the bootstrap distribution quantile. Default is 0.99.
#' @param p The support of the CUSUM statistic. Default is 1.
#' @param start.values Warm starts for the optimizers of the likelihood functions.
#' @param dampen.factor The dampen factor in the denominator of the residual process. Default is "auto".
#' @param epsilon A parameter added to ensure the boundness of the residual process. Default is  1e-5.
#' @param LOG Take the log of the residual process. Default is TRUE.
#' @param process Choose between "acd" or "hawkes" or "additive" (signal +iid noise). Default is "acd".
#' @param acd_p The p order of the ACD model. Default is 0.
#' @param acd_q The q order of the ACD model. Default is 1.
#' @param thresh2 Keep only the change-points that appear more than thresh2 M times.
#' @param num_ens Number of ensembles denoted by M in the paper. Default is 500.
#' @param pp Post-process the change-points based on the distance from the highest ranked change-points.
#' @param min_dist The minimum distance as percentage of sample size to use in the post-processing. Default is 0.005.
#' @param do.parallel Choose the number of cores for parallel computation. If 0 no parallelism is done. Default is 2.
#' @param b A parameter to control how close the random end points are to the start points. A large value will on average return shorter random intervals. If NULL all points have an equal chance to be selected (uniformly distributed). Default is NULL.
#' @examples
#' pw.acd.obj <- new("simACD")
#' pw.acd.obj@cp.loc <- seq(0.1,0.95,by=0.025)
#' pw.acd.obj@lambda_0 <- rep(c(0.5,2),1+length(pw.acd.obj@cp.loc)/2)
#' pw.acd.obj@alpha <- rep(0.2,1+length(pw.acd.obj@cp.loc))
#' pw.acd.obj@beta <- rep(0.4,1+length(pw.acd.obj@cp.loc))
#' pw.acd.obj@N <- 5000
#' pw.acd.obj <- pc_acdsim(pw.acd.obj)
#' ts.plot(pw.acd.obj@x,main="Ensemble BS");abline(v=EnBinSeg(pw.acd.obj@x)[[1]],col="red")
#' #real change-points in grey
#' abline(v=floor(pw.acd.obj@cp.loc*pw.acd.obj@N),col="grey",lty=2) 
#' ts.plot(pw.acd.obj@x,main="Standard BS");abline(v=BinSeg(pw.acd.obj@x)[[1]],col="blue")
#' #real change-points in grey
#' abline(v=floor(pw.acd.obj@cp.loc*pw.acd.obj@N),col="grey",lty=2)
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats rnorm rgeom runif rbinom
#' @useDynLib eNchange, .registration = TRUE
#' @export
#' @return Returns a list with the detected change-points and the frequency table of the ensembles across \code{M} applications.
#' @docType methods
#' @aliases EnBinSeg EnBinSeg-class EnBinSeg-methods
setGeneric(name="EnBinSeg",
           def=function(H, thresh="universal", q=0.99, p= 1,start.values=c(0.9,0.6),dampen.factor="auto",epsilon= 1e-5,LOG=TRUE,process="acd",thresh2=0.05,num_ens=500,min_dist=0.005,pp=1,do.parallel=2,b=NULL,acd_p=0,acd_q=1)
           {
             standardGeneric("EnBinSeg")
           }
)
#' @rdname EnBinSeg-methods
setMethod(f="EnBinSeg", definition = function(H, thresh="universal",q=0.99, p= 1,start.values=c(0.9,0.6),dampen.factor="auto",epsilon= 1e-5,LOG=TRUE,process="acd",thresh2=0.05,num_ens=500,min_dist=0.005,pp=1,do.parallel=2,b=NULL,acd_p=0,acd_q=1) {
    if (thresh == "universal"){
        thresh = pi_thresh(N=length(H),q=q,process=process)*log(length(H))
    } else if (thresh == "boot"){
        thresh = boot_thresh(H=H,q=q,r=100,p=p,start.values=start.values,process=process,do.parallel=do.parallel,dampen.factor=dampen.factor,epsilon= epsilon,LOG=LOG,acd_p=acd_p,acd_q=acd_q)
    }
    if (process=="acd" | process=="hawkes"){
        ens1 = BinSeg(H=H,thresh = thresh,p=p,z=NULL,start.values = start.values,dampen.factor=dampen.factor,epsilon = epsilon,LOG=LOG,process = process,acd_p=acd_p,acd_q=acd_q)
    } else{
        ens1 = BinSeg(H=H,thresh = thresh,p=p,z=H,start.values = start.values,dampen.factor=dampen.factor,epsilon = epsilon,LOG=LOG,process = process,acd_p=acd_p,acd_q=acd_q)
    }
    z=ens1[[2]]
    len_of_z = length(z)
    if (!process == "additive") min_dist=floor(min_dist*len_of_z)
    ENS_mat = c()
    rdn_mat=list(NULL)
    j=1
    for (i in 1:num_ens){
        if(!is.null(b)){
            prob=1/b
            rdn1=floor(runif(1,1,len_of_z-b))
            left_or_right=rbinom(1,1,0.5)
            if (left_or_right==1){
                left_or_right=1
            } else left_or_right = -1
            rdn2 = rgeom(1,prob)*left_or_right + rdn1
            if (rdn2 >= len_of_z) rdn2=len_of_z
            if (rdn2 <= 1) rdn2 = 1
            rdn = sort(c(rdn1,rdn2))
            ens=BinSeg(H=H,thresh = thresh,p=p,z = z[rdn[1]:rdn[2]],start.values = start.values,dampen.factor = dampen.factor,epsilon = epsilon,LOG = LOG,process = process,acd_p = acd_p,acd_q = acd_q)
            if (!is.null(ens[[1]])) ENS_mat = c(ENS_mat,ens[[1]]+rdn[1])
        } else {
            rdn=sort(sample(1:len_of_z,2))
            ens=BinSeg(H=H,thresh = thresh,p=p,z = z[rdn[1]:rdn[2]],start.values = start.values,dampen.factor = dampen.factor,epsilon = epsilon,LOG = LOG,process = process,acd_p = acd_p,acd_q = acd_q)
            if (!is.null(ens[[1]])) {
                ENS_mat = c(ENS_mat,ens[[1]]+rdn[1])
                rdn_mat[[j]]=rdn
                j=j+1
            }
        }
    }
    thresh2 = floor(thresh2*num_ens)
    out_ENS_mat = ENS_mat = table(ENS_mat)
    ENS_mat = ENS_mat[ENS_mat>thresh2]
    cp.est = list()
    cp.est[[2]] = z
    cp.est[[3]] = out_ENS_mat
    cp.est[[4]] = rdn_mat
    if (length(ENS_mat)==0){
        return(cp.est)
    } else if (pp){
        message("--post-processing")
        ENS_mat=sort(ENS_mat,decreasing = T)
        cp.est[[1]]=as.numeric(names(ENS_mat))
        len_cp_est = length(cp.est[[1]])
        if (len_cp_est==1) return(cp.est)
        else {
            i=1
            while (i < len_cp_est){
                cp.est[[1]]=cp.est[[1]][cp.est[[1]] != 0]
                if (i >= length(cp.est[[1]])){
                    cp.est[[1]]=cp.est[[1]][cp.est[[1]] != 0]
                    return(cp.est)
                } else if (max(abs(cp.est[[1]][i]-cp.est[[1]][-i]) < min_dist)==TRUE) {
                    len.temp=length(cp.est[[1]])
                    for (j in (i+1):len.temp){
                        if (abs(cp.est[[1]][i]-cp.est[[1]][j]) < min_dist) cp.est[[1]][j] = 0
                    }
                    i=1
                    if (length(cp.est[[1]])==1){
                        cp.est[[1]] = cp.est[[1]][cp.est[[1]] !=0 ]
                        return(cp.est)
                        break
                    }
                } else i=i+1
            }
            cp.est[[1]] = cp.est[[1]][cp.est[[1]]!= 0]
            return(cp.est)
        }
    } else {
        ENS_mat =sort(ENS_mat,decreasing = T)
        cp.est[[1]] = as.numeric(names(ENS_mat))
        return(cp.est)
    }
})