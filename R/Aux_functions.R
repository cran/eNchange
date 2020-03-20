BinSegm <- function(x,thresh,p=.75){
  n=length(x)
  if(n>1){
    bb=finner_prod_maxp(cusum(x),p=p)
    bb[[2]]=bb[[2]]+1
    if(abs(bb[[1]])>thresh){
      return(c(BinSegm(x[1:bb[[2]]],thresh),BinSegm(x[(bb[[2]]+1):n],thresh)))
    } else return(rep(mean(x),n))
  }else return(x)
}


pi_thresh <- function(N,q=.99,process="acd"){
  if (N<100) stop("sample size too small to obtain a reliable parameter")
  if (N>100000) N=100000
  if (q==.95){
    coef=c( 3.385e-01,  -4.919e-07 , 2.610e-12 ,  1.351e+02)
  } else if (q==.99){
    coef=c(  4.453e-01, -9.245e-07 , 4.556e-12 ,  1.573e+02)
  } else if (q==.75){
    coef=c(  2.661e-01, -7.586e-07 , 4.596e-12 ,  7.796e+01)
  } else stop("this quantile is not supported (yet)")
  return(coef%*%c(1,N,N^2,1/N))
}

dyn_dampen <- function(paramA,paramB,type="sum"){
  if (type == "sum"){
    return(max(1,min(.99,sum(paramA)+sum(paramB)))/(max(0.01,1-(sum(paramA)+sum(paramB)))))
  } else if ((type=="div")) return(max(1,min(.99,sum(paramA)/sum(paramB)))/(max(0.01,1-(sum(paramA)/sum(paramB)))))
}

sim.ACD <- function(N,lambda_0,alpha,beta,BurnIn=200,resids=NULL){
  if (is.null(BurnIn)) BurnIn=0
  x=rep(0,N+BurnIn)
  psi=rep(0,N+BurnIn)
  if (is.null(resids)){
    x[1]=rexp(1)
  } else x[1]=sample(resids,1,replace = TRUE)
  psi[1]=lambda_0+alpha*x[1]
  for (i in 2:(N+BurnIn)){
    psi[i]=lambda_0+alpha*x[i-1]+beta*psi[i-1]
    if (is.null(resids)){
      x[i]=psi[i]*rexp(1)
    } else x[i]=psi[i]*sample(resids,1,replace = TRUE)
  }
  out=list(NULL)
  out[[1]]=tail(x,N)
  out[[2]]=tail(psi,N)
  return(out)
}

acd.est <- function(H,lambda_0,alpha,beta){
  order=c(length(alpha),length(beta))
  n=length(H)
  psi_t=matrix(n,1,0)
  psi_t[1:(max(order))]=lambda_0[1]/(1-sum(c(alpha,beta)))
  for (i in (max(order)+1):n) psi_t[i]=sum(lambda_0,alpha * H[i-(1:order[1])],beta * psi_t[i-(1:order[2])])
  acd.obj=list()
  acd.obj$muHats=psi_t
  acd.obj$order=order
  acd.obj$mPara['omega'] = lambda_0
  if (order[1] !=0) acd.obj$mPara[paste("alpha",1:order[1],sep="")] = alpha
  if (order[2] !=0) acd.obj$mPara[paste("beta",1:order[2],sep="")] = beta
  return(acd.obj)
}

BinSegTree <- function(z,thresh, p=0.75,del=20){
  n = length(z)
  estimates = br.list = NULL
  criterion = thresh
  breakpoints = NULL
  f = NULL
  tree = list(matrix(0, 6, 1))
  tree[[1]][1,1] = 1
  s = tree[[1]][4,1] = 1
  e = tree[[1]][6,1] = n
  temp.r = finner_prod_maxp(cusum(z[s:e]),p=p)
  b = temp.r[[2]]
  tree[[1]][5,1] <- b
  d<-temp.r[[1]]
  
  if(abs(d)>criterion){
    tree[[1]][2,1]<-d
    breakpoints <- c(breakpoints, b)
    f<-c(f, d)  
    j<-1
    while(length(tree)==j){
      if(sum(tree[[j]][6, ]-tree[[j]][4, ]-rep(2*del, dim(tree[[j]])[2]))>0){ 
        no.parent.coeffs<-dim(tree[[j]])[2]
        no.child.coeffs<-0
        for(i in 1:no.parent.coeffs){
          if(tree[[j]][5, i]-tree[[j]][4, i] > del+1){
            s<-tree[[j]][4, i]
            e<-tree[[j]][5, i]
            temp.r = finner_prod_maxp(cusum(z[s:e]),p=p)
            b <- s+temp.r[[2]]-1
            d <- temp.r[[1]]
            
            if(abs(d)>criterion){
              if(length(tree)==j) tree<-c(tree, list(matrix(0, 6, 0)))
              no.child.coeffs<-no.child.coeffs+1
              tree[[j+1]]<-matrix(c(tree[[j+1]], matrix(0, 6, 1)), 6, no.child.coeffs)
              tree[[j+1]][1, no.child.coeffs]<-2*tree[[j]][1, i]-1
              tree[[j+1]][2, no.child.coeffs]<-d
              tree[[j+1]][4, no.child.coeffs]<-s
              tree[[j+1]][6, no.child.coeffs]<-e
              tree[[j+1]][5, no.child.coeffs]<-b
              if (sum(abs(breakpoints-b)>del)==length(breakpoints))  breakpoints<-c(breakpoints,b);
              f<-c(f, d)
            }
          }
          if(tree[[j]][6, i]-tree[[j]][5, i]> del+1){
            s<-tree[[j]][5, i]+1
            e<-tree[[j]][6, i]
            temp.r = finner_prod_maxp(cusum(z[s:e]),p=p)
            b <- s+temp.r[[2]]-1
            d <- temp.r[[1]]
            
            if(abs(d)>criterion){
              if(length(tree)==j) tree<-c(tree, list(matrix(0, 6, 0)))
              no.child.coeffs<-no.child.coeffs+1
              tree[[j+1]]<-matrix(c(tree[[j+1]], matrix(0, 6, 1)), 6, no.child.coeffs)
              tree[[j+1]][1, no.child.coeffs]<-2*tree[[j]][1, i]
              tree[[j+1]][2, no.child.coeffs]<-d 
              tree[[j+1]][4, no.child.coeffs]<-s
              tree[[j+1]][6, no.child.coeffs]<-e
              tree[[j+1]][5, no.child.coeffs]<-b
              if (sum(abs(breakpoints-b)>del)==length(breakpoints))  breakpoints<-c(breakpoints,b)
              f<-c(f, d)
              
            }
          }
        }
        
      }
      j<-j+1
    }
  }
  list(tree=tree, breakpoints=breakpoints, f=f)  
}

HR <- function(cp.real,cp.est,max.len=50){
  counter=0
  if(length(cp.real)==0) cp.real=NULL
  if(length(cp.est)==0) cp.est=NULL
  if(is.null(cp.real) & is.null(cp.est)){
    list(HR=1,HR2=1)
  } else if(is.null(cp.real) & !is.null(cp.est)){
    list(HR=0,HR2=0)
  } else if(is.null(cp.est)) {
    list(HR=0,HR2=0)
  } else if(is.na(cp.est[[1]])){
    list(HR=0,HR2=0)
  } else {
    for (i in 1:length(cp.real)){
      if(sum(abs(cp.real[i]-cp.est)<max.len)>=1) counter=counter+1
      if(i==length(cp.real)){
        HR=counter
        HR2=HR/max(length(cp.real),length(cp.est))
      }
    }
    list(HR=HR,HR2=HR2)
  }
}


