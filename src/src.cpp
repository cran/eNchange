#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector cusum(NumericVector x){
  int n = x.size();
  NumericVector out(n);
  NumericVector I_plus(n);
  NumericVector I_minus(n);
  double sum_of_x, i_a,I_plusinv,fct,npw2,n_inv;
  int iter;
  double n_a= (double)n;
  i_a=(double)0;
  I_plusinv = (double)0;
  fct=(double)0;
  n_inv=1/n_a;
  npw2=n_a*n_a;
  sum_of_x = (double)0;
  for (iter=1;iter < n;iter++) sum_of_x += x[iter];
  
  I_minus[0]= 1/sqrt(npw2-n_a)*sum_of_x;
  I_plus[0]=sqrt(1-n_inv)*x[0];
  out[0]=I_plus[0]-I_minus[0];
  
  for (iter=1;iter < n-1;iter++) {
    i_a=(double)iter;
    I_plusinv = 1/(i_a+1);
   fct = sqrt((n_a-i_a-1)*i_a*I_plusinv/(n_a-i_a));
   I_plus[iter]=x[iter]*sqrt(I_plusinv-n_inv)+I_plus[iter-1]*fct;
   I_minus[iter]=I_minus[iter-1] / fct - x[iter] / sqrt(npw2*I_plusinv-n_a);
    out[iter]= (I_plus[iter]-I_minus[iter]);
  }
  //for (iter=0;iter < n-1;iter++) out[iter]=out[iter]/pow(sum_of_x/n_a,1.0);//divide by the mean
  return(out);
}

//[[Rcpp::export]]
List finner_prod_maxp(NumericVector out,double p){
  int n=out.size();
  int max_b=0;
  int iter;
  double aux,maxipi,max_inner;
  max_inner=0;
  maxipi = -1;
  for (iter= 0;iter < n-1;iter++){
    aux = fabs(out[iter]);
    if (aux>maxipi){
      max_b = iter+1;
      maxipi=aux;
    }
  }
  max_inner =out[max_b-1];
  List ret;
  ret["max_inner"] = max_inner;
  ret["max_b"]=max_b;
  return ret;
}

//[[Rcpp::export]]
NumericVector timeChangeTrans(NumericVector t_i,NumericVector par){
  int n=t_i.size();
  NumericVector durationTrans(n);
  NumericVector A(n);
  A[0]=0;
  durationTrans[0]=t_i[0];
  
  for (int i=1;i < n;i++){
    A[i]=exp(-par[2]*(t_i[i]-t_i[i-1]))*(A[i-1]+1);
  }
  
  for (int j=1;j < n;j++){
    durationTrans[j]=par[0]*(t_i[j]-t_i[j-1])+(par[1]/par[2])*(1-exp(-par[2]*(t_i[j]-t_i[j-1])))*A[j-1];
  }
  return durationTrans;
}

//[[Rcpp::export]]
double hawkesLike(NumericVector t_i,NumericVector par){
  double logLike;
  int n=t_i.size();
  double t_n=t_i[n-1];
  NumericVector A(n);
  double Sum1=0;
  double Sum2=0;
  A[0]=0;
  
  for (int i=1;i<n;i++){
    A[i]=exp(-par[2]*(t_i[i]-t_i[i-1]))*(A[i-1]+1);
  }
  for (int j=0;j<n;j++) {
    Sum1 += (par[1]/par[2]) * (exp(-par[2]*(t_n-t_i[j]))-1);
  }
  for (int k=0;k<n;k++) {
  Sum2 += log(par[0]+par[1]*A(k));
  }
  logLike = -par[0]*t_n + Sum1 + Sum2;
  return(-logLike);
}
