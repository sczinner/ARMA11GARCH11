#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

List ARMA11GARCH11(NumericVector X, NumericVector params){
  double mu = params[0];
  double aa = params[1];
  double bb = params[2];
  double omega = params[3];
  double alpha = params[4];
  double beta = params[5];
  
  NumericVector Xhat(X.length()+1);
  NumericVector Sigma2hat(X.length()+1);
  NumericVector error(X.length());
  
  Xhat[0] = 0;
  Sigma2hat[0] = var(X);
  error[0] = X[0]-Xhat[0];
  
  for(int tt=0; tt<X.length(); tt++){
    Xhat[tt+1]=mu+aa*error[tt]+bb*X[tt];
    Sigma2hat[tt+1]=omega+alpha*pow(error[tt],2)+beta*Sigma2hat[tt];
    if((tt+1)!=X.length()){
      error[tt+1]=Xhat[tt+1]-X[tt+1];
    }
  }
  List ret; ret["Xhat"]=Xhat; ret["Sigma2hat"]=Sigma2hat; ret["error"]=error;
  return ret;
}

// [[Rcpp::export]]

NumericVector llgr(NumericVector X, NumericVector params){
  double mu = params[0];
  double aa = params[1];
  double bb = params[2];
  double omega = params[3];
  double alpha = params[4];
  double beta = params[5];
  
  NumericVector Xhat(X.length()+1);
  NumericVector Sigma2hat(X.length()+1);
  NumericVector error(X.length());
  
  NumericMatrix dSigma2hat(X.length(),6);
  NumericMatrix derror(X.length(),6);
  
  Xhat[0] = 0;
  Sigma2hat[0] = var(X);
  error[0] = X[0]-Xhat[0];
  
  NumericVector gr = (Sigma2hat[0]*2*error[0]*derror(0,_)-
    pow(error[0],2)*dSigma2hat(0,_))/
      pow(Sigma2hat[0],2)+(1/Sigma2hat[0])*dSigma2hat(0,_);
  
  for(int tt=0; tt<X.length(); tt++){
    Xhat[tt+1]=mu+aa*error[tt]+bb*X[tt];
    Sigma2hat[tt+1]=omega+alpha*pow(error[tt],2)+beta*Sigma2hat[tt];
    if((tt+1)!=X.length()){
      error[tt+1]=Xhat[tt+1]-X[tt+1];
      
      NumericVector vec1 = NumericVector::create(1,error[tt],X[tt],0,0,0);
      NumericVector vec2 = NumericVector::create(aa*derror(tt,0),aa*derror(tt,1),aa*derror(tt,2),0,0,0);
      NumericVector vec3 = NumericVector::create(alpha*2*error[tt]*derror(tt,0),
                                                 alpha*2*error[tt]*derror(tt,1),
                                                 alpha*2*error[tt]*derror(tt,2),1,pow(error[tt],2),Sigma2hat[tt]);
      
      derror(tt+1,_)=vec1+vec2;
      dSigma2hat(tt+1,_)=vec3+beta*dSigma2hat(tt,_);
      
      
      gr = gr+(Sigma2hat[tt]*2*error[tt]*derror(tt,_)-
        pow(error[tt],2)*dSigma2hat(tt,_))/
          pow(Sigma2hat[tt],2)+(1/Sigma2hat[tt])*dSigma2hat(tt,_);
    }
  }
  return gr;
}
// [[Rcpp::export]]
List optgd(NumericVector X, int nits, NumericVector initpar){
  NumericMatrix pars(nits+1,6);
  NumericVector lls(nits+1);
  pars(0,_) = initpar;
  Function ll("log.likelihood");
  lls[0] = 1;
  Function concat("c");
  for(int ii = 1; ii <nits+1; ii++){
    NumericVector gr = llgr(X,pars(ii-1,_));
    gr=clamp(-pow(10,5),gr,pow(10,5));
    NumericVector eta=concat(rep(pow(10,-10),3),pow(10,-11),rep(pow(10,-6),2));
    NumericVector parsii=pars(ii-1,_)-eta*gr;
    parsii[Range(3,5)]=pmax(parsii[ Range(3,5) ],rep(0.000000001,3));
    pars(ii,_)=parsii;
    lls[ii]=1;
  }
  
  List ret; ret["pars"]=pars; ret["lls"]=lls;
  return ret;
}




/*** R

*/
