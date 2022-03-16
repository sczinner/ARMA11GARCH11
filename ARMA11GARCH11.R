#The ARMA(1,1)-GARCH(1,1) specification. 
#Inputs are the stock returns, and the parameters, c(mu,a,b,omega,alpha,beta)
ARMA11GARCH11.R<-function(X,params){
  mu<-params[1]
  aa<-params[2]
  bb<-params[3]
  omega<-params[4]
  alpha<-params[5]
  beta<-params[6]
  
  Xhat<-numeric(length(X)+1)
  Sigma2hat<-numeric(length(X)+1)
  error<-numeric(length(X))
  Xhat[1]<-0
  Sigma2hat[1]<-var(X)#freebie
  error[1]<-X[1]-Xhat[1]
  
  for(tt in 1:length(X)){
    Xhat[tt+1]<-mu+aa*error[tt]+bb*X[tt]
    Sigma2hat[tt+1]<-omega+alpha*(error[tt])^2+beta*Sigma2hat[tt]
    if(tt!=length(X))
      error[tt+1]<-Xhat[tt+1]-X[tt+1]
  }
  return(list(Xhat=Xhat,Sigma2hat=Sigma2hat,error=error))
}

#The negative of the log-likelihood. Minimizing this function means better 
#fitting parameters.
#Inputs are the stock returns, and the parameters, c(mu,a,b,omega,alpha,beta)
log.likelihood<-function(X,params){
  X<-as.numeric(X)
  AG<-ARMA11GARCH11.R(X,params)
  ll<-sum(AG$error^2/head(AG$Sigma2hat,-1)+log(head(AG$Sigma2hat,-1)))
  return(ll)
}

##Computes the gradient of the log-likelihood.
#Inputs are the stock returns, and the parameters, c(mu,a,b,omega,alpha,beta)
ll.gr.R<-function(X,params){
  X<-as.numeric(X)
  mu<-params[1]
  aa<-params[2]
  bb<-params[3]
  omega<-params[4]
  alpha<-params[5]
  beta<-params[6]
  
  Xhat<-numeric(length(X)+1)
  Sigma2hat<-numeric(length(X)+1)
  error<-numeric(length(X))
  
  Xhat[1]<-0
  Sigma2hat[1]<-var(X)#freebie
  error[1]<-X[1]-Xhat[1]
  
  derror<-matrix(data=0,nrow=length(error),ncol=6)
  dSigma2hat<-matrix(data=0,nrow=length(error),ncol=6)
  
  gr<-(Sigma2hat[1]*2*error[1]*derror[1,]-
         error[1]^2*dSigma2hat[1,])/
    Sigma2hat[1]^2+(1/Sigma2hat[1])*dSigma2hat[1,]
  
  for(tt in 1:length(X)){
    Xhat[tt+1]<-mu+aa*error[tt]+bb*X[tt]
    Sigma2hat[tt+1]<-omega+alpha*(error[tt])^2+beta*Sigma2hat[tt]
    if(tt!=length(X)){
      derror[tt+1,]<-c(1,error[tt],X[tt],0,0,0)+c(aa*derror[tt,1:3],0,0,0)
      dSigma2hat[tt+1,]<-c(alpha*2*error[tt]*derror[tt,1:3],
                           1,error[tt]^2,Sigma2hat[tt])+beta*dSigma2hat[tt,]
      error[tt+1]<-Xhat[tt+1]-X[tt+1]
      
      gr<-gr+(Sigma2hat[tt+1]*2*error[tt+1]*derror[tt+1,]-
                error[tt+1]^2*dSigma2hat[tt+1,])/
        Sigma2hat[tt+1]^2+(1/Sigma2hat[tt+1])*dSigma2hat[tt+1,]
    }
  }
  
  # gr<-apply((head(Sigma2hat,-1)*2*error*derror-error^2*dSigma2hat)/
  #             (head(Sigma2hat,-1))^2+1/(head(Sigma2hat,-1))*dSigma2hat,2,sum)

  return(gr)
}

#Example:
#loading return data#############
X<-diff(log(EuStockMarkets[,1]))

#Optimizing numerically#############
init.par<-c(0,0,0,0.0001,0.05,0.90)
names(init.par)<-c("mu","a","b","omega","alpha","beta")
# opt<-optim(init.par,log.likelihood,X=X)
# opt.par<-opt$par
# print(opt)
# AG<-ARMA11GARCH11(X,opt.par)
# plot(time(X),as.numeric(X),type="l",ylab="X",xlab="date")
# points(time(X),head(AG$Xhat,-1)+1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")
# points(time(X),head(AG$Xhat,-1)-1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")

#Optimizing by gradient descent###########
opt.gd.R<-function(X, n.its,init.par){
  pars<-matrix(0,ncol=6,nrow=n.its+1)
  colnames(pars)<-c("mu","a","b","omega","alpha","beta")
  lls<-numeric(n.its+1)
  pars[1,]<-init.par
  #lls[1]<-log.likelihood(X,pars[1,])
  for(ii in 2:(n.its+1)){
    gr<-ll.gr.R(X,pars[ii-1,])
    gr<-pmax(-10^(5),pmin(gr,10^(5)))#clipping
    eta=c(rep(10^(-10),3),10^(-11),rep(10^(-6),2))
    pars[ii,]<-pars[ii-1,]-eta*gr
    pars[ii,4:6]<-pmax(pars[ii,4:6],0.000000001)
    #lls[ii]<-log.likelihood(X,pars[ii,])
  }
  return(list(pars=pars,lls=lls))
}

opt2<-opt.gd.R(X, 100,init.par)
pars<-opt2$pars
lls<-opt2$lls

print(tail(pars,1))
#plot(lls,type="l")
AG<-ARMA11GARCH11.R(X,tail(pars,1))
plot(time(X),as.numeric(X),type="l",ylab="X",xlab="date")
points(time(X),head(AG$Xhat,-1)+1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")
points(time(X),head(AG$Xhat,-1)-1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")

