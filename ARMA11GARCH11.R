#The ARMA(1,1)-GARCH(1,1) specification. 
#Inputs are the stock returns, and the parameters, c(mu,a,b,omega,alpha,beta)
ARMA11GARCH11<-function(X,params){
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
  AG<-ARMA11GARCH11(X,params)
  ll<-sum(AG$error^2/head(AG$Sigma2hat,-1)+log(head(AG$Sigma2hat,-1)))
  return(ll)
}

##Computes the gradient of the log-likelihood.
#Inputs are the stock returns, and the parameters, c(mu,a,b,omega,alpha,beta)
ll.gr<-function(X,params){
  X<-as.numeric(X)
  mu<-params[1]
  aa<-params[2]
  bb<-params[3]
  omega<-params[4]
  alpha<-params[5]
  beta<-params[6]
  
  AG<-ARMA11GARCH11(X,params)
  error<-AG$error
  Sigma2hat<-AG$Sigma2hat
  Xhat<-AG$Xhat
  
  derror<-matrix(data=0,nrow=length(error),ncol=6)
  dSigma2hat<-matrix(data=0,nrow=length(error),ncol=6)
  
  for(tt in 1:(nrow(derror)-1)){
    derror[tt+1,]<-c(1,error[tt],X[tt],0,0,0)+c(aa*derror[tt,1:3],0,0,0)
    dSigma2hat[tt+1,]<-c(alpha*2*error[tt]*derror[tt,1:3],
                         1,error[tt]^2,Sigma2hat[tt])+beta*dSigma2hat[tt,]
  }
  
  gr<-apply((head(Sigma2hat,-1)*2*error*derror-error^2*dSigma2hat)/
              (head(Sigma2hat,-1))^2+1/(head(Sigma2hat,-1))*dSigma2hat,2,sum)
  
  return(gr)
}

#Example:
#loading return data#############
X<-diff(log(EuStockMarkets[,1]))

#Optimizing numerically#############
init.par<-c(0,0,0,0.0001,0.05,0.90)
names(init.par)<-c("mu","a","b","omega","alpha","beta")
opt<-optim(init.par,log.likelihood,X=X)
opt.par<-opt$par
print(opt)
AG<-ARMA11GARCH11(X,opt.par)
plot(time(X),as.numeric(X),type="l",ylab="X",xlab="date")
points(time(X),head(AG$Xhat,-1)+1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")
points(time(X),head(AG$Xhat,-1)-1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")

#Optimizing by gradient descent###########
n.its<-100
pars<-matrix(0,ncol=6,nrow=n.its+1)
colnames(pars)<-c("mu","a","b","omega","alpha","beta")
lls<-numeric(n.its+1)
pars[1,]<-init.par
lls[1]<-log.likelihood(X,pars[1,])
for(ii in 2:(n.its+1)){
  gr<-ll.gr(X,pars[ii-1,])
  gr<-pmax(-10^(1),pmin(gr,10^(1)))#clipping for stability
  eta=c(rep(10^(-5),3),10^(-7),rep(10^(-3),2))
  pars[ii,]<-pars[ii-1,]-eta*gr
  pars[ii,4:6]<-pmax(pars[ii,4:6],0.000000001)
  lls[ii]<-log.likelihood(X,pars[ii,])
}
print(tail(pars,1))
plot(lls,type="l")
AG<-ARMA11GARCH11(X,tail(pars,1))
plot(time(X),as.numeric(X),type="l",ylab="X",xlab="date",main="gradient descent")
points(time(X),head(AG$Xhat,-1)+1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")
points(time(X),head(AG$Xhat,-1)-1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")


