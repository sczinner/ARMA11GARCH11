library(Rcpp)
library(rbenchmark)

sourceCpp(file="ARMA11GARCH11-cpp.cpp")

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
 return(llgr(X,params))
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
opt.gd<-function(X, n.its,init.par){
  ret<-optgd(X, n.its,init.par)
  colnames(ret[["pars"]])<-c("mu","a","b","omega","alpha","beta")
  return(ret)
}

opt2<-opt.gd(X, 100,init.par)
pars<-opt2$pars
lls<-opt2$lls

print(tail(pars,1))
plot(lls,type="l")
AG<-ARMA11GARCH11(X,tail(pars,1))
plot(time(X),as.numeric(X),type="l",ylab="X",xlab="date")
points(time(X),head(AG$Xhat,-1)+1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")
points(time(X),head(AG$Xhat,-1)-1.96*sqrt(head(AG$Sigma2hat,-1)),type="l",col="red")

source("ARMA11GARCH11.R")

benchmark("R"=opt.gd.R(X, 100,init.par),"C++"=opt.gd(X, 100,init.par),replications=5)
#Get speed up of 22x
