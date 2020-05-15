## Binomial-Poisson-Beta hybrid model code
library(readr)
library(dplyr)
library(truncnorm)
library(MASS)
library(tmvtnorm)
library(mvtnorm)
library(janitor)
library(LearnBayes)
set.seed(11)
library(tmvtnorm)
library(LearnBayes)
library(plyr)
library(dplyr)
library(janitor)
library(readr)
setwd("~/Desktop/Stat 207/Midterm/")
data<-read.csv("covid041320.txt")%>%clean_names()

ny<-data$total_cases
y<-data$deaths
cy<-data$population/10^3 #population in thousands
county<-data$county
N<-length(ny)

trace = function(colnum,theta.store){
  plot(theta.store[,colnum],type="l")
}
# view histogram of a column of theta.store
viewhist = function(colnum,theta.store,true.value=NULL){
  hist(theta.store[,colnum])
  if(!is.null(true.value))
    abline(v=true.value,lwd=2,col="red")
}

binom_pois_gamma_hier_mcmc = function(n.iter, burn.in, hypers, inits, var.tune.alpha, data){
  n.params = length(inits)
  alpha = inits[1] # mu
  beta = inits[2] # tau
  lambda = inits[3]
  a = hypers[1] # a
  b = hypers[2] # b
  theta = tail(inits,n.params-3) # theta
  ny = data$total_cases
  N = length(ny)
  y = data$deaths
  cy = data$population/10^3
  
  theta.store = matrix(NA,n.iter,n.params)
  accept = 0
  accept2 = 0
  
# modify this for alpha (mu), make another for beta (tau)
  cond.alpha <- function(alpha, beta, theta ){
    alpha <- exp(alpha)/(1+exp(alpha)) # this maps alpha to the real line
    beta <- exp(beta) # this maps tau to the real line
    sum(lchoose(ny,y)+(y+alpha*beta)*log(theta) + (ny-y+beta*(1-alpha))*log(1-theta) -theta*lambda*cy -ny*log(cy) - lfactorial(ny)) -
      N*lbeta(alpha*beta, beta*(1-alpha)) - log(alpha*(1-alpha)*(1+beta)^2) - log(alpha*(1-alpha)) - log(beta)
        # Jacobian is the last term
  }
  
  cond.theta <- function(alpha, beta, theta){
    theta <- exp(theta)/(1+exp(theta)) #exp(theta) # this maps alpha to the real line need a different transformation and Jacobian
    sum(lchoose(ny,y)+(y+alpha*beta)*log(theta) + (ny-y+beta*(1-alpha))*log(1-theta) -theta*lambda*cy -ny*log(cy) - lfactorial(ny) -
          log(theta*(1-theta))) # this last term is the jacobian
  }
  
  # main gibbs loop
  for(i in 1:n.iter){
    # gibbs for poisson param c*l
    lambda = rgamma(1,a+N,b+sum(theta*cy))
    theta.store[i,1] = lambda

    # MH for alpha
    alpha.prop = rnorm(1,log(alpha),var.alpha)
    beta.prop = rnorm(1,log(beta),var.beta) #log(beta)
    acceptance.prob = min(exp( cond.alpha(alpha.prop,beta.prop,theta) - cond.alpha(log(alpha),log(beta),theta) )  , 1) #log(beta)
    u = runif(1)
    if(u < acceptance.prob){
      accept = accept + 1
      alpha = exp(alpha.prop)
      beta = exp(beta.prop)
    }
    theta.store[i,2] = alpha# mu
    theta.store[i,3] = beta #tau
    
    # MH for theta vector
    theta.prop = rtruncnorm(N,log(theta),var.theta,a=-Inf,b=0)
    acceptance.prob = min(exp( cond.theta(alpha,beta,theta.prop) - cond.theta(alpha,beta,log(theta)) )  , 1)
    u = runif(1)
    if(u < acceptance.prob){
      accept2 = accept2 + 1
      theta = exp(theta.prop)
    }
    theta.store[i,4:n.params] = theta# mu
  }
  
  #if((accept/n.iter)<.3 | (accept/n.iter)>.5) { stop("Bad acceptance ratio",accept/n.iter)}
  
  theta.post.burn = theta.store[(burn.in+1):n.iter,]
  return(list(theta.full = theta.store, theta.post.burn = theta.post.burn, accept.ratio = accept/n.iter, accept.ratio2 = accept2/n.iter))
}


##
{
b = 1
a = 5

alpha = 0.5
beta = 1
theta = rbeta(N,1/2,1/2)
lambda = 1

var.theta = .005
var.beta = .5
var.alpha = .5
n.iter = 20000
burn.in = 10000
#  function(n.iter, burn.in, hypers, inits, var.tune.alpha, data)
theta.gibbs2 = binom_pois_gamma_hier_mcmc(n.iter, burn.in, c(a,b), c(alpha,beta,lambda,theta), var.alpha, data)
theta.gibbs2
par(mfrow=c(3,1))
plot(theta.gibbs2$theta.post.burn[,1],type="l",ylab="lambda",main="Lambda MCMC (post burn-in)")
abline(h=mean(theta.gibbs2$theta.post.burn[,1]),col="red")
plot(theta.gibbs2$theta.post.burn[,2],type="l",ylab="alpha",main="Alpha MCMC (post burn-in)")
abline(h=mean(theta.gibbs2$theta.post.burn[,2]),col="red")
plot(theta.gibbs2$theta.post.burn[,3],type="l",ylab="beta",main="Beta MCMC (post burn-in)")
abline(h=mean(theta.gibbs2$theta.post.burn[,3]),col="red")
}
####################################

#####################################
par(mfrow=c(3,1))
for(i in 1:3){
  plot(theta.gibbs2$theta.post.burn[,i],type="l")
  abline(h=mean(theta.gibbs2$theta.post.burn[,i]),col="red")
}


par(mfrow=c(3,3))
for(i in 1:9){
trace(12+i,theta.gibbs2$theta.post.burn)
}

par(mfrow=c(2,2))
for(i in 1:4){
  trace(46+i,theta.gibbs2$theta.post.burn)
}


dev.off()
par(mfrow=c(1,2))
boxplot(log(theta.gibbs2$theta.post.burn[,4:60]),use.cols = TRUE)
boxplot(theta.gibbs2$theta.post.burn[,4:60],use.cols = TRUE,
        ylab="theta",xlab="county index",
        main="MCMC theta samples (post burn-in)")
# many bugs, but it runs. are the predicted death counts reasonable?
###################################

pred.dead.5<-matrix(NA,nrow=nrow(theta.gibbs2$theta.post.burn),ncol=nrow(data))
for(i in 1:nrow(pred.dead.5)){
  for(j in 1:ncol(pred.dead.5)){
    pred.dead.5[i,j]=rbinom(n=1,size=rpois(n=1,theta.gibbs2$theta.post.burn[i,1]*theta.gibbs2$theta.post.burn[i,(j+3)]*cy[j]),
                            prob=theta.gibbs2$theta.post.burn[i,(j+3)])
  }
}
dev.off()
boxplot(pred.dead.5,use.cols = TRUE,main="Posterior Predicted Statewide Deaths", 
        xlab="county index",ylab="deaths (1k)")
hist(rowSums(pred.dead.5,na.rm=TRUE),main="Posterior Predicted Statewide Deaths", 
     xlab="1K Deaths",breaks=40)
abline(v=mean(rowSums(pred.dead.5,na.rm=TRUE)),col="red")

pred.cases<-matrix(NA,nrow=nrow(theta.gibbs2$theta.post.burn),ncol=nrow(data))
for(i in 1:nrow(pred.cases)){
  for(j in 1:ncol(pred.cases)){
    pred.cases[i,j]=rpois(n=1,theta.gibbs2$theta.post.burn[i,1]*theta.gibbs2$theta.post.burn[i,(j+3)]*cy[j])
  }
}

par(mfrow=c(2,1))
hist(rowSums(pred.cases,na.rm=TRUE),main="Posterior Predicted Statewide Cases", 
     xlab="1K Cases",breaks=40)
abline(v=mean(rowSums(pred.cases,na.rm=TRUE)),col="red")
hist(rowSums(pred.dead.5,na.rm=TRUE),main="Posterior Predicted Statewide Deaths", 
     xlab="1K Deaths",breaks=40)
abline(v=mean(rowSums(pred.dead.5,na.rm=TRUE)),col="red")
