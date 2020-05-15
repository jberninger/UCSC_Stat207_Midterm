## Poisson Gamma model code
## Stats 207
## Take Home Analysis 2
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

################
## 1 # describe variability of cases vs population
################
ggplot(data=data,aes(x=population/1000,y=log(total_cases))) + geom_point() +
  xlab("Population (1k)") + 
  ylab("log(cases)") + ggtitle("log(Cases) vs Population")
ggplot(data=data,aes(x=total_cases/(population/1000))) + geom_histogram(fill=NA) +
  xlab("Cases per 1K") + ggtitle("Infection Rate by County")# tail
ggplot(data=data,aes(x=log(total_cases/population))) + geom_histogram() # log transformation looks more symmetric

hist(ny/cy,breaks=15,xlab="Cases per 1k",main="Infection Rate")
abline(v=mean(ny/cy),col="red")
plot(cy,ny,xlab="Population (1K)",ylab="Case Count",main="Cases vs Population")
plot(cy,log(ny),xlab="Population (1K)",ylab="log(cases)",main="log(Cases) vs Population")
plot(ny,log(cy))
plot(ny,cy); cor(ny,cy) # n,c have a positive linear correlation with increasing variance

################
## 2 #
################
## Set the values of a1, b1, a2, b2 to be compatible with the assumption that, for a given county, the mean number of cases is equal to 20% of the population
## this means we want the E(lambda_i) = .2

b1<-1
a1<-2*b1
b2<-1
a2<-10*b2
x<-seq(from=0.0001,to=5,length.out=100)

# sample lambda
rgamma(n=1000,shape=rgamma(n=1000,a1,b1),rate=rgamma(n=1000,a2,b2)) %>% mean() # this is close to .2

rgamma(n=1000,shape=rgamma(n=1000,a1,b1),rate=rgamma(n=1000,a2,b2)) %>% density() %>% plot()

###################################
# I dont think I need this density:

#################################
## Grant's code for the poisson gamma model from question 15
{pois_gamma_hier_gibbs = function(n.iter, burn.in, hypers, inits, var.tune.alpha, data){
  n.params = length(inits)
  alpha = inits[1]
  beta = inits[2]
  a = hypers[1]
  b = hypers[2]
  c = hypers[3]
  d = hypers[4]
  lambda = tail(inits,n.params-2)
  ny = data$total_cases
  N = length(ny)
  y = data$deaths
  cy = data$population/10^3
  
  theta.store = matrix(NA,n.iter,n.params)
  accept = 0
  
  # full conditional for alpha ( from brunos notes )
  # cond.alpha = function(alpha,beta,lamda){
  #  alpha = exp(alpha)
  #  (a - 1)*(log(alpha))  - lgamma(alpha) - alpha*(sum(log(lambda))-N*log(beta)+b) + log(alpha)
  # }
  # full conditional for alpha from kelsey
  cond.alpha <- function(alpha, beta, cy.ly){
    alpha <- exp(alpha)
    (a - 1)*(log(alpha)) - b*alpha + sum((ny + alpha - 1)*log(cy.ly)) - N*lgamma(alpha) + log(alpha)
  }
  
  
  # main gibbs loop
  for(i in 1:n.iter){
    # gibbs for poisson param c*l
    lambda = rgamma(N,alpha+(ny/cy),beta+1); theta.store[i,(3:n.params)] = lambda
    # gibbs for beta
    beta = rgamma(1,N*alpha+c,sum(lambda)+d); theta.store[i,2] = beta
    # MH for alpha
    alpha.prop = rnorm(1,log(alpha),var.alpha)
    acceptance.prob = min(exp( cond.alpha(alpha.prop,beta,lambda) - cond.alpha(log(alpha),beta,lambda) )  , 1)
    u = runif(1)
    if(u < acceptance.prob){
      accept = accept + 1
      alpha = exp(alpha.prop)
    }
    theta.store[i,1] = alpha
  }
  
  #if((accept/n.iter)<.3 | (accept/n.iter)>.5) { stop("Bad acceptance ratio",accept/n.iter)}
  
  theta.post.burn = theta.store[(burn.in+1):n.iter,]
  return(list(theta.full = theta.store, theta.post.burn = theta.post.burn, accept.ratio = accept/n.iter))
}
# view mcmc trace of for a column of theta.store
trace = function(colnum,theta.store){
  plot(theta.store[,colnum],type="l")
}
# view histogram of a column of theta.store
viewhist = function(colnum,theta.store,true.value=NULL){
  hist(theta.store[,colnum])
  if(!is.null(true.value))
    abline(v=true.value,lwd=2,col="red")
}
}

# MAIN

# informative hyperpriors
#a = .001
#b = .001
#c = .001
#d = .001
b = 0.00001
a = 2*b
d = 0.00001
c = 10*d

alpha = 1
beta = 1
lambda = rgamma(N,alpha+(ny/cy),beta+1)
var.alpha = .15
n.iter = 20000
burn.in = 10000
theta.gibbs = pois_gamma_hier_gibbs(n.iter, burn.in, c(a,b,c,d), c(alpha,beta,lambda), var.alpha, data)

# do samples look like data?
colMeans(theta.gibbs$theta.post.burn)[3:(N+2)]*cy

# chains
par(mfrow=c(2,1))
plot(theta.gibbs$theta.post.burn[,1],type="l",
      ylab="alpha",main="Alpha MCMC (post burn-in)")
abline(h=mean(theta.gibbs$theta.post.burn[,1]),
       col="red")# alpha
plot(theta.gibbs$theta.post.burn[,2],type="l",
      ylab="beta",main="Beta MCMC (post burn-in)") # beta
abline(h=mean(theta.gibbs$theta.post.burn[,2]),
       col="red")
trace(3,theta.gibbs$theta.post.burn) # lambda1

plot(theta.gibbs$theta.post.burn[,1], 
     type="l",ylab="alpha",main="Alpha Chain (post burn-in)")
abline(h=mean(theta.gibbs$theta.post.burn[,1]),
       col="red")# alpha)

plot(theta.gibbs$theta.post.burn[,2], 
     type="l",ylab="beta",main="Beta Chain (post burn-in)")
abline(h=mean(theta.gibbs$theta.post.burn[,2]),
       col="red")# alpha)


dev.off()
par(mfrow=c(3,3))
for(i in 1:9){
plot(theta.gibbs$theta.post.burn[,12+i],type="l",main="")
}

# posterior prob of more than 450 hawks in years with route number 120
# l~Ga(alpha,beta)
post.samp.l = rgamma(n.iter-burn.in, theta.gibbs$theta.post.burn[,1],theta.gibbs$theta.post.burn[,2])
# sample of lambdas from alpha, beta

pred.store=matrix(NA,nrow(theta.gibbs$theta.post.burn),ncol=58)
for(i in 1:nrow(pred.store)){
  for(j in 1:ncol(pred.store)){
    pred.store[i,j]=rpois(n=1,cy[j]*theta.gibbs$theta.post.burn[i,(2+j)])
  }
}


# boxplots for lambdas !!!!!!!

boxplot(theta.gibbs$theta.post.burn[,3:60], use.cols=TRUE,main="Lambda Chains (post burn-in)",
        ylab="lambda",xlab="county index") # number of cases in each county
boxplot(log(theta.gibbs$theta.post.burn[,3:60]), use.cols=TRUE,main="Lambda Chains (post burn-in)",
        ylab="log(lambda)",xlab="county index")

# Posterior predictive deaths by state !!!!!!!!!!!!!
boxplot(pred.store,use.cols=TRUE,main="Posterior Predicted Cases",
        ylab="cases (1k)",xlab="county index") # correct the units

boxplot(log(pred.store),use.cols=TRUE,main="Posterior Predicted Cases",
        ylab="log(cases (1k))",xlab="county index") # correct the units

hist(rowSums(pred.store),xlab="cases (1k)",main="Posterior Estimated Statewide Cases")
abline(v=mean(rowSums(pred.store)),col="red")
boxplot(log(pred.store),use.cols=TRUE)

# posterior predictive statewide deaths !!!!!!!!!!!!!!
hist(rowSums(pred.store),main="Posterior Predicted Statewide Deaths", 
     xlab="1K Deaths",breaks=40)
abline(v=mean(rowSums(pred.store)),col="red")
mean(rowSums(pred.store))
# 6824 almost 700k predicted deaths is the posterior average?
abline(v=quantile(rowSums(pred.store),prob=0.05),col="blue")
abline(v=quantile(rowSums(pred.store),prob=0.95),col="blue")
quantile(rowSums(pred.store),probs=c(0.05,0.50,0.95))

##############################################################
## Model 3
par(mfrow=c(3,3))
trace(3,theta.gibbs$theta.post.burn) # lambda1
trace(4,theta.gibbs$theta.post.burn) # lambda1
trace(5,theta.gibbs$theta.post.burn) # lambda1
trace(6,theta.gibbs$theta.post.burn) # lambda1
trace(7,theta.gibbs$theta.post.burn) # lambda1
trace(8,theta.gibbs$theta.post.burn) # lambda1
trace(9,theta.gibbs$theta.post.burn) # lambda1
trace(10,theta.gibbs$theta.post.burn) # lambda1
trace(11,theta.gibbs$theta.post.burn) # lambda1

