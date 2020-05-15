## Binomial-Beta Hierarchical model
## Model 2 in Midterm
## Model 3 in Take Home Quiz 1

## Stat 207 Quiz 1
## model 2 recovery
library(readr)
library(tidyverse)
library(ggmap)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(janitor)
library(LearnBayes)

setwd("~/Desktop/Stat 207/")
covid <- read_csv("covid.csv") %>% mutate(County = tolower(County)) %>% clean_names()
covid1 <- covid %>% mutate(infected = floor(0.20*population))
covid2 <- data.frame(y = covid$deaths, n = covid$total_cases)
covid3 <- covid1 %>% mutate(est.infected.20 = ceiling(0.20*population))

states <- map_data("state")
cali <- subset(states, region == "california")
counties <- map_data("county")
ca_county <- subset(counties, region == "california") %>% inner_join(covid, by = c("subregion" = "county")) %>%
  mutate(death.rate = 100*round(deaths/population, 10))

ditch_the_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank()
)

ca_base <- ggplot(data = cali, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "gray")

ca_base + 
  geom_polygon(data = ca_county, aes(fill =ordered(round(100*death.rate, 2))), color = "white") +
  geom_polygon(color = "black", fill = NA) +
  ggtitle("Death Rates by County") +
  theme_bw() +
  labs(fill = "Deaths / Infected (%)") +
  ditch_the_axes

par(mfrow = c(2,1))

ggplot(data = covid1, aes(x=infected, y=deaths, color=population)) + geom_point() + 
  geom_smooth(method="lm") + ggtitle("All California Counties")

ggplot(data = covid1[2:58,], aes(x=infected, y=deaths, color=population)) + geom_point() + 
  geom_smooth(method="lm") + ggtitle("Los Angeles Removed")

ggplot(data = covid1, aes(x = 100*deaths/total_cases, fill = county)) + geom_histogram() + 
  labs(x="Deaths / Total Cases (%)")

## first model
set.seed(112358)
n = 1000
post.theta <- rbeta(n = n, shape1 = sum(covid$deaths) + 0.5, 
                    shape2 = sum(covid$total_cases) - sum(covid$deaths) + 0.5)

hist(post.theta)
mean(post.theta)
sd(post.theta)
plot(density(post.theta), main = "Model 1 : Theta Posterior Distribution", col = "blue") # make pretty

est.posterior.deaths <- matrix(NA, nrow = nrow(covid), ncol = n)
row.names(est.posterior.deaths) <- covid$county

for(i in 1:nrow(covid)){
  for(j in 1:n){
    est.posterior.deaths[i,j] = ceiling(0.20*covid1$population[i] * post.theta[j])
  }
}

## second model
betabinexch0=function(theta, data) # this computes the log posterior density of the beta-binomial exchangeable model
{
  eta=theta[,1] # this is mu in our model
  K=theta[,2] # this is tau in our model
  y=data[,1] 
  n=data[,2]
  N=length(y)
  val=0*K
  for (i in 1:N)
    val=val+lbeta(K*eta+y[i],K*(1-eta)+n[i]-y[i])
  val=val-N*lbeta(K*eta,K*(1-eta))
  val=val-2*log(1+K)-log(eta)-log(1-eta)
  return(val)
}

betabinexch=function(theta, data){
  eta = exp(theta[1])/(1 + exp(theta[1]))
  K = exp(theta[2])
  y = data[, 1]
  n = data[, 2]
  N = length(y)
  
  logf=function(y,n,K,eta)
    lbeta(K * eta + y, K * (1 - eta) + n - y)-lbeta(K * eta, K * (1 - eta))
  
  val=sum(logf(y,n,K,eta))
  
  val = val + theta[2] - 2 * log(1 + exp(theta[2]))
  return(val)
}

betabinT=function(theta,datapar){
  data=datapar$data
  tpar=datapar$par
  d=betabinexch(theta,data)-dmt(theta,mean=c(tpar$m),
                                S=tpar$var,df=tpar$df,log=TRUE)
  return(d)
}

fit <- laplace(betabinexch,array(c(-7,6),c(1,2)), data = covid2)
tpar=list(m=fit$mode,var=2*fit$var,df=4)
dataparams=list(data=covid2,par=tpar)

start=array(c(-10,10),c(1,2))
#laplace(betabinT,theta = start, mode = 10, dataparams, method = "Brent", lower = -1000, upper= 100)

start1 <- c(0.25, 0.5)# these need to be initialized?

fit1 <- optim(start, betabinT, control = list(fnscale = -1), hessian = TRUE, datapar = dataparams)

d_max <- betabinT(theta = fit1$par, data = dataparams) # this is the value of d that I think we want

tpar1=list(m=fit1$par,var=2*solve(-1*fit1$hessian), df=4) 
# convert the hessian to a covariance, stack overflow says covar is asymptotically -hessian^-1

rejectsampling=function(logf,tpar,dmax,n,data){
  theta=rmt(n,mean=c(tpar$m),S=tpar$var,df=tpar$df)
  for(i in 1:n){
    lf[i]=logf(theta[i,],data) 
  }# this was lf before
  lg=dmt(theta,mean=c(tpar$m),S=tpar$var,df=tpar$df,log=TRUE)
  prob=exp(lf-lg-dmax)
  return(theta[runif(n)<prob,])
}

n <- 10000

lf <- rep(0,n)
theta3 <- rejectsampling(betabinexch, tpar1, dmax = d_max, n, covid2)
mycontour(betabinexch,c(-4,-3,3,8), covid2, main = "Model 2 : Parameter Rejection Sampling")# change these values as well add axes labels
points(theta3[,1],theta3[,2]) # these are mu and tau


library(extraDistr)
a <- inv.logit(theta3[,1])
b <- log(theta3[,2])

par(mfrow=c(2,1))
hist(a, main = "Model 2 : Posterior Distribution of mu",xlab="mu")
hist(b, main = "Model 2 : Posterior Distribution of tau",xlab="tau")

m2.post.y<-function(theta,data){
  p<-matrix(NA,nrow=length(theta[[1]]),ncol=nrow(covid))
  for(i in 1:nrow(p)){
    for(j in 1:ncol(p)){
      p[i,j]=rbbinom(n=1,size=as.numeric(covid3[j,6]),alpha=theta[[1]],beta=theta[[2]])
      # assuming 20% in each county become infected
    }
  }
  return(p)
}

m2.post.deaths<-m2.post.y(theta=list(a,b),data=covid)
mean(m2.post.deaths)

boxplot(m2.post.deaths/1000,use.cols = TRUE,main="Posterior Predicted Statewide Deaths", 
        xlab="county index",ylab="deaths (1k)")
hist(rowSums(m2.post.deaths/1000),main="Posterior Predicted Statewide Deaths", 
     xlab="1K Deaths",breaks=40)
abline(v=mean(rowSums(m2.post.deaths)/1000),col="red")
#hist(a/(a+b))
