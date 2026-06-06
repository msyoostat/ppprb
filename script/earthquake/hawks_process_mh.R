library(parallel)
library(tictoc)
library(ggplot2)

set.seed(10) 
#1989-10-18 01:04 UTC
# ==============================
# load data
# ==============================
load("world_series_earthquake_2_2.RData")
T1=400
T2=500

events_1 <- data[data$t <= T1,3]
events_2 <- data[data$t > T1 & data$t <= T2,3]
events_3 <- data[data$t > T2 ,3]
capital_T<- max(data$t)

# Bayesian ----------------------------------------------------------------
iterations<-30000
burnin<-5000
a<-0.1
b<-0.1

mu_post<-numeric(iterations)
eta_post<-numeric(iterations)
beta_post<-numeric(iterations)

mu_post[1]<-1
eta_post[1]<-0.5
beta_post[1]<-1

mu_jump<-0.004
eta_jump<-0.004
beta_jump<-1 # 0.8


accept_ratio<-numeric(iterations)

hawkes_loglik_naive <- function(times, mu, alpha, beta, capital_T) {
  N <- length(times)
  
  loglik <- 0
  
  for (i in 1:N) {
    intensity_i <- mu
    if (i > 1) {
      intensity_i <- intensity_i +
        sum(alpha * exp(-beta * (times[i] - times[1:(i-1)])))
    }
    loglik <- loglik + log(intensity_i)
  }
  
  compensator <- mu * capital_T +
    (alpha / beta) * sum(1 - exp(-beta * (capital_T - times)))
  
  loglik - compensator
}

tic()
for(i in 2:iterations){
  mu_current<-mu_post[i-1]
  beta_current<-beta_post[i-1]
  eta_current<-eta_post[i-1]
  alpha_current<-eta_current*beta_current
  #
  mu_candidate<- rnorm(1,mu_current,
                       sd=sqrt(mu_jump))
  beta_candidate<- rnorm(1,beta_current,
                           sd=sqrt(beta_jump))
  eta_candidate<- rnorm(1,eta_current,sd=sqrt(eta_jump))
  alpha_candidate<- eta_candidate * beta_candidate
  if(mu_candidate>0 &beta_candidate>0 &eta_candidate>0 &eta_candidate<1){
      # likelihood and prior
    nom<- hawkes_loglik_naive(data$t, mu_candidate,
                              alpha_candidate,
                              beta_candidate,
                              capital_T)
    denom<-hawkes_loglik_naive(data$t, mu_current,
                               alpha_current,
                               beta_current,
                               capital_T)
    ## add prior 
    nom<- nom+ dgamma(mu_candidate,1,rate=1,log=TRUE)+
      dgamma(beta_candidate,shape=2, 
             rate=0.5,log=TRUE)+
      dbeta(eta_candidate,2,2,log=TRUE)
    denom<- denom+ dgamma(mu_current,1,rate=1,log=TRUE)+
      dgamma(beta_current,shape=2,
             rate=0.5,log=TRUE)+
      dbeta(eta_current,2,2,log=TRUE)

    # ratio
    ratio<-nom-denom 
    if(log(runif(1)) <= ratio){
      beta_post[i]<-beta_candidate
      eta_post[i]<-eta_candidate
      mu_post[i]<-mu_candidate
      accept_ratio[i]<-1
    } else{
      beta_post[i]<-beta_current
     eta_post[i]<-eta_current
      mu_post[i]<-mu_current
    }
  } else{
    beta_post[i]<-beta_current
    eta_post[i]<-eta_current
    mu_post[i]<-mu_current
  }
  if (i %% 1000 == 0) cat(i, "\n")
}
toc()
accept_ratio<-accept_ratio[-c(1:burnin)]
mean(accept_ratio)
beta_post<-beta_post[-c(1:burnin)]
eta_post<-eta_post[-c(1:burnin)]
mu_post<-mu_post[-c(1:burnin)]

save.image("hawks_process_mh.RData")
