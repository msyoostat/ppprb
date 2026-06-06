
library(parallel)
library(tictoc)
library(ggplot2)

set.seed(10) 
#1989-10-18 01:04 UTC
# ==============================
# load data
# ==============================
load("hawks_process_mh.RData")


# PPRB --------------------------------------------------------------------

## 1st stage
mu_post_1st<-numeric(iterations)
beta_post_1st<-numeric(iterations)
eta_post_1st<-numeric(iterations)

mu_post_1st[1]<-1
beta_post_1st[1]<-1
eta_post_1st[1]<-0.5

accept_ratio_1st<-numeric(iterations)

i<-2
mu_jump_1st<-mu_jump
beta_jump_1st<-beta_jump
eta_jump_1st<-eta_jump
tic()
for(i in 2:iterations){
  mu_current<-mu_post_1st[i-1]
  beta_current<-beta_post_1st[i-1]
  eta_current<-eta_post_1st[i-1]
  alpha_current<-eta_current*beta_current
  #
  mu_candidate<- rnorm(1,mu_current,
                       sd=sqrt(mu_jump_1st))
  beta_candidate<- rnorm(1,beta_current,
                         sd=sqrt(beta_jump_1st))
  eta_candidate<- rnorm(1,eta_current,sd=sqrt(eta_jump_1st))
  alpha_candidate<- eta_candidate * beta_candidate
  if(mu_candidate>0 &beta_candidate>0 &eta_candidate>0 &eta_candidate<1){
    # likelihood and prior
    nom<- hawkes_loglik_naive(events_1, mu_candidate,
                              alpha_candidate,
                              beta_candidate,
                              T1)
    denom<-hawkes_loglik_naive(events_1, mu_current,
                               alpha_current,
                               beta_current,
                               T1)
    ## add prior 
    nom<- nom+ dgamma(mu_candidate,1,1,log=TRUE)+
      dgamma(beta_candidate,shape=2, 
             rate=0.5,log=TRUE)+
      dbeta(eta_candidate,2,2,log=TRUE)
    denom<- denom+ dgamma(mu_current,1,1,log=TRUE)+
      dgamma(beta_current,shape=2,
             rate=0.5,log=TRUE)+
      dbeta(eta_current,2,2,log=TRUE)
    # ratio
    ratio<-nom-denom 
    if(log(runif(1)) <= ratio){
      beta_post_1st[i]<-beta_candidate
      eta_post_1st[i]<-eta_candidate
      mu_post_1st[i]<-mu_candidate
      accept_ratio_1st[i]<-1
    } else{
      beta_post_1st[i]<-beta_current
      eta_post_1st[i]<-eta_current
      mu_post_1st[i]<-mu_current
    }
  } else{
    beta_post_1st[i]<-beta_current
    eta_post_1st[i]<-eta_current
    mu_post_1st[i]<-mu_current
  }
}
toc()


accept_ratio_1st<-accept_ratio_1st[-c(1:burnin)]
mean(accept_ratio_1st)

beta_post_1st<-beta_post_1st[-c(1:burnin)]

eta_post_1st<-eta_post_1st[-c(1:burnin)]
mu_post_1st<-mu_post_1st[-c(1:burnin)]

# parallel computing between stages ---------------------------------------


hawkes_loglik_conditional_fast <- function(times, mu, alpha, beta, T1, 
                                           T2,capital_T) {
  
  times <- sort(times)
  
  y1 <- times[times <= T1]
  y2 <- times[times >  T1 & times <=T2]
  
  N2 <- length(y2)
  if (N2 == 0) return(0)
  
  # -------------------------
  # 1. log intensity term (O(N))
  # -------------------------
  
  loglik <- 0
  
  # 
  R <- if (length(y1) > 0) {
    sum(exp(-beta * (T1 - y1)))
  } else {
    0
  }
  
  prev_time <- T1
  
  for (i in seq_len(N2)) {
    ti <- y2[i]
    
    # decay
    R <- exp(-beta * (ti - prev_time)) * R
    
    intensity_i <- mu + alpha * R
    if (intensity_i <= 0) return(-Inf)
    
    loglik <- loglik + log(intensity_i)
    
    # 
    R <- R + 1
    prev_time <- ti
  }
  
  # -------------------------
  # 2. compensator
  # -------------------------
  
  compensator <- mu * (capital_T - T1)
  
  if (length(y1) > 0) {
    compensator <- compensator +
      (alpha / beta) * sum(
        exp(-beta * (T1 - y1)) - exp(-beta * (capital_T - y1))
      )
  }
  
  if (length(y2) > 0) {
    compensator <- compensator +
      (alpha / beta) * sum(
        1 - exp(-beta * (capital_T - y2))
      )
  }
  
  loglik - compensator
}

tic()
log_likelihood_pre <- mclapply(
  X = seq_along(mu_post_1st),
  FUN = function(i) {
    hawkes_loglik_conditional_fast(
      times=data$t,
      mu=mu_post_1st[i],
      alpha=(eta_post_1st*beta_post_1st)[i],
      beta=beta_post_1st[i],
      T1,
      T2,
      capital_T=T2
    )
  },
  mc.cores = detectCores()-1
)
log_likelihood_pre <- unlist(log_likelihood_pre)
toc()

## 2nd stage


mu_post_2nd<-numeric(iterations)
beta_post_2nd<-numeric(iterations)
eta_post_2nd<-numeric(iterations)

mu_post_2nd[1]<-mu_post_1st[100]
beta_post_2nd[1]<-beta_post_1st[100]
eta_post_2nd[1]<-eta_post_1st[100]

accept_ratio_2nd<-numeric(iterations)

i<-2

current_idx<-100
tic()
for(i in 2:iterations){
  idx<- sample(1:length(mu_post_1st),1,replace=TRUE )
  # 
  beta_candidate<- beta_post_1st[idx]
  mu_candidate<- mu_post_1st[idx]
  eta_candidate<- eta_post_1st[idx]
  # #
  beta_current<- beta_post_2nd[i-1]
  mu_current<- mu_post_2nd[i-1]
  eta_current<- eta_post_2nd[i-1]
  
  nom<- log_likelihood_pre[idx]
  denom<- log_likelihood_pre[current_idx]
  # ratio
  ratio<-nom-denom 
  if(log(runif(1)) <= ratio){
    beta_post_2nd[i]<-beta_candidate
    mu_post_2nd[i]<-mu_candidate
    eta_post_2nd[i]<-eta_candidate
    accept_ratio_2nd[i]<-1
    current_idx<-idx
  } else{
    beta_post_2nd[i]<-beta_current
    mu_post_2nd[i]<-mu_current
    eta_post_2nd[i]<-eta_current
  }
}
toc()
beta_post_2nd<-beta_post_2nd[-c(1:burnin)]
mu_post_2nd<-mu_post_2nd[-c(1:burnin)]
eta_post_2nd<-eta_post_2nd[-c(1:burnin)]
accept_ratio_2nd<-accept_ratio_2nd[-c(1:burnin)]


mean(accept_ratio_2nd)



# parallel computing between stages (2 and 3) ---------------------------------------
tic()
log_likelihood_pre <- mclapply(
  X = seq_along(mu_post_2nd),
  FUN = function(i) {
    hawkes_loglik_conditional_fast(
      times=data$t,
      mu=mu_post_2nd[i],
      alpha=(eta_post_2nd*beta_post_2nd)[i],
      beta=beta_post_2nd[i],
      T2,
      capital_T,
      capital_T=capital_T
    )
  },
  mc.cores = detectCores()-1
)
log_likelihood_pre <- unlist(log_likelihood_pre)
toc()

## 3rd stage


mu_post_3rd<-numeric(iterations)
beta_post_3rd<-numeric(iterations)
eta_post_3rd<-numeric(iterations)

mu_post_3rd[1]<-mu_post_2nd[100]
beta_post_3rd[1]<-beta_post_2nd[100]
eta_post_3rd[1]<-eta_post_2nd[100]

accept_ratio_3rd<-numeric(iterations)

i<-2

current_idx<-100
tic()
for(i in 2:iterations){
  idx<- sample(1:length(mu_post_2nd),1,replace=TRUE )
  # 
  beta_candidate<- beta_post_2nd[idx]
  mu_candidate<- mu_post_2nd[idx]
  eta_candidate<- eta_post_2nd[idx]
  # #
  beta_current<- beta_post_3rd[i-1]
  mu_current<- mu_post_3rd[i-1]
  eta_current<- eta_post_3rd[i-1]

  
  nom<- log_likelihood_pre[idx]
  denom<- log_likelihood_pre[current_idx]
  # ratio
  ratio<-nom-denom 
  if(log(runif(1)) <= ratio){
    beta_post_3rd[i]<-beta_candidate
    mu_post_3rd[i]<-mu_candidate
    eta_post_3rd[i]<-eta_candidate
    accept_ratio_3rd[i]<-1
    current_idx<-idx
  } else{
    beta_post_3rd[i]<-beta_current
    mu_post_3rd[i]<-mu_current
    eta_post_3rd[i]<-eta_current
  }
}
toc()
beta_post_3rd<-beta_post_3rd[-c(1:burnin)]
mu_post_3rd<-mu_post_3rd[-c(1:burnin)]
eta_post_3rd<-eta_post_3rd[-c(1:burnin)]
accept_ratio_3rd<-accept_ratio_3rd[-c(1:burnin)]

mean(accept_ratio_3rd)

save.image("hawks_process_pprb.RData")

