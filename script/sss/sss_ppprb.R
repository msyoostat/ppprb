library(ncdf4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(invgamma)
library(truncnorm)
library(parallel)
library("geosphere")
library(tictoc)
library(Rcpp.conditional.gp)


RNGkind("L'Ecuyer-CMRG")
set.seed(2)



load("sss_pprb.RData")


# PPPRB -------------------------------------------------------------------
## let's do ppprb


PPPRB_1st_stage<-function(data_1,
                          dist_mat_1,
                          temperature,
                          iterations,
                          burnin,
                          jump_sigma_s2,
                          jump_sigma_n2,
                          jump_phi,
                          mu_log_sigma_s,
                          sd_log_sigma_s,
                          mu_log_sigma_n,
                          sd_log_sigma_n,
                          mu_log_phi,
                          sd_log_phi
                          ){
  log_sigma_s2_post_pprb_1 <- numeric(iterations)
  log_sigma_n2_post_pprb_1 <- numeric(iterations)
  log_phi_post_pprb_1 <- numeric(iterations)
  # Initial values
  log_sigma_s2_post_pprb_1[1] <- log(0.6)
  log_sigma_n2_post_pprb_1[1] <- log(0.2)
  log_phi_post_pprb_1[1] <- log(median(dist_mat_1)/3)
  
  accept_pprb_1 <- numeric(iterations)
  log_likelihood_vec<-numeric(iterations)
  
  for(i in 2:iterations){
    
    # --- Propose new values ---
    proposal_log_sigma_s <- rnorm(1, log_sigma_s2_post_pprb_1[i-1], jump_sigma_s2)
    proposal_log_sigma_n <- rnorm(1, log_sigma_n2_post_pprb_1[i-1], jump_sigma_n2)
    proposal_log_phi     <- rnorm(1, log_phi_post_pprb_1[i-1], jump_phi)
    
    # Exponentiate to get actual parameters
    sigma_s2_prop <- exp(proposal_log_sigma_s)
    sigma_n2_prop <- exp(proposal_log_sigma_n)
    phi_prop      <- exp(proposal_log_phi)
    
    sigma_s2_curr <- exp(log_sigma_s2_post_pprb_1[i-1])
    sigma_n2_curr <- exp(log_sigma_n2_post_pprb_1[i-1])
    phi_curr      <- exp(log_phi_post_pprb_1[i-1])
    
    # --- Log posterior ---
    log_post_curr_temp<-loglik_fun(data_1, sigma_s2_curr, sigma_n2_curr, phi_curr,
                                   dist_mat_1,length(data_1)) 
    log_post_curr<- log_post_curr_temp* (1/ temperature)
    log_post_curr <- log_post_curr +
      dnorm(log_sigma_s2_post_pprb_1[i-1], mu_log_sigma_s, sd_log_sigma_s, log=TRUE) +
      dnorm(log_sigma_n2_post_pprb_1[i-1], mu_log_sigma_n, sd_log_sigma_n, log=TRUE) +
      dnorm(log_phi_post_pprb_1[i-1], mu_log_phi, sd_log_phi, log=TRUE)
    
    log_post_prop_temp<-loglik_fun(data_1, sigma_s2_prop, sigma_n2_prop, phi_prop,
                                   dist_mat_1,length(data_1))
    log_post_prop<- log_post_prop_temp* (1/temperature)
    log_post_prop <-  log_post_prop+
      dnorm(proposal_log_sigma_s, mu_log_sigma_s, sd_log_sigma_s, log=TRUE) +
      dnorm(proposal_log_sigma_n, mu_log_sigma_n, sd_log_sigma_n, log=TRUE) +
      dnorm(proposal_log_phi, mu_log_phi, sd_log_phi, log=TRUE)
    
    # --- MH accept/reject ---
    log_ratio <- log_post_prop - log_post_curr
    if(log(runif(1)) < log_ratio){
      log_sigma_s2_post_pprb_1[i] <- proposal_log_sigma_s
      log_sigma_n2_post_pprb_1[i] <- proposal_log_sigma_n
      log_phi_post_pprb_1[i]     <- proposal_log_phi
      accept_pprb_1[i] <- 1
      log_likelihood_vec[i]<-log_post_prop_temp
      
    } else {
      log_sigma_s2_post_pprb_1[i] <- log_sigma_s2_post_pprb_1[i-1]
      log_sigma_n2_post_pprb_1[i] <- log_sigma_n2_post_pprb_1[i-1]
      log_phi_post_pprb_1[i]     <- log_phi_post_pprb_1[i-1]
      log_likelihood_vec[i]<-log_post_curr_temp
    }
  }
  log_sigma_s2_post_pprb_1<-log_sigma_s2_post_pprb_1[-c(1:burnin)]
  log_sigma_n2_post_pprb_1<-log_sigma_n2_post_pprb_1[-c(1:burnin)]
  log_phi_post_pprb_1<-log_phi_post_pprb_1[-c(1:burnin)]
  accept_pprb_1<-accept_pprb_1[-c(1:burnin)]
  log_likelihood_vec<-log_likelihood_vec[-c(1:burnin)]
  
  return(list(log_sigma_s2_post_pprb_1=log_sigma_s2_post_pprb_1,
              log_sigma_n2_post_pprb_1=log_sigma_n2_post_pprb_1,
              log_phi_post_pprb_1=log_phi_post_pprb_1,
              accept_pprb_1=accept_pprb_1,
              log_likelihood_vec=log_likelihood_vec))
  
}

temperature_list<-list()
temperature_list[[1]]<-1
temperature_list[[2]]<-exp( seq(0,3,length.out=10) )[2]
temperature_list[[3]]<-exp( seq(0,3,length.out=10) )[3]
temperature_list[[4]]<-exp( seq(0,3,length.out=10) )[4]
temperature_list[[5]]<-exp( seq(0,3,length.out=10) )[5]
temperature_list[[6]]<-exp( seq(0,3,length.out=10) )[6]
temperature_list[[7]]<-exp( seq(0,3,length.out=10) )[7]
temperature_list[[8]]<-exp( seq(0,3,length.out=10) )[8]
temperature_list[[9]]<-exp( seq(0,3,length.out=10) )[9]
temperature_list[[10]]<-exp( seq(0,3,length.out=10) )[10]

jump_sigma_s2_list<-list()
jump_sigma_s2_list[[1]]<-jump_sigma_s2
jump_sigma_s2_list[[2]]<-seq(jump_sigma_s2,0.5,length.out=10)[2]
jump_sigma_s2_list[[3]]<-seq(jump_sigma_s2,0.5,length.out=10)[3]
jump_sigma_s2_list[[4]]<-seq(jump_sigma_s2,0.5,length.out=10)[4]
jump_sigma_s2_list[[5]]<-seq(jump_sigma_s2,0.5,length.out=10)[5]
jump_sigma_s2_list[[6]]<-seq(jump_sigma_s2,0.5,length.out=10)[6]
jump_sigma_s2_list[[7]]<-seq(jump_sigma_s2,0.5,length.out=10)[7]
jump_sigma_s2_list[[8]]<-seq(jump_sigma_s2,0.5,length.out=10)[8]
jump_sigma_s2_list[[9]]<-seq(jump_sigma_s2,0.5,length.out=10)[9]
jump_sigma_s2_list[[10]]<-seq(jump_sigma_s2,0.5,length.out=10)[10]

jump_sigma_n2_list<-list()
jump_sigma_n2_list[[1]]<-jump_sigma_n2
jump_sigma_n2_list[[2]]<-seq(jump_sigma_n2,0.5,length.out=10)[2]
jump_sigma_n2_list[[3]]<-seq(jump_sigma_n2,0.5,length.out=10)[3]
jump_sigma_n2_list[[4]]<-seq(jump_sigma_n2,0.5,length.out=10)[4]
jump_sigma_n2_list[[5]]<-seq(jump_sigma_n2,0.5,length.out=10)[5]
jump_sigma_n2_list[[6]]<-seq(jump_sigma_n2,0.5,length.out=10)[6]
jump_sigma_n2_list[[7]]<-seq(jump_sigma_n2,0.5,length.out=10)[7]
jump_sigma_n2_list[[8]]<-seq(jump_sigma_n2,0.5,length.out=10)[8]
jump_sigma_n2_list[[9]]<-seq(jump_sigma_n2,0.5,length.out=10)[9]
jump_sigma_n2_list[[10]]<-seq(jump_sigma_n2,0.5,length.out=10)[10]

jump_phi_list<-list()
jump_phi_list[[1]]<-jump_phi
jump_phi_list[[2]]<-seq(jump_phi,0.5,length.out=10)[2]
jump_phi_list[[3]]<-seq(jump_phi,0.5,length.out=10)[3]
jump_phi_list[[4]]<-seq(jump_phi,0.5,length.out=10)[4]
jump_phi_list[[5]]<-seq(jump_phi,0.5,length.out=10)[5]
jump_phi_list[[6]]<-seq(jump_phi,0.5,length.out=10)[6]
jump_phi_list[[7]]<-seq(jump_phi,0.5,length.out=10)[7]
jump_phi_list[[8]]<-seq(jump_phi,0.5,length.out=10)[8]
jump_phi_list[[9]]<-seq(jump_phi,0.5,length.out=10)[9]
jump_phi_list[[10]]<-seq(jump_phi,0.5,length.out=10)[10]


tic()

results <- mclapply(1:10, function(i) {
  
  PPPRB_1st_stage(
    data_1 = data_1,
    dist_mat_1 = dist_mat_1,
    temperature = temperature_list[[i]],
    iterations  = iterations,
    burnin      = burnin,
    jump_sigma_s2 = jump_sigma_s2_list[[i]],
    jump_sigma_n2  = jump_sigma_n2_list[[i]],
    jump_phi = jump_phi_list[[i]],
    mu_log_sigma_s = mu_log_sigma_s,
    sd_log_sigma_s = sd_log_sigma_s,
    mu_log_sigma_n = mu_log_sigma_n,
    sd_log_sigma_n = sd_log_sigma_n,
    mu_log_phi = mu_log_phi,
    sd_log_phi = sd_log_phi
  )
  
}, mc.cores = min(length(temperature_list), parallel::detectCores()-1),mc.set.seed = TRUE)

toc()

names(results) <- c("cold", "hot_1", "hot_2", "hot_3",
                    "hot_4","hot_5","hot_6","hot_7","hot_8","hot_9")

mean( results$cold$accept_pprb_1)
mean( results$hot_1$accept_pprb_1)
mean( results$hot_2$accept_pprb_1)
mean( results$hot_3$accept_pprb_1)
mean( results$hot_4$accept_pprb_1)
mean( results$hot_5$accept_pprb_1)
mean( results$hot_6$accept_pprb_1)
mean( results$hot_7$accept_pprb_1)
mean( results$hot_8$accept_pprb_1)
mean( results$hot_9$accept_pprb_1)


### PPPRB between stages 
### PPPRB between stages 
### PPPRB between stages 
### PPPRB between stages 

all_log_sigma_s2   <- c(
  results$cold$log_sigma_s2_post_pprb_1,
  results$hot_1$log_sigma_s2_post_pprb_1,
  results$hot_2$log_sigma_s2_post_pprb_1,
  results$hot_3$log_sigma_s2_post_pprb_1,
  results$hot_4$log_sigma_s2_post_pprb_1,
  results$hot_5$log_sigma_s2_post_pprb_1,
  results$hot_6$log_sigma_s2_post_pprb_1,
  results$hot_7$log_sigma_s2_post_pprb_1,
  results$hot_8$log_sigma_s2_post_pprb_1,
  results$hot_9$log_sigma_s2_post_pprb_1
)

all_log_sigma_n2  <- c(
  results$cold$log_sigma_n2_post_pprb_1,
  results$hot_1$log_sigma_n2_post_pprb_1,
  results$hot_2$log_sigma_n2_post_pprb_1,
  results$hot_3$log_sigma_n2_post_pprb_1,
  results$hot_4$log_sigma_n2_post_pprb_1,
  results$hot_5$log_sigma_n2_post_pprb_1,
  results$hot_6$log_sigma_n2_post_pprb_1,
  results$hot_7$log_sigma_n2_post_pprb_1,
  results$hot_8$log_sigma_n2_post_pprb_1,
  results$hot_9$log_sigma_n2_post_pprb_1
)


all_log_phi <- c(
  results$cold$log_phi_post_pprb_1,
  results$hot_1$log_phi_post_pprb_1,
  results$hot_2$log_phi_post_pprb_1,
  results$hot_3$log_phi_post_pprb_1,
  results$hot_4$log_phi_post_pprb_1,
  results$hot_5$log_phi_post_pprb_1,
  results$hot_6$log_phi_post_pprb_1,
  results$hot_7$log_phi_post_pprb_1,
  results$hot_8$log_phi_post_pprb_1,
  results$hot_9$log_phi_post_pprb_1
)


all_param_likelihood_ppprb<-cbind(all_log_sigma_s2,
                                 all_log_sigma_n2,
                                 all_log_phi)
colnames(all_param_likelihood_ppprb)<-c("log_sigma_s2","log_sigma_n2","log_phi_post")

all_param_likelihood_ppprb_uni<-unique(all_param_likelihood_ppprb)


tic()
log_likelihood_pre <- unlist(
  mclapply(
    1:nrow(all_param_likelihood_ppprb_uni),
    function(i) {
      C_full <- cov_mat_fun(
        dist = dist_mat_2_1,
        sigma_s2 = exp(all_param_likelihood_ppprb_uni[i, "log_sigma_s2"]),
        sigma_n2 = exp(all_param_likelihood_ppprb_uni[i, "log_sigma_n2"]),
        phi      = exp(all_param_likelihood_ppprb_uni[i, "log_phi_post"]),
        N        = length(data_1) + length(data_2)
      )
      
      log_likelihood_pre_fast_rcpp(
        new_data = data_2,
        old_data = data_1,
        C_full   = C_full
      )
    },
    mc.cores = parallel::detectCores() - 1,mc.set.seed=TRUE)
)

toc()

all_param_likelihood_ppprb_uni<-cbind(all_param_likelihood_ppprb_uni,log_likelihood_pre)
colnames(all_param_likelihood_ppprb_uni)[4]<-"log_likeli"
all_param_likelihood_ppprb<-cbind(all_param_likelihood_ppprb,NA)
colnames(all_param_likelihood_ppprb)[4]<-"log_likeli"

for(i in 1:nrow(all_param_likelihood_ppprb_uni)){
  idx<-which(all_param_likelihood_ppprb[,1]==all_param_likelihood_ppprb_uni[i,1]&
               all_param_likelihood_ppprb[,2]==all_param_likelihood_ppprb_uni[i,2] &
               all_param_likelihood_ppprb[,3]==all_param_likelihood_ppprb_uni[i,3])
  all_param_likelihood_ppprb[idx,4]<-all_param_likelihood_ppprb_uni[i,4]
}


log_likelihood_pre<-all_param_likelihood_ppprb[,4]

rm(all_param_likelihood_ppprb,
   all_param_likelihood_ppprb_uni
)


log_likelihood_conditional<-log_likelihood_pre
rm(log_likelihood_pre)



all_param_likeli_cond<-cbind(all_log_sigma_n2,all_log_sigma_s2,all_log_phi,log_likelihood_conditional)


colnames(all_param_likeli_cond)<-c("log_sigma_n2",
                                   "log_sigma_s2",
                                   "log_phi",
                                   "log_likeli_condi")



rm(log_likelihood_conditional)

dim(all_param_likeli_cond)

#####

all_param_likeli_1st<-cbind(all_log_sigma_n2,
                             all_log_sigma_s2,
                             all_log_phi)

                             
all_param_likeli_1st<-cbind(all_param_likeli_1st, c(results$cold$log_likelihood_vec,
                                                    results$hot_1$log_likelihood_vec,
                                                    results$hot_2$log_likelihood_vec,
                                                    results$hot_3$log_likelihood_vec,
                                                    results$hot_4$log_likelihood_vec,
                                                    results$hot_5$log_likelihood_vec,
                                                    results$hot_6$log_likelihood_vec,
                                                    results$hot_7$log_likelihood_vec,
                                                    results$hot_8$log_likelihood_vec,
                                                    results$hot_9$log_likelihood_vec))


colnames(all_param_likeli_1st)<-c("log_sigma_n2",
                                  "log_sigma_s2",
                                  "log_phi",
                                  "log_likeli_1st")

all_param_likeli_cond[,1:3]==all_param_likeli_1st[,1:3]
all_param_likeli_1_2 <-all_param_likeli_cond[,1:3]
all_param_likeli_1_2<-cbind(all_param_likeli_1_2,NA)
all_param_likeli_1_2[,4]<-all_param_likeli_cond[,4]+all_param_likeli_1st[,4]
colnames(all_param_likeli_1_2)<-c("log_sigma_n2",
                                  "log_sigma_s2",
                                  "log_phi",
                                  "log_likeli")

### PPPRB second stage 
### PPPRB second stage 
### PPPRB second stage 
### PPPRB second stage 

all_param_likeli_cond<-cbind( all_param_likeli_cond,rep(1:10,each=iterations-burnin) )

colnames(all_param_likeli_cond)[5]<-"chain"

all_param_likeli_cond<-cbind(1:nrow(all_param_likeli_cond),all_param_likeli_cond)
colnames(all_param_likeli_cond)[1]<-"idxidx"

all_param_likeli_1_2<-cbind(1:nrow(all_param_likeli_1_2),all_param_likeli_1_2)
colnames(all_param_likeli_1_2)[1]<-"idxidx"


PPPRB_2nd_stage<-function(temperature,
                          all_param_likeli_cond,
                          chain_rows,
                          current_idx,
                          chain_idx,
                          #
                          PPPRB_log_s2_1st, ## all
                          PPPRB_log_n2_1st,
                          PPPRB_log_phi_1st,
                          #
                          current_log_s2,
                          current_log_n2,
                          current_log_phi){
  rows_chain <- chain_rows[[chain_idx]]   # 
  idx <- sample(rows_chain, 1)            # 
  ## here idx is for all set of parameters
  proposal_log_s2<-PPPRB_log_s2_1st[idx]
  proposal_log_n2<-PPPRB_log_n2_1st[idx]
  proposal_log_phi<-PPPRB_log_phi_1st[idx]
  #

  nom<-all_param_likeli_cond[idx,5] * (1/temperature)
  

  denom<-all_param_likeli_cond[current_idx,5] * (1/temperature)
  
  ratio<-nom-denom 
  ## accept / reject
  if (log(runif(1)) <= ratio) {
    log_s2_final<-proposal_log_s2
    log_n2_final<-proposal_log_n2
    log_phi_final<-proposal_log_phi
    accept_or_not<-1
    idx_result<-idx
  } else {
    log_s2_final<-current_log_s2
    log_n2_final<-current_log_n2
    log_phi_final<-current_log_phi
    accept_or_not<-0
    idx_result<-current_idx
  } 
  return(list(log_s2_final=log_s2_final,
              log_n2_final=log_n2_final,
              log_phi_final=log_phi_final,
              accept_or_not=accept_or_not,
              idx_result=idx_result)
  )
}


PPPRB_swap<-function(temperature,
                     #
                     all_param_likeli_cond,
                     all_param_likeli_1st,
                     #
                     cold_idx,
                     hot_idx,
                     #
                     PPPRB_log_s2_cold,
                     PPPRB_log_s2_hot,
                     PPPRB_log_n2_cold,
                     PPPRB_log_n2_hot,
                     PPPRB_log_phi_cold,
                     PPPRB_log_phi_hot){
  

  nom_conditional_log_lik_hot<-all_param_likeli_cond[hot_idx,5]
  nom1<-nom_conditional_log_lik_hot *(1-1/temperature)
  #
  nom_log_lik_1st_hot<-all_param_likeli_1st[hot_idx,4]
  nom2<-nom_log_lik_1st_hot * (1-1/temperature)
  #
  

  denom_conditional_log_lik_cold<-all_param_likeli_cond[cold_idx,5]
  denom1<-denom_conditional_log_lik_cold* (1-1/temperature)

  denom_log_lik_1st_cold<-all_param_likeli_1st[cold_idx,4]
  denom2<- denom_log_lik_1st_cold * (1-1/temperature)
  
  #
  ratio<-nom1+nom2-(denom1+denom2)
  ## accept / reject
  if (log(runif(1)) <= ratio) {
    log_s2_cold_chain<-PPPRB_log_s2_hot
    log_n2_cold_chain<-PPPRB_log_n2_hot
    log_phi_cold_chain<-PPPRB_log_phi_hot
    #
    log_s2_hot_chain<-PPPRB_log_s2_cold
    log_n2_hot_chain<-PPPRB_log_n2_cold
    log_phi_hot_chain<-PPPRB_log_phi_cold
    #
    accept_or_not<-1
  } else {
    log_s2_cold_chain<-PPPRB_log_s2_cold
    log_n2_cold_chain<-PPPRB_log_n2_cold
    log_phi_cold_chain<-PPPRB_log_phi_cold
    #
    log_s2_hot_chain<-PPPRB_log_s2_hot
    log_n2_hot_chain<-PPPRB_log_n2_hot
    log_phi_hot_chain<-PPPRB_log_phi_hot
    #
    accept_or_not<-0

  } 
  return(list(log_s2_cold_chain=log_s2_cold_chain,
              log_n2_cold_chain=log_n2_cold_chain,
              log_phi_cold_chain=log_phi_cold_chain,
              #
              log_s2_hot_chain=log_s2_hot_chain,
              log_n2_hot_chain=log_n2_hot_chain,
              log_phi_hot_chain=log_phi_hot_chain,
              accept_or_not=accept_or_not)
  )
}

PPPRB_log_s2_2nd_cold<-numeric(iterations)
PPPRB_log_s2_2nd_hot_1<-numeric(iterations)
PPPRB_log_s2_2nd_hot_2<-numeric(iterations)
PPPRB_log_s2_2nd_hot_3<-numeric(iterations)
PPPRB_log_s2_2nd_hot_4<-numeric(iterations)
PPPRB_log_s2_2nd_hot_5<-numeric(iterations)
PPPRB_log_s2_2nd_hot_6<-numeric(iterations)
PPPRB_log_s2_2nd_hot_7<-numeric(iterations)
PPPRB_log_s2_2nd_hot_8<-numeric(iterations)
PPPRB_log_s2_2nd_hot_9<-numeric(iterations)
#
PPPRB_log_n2_2nd_cold<-numeric(iterations)
PPPRB_log_n2_2nd_hot_1<-numeric(iterations)
PPPRB_log_n2_2nd_hot_2<-numeric(iterations)
PPPRB_log_n2_2nd_hot_3<-numeric(iterations)
PPPRB_log_n2_2nd_hot_4<-numeric(iterations)
PPPRB_log_n2_2nd_hot_5<-numeric(iterations)
PPPRB_log_n2_2nd_hot_6<-numeric(iterations)
PPPRB_log_n2_2nd_hot_7<-numeric(iterations)
PPPRB_log_n2_2nd_hot_8<-numeric(iterations)
PPPRB_log_n2_2nd_hot_9<-numeric(iterations)
#
PPPRB_log_phi_2nd_cold<-numeric(iterations)
PPPRB_log_phi_2nd_hot_1<-numeric(iterations)
PPPRB_log_phi_2nd_hot_2<-numeric(iterations)
PPPRB_log_phi_2nd_hot_3<-numeric(iterations)
PPPRB_log_phi_2nd_hot_4<-numeric(iterations)
PPPRB_log_phi_2nd_hot_5<-numeric(iterations)
PPPRB_log_phi_2nd_hot_6<-numeric(iterations)
PPPRB_log_phi_2nd_hot_7<-numeric(iterations)
PPPRB_log_phi_2nd_hot_8<-numeric(iterations)
PPPRB_log_phi_2nd_hot_9<-numeric(iterations)

# acceptance ratio
accept_ratio_2nd_cold<-numeric(iterations)
accept_ratio_2nd_hot_1<-numeric(iterations)
accept_ratio_2nd_hot_2<-numeric(iterations)
accept_ratio_2nd_hot_3<-numeric(iterations)
accept_ratio_2nd_hot_4<-numeric(iterations)
accept_ratio_2nd_hot_5<-numeric(iterations)
accept_ratio_2nd_hot_6<-numeric(iterations)
accept_ratio_2nd_hot_7<-numeric(iterations)
accept_ratio_2nd_hot_8<-numeric(iterations)
accept_ratio_2nd_hot_9<-numeric(iterations)


## result from the 1st stage
PPPRB_log_s2_1st_list<-c(results$cold$log_sigma_s2_post_pprb_1,
                         results$hot_1$log_sigma_s2_post_pprb_1,
                         results$hot_2$log_sigma_s2_post_pprb_1,
                         results$hot_3$log_sigma_s2_post_pprb_1,
                         results$hot_4$log_sigma_s2_post_pprb_1,
                         results$hot_5$log_sigma_s2_post_pprb_1,
                         results$hot_6$log_sigma_s2_post_pprb_1,
                         results$hot_7$log_sigma_s2_post_pprb_1,
                         results$hot_8$log_sigma_s2_post_pprb_1,
                         results$hot_9$log_sigma_s2_post_pprb_1)

PPPRB_log_n2_1st_list<-c(results$cold$log_sigma_n2_post_pprb_1,
                         results$hot_1$log_sigma_n2_post_pprb_1,
                         results$hot_2$log_sigma_n2_post_pprb_1,
                         results$hot_3$log_sigma_n2_post_pprb_1,
                         results$hot_4$log_sigma_n2_post_pprb_1,
                         results$hot_5$log_sigma_n2_post_pprb_1,
                         results$hot_6$log_sigma_n2_post_pprb_1,
                         results$hot_7$log_sigma_n2_post_pprb_1,
                         results$hot_8$log_sigma_n2_post_pprb_1,
                         results$hot_9$log_sigma_n2_post_pprb_1)

PPPRB_log_phi_1st_list<-c(results$cold$log_phi_post_pprb_1,
                          results$hot_1$log_phi_post_pprb_1,
                          results$hot_2$log_phi_post_pprb_1,
                          results$hot_3$log_phi_post_pprb_1,
                          results$hot_4$log_phi_post_pprb_1,
                          results$hot_5$log_phi_post_pprb_1,
                          results$hot_6$log_phi_post_pprb_1,
                          results$hot_7$log_phi_post_pprb_1,
                          results$hot_8$log_phi_post_pprb_1,
                          results$hot_9$log_phi_post_pprb_1)
## initial
iterations-burnin
current_idx<-c(100,
               (iterations-burnin)+100,
               (iterations-burnin)*2+100,
               (iterations-burnin)*3+100,
               (iterations-burnin)*4+100,
               (iterations-burnin)*5+100,
               (iterations-burnin)*6+100,
               (iterations-burnin)*7+100,
               (iterations-burnin)*8+100,
               (iterations-burnin)*9+100)

PPPRB_log_s2_2nd_cold[1]<-PPPRB_log_s2_1st_list[current_idx[1]]
PPPRB_log_s2_2nd_hot_1[1]<-PPPRB_log_s2_1st_list[current_idx[2]]
PPPRB_log_s2_2nd_hot_2[1]<-PPPRB_log_s2_1st_list[current_idx[3]]
PPPRB_log_s2_2nd_hot_3[1]<-PPPRB_log_s2_1st_list[current_idx[4]]
PPPRB_log_s2_2nd_hot_4[1]<-PPPRB_log_s2_1st_list[current_idx[5]]
PPPRB_log_s2_2nd_hot_5[1]<-PPPRB_log_s2_1st_list[current_idx[6]]
PPPRB_log_s2_2nd_hot_6[1]<-PPPRB_log_s2_1st_list[current_idx[7]]
PPPRB_log_s2_2nd_hot_7[1]<-PPPRB_log_s2_1st_list[current_idx[8]]
PPPRB_log_s2_2nd_hot_8[1]<-PPPRB_log_s2_1st_list[current_idx[9]]
PPPRB_log_s2_2nd_hot_9[1]<-PPPRB_log_s2_1st_list[current_idx[10]]

## initial 
PPPRB_log_n2_2nd_cold[1]<-PPPRB_log_n2_1st_list[current_idx[1]]
PPPRB_log_n2_2nd_hot_1[1]<-PPPRB_log_n2_1st_list[current_idx[2]]
PPPRB_log_n2_2nd_hot_2[1]<-PPPRB_log_n2_1st_list[current_idx[3]]
PPPRB_log_n2_2nd_hot_3[1]<-PPPRB_log_n2_1st_list[current_idx[4]]
PPPRB_log_n2_2nd_hot_4[1]<-PPPRB_log_n2_1st_list[current_idx[5]]
PPPRB_log_n2_2nd_hot_5[1]<-PPPRB_log_n2_1st_list[current_idx[6]]
PPPRB_log_n2_2nd_hot_6[1]<-PPPRB_log_n2_1st_list[current_idx[7]]
PPPRB_log_n2_2nd_hot_7[1]<-PPPRB_log_n2_1st_list[current_idx[8]]
PPPRB_log_n2_2nd_hot_8[1]<-PPPRB_log_n2_1st_list[current_idx[9]]
PPPRB_log_n2_2nd_hot_9[1]<-PPPRB_log_n2_1st_list[current_idx[10]]

## initial 
PPPRB_log_phi_2nd_cold[1]<-PPPRB_log_phi_1st_list[current_idx[[1]]]
PPPRB_log_phi_2nd_hot_1[1]<-PPPRB_log_phi_1st_list[current_idx[[2]]]
PPPRB_log_phi_2nd_hot_2[1]<-PPPRB_log_phi_1st_list[current_idx[[3]]]
PPPRB_log_phi_2nd_hot_3[1]<-PPPRB_log_phi_1st_list[current_idx[[4]]]
PPPRB_log_phi_2nd_hot_4[1]<-PPPRB_log_phi_1st_list[current_idx[[5]]]
PPPRB_log_phi_2nd_hot_5[1]<-PPPRB_log_phi_1st_list[current_idx[[6]]]
PPPRB_log_phi_2nd_hot_6[1]<-PPPRB_log_phi_1st_list[current_idx[[7]]]
PPPRB_log_phi_2nd_hot_7[1]<-PPPRB_log_phi_1st_list[current_idx[[8]]]
PPPRB_log_phi_2nd_hot_8[1]<-PPPRB_log_phi_1st_list[current_idx[[9]]]
PPPRB_log_phi_2nd_hot_9[1]<-PPPRB_log_phi_1st_list[current_idx[[10]]]



j<-2

accept_ratio_swap<-numeric(iterations)

chain_rows <- split(
  seq_len(nrow(all_param_likeli_cond)),
  all_param_likeli_cond[, "chain"]   # 
)


tic()
j<-2
for(j in 2: iterations){
  current_log_s2 <- c(
    PPPRB_log_s2_2nd_cold[j-1],
    PPPRB_log_s2_2nd_hot_1[j-1],
    PPPRB_log_s2_2nd_hot_2[j-1],
    PPPRB_log_s2_2nd_hot_3[j-1],
    PPPRB_log_s2_2nd_hot_4[j-1],
    PPPRB_log_s2_2nd_hot_5[j-1],
    PPPRB_log_s2_2nd_hot_6[j-1],
    PPPRB_log_s2_2nd_hot_7[j-1],
    PPPRB_log_s2_2nd_hot_8[j-1],
    PPPRB_log_s2_2nd_hot_9[j-1]
  )
  
  current_log_n2 <- c(
    PPPRB_log_n2_2nd_cold[j-1],
    PPPRB_log_n2_2nd_hot_1[j-1],
    PPPRB_log_n2_2nd_hot_2[j-1],
    PPPRB_log_n2_2nd_hot_3[j-1],
    PPPRB_log_n2_2nd_hot_4[j-1],
    PPPRB_log_n2_2nd_hot_5[j-1],
    PPPRB_log_n2_2nd_hot_6[j-1],
    PPPRB_log_n2_2nd_hot_7[j-1],
    PPPRB_log_n2_2nd_hot_8[j-1],
    PPPRB_log_n2_2nd_hot_9[j-1]
  )
  
  current_log_phi <- c(
    PPPRB_log_phi_2nd_cold[j-1],
    PPPRB_log_phi_2nd_hot_1[j-1],
    PPPRB_log_phi_2nd_hot_2[j-1],
    PPPRB_log_phi_2nd_hot_3[j-1],
    PPPRB_log_phi_2nd_hot_4[j-1],
    PPPRB_log_phi_2nd_hot_5[j-1],
    PPPRB_log_phi_2nd_hot_6[j-1],
    PPPRB_log_phi_2nd_hot_7[j-1],
    PPPRB_log_phi_2nd_hot_8[j-1],
    PPPRB_log_phi_2nd_hot_9[j-1]
  )

  results_ppprb <- lapply(
    1:10,
    function(i) {
      PPPRB_2nd_stage(
        temperature = temperature_list[[i]],
        all_param_likeli_cond = all_param_likeli_cond,
        chain_rows=chain_rows,
        current_idx=current_idx[i],
        chain_idx=c(1:10)[i],
        #
        PPPRB_log_s2_1st = PPPRB_log_s2_1st_list,
        PPPRB_log_n2_1st  = PPPRB_log_n2_1st_list,
        PPPRB_log_phi_1st   = PPPRB_log_phi_1st_list,
        #
        current_log_s2  = current_log_s2[i],
        current_log_n2    = current_log_n2[i],
        current_log_phi   = current_log_phi[i]
      )
    }
  )
  names(results_ppprb) <- c("cold", "hot_1", "hot_2", "hot_3","hot_4","hot_5",
                            "hot_6","hot_7","hot_8","hot_9")
  #
  log_s2_final    <- sapply(results_ppprb, `[[`, "log_s2_final")
  log_n2_final  <- sapply(results_ppprb, `[[`, "log_n2_final")
  log_phi_final   <- sapply(results_ppprb, `[[`, "log_phi_final")
  accept      <- sapply(results_ppprb, `[[`, "accept_or_not")
  idx_final      <- sapply(results_ppprb, `[[`, "idx_result")
  #
  PPPRB_log_s2_2nd_cold[j]    <- log_s2_final["cold"]
  PPPRB_log_s2_2nd_hot_1[j]   <- log_s2_final["hot_1"]
  PPPRB_log_s2_2nd_hot_2[j]   <- log_s2_final["hot_2"]
  PPPRB_log_s2_2nd_hot_3[j]   <- log_s2_final["hot_3"]
  PPPRB_log_s2_2nd_hot_4[j]   <- log_s2_final["hot_4"]
  PPPRB_log_s2_2nd_hot_5[j]   <- log_s2_final["hot_5"]
  PPPRB_log_s2_2nd_hot_6[j]   <- log_s2_final["hot_6"]
  PPPRB_log_s2_2nd_hot_7[j]   <- log_s2_final["hot_7"]
  PPPRB_log_s2_2nd_hot_8[j]   <- log_s2_final["hot_8"]
  PPPRB_log_s2_2nd_hot_9[j]   <- log_s2_final["hot_9"]
  
  PPPRB_log_n2_2nd_cold[j]  <- log_n2_final["cold"]
  PPPRB_log_n2_2nd_hot_1[j] <- log_n2_final["hot_1"]
  PPPRB_log_n2_2nd_hot_2[j] <- log_n2_final["hot_2"]
  PPPRB_log_n2_2nd_hot_3[j] <- log_n2_final["hot_3"]
  PPPRB_log_n2_2nd_hot_4[j] <- log_n2_final["hot_4"]
  PPPRB_log_n2_2nd_hot_5[j] <- log_n2_final["hot_5"]
  PPPRB_log_n2_2nd_hot_6[j] <- log_n2_final["hot_6"]
  PPPRB_log_n2_2nd_hot_7[j] <- log_n2_final["hot_7"]
  PPPRB_log_n2_2nd_hot_8[j] <- log_n2_final["hot_8"]
  PPPRB_log_n2_2nd_hot_9[j] <- log_n2_final["hot_9"]
  
  PPPRB_log_phi_2nd_cold[j]   <- log_phi_final["cold"]
  PPPRB_log_phi_2nd_hot_1[j]  <- log_phi_final["hot_1"]
  PPPRB_log_phi_2nd_hot_2[j]  <- log_phi_final["hot_2"]
  PPPRB_log_phi_2nd_hot_3[j]  <- log_phi_final["hot_3"]
  PPPRB_log_phi_2nd_hot_4[j]  <- log_phi_final["hot_4"]
  PPPRB_log_phi_2nd_hot_5[j]  <- log_phi_final["hot_5"]
  PPPRB_log_phi_2nd_hot_6[j]  <- log_phi_final["hot_6"]
  PPPRB_log_phi_2nd_hot_7[j]  <- log_phi_final["hot_7"]
  PPPRB_log_phi_2nd_hot_8[j]  <- log_phi_final["hot_8"]
  PPPRB_log_phi_2nd_hot_9[j]  <- log_phi_final["hot_9"]
  
  accept_ratio_2nd_cold[j]  <- accept["cold"]
  accept_ratio_2nd_hot_1[j] <- accept["hot_1"]
  accept_ratio_2nd_hot_2[j] <- accept["hot_2"]
  accept_ratio_2nd_hot_3[j] <- accept["hot_3"]
  accept_ratio_2nd_hot_4[j] <- accept["hot_4"]
  accept_ratio_2nd_hot_5[j] <- accept["hot_5"]
  accept_ratio_2nd_hot_6[j] <- accept["hot_6"]
  accept_ratio_2nd_hot_7[j] <- accept["hot_7"]
  accept_ratio_2nd_hot_8[j] <- accept["hot_8"]
  accept_ratio_2nd_hot_9[j] <- accept["hot_9"]
  ## current_idx_update
  current_idx<-as.numeric(idx_final)
  
  #
  if (j %% 1000 == 0) cat(j, "\n")
  # swapping
  # swapping
  # swapping
  idx<-sample(2:10,1,replace=FALSE)
  if(idx==2){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_1[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_1[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_1[j]
  } else if(idx==3){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_2[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_2[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_2[j]
  } else if(idx==4){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_3[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_3[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_3[j]

  } else if(idx==5){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_4[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_4[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_4[j]
  } else if(idx==6){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_5[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_5[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_5[j]
  } else if(idx==7){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_6[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_6[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_6[j]
  } else if(idx==8){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_7[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_7[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_7[j]
  } else if(idx==9){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_8[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_8[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_8[j]
  } else if(idx==10){
    PPPRB_log_s2_hot<- PPPRB_log_s2_2nd_hot_9[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_2nd_hot_9[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_2nd_hot_9[j]
  } 

  swap_result<-PPPRB_swap(temperature = temperature_list[[idx]],
                          #
                          all_param_likeli_cond=all_param_likeli_cond,
                          all_param_likeli_1st=all_param_likeli_1st,
                          #
                          cold_idx=current_idx[1],
                          hot_idx=current_idx[idx],
                          
                          PPPRB_log_s2_cold=PPPRB_log_s2_2nd_cold[j],
                          PPPRB_log_s2_hot=PPPRB_log_s2_hot,
                          #
                          PPPRB_log_n2_cold=PPPRB_log_n2_2nd_cold[j],
                          PPPRB_log_n2_hot=PPPRB_log_n2_hot,
                          #
                          PPPRB_log_phi_cold=PPPRB_log_phi_2nd_cold[j],
                          PPPRB_log_phi_hot=PPPRB_log_phi_hot)
  if(swap_result$accept_or_not==1){
    PPPRB_log_s2_2nd_cold[j] <-swap_result$log_s2_cold_chain
    PPPRB_log_n2_2nd_cold[j]<-swap_result$log_n2_cold_chain
    PPPRB_log_phi_2nd_cold[j]<-swap_result$log_phi_cold_chain
    accept_ratio_swap[j]<-1
    current_idx_cold<-current_idx[1]
    current_idx_hot<-current_idx[idx]
    current_idx[1]<-current_idx_hot
    #
    if(idx==2){
      PPPRB_log_s2_2nd_hot_1[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_1[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_1[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold

    } else if(idx==3){
      PPPRB_log_s2_2nd_hot_2[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_2[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_2[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==4){
      PPPRB_log_s2_2nd_hot_3[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_3[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_3[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==5){
      PPPRB_log_s2_2nd_hot_4[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_4[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_4[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==6){
      PPPRB_log_s2_2nd_hot_5[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_5[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_5[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==7){
      PPPRB_log_s2_2nd_hot_6[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_6[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_6[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==8){
      PPPRB_log_s2_2nd_hot_7[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_7[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_7[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==9){
      PPPRB_log_s2_2nd_hot_8[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_8[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_8[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==10){
      PPPRB_log_s2_2nd_hot_9[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_2nd_hot_9[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_2nd_hot_9[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    }
  }
}

toc()


mean(accept_ratio_swap)
mean( accept_ratio_2nd_cold )
mean( accept_ratio_2nd_hot_1 )
mean( accept_ratio_2nd_hot_2 )
mean( accept_ratio_2nd_hot_3 )
mean( accept_ratio_2nd_hot_4 )
mean( accept_ratio_2nd_hot_5 )
mean( accept_ratio_2nd_hot_6 )
mean( accept_ratio_2nd_hot_7 )
mean( accept_ratio_2nd_hot_8 )
mean( accept_ratio_2nd_hot_9 )

PPPRB_log_s2_2nd_cold<-PPPRB_log_s2_2nd_cold[-c(1:burnin)]
PPPRB_log_n2_2nd_cold<-PPPRB_log_n2_2nd_cold[-c(1:burnin)]
PPPRB_log_phi_2nd_cold<-PPPRB_log_phi_2nd_cold[-c(1:burnin)]




## hot chains burnin remove
PPPRB_log_s2_2nd_hot_1<-PPPRB_log_s2_2nd_hot_1[-c(1:burnin)]
PPPRB_log_s2_2nd_hot_2<-PPPRB_log_s2_2nd_hot_2[-c(1:burnin)]
PPPRB_log_s2_2nd_hot_3<-PPPRB_log_s2_2nd_hot_3[-c(1:burnin)]
PPPRB_log_s2_2nd_hot_4<-PPPRB_log_s2_2nd_hot_4[-c(1:burnin)]
PPPRB_log_s2_2nd_hot_5<-PPPRB_log_s2_2nd_hot_5[-c(1:burnin)]
PPPRB_log_s2_2nd_hot_6<-PPPRB_log_s2_2nd_hot_6[-c(1:burnin)]
PPPRB_log_s2_2nd_hot_7<-PPPRB_log_s2_2nd_hot_7[-c(1:burnin)]
PPPRB_log_s2_2nd_hot_8<-PPPRB_log_s2_2nd_hot_8[-c(1:burnin)]
PPPRB_log_s2_2nd_hot_9<-PPPRB_log_s2_2nd_hot_9[-c(1:burnin)]


PPPRB_log_n2_2nd_hot_1<-PPPRB_log_n2_2nd_hot_1[-c(1:burnin)]
PPPRB_log_n2_2nd_hot_2<-PPPRB_log_n2_2nd_hot_2[-c(1:burnin)]
PPPRB_log_n2_2nd_hot_3<-PPPRB_log_n2_2nd_hot_3[-c(1:burnin)]
PPPRB_log_n2_2nd_hot_4<-PPPRB_log_n2_2nd_hot_4[-c(1:burnin)]
PPPRB_log_n2_2nd_hot_5<-PPPRB_log_n2_2nd_hot_5[-c(1:burnin)]
PPPRB_log_n2_2nd_hot_6<-PPPRB_log_n2_2nd_hot_6[-c(1:burnin)]
PPPRB_log_n2_2nd_hot_7<-PPPRB_log_n2_2nd_hot_7[-c(1:burnin)]
PPPRB_log_n2_2nd_hot_8<-PPPRB_log_n2_2nd_hot_8[-c(1:burnin)]
PPPRB_log_n2_2nd_hot_9<-PPPRB_log_n2_2nd_hot_9[-c(1:burnin)]

PPPRB_log_phi_2nd_hot_1<-PPPRB_log_phi_2nd_hot_1[-c(1:burnin)]
PPPRB_log_phi_2nd_hot_2<-PPPRB_log_phi_2nd_hot_2[-c(1:burnin)]
PPPRB_log_phi_2nd_hot_3<-PPPRB_log_phi_2nd_hot_3[-c(1:burnin)]
PPPRB_log_phi_2nd_hot_4<-PPPRB_log_phi_2nd_hot_4[-c(1:burnin)]
PPPRB_log_phi_2nd_hot_5<-PPPRB_log_phi_2nd_hot_5[-c(1:burnin)]
PPPRB_log_phi_2nd_hot_6<-PPPRB_log_phi_2nd_hot_6[-c(1:burnin)]
PPPRB_log_phi_2nd_hot_7<-PPPRB_log_phi_2nd_hot_7[-c(1:burnin)]
PPPRB_log_phi_2nd_hot_8<-PPPRB_log_phi_2nd_hot_8[-c(1:burnin)]
PPPRB_log_phi_2nd_hot_9<-PPPRB_log_phi_2nd_hot_9[-c(1:burnin)]

rm(all_log_phi,all_log_sigma_n2,all_log_sigma_s2,current_idx,current_idx_cold,
   current_idx_hot,current_log_n2,current_log_phi,current_log_s2)

rm(idx,idx_final)
rm(log_n2_final,log_phi_final)
rm(log_s2_final)
rm(PPPRB_log_n2_hot,PPPRB_log_phi_hot)
rm(PPPRB_log_s2_hot)
# ppprb 3rd stage ---------------------------------------------------------
all_param_likeli_1_2
### PPPRB between stages 
### PPPRB between stages 
### PPPRB between stages 
### PPPRB between stages 

all_log_sigma_s2   <- c(
  PPPRB_log_s2_2nd_cold,
  PPPRB_log_s2_2nd_hot_1,
  PPPRB_log_s2_2nd_hot_2,
  PPPRB_log_s2_2nd_hot_3,
  PPPRB_log_s2_2nd_hot_4,
  PPPRB_log_s2_2nd_hot_5,
  PPPRB_log_s2_2nd_hot_6,
  PPPRB_log_s2_2nd_hot_7,
  PPPRB_log_s2_2nd_hot_8,
  PPPRB_log_s2_2nd_hot_9
)

all_log_sigma_n2  <- c(
  PPPRB_log_n2_2nd_cold,
  PPPRB_log_n2_2nd_hot_1,
  PPPRB_log_n2_2nd_hot_2,
  PPPRB_log_n2_2nd_hot_3,
  PPPRB_log_n2_2nd_hot_4,
  PPPRB_log_n2_2nd_hot_5,
  PPPRB_log_n2_2nd_hot_6,
  PPPRB_log_n2_2nd_hot_7,
  PPPRB_log_n2_2nd_hot_8,
  PPPRB_log_n2_2nd_hot_9
)


all_log_phi <- c(
  PPPRB_log_phi_2nd_cold,
  PPPRB_log_phi_2nd_hot_1,
  PPPRB_log_phi_2nd_hot_2,
  PPPRB_log_phi_2nd_hot_3,
  PPPRB_log_phi_2nd_hot_4,
  PPPRB_log_phi_2nd_hot_5,
  PPPRB_log_phi_2nd_hot_6,
  PPPRB_log_phi_2nd_hot_7,
  PPPRB_log_phi_2nd_hot_8,
  PPPRB_log_phi_2nd_hot_9
)



all_param_likelihood_ppprb<-cbind(all_log_sigma_s2,
                                  all_log_sigma_n2,
                                  all_log_phi)
colnames(all_param_likelihood_ppprb)<-c("log_sigma_s2","log_sigma_n2","log_phi_post")

all_param_likelihood_ppprb_uni<-unique(all_param_likelihood_ppprb)





tic()
log_likelihood_pre <- unlist(
  mclapply(
    1:nrow(all_param_likelihood_ppprb_uni),
    function(i) {
      # 
      C_full <- cov_mat_fun(
        dist = dist_mat_3_21,
        sigma_s2 = exp(all_param_likelihood_ppprb_uni[i, 1]),
        sigma_n2 = exp(all_param_likelihood_ppprb_uni[i, 2]),
        phi      = exp(all_param_likelihood_ppprb_uni[i, 3]),
        N        = length(data_1) + length(data_2) + length(data_3)
      )
      
      # 
      log_likelihood_pre_fast_rcpp(
        new_data = data_3,
        old_data = c(data_2, data_1),
        C_full   = C_full
      )
    },
    mc.cores =  parallel::detectCores() - 1,
    mc.set.seed = TRUE
  )
)
toc()


all_param_likelihood_ppprb_uni<-cbind(all_param_likelihood_ppprb_uni,log_likelihood_pre)
colnames(all_param_likelihood_ppprb_uni)[4]<-"log_likeli"
all_param_likelihood_ppprb<-cbind(all_param_likelihood_ppprb,NA)
colnames(all_param_likelihood_ppprb)[4]<-"log_likeli"


for(i in 1:nrow(all_param_likelihood_ppprb_uni)){
  idx<-which(all_param_likelihood_ppprb[,1]==all_param_likelihood_ppprb_uni[i,1]&
               all_param_likelihood_ppprb[,2]==all_param_likelihood_ppprb_uni[i,2] &
               all_param_likelihood_ppprb[,3]==all_param_likelihood_ppprb_uni[i,3])
  all_param_likelihood_ppprb[idx,4]<-all_param_likelihood_ppprb_uni[i,4]
}


log_likelihood_pre<-all_param_likelihood_ppprb[,4]
rm(all_param_likelihood_ppprb,
   all_param_likelihood_ppprb_uni
)



log_likelihood_pre

log_likelihood_conditional_2<-log_likelihood_pre
rm(log_likelihood_pre)


all_param_likeli_cond_2nd<-cbind(all_log_sigma_n2,all_log_sigma_s2,all_log_phi,log_likelihood_conditional_2)

head(all_param_likeli_cond_2nd)
colnames(all_param_likeli_cond_2nd)<-c("log_sigma_n2",
                                   "log_sigma_s2",
                                   "log_phi",
                                   "log_likeli_condi")
rm(log_likelihood_conditional_2)


#####


all_param_likeli_2nd<-cbind(all_log_sigma_n2,
                            all_log_sigma_s2,
                            all_log_phi)
colnames(all_param_likeli_2nd)<-c("log_sigma_n2",
                                  "log_sigma_s2",
                                  "log_phi")
all_param_likeli_2nd<-cbind(all_param_likeli_2nd,NA)
colnames(all_param_likeli_2nd)[4]<-"log_likeli"

all_param_likeli_1_2<-all_param_likeli_1_2[,-c(1)]

all_param_likeli_1_2_uni<-unique(all_param_likeli_1_2)

for(i in 1:nrow(all_param_likeli_1_2_uni)){
  idx<-which(all_param_likeli_2nd[,1]==all_param_likeli_1_2_uni[i,1]&
               all_param_likeli_2nd[,2]==all_param_likeli_1_2_uni[i,2] &
               all_param_likeli_2nd[,3]==all_param_likeli_1_2_uni[i,3])
  all_param_likeli_2nd[idx,4]<-all_param_likeli_1_2_uni[i,4]
}

tail(all_param_likeli_2nd)
dist_mat_21 <- distm(rbind(locations_2,locations_1), fun = distCosine) / 1000


loglik_fun(c(data_2,data_1),
           exp(all_param_likeli_2nd[150000,2]),
           exp(all_param_likeli_2nd[150000,1]),
           exp(all_param_likeli_2nd[150000,3]),
           dist_mat_21,
           length(c(data_2,data_1))
)
all_param_likeli_cond_2nd

sum( all_param_likeli_2nd[,c(1)]==all_param_likeli_cond_2nd[,c(1)] )
### PPPRB third stage 
### PPPRB third stage 
### PPPRB third stage 
### PPPRB third stage 

all_param_likeli_cond_2nd<-cbind( all_param_likeli_cond_2nd,rep(1:10,each=iterations-burnin) )

colnames(all_param_likeli_cond_2nd)[5]<-"chain"
all_param_likeli_cond_2nd<-cbind(1:nrow(all_param_likeli_cond_2nd),all_param_likeli_cond_2nd)
colnames(all_param_likeli_cond_2nd)[1]<-"idxidx"

PPPRB_log_s2_3rd_cold<-numeric(iterations)
PPPRB_log_s2_3rd_hot_1<-numeric(iterations)
PPPRB_log_s2_3rd_hot_2<-numeric(iterations)
PPPRB_log_s2_3rd_hot_3<-numeric(iterations)
PPPRB_log_s2_3rd_hot_4<-numeric(iterations)
PPPRB_log_s2_3rd_hot_5<-numeric(iterations)
PPPRB_log_s2_3rd_hot_6<-numeric(iterations)
PPPRB_log_s2_3rd_hot_7<-numeric(iterations)
PPPRB_log_s2_3rd_hot_8<-numeric(iterations)
PPPRB_log_s2_3rd_hot_9<-numeric(iterations)
#
PPPRB_log_n2_3rd_cold<-numeric(iterations)
PPPRB_log_n2_3rd_hot_1<-numeric(iterations)
PPPRB_log_n2_3rd_hot_2<-numeric(iterations)
PPPRB_log_n2_3rd_hot_3<-numeric(iterations)
PPPRB_log_n2_3rd_hot_4<-numeric(iterations)
PPPRB_log_n2_3rd_hot_5<-numeric(iterations)
PPPRB_log_n2_3rd_hot_6<-numeric(iterations)
PPPRB_log_n2_3rd_hot_7<-numeric(iterations)
PPPRB_log_n2_3rd_hot_8<-numeric(iterations)
PPPRB_log_n2_3rd_hot_9<-numeric(iterations)
#
PPPRB_log_phi_3rd_cold<-numeric(iterations)
PPPRB_log_phi_3rd_hot_1<-numeric(iterations)
PPPRB_log_phi_3rd_hot_2<-numeric(iterations)
PPPRB_log_phi_3rd_hot_3<-numeric(iterations)
PPPRB_log_phi_3rd_hot_4<-numeric(iterations)
PPPRB_log_phi_3rd_hot_5<-numeric(iterations)
PPPRB_log_phi_3rd_hot_6<-numeric(iterations)
PPPRB_log_phi_3rd_hot_7<-numeric(iterations)
PPPRB_log_phi_3rd_hot_8<-numeric(iterations)
PPPRB_log_phi_3rd_hot_9<-numeric(iterations)

# acceptance ratio
accept_ratio_3rd_cold<-numeric(iterations)
accept_ratio_3rd_hot_1<-numeric(iterations)
accept_ratio_3rd_hot_2<-numeric(iterations)
accept_ratio_3rd_hot_3<-numeric(iterations)
accept_ratio_3rd_hot_4<-numeric(iterations)
accept_ratio_3rd_hot_5<-numeric(iterations)
accept_ratio_3rd_hot_6<-numeric(iterations)
accept_ratio_3rd_hot_7<-numeric(iterations)
accept_ratio_3rd_hot_8<-numeric(iterations)
accept_ratio_3rd_hot_9<-numeric(iterations)


## result from the 2nd stage
PPPRB_log_s2_2nd_list<-c(PPPRB_log_s2_2nd_cold,
                         PPPRB_log_s2_2nd_hot_1,
                         PPPRB_log_s2_2nd_hot_2,
                         PPPRB_log_s2_2nd_hot_3,
                         PPPRB_log_s2_2nd_hot_4,
                         PPPRB_log_s2_2nd_hot_5,
                         PPPRB_log_s2_2nd_hot_6,
                         PPPRB_log_s2_2nd_hot_7,
                         PPPRB_log_s2_2nd_hot_8,
                         PPPRB_log_s2_2nd_hot_9)

PPPRB_log_n2_2nd_list<-c(PPPRB_log_n2_2nd_cold,
                         PPPRB_log_n2_2nd_hot_1,
                         PPPRB_log_n2_2nd_hot_2,
                         PPPRB_log_n2_2nd_hot_3,
                         PPPRB_log_n2_2nd_hot_4,
                         PPPRB_log_n2_2nd_hot_5,
                         PPPRB_log_n2_2nd_hot_6,
                         PPPRB_log_n2_2nd_hot_7,
                         PPPRB_log_n2_2nd_hot_8,
                         PPPRB_log_n2_2nd_hot_9)

PPPRB_log_phi_2nd_list<-c(PPPRB_log_phi_2nd_cold,
                          PPPRB_log_phi_2nd_hot_1,
                          PPPRB_log_phi_2nd_hot_2,
                          PPPRB_log_phi_2nd_hot_3,
                          PPPRB_log_phi_2nd_hot_4,
                          PPPRB_log_phi_2nd_hot_5,
                          PPPRB_log_phi_2nd_hot_6,
                          PPPRB_log_phi_2nd_hot_7,
                          PPPRB_log_phi_2nd_hot_8,
                          PPPRB_log_phi_2nd_hot_9)
## initial
iterations-burnin
current_idx<-c(100,
               (iterations-burnin)+100,
               (iterations-burnin)*2+100,
               (iterations-burnin)*3+100,
               (iterations-burnin)*4+100,
               (iterations-burnin)*5+100,
               (iterations-burnin)*6+100,
               (iterations-burnin)*7+100,
               (iterations-burnin)*8+100,
               (iterations-burnin)*9+100
)

PPPRB_log_s2_3rd_cold[1]<-PPPRB_log_s2_2nd_list[current_idx[1]]
PPPRB_log_s2_3rd_hot_1[1]<-PPPRB_log_s2_2nd_list[current_idx[2]]
PPPRB_log_s2_3rd_hot_2[1]<-PPPRB_log_s2_2nd_list[current_idx[3]]
PPPRB_log_s2_3rd_hot_3[1]<-PPPRB_log_s2_2nd_list[current_idx[4]]
PPPRB_log_s2_3rd_hot_4[1]<-PPPRB_log_s2_2nd_list[current_idx[5]]
PPPRB_log_s2_3rd_hot_5[1]<-PPPRB_log_s2_2nd_list[current_idx[6]]
PPPRB_log_s2_3rd_hot_6[1]<-PPPRB_log_s2_2nd_list[current_idx[7]]
PPPRB_log_s2_3rd_hot_7[1]<-PPPRB_log_s2_2nd_list[current_idx[8]]
PPPRB_log_s2_3rd_hot_8[1]<-PPPRB_log_s2_2nd_list[current_idx[9]]
PPPRB_log_s2_3rd_hot_9[1]<-PPPRB_log_s2_2nd_list[current_idx[10]]


## initial 
PPPRB_log_n2_3rd_cold[1]<-PPPRB_log_n2_2nd_list[current_idx[1]]
PPPRB_log_n2_3rd_hot_1[1]<-PPPRB_log_n2_2nd_list[current_idx[2]]
PPPRB_log_n2_3rd_hot_2[1]<-PPPRB_log_n2_2nd_list[current_idx[3]]
PPPRB_log_n2_3rd_hot_3[1]<-PPPRB_log_n2_2nd_list[current_idx[4]]
PPPRB_log_n2_3rd_hot_4[1]<-PPPRB_log_n2_2nd_list[current_idx[5]]
PPPRB_log_n2_3rd_hot_5[1]<-PPPRB_log_n2_2nd_list[current_idx[6]]
PPPRB_log_n2_3rd_hot_6[1]<-PPPRB_log_n2_2nd_list[current_idx[7]]
PPPRB_log_n2_3rd_hot_7[1]<-PPPRB_log_n2_2nd_list[current_idx[8]]
PPPRB_log_n2_3rd_hot_8[1]<-PPPRB_log_n2_2nd_list[current_idx[9]]
PPPRB_log_n2_3rd_hot_9[1]<-PPPRB_log_n2_2nd_list[current_idx[10]]

## initial 
PPPRB_log_phi_3rd_cold[1]<-PPPRB_log_phi_2nd_list[current_idx[[1]]]
PPPRB_log_phi_3rd_hot_1[1]<-PPPRB_log_phi_2nd_list[current_idx[[2]]]
PPPRB_log_phi_3rd_hot_2[1]<-PPPRB_log_phi_2nd_list[current_idx[[3]]]
PPPRB_log_phi_3rd_hot_3[1]<-PPPRB_log_phi_2nd_list[current_idx[[4]]]
PPPRB_log_phi_3rd_hot_4[1]<-PPPRB_log_phi_2nd_list[current_idx[[5]]]
PPPRB_log_phi_3rd_hot_5[1]<-PPPRB_log_phi_2nd_list[current_idx[[6]]]
PPPRB_log_phi_3rd_hot_6[1]<-PPPRB_log_phi_2nd_list[current_idx[[7]]]
PPPRB_log_phi_3rd_hot_7[1]<-PPPRB_log_phi_2nd_list[current_idx[[8]]]
PPPRB_log_phi_3rd_hot_8[1]<-PPPRB_log_phi_2nd_list[current_idx[[9]]]
PPPRB_log_phi_3rd_hot_9[1]<-PPPRB_log_phi_2nd_list[current_idx[[10]]]



j<-2


accept_ratio_swap_3rd<-numeric(iterations)


all_param_likeli_cond_2nd[,c(2:4)]==all_param_likeli_2nd[,c(1:3)]

chain_rows <- split(
  seq_len(nrow(all_param_likeli_cond_2nd)),
  all_param_likeli_cond_2nd[, "chain"]   # 
)


tic()
j<-2
for(j in 2: iterations){
  current_log_s2 <- c(
    PPPRB_log_s2_3rd_cold[j-1],
    PPPRB_log_s2_3rd_hot_1[j-1],
    PPPRB_log_s2_3rd_hot_2[j-1],
    PPPRB_log_s2_3rd_hot_3[j-1],
    PPPRB_log_s2_3rd_hot_4[j-1],
    PPPRB_log_s2_3rd_hot_5[j-1],
    PPPRB_log_s2_3rd_hot_6[j-1],
    PPPRB_log_s2_3rd_hot_7[j-1],
    PPPRB_log_s2_3rd_hot_8[j-1],
    PPPRB_log_s2_3rd_hot_9[j-1]
  )
  
  current_log_n2 <- c(
    PPPRB_log_n2_3rd_cold[j-1],
    PPPRB_log_n2_3rd_hot_1[j-1],
    PPPRB_log_n2_3rd_hot_2[j-1],
    PPPRB_log_n2_3rd_hot_3[j-1],
    PPPRB_log_n2_3rd_hot_4[j-1],
    PPPRB_log_n2_3rd_hot_5[j-1],
    PPPRB_log_n2_3rd_hot_6[j-1],
    PPPRB_log_n2_3rd_hot_7[j-1],
    PPPRB_log_n2_3rd_hot_8[j-1],
    PPPRB_log_n2_3rd_hot_9[j-1]
    
  )
  
  current_log_phi <- c(
    PPPRB_log_phi_3rd_cold[j-1],
    PPPRB_log_phi_3rd_hot_1[j-1],
    PPPRB_log_phi_3rd_hot_2[j-1],
    PPPRB_log_phi_3rd_hot_3[j-1],
    PPPRB_log_phi_3rd_hot_4[j-1],
    PPPRB_log_phi_3rd_hot_5[j-1],
    PPPRB_log_phi_3rd_hot_6[j-1],
    PPPRB_log_phi_3rd_hot_7[j-1],
    PPPRB_log_phi_3rd_hot_8[j-1],
    PPPRB_log_phi_3rd_hot_9[j-1]
  )
  
  results_ppprb <- lapply(
    1:10,
    function(i) {
      PPPRB_2nd_stage(
        temperature = temperature_list[[i]],
        all_param_likeli_cond = all_param_likeli_cond_2nd,
        chain_rows=chain_rows,
        current_idx=current_idx[i],
        chain_idx=c(1:10)[i],
        #
        PPPRB_log_s2_1st = PPPRB_log_s2_2nd_list,
        PPPRB_log_n2_1st  = PPPRB_log_n2_2nd_list,
        PPPRB_log_phi_1st   = PPPRB_log_phi_2nd_list,
        #
        current_log_s2  = current_log_s2[i],
        current_log_n2    = current_log_n2[i],
        current_log_phi   = current_log_phi[i]
      )
    }
  )
  names(results_ppprb) <- c("cold", "hot_1", "hot_2", "hot_3","hot_4","hot_5",
                            "hot_6","hot_7","hot_8","hot_9")
  #
  log_s2_final    <- sapply(results_ppprb, `[[`, "log_s2_final")
  log_n2_final  <- sapply(results_ppprb, `[[`, "log_n2_final")
  log_phi_final   <- sapply(results_ppprb, `[[`, "log_phi_final")
  accept      <- sapply(results_ppprb, `[[`, "accept_or_not")
  idx_final      <- sapply(results_ppprb, `[[`, "idx_result")
  #
  PPPRB_log_s2_3rd_cold[j]    <- log_s2_final["cold"]
  PPPRB_log_s2_3rd_hot_1[j]   <- log_s2_final["hot_1"]
  PPPRB_log_s2_3rd_hot_2[j]   <- log_s2_final["hot_2"]
  PPPRB_log_s2_3rd_hot_3[j]   <- log_s2_final["hot_3"]
  PPPRB_log_s2_3rd_hot_4[j]   <- log_s2_final["hot_4"]
  PPPRB_log_s2_3rd_hot_5[j]   <- log_s2_final["hot_5"]
  PPPRB_log_s2_3rd_hot_6[j]   <- log_s2_final["hot_6"]
  PPPRB_log_s2_3rd_hot_7[j]   <- log_s2_final["hot_7"]
  PPPRB_log_s2_3rd_hot_8[j]   <- log_s2_final["hot_8"]
  PPPRB_log_s2_3rd_hot_9[j]   <- log_s2_final["hot_9"]
  
  PPPRB_log_n2_3rd_cold[j]  <- log_n2_final["cold"]
  PPPRB_log_n2_3rd_hot_1[j] <- log_n2_final["hot_1"]
  PPPRB_log_n2_3rd_hot_2[j] <- log_n2_final["hot_2"]
  PPPRB_log_n2_3rd_hot_3[j] <- log_n2_final["hot_3"]
  PPPRB_log_n2_3rd_hot_4[j] <- log_n2_final["hot_4"]
  PPPRB_log_n2_3rd_hot_5[j] <- log_n2_final["hot_5"]
  PPPRB_log_n2_3rd_hot_6[j] <- log_n2_final["hot_6"]
  PPPRB_log_n2_3rd_hot_7[j] <- log_n2_final["hot_7"]
  PPPRB_log_n2_3rd_hot_8[j] <- log_n2_final["hot_8"]
  PPPRB_log_n2_3rd_hot_9[j] <- log_n2_final["hot_9"]
  
  PPPRB_log_phi_3rd_cold[j]   <- log_phi_final["cold"]
  PPPRB_log_phi_3rd_hot_1[j]  <- log_phi_final["hot_1"]
  PPPRB_log_phi_3rd_hot_2[j]  <- log_phi_final["hot_2"]
  PPPRB_log_phi_3rd_hot_3[j]  <- log_phi_final["hot_3"]
  PPPRB_log_phi_3rd_hot_4[j]  <- log_phi_final["hot_4"]
  PPPRB_log_phi_3rd_hot_5[j]  <- log_phi_final["hot_5"]
  PPPRB_log_phi_3rd_hot_6[j]  <- log_phi_final["hot_6"]
  PPPRB_log_phi_3rd_hot_7[j]  <- log_phi_final["hot_7"]
  PPPRB_log_phi_3rd_hot_8[j]  <- log_phi_final["hot_8"]
  PPPRB_log_phi_3rd_hot_9[j]  <- log_phi_final["hot_9"]
  
  accept_ratio_3rd_cold[j]  <- accept["cold"]
  accept_ratio_3rd_hot_1[j] <- accept["hot_1"]
  accept_ratio_3rd_hot_2[j] <- accept["hot_2"]
  accept_ratio_3rd_hot_3[j] <- accept["hot_3"]
  accept_ratio_3rd_hot_4[j] <- accept["hot_4"]
  accept_ratio_3rd_hot_5[j] <- accept["hot_5"]
  accept_ratio_3rd_hot_6[j] <- accept["hot_6"]
  accept_ratio_3rd_hot_7[j] <- accept["hot_7"]
  accept_ratio_3rd_hot_8[j] <- accept["hot_8"]
  accept_ratio_3rd_hot_9[j] <- accept["hot_9"]
  ## current_idx_update
  current_idx<-as.numeric(idx_final)
  
  #
  if (j %% 1000 == 0) cat(j, "\n")
  # swapping
  # swapping
  # swapping
  idx<-sample(2:10,1,replace=FALSE)
  if(idx==2){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_1[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_1[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_1[j]
  } else if(idx==3){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_2[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_2[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_2[j]
  } else if(idx==4){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_3[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_3[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_3[j]
    
  } else if(idx==5){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_4[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_4[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_4[j]
  } else if(idx==6){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_5[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_5[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_5[j]
  } else if(idx==7){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_6[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_6[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_6[j]
  } else if(idx==8){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_7[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_7[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_7[j]
  } else if(idx==9){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_8[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_8[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_8[j]
  } else if(idx==10){
    PPPRB_log_s2_hot<- PPPRB_log_s2_3rd_hot_9[j]
    PPPRB_log_n2_hot<- PPPRB_log_n2_3rd_hot_9[j]
    PPPRB_log_phi_hot<- PPPRB_log_phi_3rd_hot_9[j]
  }
  
  swap_result<-PPPRB_swap(temperature = temperature_list[[idx]],
                          #
                          all_param_likeli_cond=all_param_likeli_cond_2nd,
                          all_param_likeli_1st=all_param_likeli_2nd,
                          #
                          cold_idx=current_idx[1],
                          hot_idx=current_idx[idx],
                          
                          PPPRB_log_s2_cold=PPPRB_log_s2_3rd_cold[j],
                          PPPRB_log_s2_hot=PPPRB_log_s2_hot,
                          #
                          PPPRB_log_n2_cold=PPPRB_log_n2_3rd_cold[j],
                          PPPRB_log_n2_hot=PPPRB_log_n2_hot,
                          #
                          PPPRB_log_phi_cold=PPPRB_log_phi_3rd_cold[j],
                          PPPRB_log_phi_hot=PPPRB_log_phi_hot)
  if(swap_result$accept_or_not==1){
    PPPRB_log_s2_3rd_cold[j] <-swap_result$log_s2_cold_chain
    PPPRB_log_n2_3rd_cold[j]<-swap_result$log_n2_cold_chain
    PPPRB_log_phi_3rd_cold[j]<-swap_result$log_phi_cold_chain
    accept_ratio_swap_3rd[j]<-1
    current_idx_cold<-current_idx[1]
    current_idx_hot<-current_idx[idx]
    current_idx[1]<-current_idx_hot
    #
    if(idx==2){
      PPPRB_log_s2_3rd_hot_1[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_1[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_1[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==3){
      PPPRB_log_s2_3rd_hot_2[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_2[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_2[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==4){
      PPPRB_log_s2_3rd_hot_3[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_3[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_3[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==5){
      PPPRB_log_s2_3rd_hot_4[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_4[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_4[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==6){
      PPPRB_log_s2_3rd_hot_5[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_5[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_5[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==7){
      PPPRB_log_s2_3rd_hot_6[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_6[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_6[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==8){
      PPPRB_log_s2_3rd_hot_7[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_7[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_7[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==9){
      PPPRB_log_s2_3rd_hot_8[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_8[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_8[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==10){
      PPPRB_log_s2_3rd_hot_9[j]<-swap_result$log_s2_hot_chain
      PPPRB_log_n2_3rd_hot_9[j]<-swap_result$log_n2_hot_chain
      PPPRB_log_phi_3rd_hot_9[j]<-swap_result$log_phi_hot_chain
      current_idx[idx]<-current_idx_cold
    }
  }
}

toc()


mean(accept_ratio_swap_3rd)
mean( accept_ratio_3rd_cold )
mean( accept_ratio_3rd_hot_1 )
mean( accept_ratio_3rd_hot_2 )
mean( accept_ratio_3rd_hot_3 )
mean( accept_ratio_3rd_hot_4 )
mean( accept_ratio_3rd_hot_5 )
mean( accept_ratio_3rd_hot_6 )
mean( accept_ratio_3rd_hot_7 )
mean( accept_ratio_3rd_hot_8 )
mean( accept_ratio_3rd_hot_9 )


PPPRB_log_s2_3rd_cold<-PPPRB_log_s2_3rd_cold[-c(1:burnin)]
PPPRB_log_n2_3rd_cold<-PPPRB_log_n2_3rd_cold[-c(1:burnin)]
PPPRB_log_phi_3rd_cold<-PPPRB_log_phi_3rd_cold[-c(1:burnin)]



#save.image("sss_ppprb.RData")


library(posterior)
draws_arr <- array(NA, dim = c(length(log_sigma_n2_post), 2, 1))
draws_arr[, 1, 1] <- log_sigma_n2_post
draws_arr[, 2, 1] <- PPPRB_log_n2_3rd_cold
dimnames(draws_arr) <- list(NULL, NULL, c("log_sigma_n2"))
d <- as_draws_array(draws_arr)
summarise_draws(d)


draws_arr <- array(NA, dim = c(length(log_sigma_s2_post), 2, 1))
draws_arr[, 1, 1] <- log_sigma_s2_post
draws_arr[, 2, 1] <- PPPRB_log_s2_3rd_cold
dimnames(draws_arr) <- list(NULL, NULL, c("log_sigma_s2"))
d <- as_draws_array(draws_arr)
summarise_draws(d)



draws_arr <- array(NA, dim = c(length(log_phi_post), 2, 1))
draws_arr[, 1, 1] <- log_phi_post
draws_arr[, 2, 1] <- PPPRB_log_phi_3rd_cold
dimnames(draws_arr) <- list(NULL, NULL, c("log_phi"))
d <- as_draws_array(draws_arr)
summarise_draws(d)

library(coda)

combined <- mcmc.list(as.mcmc(log_phi_post), as.mcmc(PPPRB_log_phi_3rd_cold))
gelman.diag(combined)

combined <- mcmc.list(as.mcmc(log_sigma_n2_post), as.mcmc(PPPRB_log_n2_3rd_cold))
gelman.diag(combined)

combined <- mcmc.list(as.mcmc(log_sigma_s2_post), as.mcmc(PPPRB_log_s2_3rd_cold))
gelman.diag(combined)


### PPRB


combined <- mcmc.list(as.mcmc(log_phi_post), as.mcmc(log_phi_post_pprb_3))
gelman.diag(combined)

combined <- mcmc.list(as.mcmc(log_sigma_n2_post), as.mcmc(log_sigma_n2_post_pprb_3))
gelman.diag(combined)

combined <- mcmc.list(as.mcmc(log_sigma_s2_post), as.mcmc(log_sigma_s2_post_pprb_3))
gelman.diag(combined)


draws_arr <- array(NA, dim = c(length(log_sigma_n2_post), 2, 1))
draws_arr[, 1, 1] <- log_sigma_n2_post
draws_arr[, 2, 1] <- log_sigma_n2_post_pprb_3
dimnames(draws_arr) <- list(NULL, NULL, c("log_sigma_n2"))
d <- as_draws_array(draws_arr)
summarise_draws(d)


draws_arr <- array(NA, dim = c(length(log_sigma_s2_post), 2, 1))
draws_arr[, 1, 1] <- log_sigma_s2_post
draws_arr[, 2, 1] <- log_sigma_s2_post_pprb_3
dimnames(draws_arr) <- list(NULL, NULL, c("log_sigma_s2"))
d <- as_draws_array(draws_arr)
summarise_draws(d)



draws_arr <- array(NA, dim = c(length(log_phi_post), 2, 1))
draws_arr[, 1, 1] <- log_phi_post
draws_arr[, 2, 1] <- log_phi_post_pprb_3
dimnames(draws_arr) <- list(NULL, NULL, c("log_phi"))
d <- as_draws_array(draws_arr)
summarise_draws(d)

# plots for the paper -----------------------------------------------------

summary_stat <- function(x) {
  c(
    mean = mean(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975)
  )
}
library(dplyr)
library(tidyr)


log_phi_df <- bind_rows(
  data.frame(stage = "PP-RB 1st",  parameter = "log(phi)",  as.list(summary_stat(log_phi_post_pprb_1))),
  data.frame(stage = "PP-RB 2nd",  parameter = "log(phi)",  as.list(summary_stat(log_phi_post_pprb_2))),
  data.frame(stage = "PP-RB 3rd",  parameter = "log(phi)",  as.list(summary_stat(log_phi_post_pprb_3))),
  data.frame(stage = "PPPRB 1st",  parameter = "log(phi)",  as.list(summary_stat(results$cold$log_phi_post_pprb_1))),
  data.frame(stage = "PPPRB 2nd",  parameter = "log(phi)",  as.list(summary_stat(PPPRB_log_phi_2nd_cold))),
  data.frame(stage = "PPPRB 3rd",  parameter = "log(phi)",  as.list(summary_stat(PPPRB_log_phi_3rd_cold))),
  data.frame(stage = "Full",       parameter = "log(phi)",  as.list(summary_stat(log_phi_post)))
)


log_sigma_n2_df <- bind_rows(
  data.frame(stage = "PP-RB 1st",  parameter = "beta",  as.list(summary_stat(log_sigma_n2_post_pprb_1))),
  data.frame(stage = "PP-RB 2nd",  parameter = "beta",  as.list(summary_stat(log_sigma_n2_post_pprb_2))),
  data.frame(stage = "PP-RB 3rd",  parameter = "beta",  as.list(summary_stat(log_sigma_n2_post_pprb_3))),
  data.frame(stage = "PPPRB 1st",  parameter = "beta",  as.list(summary_stat(results$cold$log_sigma_n2_post_pprb_1))),
  data.frame(stage = "PPPRB 2nd",  parameter = "beta",  as.list(summary_stat(PPPRB_log_n2_2nd_cold))),
  data.frame(stage = "PPPRB 3rd",  parameter = "log(phi)",  as.list(summary_stat(PPPRB_log_n2_3rd_cold))),
  data.frame(stage = "Full",       parameter = "beta",  as.list(summary_stat(log_sigma_n2_post)))
)

log_sigma_s2_df <- bind_rows(
  data.frame(stage = "PP-RB 1st",  parameter = "beta",  as.list(summary_stat(log_sigma_s2_post_pprb_1))),
  data.frame(stage = "PP-RB 2nd",  parameter = "beta",  as.list(summary_stat(log_sigma_s2_post_pprb_2))),
  data.frame(stage = "PP-RB 3rd",  parameter = "beta",  as.list(summary_stat(log_sigma_s2_post_pprb_3))),
  data.frame(stage = "PPPRB 1st",  parameter = "beta",  as.list(summary_stat(results$cold$log_sigma_s2_post_pprb_1))),
  data.frame(stage = "PPPRB 2nd",  parameter = "beta",  as.list(summary_stat(PPPRB_log_s2_2nd_cold))),
  data.frame(stage = "PPPRB 3rd",  parameter = "log(phi)",  as.list(summary_stat(PPPRB_log_s2_3rd_cold))),
  data.frame(stage = "Full",       parameter = "beta",  as.list(summary_stat(log_sigma_s2_post)))
)

colnames(log_phi_df)[c(4:5)]<-c("lower","upper")
colnames(log_sigma_n2_df)[c(4:5)]<-c("lower","upper")
colnames(log_sigma_s2_df)[c(4:5)]<-c("lower","upper")


# 1.
p_log_phi <- ggplot(log_phi_df, aes(x = stage, y = mean, color = stage)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.8) +
  labs(title = expression(log(phi)), x = NULL, y = "Value")+
  scale_color_manual(values = c("black", "red", "orange", "green", "blue", "purple", "cyan")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. 
p_log_sigma_n2 <- ggplot(log_sigma_n2_df, aes(x = stage, y = mean, color = stage)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.8) +
  labs(title = expression(log(sigma[n]^2)), x = NULL, y = "Value")+
  scale_color_manual(values = c("black", "red", "orange", "green", "blue", "purple", "cyan")) +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1))

# 3.
p_log_sigma_s2 <- ggplot(log_sigma_s2_df, aes(x = stage, y = mean, color = stage)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.8) +
  labs(title = expression(log(sigma[s]^2)), x = NULL, y = "Value")+
  scale_color_manual(values = c("black", "red", "orange", "green", "blue", "purple", "cyan")) +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1))

library(patchwork)
common_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  aspect.ratio = 1
)

#fig1_2 <- fig1 + common_theme
p_log_phi_2 <- p_log_phi + common_theme
p_log_sigma_n2_2 <- p_log_sigma_n2 + common_theme
p_log_sigma_s2_2 <- p_log_sigma_s2 + common_theme
fig1_2<-fig1+common_theme
p_log_phi_2+p_log_sigma_n2_2+p_log_sigma_s2_2
fig1_2+p_log_phi_2+p_log_sigma_n2_2+p_log_sigma_s2_2
# 1200 auto
fig1 /
  (p_log_phi_2 + p_log_sigma_n2_2 + p_log_sigma_s2_2)

p_log_phi_2  <- p_log_phi_2  + theme(legend.position = "none")+ scale_color_brewer(palette = "Dark2")
p_log_sigma_n2_2 <- p_log_sigma_n2_2 + theme(legend.position = "none")+ scale_color_brewer(palette = "Dark2")
p_log_sigma_s2_2 <- p_log_sigma_s2_2 + theme(legend.position = "none")+ scale_color_brewer(palette = "Dark2")

fig1 /
  (p_log_phi_2 + p_log_sigma_n2_2 + p_log_sigma_s2_2) +
  plot_layout(guides = "collect",heights = c(1, 1)) &theme(legend.position = "right")

# ESS ---------------------------------------------------------------------



library(coda)
single_log_sigma_s2<-mcmc(log_sigma_s2_post)
single_log_sigma_n2<-mcmc(log_sigma_n2_post)
single_log_phi<-mcmc(log_phi_post)


pprb_log_sigma_s2<-mcmc(log_sigma_s2_post_pprb_3)
pprb_log_sigma_n2<-mcmc(log_sigma_n2_post_pprb_3)
pprb_log_phi<-mcmc(log_phi_post_pprb_3)

ppprb_log_sigma_s2<-mcmc(PPPRB_log_s2_3rd_cold)
ppprb_log_sigma_n2<-mcmc(PPPRB_log_n2_3rd_cold)
ppprb_log_phi<-mcmc(PPPRB_log_phi_3rd_cold)


mcmc_list <- mcmc.list(single_log_sigma_s2,
                       single_log_sigma_n2,
                       single_log_phi,
                       #
                       pprb_log_sigma_s2,
                       pprb_log_sigma_n2,
                       pprb_log_phi,
                       #
                       ppprb_log_sigma_s2,
                       ppprb_log_sigma_n2,
                       ppprb_log_phi)

ess_sigma_each <- sapply(mcmc_list, effectiveSize)

names(ess_sigma_each)<-c("single_log_sigma_s2",
                         "single_log_sigma_n2",
                         "single_log_phi",
                         #
                         "pprb_log_sigma_s2",
                         "pprb_log_sigma_n2",
                         "pprb_log_phi",
                         #
                         "ppprb_log_sigma_s2",
                         "ppprb_log_sigma_n2",
                         "ppprb_log_phi")
                         
round(ess_sigma_each,1)
round(2685.367/ess_sigma_each[1:3],1)
round(309.98/ess_sigma_each[4:6],1)
round(828.758/ess_sigma_each[7:9],1)


