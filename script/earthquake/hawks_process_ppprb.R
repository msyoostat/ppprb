library(parallel)
library(tictoc)
library(ggplot2)

#1989-10-18 01:04 UTC
# ==============================
# load data
# ==============================
load("hawks_process_pprb.RData")

RNGkind("L'Ecuyer-CMRG")
set.seed(10) 

# ppprb -------------------------------------------------------------------


PPPRB_1st_stage<-function(temperature, 
                          iterations,
                          burnin,
                          beta_jump_1st,
                          eta_jump_1st,
                          mu_jump_1st,
                          events_1,
                          capital_T,
                          T1){
  
  post_beta_tempered_1st<-numeric(iterations)
  post_eta_tempered_1st<-numeric(iterations)
  post_mu_tempered_1st<-numeric(iterations)
  accept_ratio_1st_tempered<-numeric(iterations)
  #
  log_likelihood_vec<-numeric(iterations)
  #
  post_mu_tempered_1st[1]<-1
  post_beta_tempered_1st[1]<-10
  post_eta_tempered_1st[1]<-0.5
  i<-2
  for(i in 2: iterations){
    mu_current<-post_mu_tempered_1st[i-1]
    beta_current<-post_beta_tempered_1st[i-1]
    eta_current<-post_eta_tempered_1st[i-1]
    
    alpha_current<-eta_current*beta_current
    #
    mu_candidate<- rnorm(1,mu_current,
                         sd=sqrt(mu_jump_1st))
    beta_candidate<- rnorm(1,beta_current,
                           sd=sqrt(beta_jump_1st))
    eta_candidate<- rnorm(1,eta_current,sd=sqrt(eta_jump_1st))
    alpha_candidate<- eta_candidate * beta_candidate
    ## compute below before if statement
    denom_temp<-hawkes_loglik_naive(events_1, mu_current,
                                    alpha_current,
                                    beta_current,
                                    T1)
    
    if(mu_candidate>0 &beta_candidate>0 &eta_candidate>0 &eta_candidate<1){
      # likelihood and prior
      nom_temp<-hawkes_loglik_naive(events_1, mu_candidate,
                                    alpha_candidate,
                                    beta_candidate,
                                    T1) 
      nom<-nom_temp * (1/temperature )
      #
      denom<- denom_temp* (1/temperature )
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
        post_beta_tempered_1st[i]<-beta_candidate
        post_eta_tempered_1st[i]<-eta_candidate
        post_mu_tempered_1st[i]<-mu_candidate
        accept_ratio_1st_tempered[i]<-1
        #
        log_likelihood_vec[i]<-nom_temp
      } else{
        post_beta_tempered_1st[i]<-beta_current
        post_eta_tempered_1st[i]<-eta_current
        post_mu_tempered_1st[i]<-mu_current
        #
        log_likelihood_vec[i]<-denom_temp
        
      }
    } else{
      post_beta_tempered_1st[i]<-beta_current
      post_eta_tempered_1st[i]<-eta_current
      post_mu_tempered_1st[i]<-mu_current
      #
      log_likelihood_vec[i]<-denom_temp
      
    }
  }

  post_beta_tempered_1st<-post_beta_tempered_1st[-c(1:burnin)]
  post_eta_tempered_1st<-post_eta_tempered_1st[-c(1:burnin)]
  post_mu_tempered_1st<-post_mu_tempered_1st[-c(1:burnin)]
  accept_ratio_1st_tempered<-accept_ratio_1st_tempered[-c(1:burnin)]
  log_likelihood_vec<-log_likelihood_vec[-c(1:burnin)]
  return(list(post_beta_tempered_1st=post_beta_tempered_1st,
              post_eta_tempered_1st=post_eta_tempered_1st,
              post_mu_tempered_1st=post_mu_tempered_1st,
              accept_ratio_1st_tempered=accept_ratio_1st_tempered,
              log_likelihood_vec=log_likelihood_vec))
}
iterations
burnin


temperature_list<-list()
temperature_list[[1]]<-1
temperature_list[[2]]<-exp( seq(0,2,length.out=10) )[2]
temperature_list[[3]]<-exp( seq(0,2,length.out=10) )[3]
temperature_list[[4]]<-exp( seq(0,2,length.out=10) )[4]
temperature_list[[5]]<-exp( seq(0,2,length.out=10) )[5]
temperature_list[[6]]<-exp( seq(0,2,length.out=10) )[6]
temperature_list[[7]]<-exp( seq(0,2,length.out=10) )[7]
temperature_list[[8]]<-exp( seq(0,2,length.out=10) )[8]
temperature_list[[9]]<-exp( seq(0,2,length.out=10) )[9]
temperature_list[[10]]<-exp( seq(0,2,length.out=10) )[10]


jump_beta_list<-list()
jump_beta_list[[1]]<-beta_jump
jump_beta_list[[2]]<-seq(beta_jump,3,length.out=10)[2]
jump_beta_list[[3]]<-seq(beta_jump,3,length.out=10)[3]
jump_beta_list[[4]]<-seq(beta_jump,3,length.out=10)[4]
jump_beta_list[[5]]<-seq(beta_jump,3,length.out=10)[5]
jump_beta_list[[6]]<-seq(beta_jump,3,length.out=10)[6]
jump_beta_list[[7]]<-seq(beta_jump,3,length.out=10)[7]
jump_beta_list[[8]]<-seq(beta_jump,3,length.out=10)[8]
jump_beta_list[[9]]<-seq(beta_jump,3,length.out=10)[9]
jump_beta_list[[10]]<-seq(beta_jump,3,length.out=10)[10]

jump_eta_list<-list()
jump_eta_list[[1]]<-eta_jump
jump_eta_list[[2]]<-seq(eta_jump,0.02,length.out=10)[2]
jump_eta_list[[3]]<-seq(eta_jump,0.02,length.out=10)[3]
jump_eta_list[[4]]<-seq(eta_jump,0.02,length.out=10)[4]
jump_eta_list[[5]]<-seq(eta_jump,0.02,length.out=10)[5]
jump_eta_list[[6]]<-seq(eta_jump,0.02,length.out=10)[6]
jump_eta_list[[7]]<-seq(eta_jump,0.02,length.out=10)[7]
jump_eta_list[[8]]<-seq(eta_jump,0.02,length.out=10)[8]
jump_eta_list[[9]]<-seq(eta_jump,0.02,length.out=10)[9]
jump_eta_list[[10]]<-seq(eta_jump,0.02,length.out=10)[10]


jump_mu_list<-list()
jump_mu_list[[1]]<-mu_jump
jump_mu_list[[2]]<-seq(mu_jump,0.020,length.out=10)[2]
jump_mu_list[[3]]<-seq(mu_jump,0.020,length.out=10)[3]
jump_mu_list[[4]]<-seq(mu_jump,0.020,length.out=10)[4]
jump_mu_list[[5]]<-seq(mu_jump,0.020,length.out=10)[5]
jump_mu_list[[6]]<-seq(mu_jump,0.020,length.out=10)[6]
jump_mu_list[[7]]<-seq(mu_jump,0.020,length.out=10)[7]
jump_mu_list[[8]]<-seq(mu_jump,0.020,length.out=10)[8]
jump_mu_list[[9]]<-seq(mu_jump,0.020,length.out=10)[9]
jump_mu_list[[10]]<-seq(mu_jump,0.020,length.out=10)[10]



tic()

results <- mclapply(
  X = 1:10,
  FUN = function(i) {
    PPPRB_1st_stage(
      temperature = temperature_list[[i]],
      iterations  = iterations,
      burnin      = burnin,
      beta_jump_1st = jump_beta_list[[i]],
      eta_jump_1st  = jump_eta_list[[i]],
      mu_jump_1st   = jump_mu_list[[i]],
      events_1   = events_1,
      capital_T = capital_T,
      T1        = T1
    )
  },
  mc.cores = min(length(temperature_list), detectCores()-1)
)

toc()

names(results) <- c("cold", "hot_1", "hot_2", "hot_3","hot_4","hot_5","hot_6",
                    "hot_7","hot_8","hot_9")

mean( results$cold$accept_ratio_1st_tempered)
mean( results$hot_1$accept_ratio_1st_tempered)
mean( results$hot_2$accept_ratio_1st_tempered)
mean( results$hot_3$accept_ratio_1st_tempered)
mean( results$hot_4$accept_ratio_1st_tempered)
mean( results$hot_5$accept_ratio_1st_tempered)
mean( results$hot_6$accept_ratio_1st_tempered)
mean( results$hot_7$accept_ratio_1st_tempered)
mean( results$hot_8$accept_ratio_1st_tempered)
mean( results$hot_9$accept_ratio_1st_tempered)



all_mu   <- c(
  results$cold$post_mu_tempered_1st,
  results$hot_1$post_mu_tempered_1st,
  results$hot_2$post_mu_tempered_1st,
  results$hot_3$post_mu_tempered_1st,
  results$hot_4$post_mu_tempered_1st,
  results$hot_5$post_mu_tempered_1st,
  results$hot_6$post_mu_tempered_1st,
  results$hot_7$post_mu_tempered_1st,
  results$hot_8$post_mu_tempered_1st,
  results$hot_9$post_mu_tempered_1st
)

all_eta  <- c(
  results$cold$post_eta_tempered_1st,
  results$hot_1$post_eta_tempered_1st,
  results$hot_2$post_eta_tempered_1st,
  results$hot_3$post_eta_tempered_1st,
  results$hot_4$post_eta_tempered_1st,
  results$hot_5$post_eta_tempered_1st,
  results$hot_6$post_eta_tempered_1st,
  results$hot_7$post_eta_tempered_1st,
  results$hot_8$post_eta_tempered_1st,
  results$hot_9$post_eta_tempered_1st
)

all_beta <- c(
  results$cold$post_beta_tempered_1st,
  results$hot_1$post_beta_tempered_1st,
  results$hot_2$post_beta_tempered_1st,
  results$hot_3$post_beta_tempered_1st,
  results$hot_4$post_beta_tempered_1st,
  results$hot_5$post_beta_tempered_1st,
  results$hot_6$post_beta_tempered_1st,
  results$hot_7$post_beta_tempered_1st,
  results$hot_8$post_beta_tempered_1st,
  results$hot_9$post_beta_tempered_1st
)

tic()


hawkes_loglik_conditional_all <- mclapply(
  seq_along(all_mu),
  function(i) {
    hawkes_loglik_conditional_fast(
      data$t,
      all_mu[i],
      all_eta[i] * all_beta[i],
      all_beta[i],
      T1,
      T2,
      T2
    )
  },
  mc.cores = detectCores()-1
)
toc()

hawkes_loglik_conditional_all <- unlist(hawkes_loglik_conditional_all)

all_param_likeli_cond<-cbind(all_mu,
                             all_eta,
                             all_beta,
                             hawkes_loglik_conditional_all)

colnames(all_param_likeli_cond)<-c("mu",
                                   "eta",
                                   "beta",
                                   "log_likeli_cond")

rm(hawkes_loglik_conditional_all)

all_param_likeli_1st<-cbind(all_mu,
                            all_eta,
                            all_beta)

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

colnames(all_param_likeli_1st)<-c("mu",
                                  "eta",
                                  "beta",
                                  "log_likeli_1st")

all_param_likeli_cond[,1:3]==all_param_likeli_1st[,1:3]
all_param_likeli_1_2 <-all_param_likeli_cond[,1:3]
all_param_likeli_1_2<-cbind(all_param_likeli_1_2,NA)
all_param_likeli_1_2[,4]<-all_param_likeli_cond[,4]+all_param_likeli_1st[,4]
colnames(all_param_likeli_1_2)<-c("mu",
                                  "eta",
                                  "beta",
                                  "log_likeli")

### PPPRB second stage 
### PPPRB second stage 
### PPPRB second stage 
### PPPRB second stage 

all_param_likeli_cond<-cbind( all_param_likeli_cond,rep(1:10,each=iterations-burnin) )
colnames(all_param_likeli_cond)[5]<-"chain"
all_param_likeli_cond<-cbind(1:nrow(all_param_likeli_cond),all_param_likeli_cond)
colnames(all_param_likeli_cond)[1]<-"idxidx"


PPPRB_2nd_stage<-function(temperature,
                          all_param_likeli_cond,
                          chain_rows,
                          current_idx,
                          chain_idx,
                          #
                          PPPRB_beta_1st,
                          PPPRB_eta_1st,
                          PPPRB_mu_1st,
                          #
                          current_beta,
                          current_mu,
                          current_eta){
  
  rows_chain <- chain_rows[[chain_idx]]   # 
  idx <- sample(rows_chain, 1)            # 
  
  proposal_beta<-PPPRB_beta_1st[idx]
  proposal_eta<-PPPRB_eta_1st[idx]
  proposal_mu<-PPPRB_mu_1st[idx]
  #

  nom<-all_param_likeli_cond[idx,5] * (1/temperature)
  

  denom<-all_param_likeli_cond[current_idx,5] * (1/temperature)

  ratio<-nom-denom 
  ## accept / reject
  if (log(runif(1)) <= ratio) {
    beta_final<-proposal_beta
    eta_final<-proposal_eta
    mu_final<-proposal_mu
    accept_or_not<-1
    idx_result<-idx
  } else {
    beta_final<-current_beta
    eta_final<-current_eta
    mu_final<-current_mu
    accept_or_not<-0
    idx_result<-current_idx
  } 
  return(list(beta_final=beta_final,
              eta_final=eta_final,
              mu_final=mu_final,
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
                     PPPRB_beta_cold,
                     PPPRB_beta_hot,
                     PPPRB_eta_cold,
                     PPPRB_eta_hot,
                     PPPRB_mu_cold,
                     PPPRB_mu_hot){
  

  # cold term
  nom_conditional_log_lik_hot<-all_param_likeli_cond[hot_idx,5]
  nom1<-nom_conditional_log_lik_hot *(1-1/temperature)

  nom_log_lik_1st_hot<-all_param_likeli_1st[hot_idx,4]
  nom2<-nom_log_lik_1st_hot * (1-1/temperature)
  
  denom_conditional_log_lik_cold<-all_param_likeli_cond[cold_idx,5]
  denom1<-denom_conditional_log_lik_cold* (1-1/temperature)
  
  denom_log_lik_1st_cold<-all_param_likeli_1st[cold_idx,4]
  denom2<- denom_log_lik_1st_cold * (1-1/temperature)
  
  #
  ratio<-nom1+nom2-(denom1+denom2)
  ## accept / reject
  if (log(runif(1)) <= ratio) {
    mu_cold_chain<-PPPRB_mu_hot
    beta_cold_chain<-PPPRB_beta_hot
    eta_cold_chain<-PPPRB_eta_hot
    #
    mu_hot_chain<-PPPRB_mu_cold
    beta_hot_chain<-PPPRB_beta_cold
    eta_hot_chain<-PPPRB_eta_cold
    #
    accept_or_not<-1
  } else {
    mu_cold_chain<-PPPRB_mu_cold
    beta_cold_chain<-PPPRB_beta_cold
    eta_cold_chain<-PPPRB_eta_cold
    #
    mu_hot_chain<-PPPRB_mu_hot
    beta_hot_chain<-PPPRB_beta_hot
    eta_hot_chain<-PPPRB_eta_hot
    #
    accept_or_not<-0
  } 
  return(list(mu_cold_chain=mu_cold_chain,
              beta_cold_chain=beta_cold_chain,
              eta_cold_chain=eta_cold_chain,
              #
              mu_hot_chain=mu_hot_chain,
              beta_hot_chain=beta_hot_chain,
              eta_hot_chain=eta_hot_chain,
              accept_or_not=accept_or_not)
  )
}

PPPRB_mu_2nd_cold<-numeric(iterations)
PPPRB_mu_2nd_hot_1<-numeric(iterations)
PPPRB_mu_2nd_hot_2<-numeric(iterations)
PPPRB_mu_2nd_hot_3<-numeric(iterations)
PPPRB_mu_2nd_hot_4<-numeric(iterations)
PPPRB_mu_2nd_hot_5<-numeric(iterations)
PPPRB_mu_2nd_hot_6<-numeric(iterations)
PPPRB_mu_2nd_hot_7<-numeric(iterations)
PPPRB_mu_2nd_hot_8<-numeric(iterations)
PPPRB_mu_2nd_hot_9<-numeric(iterations)

#
PPPRB_beta_2nd_cold<-numeric(iterations)
PPPRB_beta_2nd_hot_1<-numeric(iterations)
PPPRB_beta_2nd_hot_2<-numeric(iterations)
PPPRB_beta_2nd_hot_3<-numeric(iterations)
PPPRB_beta_2nd_hot_4<-numeric(iterations)
PPPRB_beta_2nd_hot_5<-numeric(iterations)
PPPRB_beta_2nd_hot_6<-numeric(iterations)
PPPRB_beta_2nd_hot_7<-numeric(iterations)
PPPRB_beta_2nd_hot_8<-numeric(iterations)
PPPRB_beta_2nd_hot_9<-numeric(iterations)
#
PPPRB_eta_2nd_cold<-numeric(iterations)
PPPRB_eta_2nd_hot_1<-numeric(iterations)
PPPRB_eta_2nd_hot_2<-numeric(iterations)
PPPRB_eta_2nd_hot_3<-numeric(iterations)
PPPRB_eta_2nd_hot_4<-numeric(iterations)
PPPRB_eta_2nd_hot_5<-numeric(iterations)
PPPRB_eta_2nd_hot_6<-numeric(iterations)
PPPRB_eta_2nd_hot_7<-numeric(iterations)
PPPRB_eta_2nd_hot_8<-numeric(iterations)
PPPRB_eta_2nd_hot_9<-numeric(iterations)

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
PPPRB_mu_1st_list<-c(results$cold$post_mu_tempered_1st,
                     results$hot_1$post_mu_tempered_1st,
                     results$hot_2$post_mu_tempered_1st,
                     results$hot_3$post_mu_tempered_1st,
                     results$hot_4$post_mu_tempered_1st,
                     results$hot_5$post_mu_tempered_1st,
                     results$hot_6$post_mu_tempered_1st,
                     results$hot_7$post_mu_tempered_1st,
                     results$hot_8$post_mu_tempered_1st,
                     results$hot_9$post_mu_tempered_1st)

PPPRB_beta_1st_list<-c(results$cold$post_beta_tempered_1st,
                       results$hot_1$post_beta_tempered_1st,
                       results$hot_2$post_beta_tempered_1st,
                       results$hot_3$post_beta_tempered_1st,
                       results$hot_4$post_beta_tempered_1st,
                       results$hot_5$post_beta_tempered_1st,
                       results$hot_6$post_beta_tempered_1st,
                       results$hot_7$post_beta_tempered_1st,
                       results$hot_8$post_beta_tempered_1st,
                       results$hot_9$post_beta_tempered_1st)


PPPRB_eta_1st_list<-c(results$cold$post_eta_tempered_1st,
                      results$hot_1$post_eta_tempered_1st,
                      results$hot_2$post_eta_tempered_1st,
                      results$hot_3$post_eta_tempered_1st,
                      results$hot_4$post_eta_tempered_1st,
                      results$hot_5$post_eta_tempered_1st,
                      results$hot_6$post_eta_tempered_1st,
                      results$hot_7$post_eta_tempered_1st,
                      results$hot_8$post_eta_tempered_1st,
                      results$hot_9$post_eta_tempered_1st)


## initial

current_idx<-c(100,
               (iterations-burnin)+100,
               2*(iterations-burnin)+100,
               3*(iterations-burnin)+100,
               4*(iterations-burnin)+100,
               5*(iterations-burnin)+100,
               6*(iterations-burnin)+100,
               7*(iterations-burnin)+100,
               8*(iterations-burnin)+100,
               9*(iterations-burnin)+100)


PPPRB_mu_2nd_cold[1]<-PPPRB_mu_1st_list[current_idx[1]]
PPPRB_mu_2nd_hot_1[1]<-PPPRB_mu_1st_list[current_idx[2]]
PPPRB_mu_2nd_hot_2[1]<-PPPRB_mu_1st_list[current_idx[3]]
PPPRB_mu_2nd_hot_3[1]<-PPPRB_mu_1st_list[current_idx[4]]
PPPRB_mu_2nd_hot_4[1]<-PPPRB_mu_1st_list[current_idx[5]]
PPPRB_mu_2nd_hot_5[1]<-PPPRB_mu_1st_list[current_idx[6]]
PPPRB_mu_2nd_hot_6[1]<-PPPRB_mu_1st_list[current_idx[7]]
PPPRB_mu_2nd_hot_7[1]<-PPPRB_mu_1st_list[current_idx[8]]
PPPRB_mu_2nd_hot_8[1]<-PPPRB_mu_1st_list[current_idx[9]]
PPPRB_mu_2nd_hot_9[1]<-PPPRB_mu_1st_list[current_idx[10]]

## initial 
PPPRB_beta_2nd_cold[1]<-PPPRB_beta_1st_list[current_idx[1]]
PPPRB_beta_2nd_hot_1[1]<-PPPRB_beta_1st_list[current_idx[2]]
PPPRB_beta_2nd_hot_2[1]<-PPPRB_beta_1st_list[current_idx[3]]
PPPRB_beta_2nd_hot_3[1]<-PPPRB_beta_1st_list[current_idx[4]]
PPPRB_beta_2nd_hot_4[1]<-PPPRB_beta_1st_list[current_idx[5]]
PPPRB_beta_2nd_hot_5[1]<-PPPRB_beta_1st_list[current_idx[6]]
PPPRB_beta_2nd_hot_6[1]<-PPPRB_beta_1st_list[current_idx[7]]
PPPRB_beta_2nd_hot_7[1]<-PPPRB_beta_1st_list[current_idx[8]]
PPPRB_beta_2nd_hot_8[1]<-PPPRB_beta_1st_list[current_idx[9]]
PPPRB_beta_2nd_hot_9[1]<-PPPRB_beta_1st_list[current_idx[10]]

## initial 
PPPRB_eta_2nd_cold[1]<-PPPRB_eta_1st_list[current_idx[[1]]]
PPPRB_eta_2nd_hot_1[1]<-PPPRB_eta_1st_list[current_idx[[2]]]
PPPRB_eta_2nd_hot_2[1]<-PPPRB_eta_1st_list[current_idx[[3]]]
PPPRB_eta_2nd_hot_3[1]<-PPPRB_eta_1st_list[current_idx[[4]]]
PPPRB_eta_2nd_hot_4[1]<-PPPRB_eta_1st_list[current_idx[[5]]]
PPPRB_eta_2nd_hot_5[1]<-PPPRB_eta_1st_list[current_idx[[6]]]
PPPRB_eta_2nd_hot_6[1]<-PPPRB_eta_1st_list[current_idx[[7]]]
PPPRB_eta_2nd_hot_7[1]<-PPPRB_eta_1st_list[current_idx[[8]]]
PPPRB_eta_2nd_hot_8[1]<-PPPRB_eta_1st_list[current_idx[[9]]]
PPPRB_eta_2nd_hot_9[1]<-PPPRB_eta_1st_list[current_idx[[10]]]


# acceptance ratio for swapping
accept_ratio_swap<-numeric(iterations)

j<-2

chain_rows <- split(
  seq_len(nrow(all_param_likeli_cond)),
  all_param_likeli_cond[, "chain"]   #
)


tic()

for(j in 2: iterations){

  current_mu <- c(
    PPPRB_mu_2nd_cold[j-1],
    PPPRB_mu_2nd_hot_1[j-1],
    PPPRB_mu_2nd_hot_2[j-1],
    PPPRB_mu_2nd_hot_3[j-1],
    PPPRB_mu_2nd_hot_4[j-1],
    PPPRB_mu_2nd_hot_5[j-1],
    PPPRB_mu_2nd_hot_6[j-1],
    PPPRB_mu_2nd_hot_7[j-1],
    PPPRB_mu_2nd_hot_8[j-1],
    PPPRB_mu_2nd_hot_9[j-1]
  )
  
  current_beta <- c(
    PPPRB_beta_2nd_cold[j-1],
    PPPRB_beta_2nd_hot_1[j-1],
    PPPRB_beta_2nd_hot_2[j-1],
    PPPRB_beta_2nd_hot_3[j-1],
    PPPRB_beta_2nd_hot_4[j-1],
    PPPRB_beta_2nd_hot_5[j-1],
    PPPRB_beta_2nd_hot_6[j-1],
    PPPRB_beta_2nd_hot_7[j-1],
    PPPRB_beta_2nd_hot_8[j-1],
    PPPRB_beta_2nd_hot_9[j-1]
  )
  
  current_eta <- c(
    PPPRB_eta_2nd_cold[j-1],
    PPPRB_eta_2nd_hot_1[j-1],
    PPPRB_eta_2nd_hot_2[j-1],
    PPPRB_eta_2nd_hot_3[j-1],
    PPPRB_eta_2nd_hot_4[j-1],
    PPPRB_eta_2nd_hot_5[j-1],
    PPPRB_eta_2nd_hot_6[j-1],
    PPPRB_eta_2nd_hot_7[j-1],
    PPPRB_eta_2nd_hot_8[j-1],
    PPPRB_eta_2nd_hot_9[j-1]
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
        PPPRB_beta_1st = PPPRB_beta_1st_list,
        PPPRB_eta_1st  = PPPRB_eta_1st_list,
        PPPRB_mu_1st   = PPPRB_mu_1st_list,
        #
        current_beta  = current_beta[i],
        current_mu    = current_mu[i],
        current_eta   = current_eta[i]
      )
    }
  )

  names(results_ppprb) <- c("cold", "hot_1", "hot_2", "hot_3","hot_4",
                            "hot_5","hot_6","hot_7","hot_8","hot_9")
  #
  mu_final    <- sapply(results_ppprb, `[[`, "mu_final")
  beta_final  <- sapply(results_ppprb, `[[`, "beta_final")
  eta_final   <- sapply(results_ppprb, `[[`, "eta_final")
  accept      <- sapply(results_ppprb, `[[`, "accept_or_not")
  idx_final     <- sapply(results_ppprb, `[[`, "idx_result")
  #
  PPPRB_mu_2nd_cold[j]    <- mu_final["cold"]
  PPPRB_mu_2nd_hot_1[j]   <- mu_final["hot_1"]
  PPPRB_mu_2nd_hot_2[j]   <- mu_final["hot_2"]
  PPPRB_mu_2nd_hot_3[j]   <- mu_final["hot_3"]
  PPPRB_mu_2nd_hot_4[j]   <- mu_final["hot_4"]
  PPPRB_mu_2nd_hot_5[j]   <- mu_final["hot_5"]
  PPPRB_mu_2nd_hot_6[j]   <- mu_final["hot_6"]
  PPPRB_mu_2nd_hot_7[j]   <- mu_final["hot_7"]
  PPPRB_mu_2nd_hot_8[j]   <- mu_final["hot_8"]
  PPPRB_mu_2nd_hot_9[j]   <- mu_final["hot_9"]
  
  PPPRB_beta_2nd_cold[j]  <- beta_final["cold"]
  PPPRB_beta_2nd_hot_1[j] <- beta_final["hot_1"]
  PPPRB_beta_2nd_hot_2[j] <- beta_final["hot_2"]
  PPPRB_beta_2nd_hot_3[j] <- beta_final["hot_3"]
  PPPRB_beta_2nd_hot_4[j] <- beta_final["hot_4"]
  PPPRB_beta_2nd_hot_5[j] <- beta_final["hot_5"]
  PPPRB_beta_2nd_hot_6[j] <- beta_final["hot_6"]
  PPPRB_beta_2nd_hot_7[j] <- beta_final["hot_7"]
  PPPRB_beta_2nd_hot_8[j] <- beta_final["hot_8"]
  PPPRB_beta_2nd_hot_9[j] <- beta_final["hot_9"]
  
  PPPRB_eta_2nd_cold[j]   <- eta_final["cold"]
  PPPRB_eta_2nd_hot_1[j]  <- eta_final["hot_1"]
  PPPRB_eta_2nd_hot_2[j]  <- eta_final["hot_2"]
  PPPRB_eta_2nd_hot_3[j]  <- eta_final["hot_3"]
  PPPRB_eta_2nd_hot_4[j]  <- eta_final["hot_4"]
  PPPRB_eta_2nd_hot_5[j]  <- eta_final["hot_5"]
  PPPRB_eta_2nd_hot_6[j]  <- eta_final["hot_6"]
  PPPRB_eta_2nd_hot_7[j]  <- eta_final["hot_7"]
  PPPRB_eta_2nd_hot_8[j]  <- eta_final["hot_8"]
  PPPRB_eta_2nd_hot_9[j]  <- eta_final["hot_9"]
  
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
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_1[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_1[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_1[j]
  } else if(idx==3){
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_2[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_2[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_2[j]
  } else if(idx==4){
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_3[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_3[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_3[j]
  } else if(idx==5){
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_4[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_4[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_4[j]
  } else if(idx==6){
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_5[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_5[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_5[j]
  } else if(idx==7){
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_6[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_6[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_6[j]
  } else if(idx==8){
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_7[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_7[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_7[j]
  } else if(idx==9){
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_8[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_8[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_8[j]
  } else if(idx==10){
    PPPRB_beta_hot<- PPPRB_beta_2nd_hot_9[j]
    PPPRB_eta_hot<- PPPRB_eta_2nd_hot_9[j]
    PPPRB_mu_hot<- PPPRB_mu_2nd_hot_9[j]
  }
  swap_result<-PPPRB_swap(temperature = temperature_list[[idx]],
                          #
                          all_param_likeli_cond=all_param_likeli_cond,
                          all_param_likeli_1st=all_param_likeli_1st,
                          #
                          cold_idx=current_idx[1],
                          hot_idx=current_idx[idx],
                          #
                          PPPRB_beta_cold=PPPRB_beta_2nd_cold[j],
                          PPPRB_beta_hot=PPPRB_beta_hot,
                          PPPRB_eta_cold=PPPRB_eta_2nd_cold[j],
                          PPPRB_eta_hot=PPPRB_eta_hot,
                          PPPRB_mu_cold=PPPRB_mu_2nd_cold[j],
                          PPPRB_mu_hot=PPPRB_mu_hot)
  
  # current_idx

  if(swap_result$accept_or_not==1){
    PPPRB_mu_2nd_cold[j] <-swap_result$mu_cold_chain
    PPPRB_beta_2nd_cold[j]<-swap_result$beta_cold_chain
    PPPRB_eta_2nd_cold[j]<-swap_result$eta_cold_chain
    accept_ratio_swap[j]<-1
    current_idx_cold<-current_idx[1]
    current_idx_hot<-current_idx[idx]
    current_idx[1]<-current_idx_hot
    #
    if(idx==2){
      PPPRB_mu_2nd_hot_1[j]<-swap_result$mu_hot_chain
      PPPRB_beta_2nd_hot_1[j]<-swap_result$beta_hot_chain
      PPPRB_eta_2nd_hot_1[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==3){
      PPPRB_mu_2nd_hot_2[j]<-swap_result$mu_hot_chain
      PPPRB_beta_2nd_hot_2[j]<-swap_result$beta_hot_chain
      PPPRB_eta_2nd_hot_2[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==4){
      PPPRB_mu_2nd_hot_3[j]<-swap_result$mu_hot_chain
      PPPRB_beta_2nd_hot_3[j]<-swap_result$beta_hot_chain
      PPPRB_eta_2nd_hot_3[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==5){
      PPPRB_mu_2nd_hot_4[j]<-swap_result$mu_hot_chain
      PPPRB_beta_2nd_hot_4[j]<-swap_result$beta_hot_chain
      PPPRB_eta_2nd_hot_4[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
     }else if(idx==6){
      PPPRB_mu_2nd_hot_5[j]<-swap_result$mu_hot_chain
      PPPRB_beta_2nd_hot_5[j]<-swap_result$beta_hot_chain
      PPPRB_eta_2nd_hot_5[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
     } else if(idx==7){
       PPPRB_mu_2nd_hot_6[j]<-swap_result$mu_hot_chain
       PPPRB_beta_2nd_hot_6[j]<-swap_result$beta_hot_chain
       PPPRB_eta_2nd_hot_6[j]<-swap_result$eta_hot_chain
       current_idx[idx]<-current_idx_cold
     } else if(idx==8){
       PPPRB_mu_2nd_hot_7[j]<-swap_result$mu_hot_chain
       PPPRB_beta_2nd_hot_7[j]<-swap_result$beta_hot_chain
       PPPRB_eta_2nd_hot_7[j]<-swap_result$eta_hot_chain
       current_idx[idx]<-current_idx_cold
     } else if(idx==9){
       PPPRB_mu_2nd_hot_8[j]<-swap_result$mu_hot_chain
       PPPRB_beta_2nd_hot_8[j]<-swap_result$beta_hot_chain
       PPPRB_eta_2nd_hot_8[j]<-swap_result$eta_hot_chain
       current_idx[idx]<-current_idx_cold
     } else if(idx==10){
       PPPRB_mu_2nd_hot_9[j]<-swap_result$mu_hot_chain
       PPPRB_beta_2nd_hot_9[j]<-swap_result$beta_hot_chain
       PPPRB_eta_2nd_hot_9[j]<-swap_result$eta_hot_chain
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



PPPRB_mu_2nd_cold<-PPPRB_mu_2nd_cold[-c(1:burnin)]
PPPRB_eta_2nd_cold<-PPPRB_eta_2nd_cold[-c(1:burnin)]
PPPRB_beta_2nd_cold<-PPPRB_beta_2nd_cold[-c(1:burnin)]


## hot chains burnin remove
PPPRB_beta_2nd_hot_1<-PPPRB_beta_2nd_hot_1[-c(1:burnin)]
PPPRB_beta_2nd_hot_2<-PPPRB_beta_2nd_hot_2[-c(1:burnin)]
PPPRB_beta_2nd_hot_3<-PPPRB_beta_2nd_hot_3[-c(1:burnin)]
PPPRB_beta_2nd_hot_4<-PPPRB_beta_2nd_hot_4[-c(1:burnin)]
PPPRB_beta_2nd_hot_5<-PPPRB_beta_2nd_hot_5[-c(1:burnin)]
PPPRB_beta_2nd_hot_6<-PPPRB_beta_2nd_hot_6[-c(1:burnin)]
PPPRB_beta_2nd_hot_7<-PPPRB_beta_2nd_hot_7[-c(1:burnin)]
PPPRB_beta_2nd_hot_8<-PPPRB_beta_2nd_hot_8[-c(1:burnin)]
PPPRB_beta_2nd_hot_9<-PPPRB_beta_2nd_hot_9[-c(1:burnin)]

PPPRB_eta_2nd_hot_1<-PPPRB_eta_2nd_hot_1[-c(1:burnin)]
PPPRB_eta_2nd_hot_2<-PPPRB_eta_2nd_hot_2[-c(1:burnin)]
PPPRB_eta_2nd_hot_3<-PPPRB_eta_2nd_hot_3[-c(1:burnin)]
PPPRB_eta_2nd_hot_4<-PPPRB_eta_2nd_hot_4[-c(1:burnin)]
PPPRB_eta_2nd_hot_5<-PPPRB_eta_2nd_hot_5[-c(1:burnin)]
PPPRB_eta_2nd_hot_6<-PPPRB_eta_2nd_hot_6[-c(1:burnin)]
PPPRB_eta_2nd_hot_7<-PPPRB_eta_2nd_hot_7[-c(1:burnin)]
PPPRB_eta_2nd_hot_8<-PPPRB_eta_2nd_hot_8[-c(1:burnin)]
PPPRB_eta_2nd_hot_9<-PPPRB_eta_2nd_hot_9[-c(1:burnin)]


PPPRB_mu_2nd_hot_1<-PPPRB_mu_2nd_hot_1[-c(1:burnin)]
PPPRB_mu_2nd_hot_2<-PPPRB_mu_2nd_hot_2[-c(1:burnin)]
PPPRB_mu_2nd_hot_3<-PPPRB_mu_2nd_hot_3[-c(1:burnin)]
PPPRB_mu_2nd_hot_4<-PPPRB_mu_2nd_hot_4[-c(1:burnin)]
PPPRB_mu_2nd_hot_5<-PPPRB_mu_2nd_hot_5[-c(1:burnin)]
PPPRB_mu_2nd_hot_6<-PPPRB_mu_2nd_hot_6[-c(1:burnin)]
PPPRB_mu_2nd_hot_7<-PPPRB_mu_2nd_hot_7[-c(1:burnin)]
PPPRB_mu_2nd_hot_8<-PPPRB_mu_2nd_hot_8[-c(1:burnin)]
PPPRB_mu_2nd_hot_9<-PPPRB_mu_2nd_hot_9[-c(1:burnin)]


# PPPRB :third stage ------------------------------------------------------


all_mu   <- c(
  PPPRB_mu_2nd_cold,
  PPPRB_mu_2nd_hot_1,
  PPPRB_mu_2nd_hot_2,
  PPPRB_mu_2nd_hot_3,
  PPPRB_mu_2nd_hot_4,
  PPPRB_mu_2nd_hot_5,
  PPPRB_mu_2nd_hot_6,
  PPPRB_mu_2nd_hot_7,
  PPPRB_mu_2nd_hot_8,
  PPPRB_mu_2nd_hot_9
)


all_eta  <- c(
  PPPRB_eta_2nd_cold,
  PPPRB_eta_2nd_hot_1,
  PPPRB_eta_2nd_hot_2,
  PPPRB_eta_2nd_hot_3,
  PPPRB_eta_2nd_hot_4,
  PPPRB_eta_2nd_hot_5,
  PPPRB_eta_2nd_hot_6,
  PPPRB_eta_2nd_hot_7,
  PPPRB_eta_2nd_hot_8,
  PPPRB_eta_2nd_hot_9
)

all_beta <- c(
  PPPRB_beta_2nd_cold,
  PPPRB_beta_2nd_hot_1,
  PPPRB_beta_2nd_hot_2,
  PPPRB_beta_2nd_hot_3,
  PPPRB_beta_2nd_hot_4,
  PPPRB_beta_2nd_hot_5,
  PPPRB_beta_2nd_hot_6,
  PPPRB_beta_2nd_hot_7,
  PPPRB_beta_2nd_hot_8,
  PPPRB_beta_2nd_hot_9
)

tic()
hawkes_loglik_conditional_all <- mclapply(
  seq_along(all_mu),
  function(i) {
    hawkes_loglik_conditional_fast(
      data$t,
      all_mu[i],
      all_eta[i] * all_beta[i],
      all_beta[i],
      T2,
      capital_T,
      capital_T
    )
  },
  mc.cores = detectCores()-1
)
toc()

hawkes_loglik_conditional_all <- unlist(hawkes_loglik_conditional_all)

all_param_likeli_cond_2nd<-cbind(all_mu,
                                 all_eta,
                                 all_beta,
                                 hawkes_loglik_conditional_all)

colnames(all_param_likeli_cond_2nd)<-c("mu",
                                       "eta",
                                       "beta",
                                       "log_likeli_cond")

rm(hawkes_loglik_conditional_all)

all_param_likeli_2nd<-cbind(all_mu,
                            all_eta,
                            all_beta)

colnames(all_param_likeli_2nd)<-c("mu","eta","beta")
all_param_likeli_2nd<-cbind(all_param_likeli_2nd,NA)
colnames(all_param_likeli_2nd)[4]<-"log_likeli"

all_param_likeli_1_2_uni<-unique(all_param_likeli_1_2)
for(i in 1:nrow(all_param_likeli_1_2_uni)){
  idx<-which(all_param_likeli_2nd[,1]==all_param_likeli_1_2_uni[i,1]&
               all_param_likeli_2nd[,2]==all_param_likeli_1_2_uni[i,2] &
               all_param_likeli_2nd[,3]==all_param_likeli_1_2_uni[i,3])
  all_param_likeli_2nd[idx,4]<-all_param_likeli_1_2_uni[i,4]
}

sum( all_param_likeli_2nd[,c(1:3)]==all_param_likeli_cond_2nd[,c(1:3)] )

### PPPRB third stage 
### PPPRB third stage 
### PPPRB third stage 
### PPPRB third stage 
all_param_likeli_cond_2nd<-cbind( all_param_likeli_cond_2nd,rep(1:10,each=iterations-burnin) )
colnames(all_param_likeli_cond_2nd)[5]<-"chain"
all_param_likeli_cond_2nd<-cbind(1:nrow(all_param_likeli_cond_2nd),all_param_likeli_cond_2nd)
colnames(all_param_likeli_cond_2nd)[1]<-"idxidx"

PPPRB_mu_3rd_cold<-numeric(iterations)
PPPRB_mu_3rd_hot_1<-numeric(iterations)
PPPRB_mu_3rd_hot_2<-numeric(iterations)
PPPRB_mu_3rd_hot_3<-numeric(iterations)
PPPRB_mu_3rd_hot_4<-numeric(iterations)
PPPRB_mu_3rd_hot_5<-numeric(iterations)
PPPRB_mu_3rd_hot_6<-numeric(iterations)
PPPRB_mu_3rd_hot_7<-numeric(iterations)
PPPRB_mu_3rd_hot_8<-numeric(iterations)
PPPRB_mu_3rd_hot_9<-numeric(iterations)

#
PPPRB_beta_3rd_cold<-numeric(iterations)
PPPRB_beta_3rd_hot_1<-numeric(iterations)
PPPRB_beta_3rd_hot_2<-numeric(iterations)
PPPRB_beta_3rd_hot_3<-numeric(iterations)
PPPRB_beta_3rd_hot_4<-numeric(iterations)
PPPRB_beta_3rd_hot_5<-numeric(iterations)
PPPRB_beta_3rd_hot_6<-numeric(iterations)
PPPRB_beta_3rd_hot_7<-numeric(iterations)
PPPRB_beta_3rd_hot_8<-numeric(iterations)
PPPRB_beta_3rd_hot_9<-numeric(iterations)
#
PPPRB_eta_3rd_cold<-numeric(iterations)
PPPRB_eta_3rd_hot_1<-numeric(iterations)
PPPRB_eta_3rd_hot_2<-numeric(iterations)
PPPRB_eta_3rd_hot_3<-numeric(iterations)
PPPRB_eta_3rd_hot_4<-numeric(iterations)
PPPRB_eta_3rd_hot_5<-numeric(iterations)
PPPRB_eta_3rd_hot_6<-numeric(iterations)
PPPRB_eta_3rd_hot_7<-numeric(iterations)
PPPRB_eta_3rd_hot_8<-numeric(iterations)
PPPRB_eta_3rd_hot_9<-numeric(iterations)

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


## result from the 1st stage
PPPRB_mu_2nd_list<-c(  PPPRB_mu_2nd_cold,
                       PPPRB_mu_2nd_hot_1,
                       PPPRB_mu_2nd_hot_2,
                       PPPRB_mu_2nd_hot_3,
                       PPPRB_mu_2nd_hot_4,
                       PPPRB_mu_2nd_hot_5,
                       PPPRB_mu_2nd_hot_6,
                       PPPRB_mu_2nd_hot_7,
                       PPPRB_mu_2nd_hot_8,
                       PPPRB_mu_2nd_hot_9 )

PPPRB_beta_2nd_list<-c(PPPRB_beta_2nd_cold,
                       PPPRB_beta_2nd_hot_1,
                       PPPRB_beta_2nd_hot_2,
                       PPPRB_beta_2nd_hot_3,
                       PPPRB_beta_2nd_hot_4,
                       PPPRB_beta_2nd_hot_5,
                       PPPRB_beta_2nd_hot_6,
                       PPPRB_beta_2nd_hot_7,
                       PPPRB_beta_2nd_hot_8,
                       PPPRB_beta_2nd_hot_9
)



PPPRB_eta_2nd_list<-c(PPPRB_eta_2nd_cold,
                      PPPRB_eta_2nd_hot_1,
                      PPPRB_eta_2nd_hot_2,
                      PPPRB_eta_2nd_hot_3,
                      PPPRB_eta_2nd_hot_4,
                      PPPRB_eta_2nd_hot_5,
                      PPPRB_eta_2nd_hot_6,
                      PPPRB_eta_2nd_hot_7,
                      PPPRB_eta_2nd_hot_8,
                      PPPRB_eta_2nd_hot_9
)





## initial

current_idx<-c(100,
               (iterations-burnin)+100,
               2*(iterations-burnin)+100,
               3*(iterations-burnin)+100,
               4*(iterations-burnin)+100,
               5*(iterations-burnin)+100,
               6*(iterations-burnin)+100,
               7*(iterations-burnin)+100,
               8*(iterations-burnin)+100,
               9*(iterations-burnin)+100)

## initial 
PPPRB_mu_3rd_cold[1]<-PPPRB_mu_2nd_list[current_idx[1]]
PPPRB_mu_3rd_hot_1[1]<-PPPRB_mu_2nd_list[current_idx[2]]
PPPRB_mu_3rd_hot_2[1]<-PPPRB_mu_2nd_list[current_idx[3]]
PPPRB_mu_3rd_hot_3[1]<-PPPRB_mu_2nd_list[current_idx[4]]
PPPRB_mu_3rd_hot_4[1]<-PPPRB_mu_2nd_list[current_idx[5]]
PPPRB_mu_3rd_hot_5[1]<-PPPRB_mu_2nd_list[current_idx[6]]
PPPRB_mu_3rd_hot_6[1]<-PPPRB_mu_2nd_list[current_idx[7]]
PPPRB_mu_3rd_hot_7[1]<-PPPRB_mu_2nd_list[current_idx[8]]
PPPRB_mu_3rd_hot_8[1]<-PPPRB_mu_2nd_list[current_idx[9]]
PPPRB_mu_3rd_hot_9[1]<-PPPRB_mu_2nd_list[current_idx[10]]

## initial 
PPPRB_beta_3rd_cold[1]<-PPPRB_beta_2nd_list[current_idx[1]]
PPPRB_beta_3rd_hot_1[1]<-PPPRB_beta_2nd_list[current_idx[2]]
PPPRB_beta_3rd_hot_2[1]<-PPPRB_beta_2nd_list[current_idx[3]]
PPPRB_beta_3rd_hot_3[1]<-PPPRB_beta_2nd_list[current_idx[4]]
PPPRB_beta_3rd_hot_4[1]<-PPPRB_beta_2nd_list[current_idx[5]]
PPPRB_beta_3rd_hot_5[1]<-PPPRB_beta_2nd_list[current_idx[6]]
PPPRB_beta_3rd_hot_6[1]<-PPPRB_beta_2nd_list[current_idx[7]]
PPPRB_beta_3rd_hot_7[1]<-PPPRB_beta_2nd_list[current_idx[8]]
PPPRB_beta_3rd_hot_8[1]<-PPPRB_beta_2nd_list[current_idx[9]]
PPPRB_beta_3rd_hot_9[1]<-PPPRB_beta_2nd_list[current_idx[10]]

## initial 
PPPRB_eta_3rd_cold[1]<-PPPRB_eta_2nd_list[current_idx[[1]]]
PPPRB_eta_3rd_hot_1[1]<-PPPRB_eta_2nd_list[current_idx[[2]]]
PPPRB_eta_3rd_hot_2[1]<-PPPRB_eta_2nd_list[current_idx[[3]]]
PPPRB_eta_3rd_hot_3[1]<-PPPRB_eta_2nd_list[current_idx[[4]]]
PPPRB_eta_3rd_hot_4[1]<-PPPRB_eta_2nd_list[current_idx[[5]]]
PPPRB_eta_3rd_hot_5[1]<-PPPRB_eta_2nd_list[current_idx[[6]]]
PPPRB_eta_3rd_hot_6[1]<-PPPRB_eta_2nd_list[current_idx[[7]]]
PPPRB_eta_3rd_hot_7[1]<-PPPRB_eta_2nd_list[current_idx[[8]]]
PPPRB_eta_3rd_hot_8[1]<-PPPRB_eta_2nd_list[current_idx[[9]]]
PPPRB_eta_3rd_hot_9[1]<-PPPRB_eta_2nd_list[current_idx[[10]]]


# acceptance ratio for swapping
accept_ratio_swap_3rd<-numeric(iterations)

chain_rows <- split(
  seq_len(nrow(all_param_likeli_cond_2nd)),
  all_param_likeli_cond_2nd[, "chain"]   # 
)

tic()

for(j in 2: iterations){

  current_mu <- c(
    PPPRB_mu_3rd_cold[j-1],
    PPPRB_mu_3rd_hot_1[j-1],
    PPPRB_mu_3rd_hot_2[j-1],
    PPPRB_mu_3rd_hot_3[j-1],
    PPPRB_mu_3rd_hot_4[j-1],
    PPPRB_mu_3rd_hot_5[j-1],
    PPPRB_mu_3rd_hot_6[j-1],
    PPPRB_mu_3rd_hot_7[j-1],
    PPPRB_mu_3rd_hot_8[j-1],
    PPPRB_mu_3rd_hot_9[j-1]
  )
  
  current_beta <- c(
    PPPRB_beta_3rd_cold[j-1],
    PPPRB_beta_3rd_hot_1[j-1],
    PPPRB_beta_3rd_hot_2[j-1],
    PPPRB_beta_3rd_hot_3[j-1],
    PPPRB_beta_3rd_hot_4[j-1],
    PPPRB_beta_3rd_hot_5[j-1],
    PPPRB_beta_3rd_hot_6[j-1],
    PPPRB_beta_3rd_hot_7[j-1],
    PPPRB_beta_3rd_hot_8[j-1],
    PPPRB_beta_3rd_hot_9[j-1]
  )
  
  current_eta <- c(
    PPPRB_eta_3rd_cold[j-1],
    PPPRB_eta_3rd_hot_1[j-1],
    PPPRB_eta_3rd_hot_2[j-1],
    PPPRB_eta_3rd_hot_3[j-1],
    PPPRB_eta_3rd_hot_4[j-1],
    PPPRB_eta_3rd_hot_5[j-1],
    PPPRB_eta_3rd_hot_6[j-1],
    PPPRB_eta_3rd_hot_7[j-1],
    PPPRB_eta_3rd_hot_8[j-1],
    PPPRB_eta_3rd_hot_9[j-1]
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
        PPPRB_beta_1st = PPPRB_beta_2nd_list,
        PPPRB_eta_1st  = PPPRB_eta_2nd_list,
        PPPRB_mu_1st   = PPPRB_mu_2nd_list,
        #
        current_beta  = current_beta[i],
        current_mu    = current_mu[i],
        current_eta   = current_eta[i]
      )
    }
  )
  
  names(results_ppprb) <- c("cold", "hot_1", "hot_2", "hot_3","hot_4",
                            "hot_5","hot_6","hot_7","hot_8","hot_9")
  #
  mu_final    <- sapply(results_ppprb, `[[`, "mu_final")
  beta_final  <- sapply(results_ppprb, `[[`, "beta_final")
  eta_final   <- sapply(results_ppprb, `[[`, "eta_final")
  accept      <- sapply(results_ppprb, `[[`, "accept_or_not")
  idx_final     <- sapply(results_ppprb, `[[`, "idx_result")
  #
  PPPRB_mu_3rd_cold[j]    <- mu_final["cold"]
  PPPRB_mu_3rd_hot_1[j]   <- mu_final["hot_1"]
  PPPRB_mu_3rd_hot_2[j]   <- mu_final["hot_2"]
  PPPRB_mu_3rd_hot_3[j]   <- mu_final["hot_3"]
  PPPRB_mu_3rd_hot_4[j]   <- mu_final["hot_4"]
  PPPRB_mu_3rd_hot_5[j]   <- mu_final["hot_5"]
  PPPRB_mu_3rd_hot_6[j]   <- mu_final["hot_6"]
  PPPRB_mu_3rd_hot_7[j]   <- mu_final["hot_7"]
  PPPRB_mu_3rd_hot_8[j]   <- mu_final["hot_8"]
  PPPRB_mu_3rd_hot_9[j]   <- mu_final["hot_9"]
  
  PPPRB_beta_3rd_cold[j]  <- beta_final["cold"]
  PPPRB_beta_3rd_hot_1[j] <- beta_final["hot_1"]
  PPPRB_beta_3rd_hot_2[j] <- beta_final["hot_2"]
  PPPRB_beta_3rd_hot_3[j] <- beta_final["hot_3"]
  PPPRB_beta_3rd_hot_4[j] <- beta_final["hot_4"]
  PPPRB_beta_3rd_hot_5[j] <- beta_final["hot_5"]
  PPPRB_beta_3rd_hot_6[j] <- beta_final["hot_6"]
  PPPRB_beta_3rd_hot_7[j] <- beta_final["hot_7"]
  PPPRB_beta_3rd_hot_8[j] <- beta_final["hot_8"]
  PPPRB_beta_3rd_hot_9[j] <- beta_final["hot_9"]
  
  PPPRB_eta_3rd_cold[j]   <- eta_final["cold"]
  PPPRB_eta_3rd_hot_1[j]  <- eta_final["hot_1"]
  PPPRB_eta_3rd_hot_2[j]  <- eta_final["hot_2"]
  PPPRB_eta_3rd_hot_3[j]  <- eta_final["hot_3"]
  PPPRB_eta_3rd_hot_4[j]  <- eta_final["hot_4"]
  PPPRB_eta_3rd_hot_5[j]  <- eta_final["hot_5"]
  PPPRB_eta_3rd_hot_6[j]  <- eta_final["hot_6"]
  PPPRB_eta_3rd_hot_7[j]  <- eta_final["hot_7"]
  PPPRB_eta_3rd_hot_8[j]  <- eta_final["hot_8"]
  PPPRB_eta_3rd_hot_9[j]  <- eta_final["hot_9"]
  
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
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_1[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_1[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_1[j]
  } else if(idx==3){
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_2[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_2[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_2[j]
  } else if(idx==4){
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_3[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_3[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_3[j]
  } else if(idx==5){
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_4[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_4[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_4[j]
  } else if(idx==6){
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_5[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_5[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_5[j]
  } else if(idx==7){
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_6[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_6[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_6[j]
  } else if(idx==8){
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_7[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_7[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_7[j]
  } else if(idx==9){
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_8[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_8[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_8[j]
  } else if(idx==10){
    PPPRB_beta_hot<- PPPRB_beta_3rd_hot_9[j]
    PPPRB_eta_hot<- PPPRB_eta_3rd_hot_9[j]
    PPPRB_mu_hot<- PPPRB_mu_3rd_hot_9[j]
  }
  swap_result<-PPPRB_swap(temperature = temperature_list[[idx]],
                          #
                          all_param_likeli_cond=all_param_likeli_cond_2nd,
                          all_param_likeli_1st=all_param_likeli_2nd,
                          #
                          cold_idx=current_idx[1],
                          hot_idx=current_idx[idx],
                          #
                          PPPRB_beta_cold=PPPRB_beta_3rd_cold[j],
                          PPPRB_beta_hot=PPPRB_beta_hot,
                          PPPRB_eta_cold=PPPRB_eta_3rd_cold[j],
                          PPPRB_eta_hot=PPPRB_eta_hot,
                          PPPRB_mu_cold=PPPRB_mu_3rd_cold[j],
                          PPPRB_mu_hot=PPPRB_mu_hot)
  
  # current_idx

  if(swap_result$accept_or_not==1){
    PPPRB_mu_3rd_cold[j] <-swap_result$mu_cold_chain
    PPPRB_beta_3rd_cold[j]<-swap_result$beta_cold_chain
    PPPRB_eta_3rd_cold[j]<-swap_result$eta_cold_chain
    accept_ratio_swap_3rd[j]<-1
    current_idx_cold<-current_idx[1]
    current_idx_hot<-current_idx[idx]
    current_idx[1]<-current_idx_hot
    #
    if(idx==2){
      PPPRB_mu_3rd_hot_1[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_1[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_1[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==3){
      PPPRB_mu_3rd_hot_2[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_2[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_2[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==4){
      PPPRB_mu_3rd_hot_3[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_3[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_3[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
      
    } else if(idx==5){
      PPPRB_mu_3rd_hot_4[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_4[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_4[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
    }else if(idx==6){
      PPPRB_mu_3rd_hot_5[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_5[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_5[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==7){
      PPPRB_mu_3rd_hot_6[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_6[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_6[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==8){
      PPPRB_mu_3rd_hot_7[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_7[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_7[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==9){
      PPPRB_mu_3rd_hot_8[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_8[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_8[j]<-swap_result$eta_hot_chain
      current_idx[idx]<-current_idx_cold
    } else if(idx==10){
      PPPRB_mu_3rd_hot_9[j]<-swap_result$mu_hot_chain
      PPPRB_beta_3rd_hot_9[j]<-swap_result$beta_hot_chain
      PPPRB_eta_3rd_hot_9[j]<-swap_result$eta_hot_chain
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


mu_post<-mu_post[-c(1:burnin)]
beta_post<-beta_post[-c(1:burnin)]
eta_post<-eta_post[-c(1:burnin)]

PPPRB_mu_3rd_cold<-PPPRB_mu_3rd_cold[-c(1:burnin)]
PPPRB_eta_3rd_cold<-PPPRB_eta_3rd_cold[-c(1:burnin)]
PPPRB_beta_3rd_cold<-PPPRB_beta_3rd_cold[-c(1:burnin)]

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


mu_df <- bind_rows(
  data.frame(stage = "PP-RB 1st",  parameter = "mu",  as.list(summary_stat(mu_post_1st))),
  data.frame(stage = "PP-RB 2nd",  parameter = "mu",  as.list(summary_stat(mu_post_2nd))),
  data.frame(stage = "PP-RB 3rd",  parameter = "mu",  as.list(summary_stat(mu_post_3rd))),
  data.frame(stage = "PPPRB 1st",  parameter = "mu",  as.list(summary_stat(results$cold$post_mu_tempered_1st))),
  data.frame(stage = "PPPRB 2nd",  parameter = "mu",  as.list(summary_stat(PPPRB_mu_2nd_cold))),
  data.frame(stage = "PPPRB 3rd",  parameter = "mu",  as.list(summary_stat(PPPRB_mu_3rd_cold))),
  data.frame(stage = "Full",       parameter = "mu",  as.list(summary_stat(mu_post)))
)

beta_df <- bind_rows(
  data.frame(stage = "PP-RB 1st",  parameter = "beta",  as.list(summary_stat(beta_post_1st))),
  data.frame(stage = "PP-RB 2nd",  parameter = "beta",  as.list(summary_stat(beta_post_2nd))),
  data.frame(stage = "PP-RB 3rd",  parameter = "beta",  as.list(summary_stat(beta_post_3rd))),
  data.frame(stage = "PPPRB 1st",  parameter = "beta",  as.list(summary_stat(results$cold$post_beta_tempered_1st))),
  data.frame(stage = "PPPRB 2nd",  parameter = "beta",  as.list(summary_stat(PPPRB_beta_2nd_cold))),
  data.frame(stage = "PPPRB 3rd",  parameter = "beta",  as.list(summary_stat(PPPRB_beta_3rd_cold))),
  data.frame(stage = "Full",       parameter = "beta",  as.list(summary_stat(beta_post)))
)

eta_df <- bind_rows(
  data.frame(stage = "PP-RB 1st",  parameter = "eta",  as.list(summary_stat(eta_post_1st))),
  data.frame(stage = "PP-RB 2nd",  parameter = "eta",  as.list(summary_stat(eta_post_2nd))),
  data.frame(stage = "PP-RB 3rd",  parameter = "eta",  as.list(summary_stat(eta_post_3rd))),
  data.frame(stage = "PPPRB 1st",  parameter = "eta",  as.list(summary_stat(results$cold$post_eta_tempered_1st))),
  data.frame(stage = "PPPRB 2nd",  parameter = "eta",  as.list(summary_stat(PPPRB_eta_2nd_cold))),
  data.frame(stage = "PPPRB 3rd",  parameter = "eta",  as.list(summary_stat(PPPRB_eta_3rd_cold))),
  data.frame(stage = "Full",       parameter = "eta",  as.list(summary_stat(eta_post)))
)
colnames(eta_df)[c(4:5)]<-c("lower","upper")
colnames(beta_df)[c(4:5)]<-c("lower","upper")
colnames(mu_df)[c(4:5)]<-c("lower","upper")


# 1. Beta Plot
p_beta <- ggplot(beta_df, aes(x = stage, y = mean, color = stage)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.8) +
  labs(title = expression(paste(beta)), x = NULL, y = "Value") +
  scale_color_manual(values = c("black", "red", "orange", "green", "blue", "purple", "cyan")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Eta Plot 
p_eta <- ggplot(eta_df, aes(x = stage, y = mean, color = stage)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.8) +
  labs(title = expression(paste(eta)), x = NULL, y = "Value") +
  scale_color_manual(values = c("black", "red", "orange", "green", "blue", "purple", "cyan")) + 
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Mu Plot 
p_mu <- ggplot(mu_df, aes(x = stage, y = mean, color = stage)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.8) +
  labs(title = expression(paste(mu)), x = NULL, y = "Value") +
  scale_color_manual(values = c("black", "red", "orange", "green", "blue", "purple", "cyan")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # 마지막 그림에만 범례나 축 설정을 유지할 수 있
  
fig1

fig1<-ggplot(data, aes(x = time, y = mag)) +
  geom_point(color = "black", shape = 16) +
  geom_vline(
    xintercept = as.POSIXct("1990-02-06 18:00"),
    linetype = "dashed",
    color = "#D55E00"   # vermillion
  ) +
  geom_vline(
    xintercept = as.POSIXct("1990-05-17 18:50"),
    linetype = "dashed",
    color = "#D55E00"   # vermillion
  ) +
  geom_vline(
    xintercept = as.POSIXct("1989-10-18 01:04"),
    linetype = "dashed",
    color = "#0072B2"
  ) +
  labs(
    x = "Date",
    y = "Magnitude",
    title = "Earthquakes > 2.5 Near 1989 Loma Prieta"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),  # 
    axis.title.x = element_text(size = 14),             # 
    axis.title.y = element_text(size = 14)              # 
  )
library(patchwork)
common_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  aspect.ratio = 1
)
fig1_2 <- fig1 + common_theme
p_beta_2 <- p_beta + common_theme
p_eta_2 <- p_eta + common_theme
p_mu_2 <- p_mu + common_theme

fig1_2+p_beta_2+p_eta_2+p_mu_2

p_beta_2  <- p_beta_2  + theme(legend.position = "none")+ scale_color_brewer(palette = "Dark2")
p_eta_2 <- p_eta_2 + theme(legend.position = "none")+ scale_color_brewer(palette = "Dark2")
p_mu_2 <- p_mu_2 + theme(legend.position = "none")+ scale_color_brewer(palette = "Dark2")

fig1 /
  (p_beta_2 + p_eta_2 + p_mu_2) +
  plot_layout(guides = "collect",heights = c(1, 1)) &theme(legend.position = "right")
# 1200 auto


library(posterior)
draws_arr <- array(NA, dim = c(length(mu_post), 2, 1))
draws_arr[, 1, 1] <- mu_post
draws_arr[, 2, 1] <- PPPRB_mu_3rd_cold
dimnames(draws_arr) <- list(NULL, NULL, c("mu"))
d <- as_draws_array(draws_arr)
summarise_draws(d)


draws_arr <- array(NA, dim = c(length(beta_post), 2, 1))
draws_arr[, 1, 1] <- beta_post
draws_arr[, 2, 1] <- PPPRB_beta_3rd_cold
dimnames(draws_arr) <- list(NULL, NULL, c("beta"))
d <- as_draws_array(draws_arr)
summarise_draws(d)



draws_arr <- array(NA, dim = c(length(eta_post), 2, 1))
draws_arr[, 1, 1] <- eta_post
draws_arr[, 2, 1] <- PPPRB_eta_3rd_cold
dimnames(draws_arr) <- list(NULL, NULL, c("eta"))
d <- as_draws_array(draws_arr)
summarise_draws(d)

save.image("hawks_process_ppprb.RData")

library(coda)
combined <- mcmc.list(as.mcmc(mu_post), as.mcmc(PPPRB_mu_3rd_cold))
gelman.diag(combined)

combined <- mcmc.list(as.mcmc(beta_post), as.mcmc(PPPRB_beta_3rd_cold))
gelman.diag(combined)

combined <- mcmc.list(as.mcmc(eta_post), as.mcmc(PPPRB_eta_3rd_cold))
gelman.diag(combined)


### PPRB and MH

combined <- mcmc.list(as.mcmc(mu_post), as.mcmc(mu_post_3rd))
gelman.diag(combined)

combined <- mcmc.list(as.mcmc(beta_post), as.mcmc(beta_post_3rd))
gelman.diag(combined)

combined <- mcmc.list(as.mcmc(eta_post), as.mcmc(eta_post_3rd))
gelman.diag(combined)


draws_arr <- array(NA, dim = c(length(mu_post), 2, 1))
draws_arr[, 1, 1] <- mu_post
draws_arr[, 2, 1] <- mu_post_3rd
dimnames(draws_arr) <- list(NULL, NULL, c("mu"))
d <- as_draws_array(draws_arr)
summarise_draws(d)


draws_arr <- array(NA, dim = c(length(beta_post), 2, 1))
draws_arr[, 1, 1] <- beta_post
draws_arr[, 2, 1] <- beta_post_3rd
dimnames(draws_arr) <- list(NULL, NULL, c("beta"))
d <- as_draws_array(draws_arr)
summarise_draws(d)



draws_arr <- array(NA, dim = c(length(eta_post), 2, 1))
draws_arr[, 1, 1] <- eta_post
draws_arr[, 2, 1] <- eta_post_3rd
dimnames(draws_arr) <- list(NULL, NULL, c("eta"))
d <- as_draws_array(draws_arr)
summarise_draws(d)

# ESS and so on ----------------------------------------------------------
library(coda)
single_mu <- mcmc(mu_post)
single_beta <- mcmc(beta_post)
single_eta <- mcmc(eta_post)


second_mu_pprb <- mcmc(mu_post_3rd)
second_beta_pprb <- mcmc(beta_post_3rd)
second_eta_pprb <- mcmc(eta_post_3rd)

second_mu_cold_ppprb <- mcmc(PPPRB_mu_3rd_cold)
second_eta_cold_ppprb <- mcmc(PPPRB_eta_3rd_cold)
second_beta_cold_ppprb <- mcmc(PPPRB_beta_3rd_cold)

plot(single_mu)
plot(second_mu_pprb)
plot(second_mu_cold_ppprb)

mcmc_list <- mcmc.list(single_mu,
                       single_beta,
                       single_eta,
                       second_mu_pprb,
                       second_beta_pprb,
                       second_eta_pprb,
                       second_mu_cold_ppprb,
                       second_beta_cold_ppprb,
                       second_eta_cold_ppprb
)

length( unique(second_beta_cold_ppprb) )
length( unique(single_beta) )


ess_theta_each <- sapply(mcmc_list, effectiveSize)

names(ess_theta_each)<-c("single_mu",
                         "single_beta",
                         "single_eta",
                         "mu_pprb_3rd",
                         "beta_pprb_3rd",
                         "eta_pprb_3rd",
                         "mu_ppprb_cold_3rd",
                         "beta_ppprb_cold_3rd",
                         "eta_ppprb_cold_3rd")

round(ess_theta_each,1)

round(ess_theta_each[1:3]/806.2 ,1)

round(ess_theta_each[4:6]/344.3,1)

round(ess_theta_each[7:9]/397.5,1)

