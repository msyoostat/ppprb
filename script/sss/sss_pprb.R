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

load("sss_mh.RData")

#RNGkind("L'Ecuyer-CMRG")
set.seed(100) #1500000 

# PPRB --------------------------------------------------------------------
### partition the data
jump_sigma_s2 <- 0.07
jump_sigma_n2 <- 0.05
jump_phi <- 0.15


parts <- split(1:N, sample(rep(1:3, length.out = N)))

idx_1<-parts$`1`
idx_2<-parts$`2`
idx_3<-parts$`3`

data_1<-data_scaled[idx_1]
data_2<-data_scaled[idx_2]
data_3<-data_scaled[idx_3]

locations_1<-locations[idx_1,]
locations_2<-locations[idx_2,]
locations_3<-locations[idx_3,]

dist_mat_1 <- distm(locations_1, fun = distCosine) / 1000

temp<-data.frame(lon=(rbind(locations_1,locations_2,locations_3) )[,1],
                 lat=(rbind(locations_1,locations_2,locations_3) )[,2],
                 sss=c(data_1,rep(NA,(length(idx_2)) ),rep(NA,(length(idx_3)) ) ))

ggplot()+
  geom_tile(data=temp,aes(x=lon,y=lat,fill=sss))+
  scale_fill_viridis()+
  ggtitle("partition_1")

temp<-data.frame(lon=(rbind(locations_1,locations_2,locations_3) )[,1],
                 lat=(rbind(locations_1,locations_2,locations_3) )[,2],
                 sss=c(rep(NA,length(idx_1) ),data_2,rep(NA,length(idx_3) ) ))

ggplot()+
  geom_tile(data=temp,aes(x=lon,y=lat,fill=sss))+
  scale_fill_viridis()+
  ggtitle("partition_2")

rm(temp,parts,idx_1,idx_2)
rm(idx_3)

log_sigma_s2_post_pprb_1 <- numeric(iterations)
log_sigma_n2_post_pprb_1 <- numeric(iterations)
log_phi_post_pprb_1 <- numeric(iterations)

# Initial values
log_sigma_s2_post_pprb_1[1] <- log(0.6)
log_sigma_n2_post_pprb_1[1] <- log(0.2)
log_phi_post_pprb_1[1] <- log(median(dist_mat)/3)

accept_pprb_1 <- numeric(iterations)





# ---------------- MCMC -----------------
tic()
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
  log_post_curr <- loglik_fun(data_1, sigma_s2_curr, sigma_n2_curr, phi_curr,
                              dist_mat_1,length(data_1)) +
    dnorm(log_sigma_s2_post_pprb_1[i-1], mu_log_sigma_s, sd_log_sigma_s, log=TRUE) +
    dnorm(log_sigma_n2_post_pprb_1[i-1], mu_log_sigma_n, sd_log_sigma_n, log=TRUE) +
    dnorm(log_phi_post_pprb_1[i-1], mu_log_phi, sd_log_phi, log=TRUE)
  
  log_post_prop <- loglik_fun(data_1, sigma_s2_prop, sigma_n2_prop, phi_prop,
                              dist_mat_1,length(data_1)) +
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
  } else {
    log_sigma_s2_post_pprb_1[i] <- log_sigma_s2_post_pprb_1[i-1]
    log_sigma_n2_post_pprb_1[i] <- log_sigma_n2_post_pprb_1[i-1]
    log_phi_post_pprb_1[i]     <- log_phi_post_pprb_1[i-1]
  }
  
  if(i %% 500 == 0){
    cat("iter:", i, "accept rate:", mean(accept_pprb_1[(i-499):i]), "\n")
  }
}
toc()

accept_pprb_1 <- accept_pprb_1[-c(1:burnin)]
mean(accept_pprb_1)

log_sigma_s2_post_pprb_1<-log_sigma_s2_post_pprb_1[-c(1:burnin)]
log_sigma_n2_post_pprb_1<-log_sigma_n2_post_pprb_1[-c(1:burnin)]
log_phi_post_pprb_1<-log_phi_post_pprb_1[-c(1:burnin)]



rm(i,log_post_curr,log_post_prop,log_ratio)
rm(phi_curr,phi_prop,proposal_log_phi,proposal_log_sigma_n,proposal_log_sigma_s)
rm(sigma_n2_curr,sigma_n2_prop,sigma_s2_curr,sigma_s2_prop)


log_likelihood_pre_fast <- function(new_data,
                                    old_data,
                                    dist_new_old,
                                    log_sigma_n2,
                                    log_sigma_s2,
                                    log_phi){
  
  # ---- transform parameters ----
  sigma_s2 <- exp(log_sigma_s2)
  sigma_n2 <- exp(log_sigma_n2)
  phi      <- exp(log_phi)
  
  N1 <- length(new_data)
  N2 <- length(old_data)
  N  <- N1 + N2
  
  # ---- full covariance ----
  C_full <- cov_mat_fun(dist = dist_new_old,
                        sigma_s2 = sigma_s2,
                        sigma_n2 = sigma_n2,
                        phi = phi,
                        N = N)
  
  # ---- old covariance ----
  C_old <- C_full[(N1+1):N, (N1+1):N]
  
  # ---- joint log-density ----
  y_full <- c(new_data, old_data)
  
  L_full <- chol(C_full)
  alpha_full <- backsolve(L_full,
                          forwardsolve(t(L_full), y_full))
  
  logdet_full <- 2 * sum(log(diag(L_full)))
  
  loglik_full <- -0.5 * (
    crossprod(y_full, alpha_full) +
      logdet_full +
      N * log(2*pi)
  )
  
  # ---- old log-density ----
  L_old <- chol(C_old)
  alpha_old <- backsolve(L_old,
                         forwardsolve(t(L_old), old_data))
  
  logdet_old <- 2 * sum(log(diag(L_old)))
  
  loglik_old <- -0.5 * (
    crossprod(old_data, alpha_old) +
      logdet_old +
      N2 * log(2*pi)
  )
  
  # ---- conditional log-likelihood ----
  return(as.numeric(loglik_full - loglik_old))
}


dist_mat_2_1 <- distm(rbind(locations_2,
                            locations_1), fun = distCosine) / 1000


library(Rcpp.conditional.gp)

all_param_likelihood_pprb<-cbind(log_sigma_s2_post_pprb_1,
                                 log_sigma_n2_post_pprb_1,
                                 log_phi_post_pprb_1)
colnames(all_param_likelihood_pprb)<-c("log_sigma_s2","log_sigma_n2","log_phi_post")

all_param_likelihood_pprb_uni<-unique(all_param_likelihood_pprb)


tic()
log_likelihood_pre2 <- unlist(
  mclapply(
    1:nrow(all_param_likelihood_pprb_uni),
    function(i) {
      
      C_full <- cov_mat_fun(
        dist = dist_mat_2_1,
        sigma_s2 = exp(all_param_likelihood_pprb_uni[i,1]),
        sigma_n2 = exp(all_param_likelihood_pprb_uni[i,2]),
        phi      = exp(all_param_likelihood_pprb_uni[i,3]),
        N = length(data_1) + length(data_2)
      )
      
      log_likelihood_pre_fast_rcpp(
        new_data = data_2,
        old_data = data_1,
        C_full   = C_full
      )
    },
    mc.cores = parallel::detectCores() - 1,
    mc.set.seed=TRUE
  )
)
toc()

all_param_likelihood_pprb_uni<-cbind(all_param_likelihood_pprb_uni,log_likelihood_pre2)
colnames(all_param_likelihood_pprb_uni)[4]<-"log_likeli"
all_param_likelihood_pprb<-cbind(all_param_likelihood_pprb,NA)
colnames(all_param_likelihood_pprb)[4]<-"log_likeli"

for(i in 1:nrow(all_param_likelihood_pprb_uni)){
  idx<-which(all_param_likelihood_pprb[,1]==all_param_likelihood_pprb_uni[i,1]&
               all_param_likelihood_pprb[,2]==all_param_likelihood_pprb_uni[i,2] &
               all_param_likelihood_pprb[,3]==all_param_likelihood_pprb_uni[i,3])
  all_param_likelihood_pprb[idx,4]<-all_param_likelihood_pprb_uni[i,4]
}

log_likelihood_pre<-all_param_likelihood_pprb[,4]
rm(all_param_likelihood_pprb,
   all_param_likelihood_pprb_uni
   )


gc()

### now let's move onto stage 2
### now let's move onto stage 2
### now let's move onto stage 2

log_sigma_n2_post_pprb_2<-numeric(iterations)
log_sigma_s2_post_pprb_2<-numeric(iterations)
log_phi_post_pprb_2<-numeric(iterations)

log_sigma_n2_post_pprb_2[1]<-log_sigma_n2_post_pprb_1[200]
log_sigma_s2_post_pprb_2[1]<-log_sigma_s2_post_pprb_1[200]
log_phi_post_pprb_2[1]<-log_phi_post_pprb_1[200]

accept_pprb_2<-numeric(iterations)


i<-2

current_idx<-200
tic()
for(i in 2:iterations){
  idx<- sample(1:length(log_sigma_n2_post_pprb_1),1,replace=TRUE )
  # 
  proposal_log_sigma_s <- log_sigma_s2_post_pprb_1[idx]
  proposal_log_sigma_n <- log_sigma_n2_post_pprb_1[idx]
  proposal_log_phi     <- log_phi_post_pprb_1[idx]
  
  nom<- log_likelihood_pre[idx]
  denom<- log_likelihood_pre[current_idx]
  # ratio
  ratio<-nom-denom 
  if(log(runif(1)) <= ratio){
    log_sigma_s2_post_pprb_2[i]<-proposal_log_sigma_s
    log_sigma_n2_post_pprb_2[i]<-proposal_log_sigma_n
    log_phi_post_pprb_2[i]<-proposal_log_phi
    accept_pprb_2[i]<-1
    current_idx<-idx
  } else{
    log_sigma_s2_post_pprb_2[i]<-log_sigma_s2_post_pprb_2[i-1]
    log_sigma_n2_post_pprb_2[i]<-log_sigma_n2_post_pprb_2[i-1]
    log_phi_post_pprb_2[i]<-log_phi_post_pprb_2[i-1]
  }
}
toc()
mean(accept_pprb_2)

log_sigma_s2_post_pprb_2<-log_sigma_s2_post_pprb_2[-c(1:burnin)]
log_sigma_n2_post_pprb_2<-log_sigma_n2_post_pprb_2[-c(1:burnin)]
log_phi_post_pprb_2<-log_phi_post_pprb_2[-c(1:burnin)]
accept_pprb_2<-accept_pprb_2[-c(1:burnin)]

rm(current_idx,denom,i,idx,log_likelihood_pre,log_likelihood_pre2)
rm(nom,proposal_log_phi,proposal_log_sigma_n,proposal_log_sigma_s)
rm(ratio)

# third stage -------------------------------------------------------------

dist_mat_3_21 <- distm(rbind(locations_3,
                             locations_2,
                             locations_1), fun = distCosine) / 1000


all_param_likelihood_pprb<-cbind(log_sigma_s2_post_pprb_2,
                                 log_sigma_n2_post_pprb_2,
                                 log_phi_post_pprb_2)
colnames(all_param_likelihood_pprb)<-c("log_sigma_s2","log_sigma_n2","log_phi_post")

all_param_likelihood_pprb_uni<-unique(all_param_likelihood_pprb)

tic()
log_likelihood_pre2 <- unlist(
  mclapply(1:nrow(all_param_likelihood_pprb_uni), function(i) {
    
    C_full <- cov_mat_fun(
      dist = dist_mat_3_21,
      sigma_s2 = exp(all_param_likelihood_pprb_uni[i, 1]),
      sigma_n2 = exp(all_param_likelihood_pprb_uni[i, 2]),
      phi      = exp(all_param_likelihood_pprb_uni[i, 3]),
      N = length(data_1) + length(data_2) + length(data_3)
    )
    
    log_likelihood_pre_fast_rcpp(
      new_data = data_3,
      old_data = c(data_2, data_1),
      C_full = C_full
    )
    
  }, mc.cores = parallel::detectCores()-1, mc.set.seed = TRUE)
)
toc()

all_param_likelihood_pprb_uni<-cbind(all_param_likelihood_pprb_uni,log_likelihood_pre2)
colnames(all_param_likelihood_pprb_uni)[4]<-"log_likeli"
all_param_likelihood_pprb<-cbind(all_param_likelihood_pprb,NA)
colnames(all_param_likelihood_pprb)[4]<-"log_likeli"

for(i in 1:nrow(all_param_likelihood_pprb_uni)){
  idx<-which(all_param_likelihood_pprb[,1]==all_param_likelihood_pprb_uni[i,1]&
               all_param_likelihood_pprb[,2]==all_param_likelihood_pprb_uni[i,2] &
               all_param_likelihood_pprb[,3]==all_param_likelihood_pprb_uni[i,3])
  all_param_likelihood_pprb[idx,4]<-all_param_likelihood_pprb_uni[i,4]
}


log_likelihood_pre<-all_param_likelihood_pprb[,4]
rm(all_param_likelihood_pprb,
   all_param_likelihood_pprb_uni
)


gc()

### now let's move onto stage 3
### now let's move onto stage 3
### now let's move onto stage 3

log_sigma_n2_post_pprb_3<-numeric(iterations)
log_sigma_s2_post_pprb_3<-numeric(iterations)
log_phi_post_pprb_3<-numeric(iterations)

log_sigma_n2_post_pprb_3[1]<-log_sigma_n2_post_pprb_2[200]
log_sigma_s2_post_pprb_3[1]<-log_sigma_s2_post_pprb_2[200]
log_phi_post_pprb_3[1]<-log_phi_post_pprb_2[200]

accept_pprb_3<-numeric(iterations)
i<-2

current_idx<-200
tic()
for(i in 2:iterations){
  idx<- sample(1:length(log_sigma_n2_post_pprb_2),1,replace=TRUE )
  # 
  proposal_log_sigma_s <- log_sigma_s2_post_pprb_2[idx]
  proposal_log_sigma_n <- log_sigma_n2_post_pprb_2[idx]
  proposal_log_phi     <- log_phi_post_pprb_2[idx]
  
  nom<- log_likelihood_pre[idx]
  denom<- log_likelihood_pre[current_idx]
  # ratio
  ratio<-nom-denom 
  if(log(runif(1)) <= ratio){
    log_sigma_s2_post_pprb_3[i]<-proposal_log_sigma_s
    log_sigma_n2_post_pprb_3[i]<-proposal_log_sigma_n
    log_phi_post_pprb_3[i]<-proposal_log_phi
    accept_pprb_3[i]<-1
    current_idx<-idx
  } else{
    log_sigma_s2_post_pprb_3[i]<-log_sigma_s2_post_pprb_3[i-1]
    log_sigma_n2_post_pprb_3[i]<-log_sigma_n2_post_pprb_3[i-1]
    log_phi_post_pprb_3[i]<-log_phi_post_pprb_3[i-1]
  }
}
toc()
accept_pprb_3<-accept_pprb_3[-c(1:burnin)]

mean(accept_pprb_3)

log_sigma_s2_post_pprb_3<-log_sigma_s2_post_pprb_3[-c(1:burnin)]
log_sigma_n2_post_pprb_3<-log_sigma_n2_post_pprb_3[-c(1:burnin)]
log_phi_post_pprb_3<-log_phi_post_pprb_3[-c(1:burnin)]

rm(current_idx,denom,i,idx,log_likelihood_pre,log_likelihood_pre2)
rm(nom,proposal_log_phi,proposal_log_sigma_n,proposal_log_sigma_s,ratio)

#save.image("sss_pprb.RData")



