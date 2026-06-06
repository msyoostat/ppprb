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
set.seed(2)

# Load data
load("/Users/myungsooyoo/Desktop/my_stuff/research/post_doc_2/gp/sss_7.RData")
## plot

common_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12)
)

temptemp<-cbind.data.frame(lon=locations[,1],lat=locations[,2],sss=data_mat)
fig1<-ggplot() +
  geom_raster(data = temptemp, aes(x = lon, y = lat, fill = sss)) +
  scale_fill_viridis_c(option = "plasma", name = "SSS") +
  # geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
  #              fill = NA, color = "black", size = 0.5) +
  coord_quickmap(xlim = c(-41, -25), ylim = c(21, 28)) +
  labs(x = "Longitude", y = "Latitude", title = "Sea Surface Salinity (SSS)") +
  theme_minimal() +
  common_theme
rm(temptemp)

full_data <- data_mat
full_locations <- locations

data<-full_data
# Scale data
data_mean <- mean(data)
data_sd   <- sd(data)
data_scaled <- (data - data_mean) / data_sd

# Distance matrix
dist_mat <- distm(locations, fun = distCosine) / 1000
N <- nrow(locations)

# ---------------- Covariance matrix -----------------
cov_mat_fun <- function(dist, sigma_s2, sigma_n2, phi, N){
  sigma_s2 * (1 + dist/phi) * exp(-dist/phi) + sigma_n2 * diag(N)
}

# ---------------- Log-likelihood -----------------
loglik_fun <- function(y, sigma_s2, sigma_n2, phi,dist_mat,N){
  C_mat <- cov_mat_fun(dist_mat, sigma_s2, sigma_n2, phi, N)
  Rchol <- chol(C_mat)
  z <- forwardsolve(t(Rchol), y)
  quad <- sum(z^2)
  logdet <- 2 * sum(log(diag(Rchol)))
  val <- -0.5 * (N * log(2*pi) + logdet + quad)
  return(val)
}

# ---------------- MCMC settings -----------------
iterations <- 30000
burnin <- 5000
log_sigma_s2_post <- numeric(iterations)
log_sigma_n2_post <- numeric(iterations)
log_phi_post <- numeric(iterations)

# Initial values
log_sigma_s2_post[1] <- log(0.6)
log_sigma_n2_post[1] <- log(0.05)
log_phi_post[1] <- log(median(dist_mat)/3)

jump_sigma_s2 <- 0.06 #
jump_sigma_n2 <- 0.04 #
jump_phi <- 0.1 #

accept <- numeric(iterations)

# Priors
mu_log_sigma_s <- log(0.6) # 
sd_log_sigma_s <- 1
mu_log_sigma_n <- log(0.05) #
sd_log_sigma_n <- 1
mu_log_phi     <- log(median(dist_mat)/5)
sd_log_phi <- 1


# ---------------- MCMC -----------------
tic()
for(i in 2:iterations){
  
  # --- Propose new values ---
  proposal_log_sigma_s <- rnorm(1, log_sigma_s2_post[i-1], jump_sigma_s2)
  proposal_log_sigma_n <- rnorm(1, log_sigma_n2_post[i-1], jump_sigma_n2)
  proposal_log_phi     <- rnorm(1, log_phi_post[i-1], jump_phi)
  
  # Exponentiate to get actual parameters
  sigma_s2_prop <- exp(proposal_log_sigma_s)
  sigma_n2_prop <- exp(proposal_log_sigma_n)
  phi_prop      <- exp(proposal_log_phi)
  
  sigma_s2_curr <- exp(log_sigma_s2_post[i-1])
  sigma_n2_curr <- exp(log_sigma_n2_post[i-1])
  phi_curr      <- exp(log_phi_post[i-1])
  
  # --- Log posterior ---
  log_post_curr <- loglik_fun(data_scaled, sigma_s2_curr, sigma_n2_curr, phi_curr,
                              dist_mat,N) +
    dnorm(log_sigma_s2_post[i-1], mu_log_sigma_s, sd_log_sigma_s, log=TRUE) +
    dnorm(log_sigma_n2_post[i-1], mu_log_sigma_n, sd_log_sigma_n, log=TRUE) +
    dnorm(log_phi_post[i-1], mu_log_phi, sd_log_phi, log=TRUE)
  
  log_post_prop <- loglik_fun(data_scaled, sigma_s2_prop, sigma_n2_prop, phi_prop,
                              dist_mat,N) +
    dnorm(proposal_log_sigma_s, mu_log_sigma_s, sd_log_sigma_s, log=TRUE) +
    dnorm(proposal_log_sigma_n, mu_log_sigma_n, sd_log_sigma_n, log=TRUE) +
    dnorm(proposal_log_phi, mu_log_phi, sd_log_phi, log=TRUE)
  
  # --- MH accept/reject ---
  log_ratio <- log_post_prop - log_post_curr
  if(log(runif(1)) < log_ratio){
    log_sigma_s2_post[i] <- proposal_log_sigma_s
    log_sigma_n2_post[i] <- proposal_log_sigma_n
    log_phi_post[i]     <- proposal_log_phi
    accept[i] <- 1
  } else {
    log_sigma_s2_post[i] <- log_sigma_s2_post[i-1]
    log_sigma_n2_post[i] <- log_sigma_n2_post[i-1]
    log_phi_post[i]     <- log_phi_post[i-1]
  }
  if(i %% 500 == 0){
    cat("iter:", i, "accept rate:", mean(accept[(i-499):i]), "\n")
  }
}
toc()
accept <- accept[-c(1:burnin)]
mean(accept)

log_sigma_s2_post<-log_sigma_s2_post[-c(1:burnin)]
log_sigma_n2_post<-log_sigma_n2_post[-c(1:burnin)]
log_phi_post<-log_phi_post[-c(1:burnin)]



rm(i,log_post_curr,log_post_prop,log_ratio)
rm(phi_curr,phi_hat,phi_prop,proposal_log_phi,proposal_log_sigma_n,proposal_log_sigma_s)
rm(sigma_n2_curr,sigma_n2_hat,sigma_n2_prop,sigma_s2_curr,sigma_s2_hat,sigma_s2_prop)

#save.image("sss_mh.RData")
