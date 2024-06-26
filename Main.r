rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(purrr)
library(knitr)
library(kableExtra)
library(tinytex)
library(tidyverse)
library(stats4)
set.seed(123)

# Data load ind
skip = 5

#data files
df <- read.csv("data/Syn_Greenland.csv",sep=",",skip = skip)
df_in <- read.csv("data/China_cave.csv",sep=",",skip = skip)

# Renaming column names for ease of reference
df <- rename(df, age = EDC3.Age..kyr., d18O =  GLT_syn.δ18O....)
df_in <- rename(df_in, age = Age..ka.BP., d18O = δ18O.carb....VPDB.)

# Data processing for duplicated data
preprocess_data <- function(data) {
  
  return_df <- data %>%
    
    mutate(rolling_mean = rollmean(d18O, k = 625, fill = NA, align = "right")) %>%
    
    mutate(d18O = d18O - rolling_mean) %>%
    
    select(-rolling_mean) %>%
    
    arrange(desc(age))
  
  return(data.frame(return_df))
}
# Data processing
df_pro <- drop_na(preprocess_data(df))

# Plotting the data after and before data processing
ggplot(data = df, aes(x = age, y = d18O)) +
  geom_line() +
  xlab("Time before present (kyrs)") +
  ylab(expression(delta^{18} * O ~ (permil)))


ggplot(data = df_pro, aes(x = age, y = d18O)) +
  geom_line() +
  xlab("Time before present (kyrs)") +
  ylab(expression(deltaˆ18 ~ O ~ (permil)))


# Define the SDE parameters and functions
F_drift <- function(X, beta1,beta2,beta4){
  return(beta4*X^3 - beta2*X - beta1)
}

# visualize drift function
plot(seq(-2,2,0.01),F_drift(seq(-2,2,0.01), 1,  1, 0.5), type='l')

F_jump <- function(x, x_val, alpha) {
  jump <- 1 + exp((x_val-x)/alpha)
  return(jump)
}

# Parameters
dt <- 0.05
T_start <- df_pro[length(df_pro[,1]), 1]
T_end <- df_pro[1, 1]
x0 <- 0
num_steps <- (T_end - T_start) / 0.05
ylim <- c(-6, 6)

# y-axis limits for plotting
T_start <- df_pro[length(df_pro[,1]),1]
T_end <- df_pro[1,1]
x0=0

# Euler-Maruyama negative log likelihood for jump-diffusion model
EM_nll_jump <- function(theta, X, dt) {
  N <- length(X) - 1
  beta1<-theta[1]
  beta2<-theta[2]
  beta4<-theta[3]
  sigma<-theta[4]
  lambda <- theta[5]
  alpha <- theta[6]
  x_val <- theta[7]
  # Initialize negative log likelihood
  nll <- 0
  
  for (n in 1:N) {
    # Drift function
    F <- F_drift(X[n], beta1, beta2, beta4)
    alpha_1 <- F_jump(X[n],x_val,alpha)
    
    # Diffusion matrix
    SigmaSigma <- sigma^2
    
    # Update negative log likelihood
    
    for (i in 1:10){
      nll <- nll + 0.5 * log(2 * pi * SigmaSigma * dt) +
        0.5 * (X[n + 1] - X[n] - F * dt - alpha_1 * i) *
        1/(SigmaSigma * dt) * (X[n + 1] - X[n] - F * dt - alpha_1 * i) - 
        log((lambda * dt)^(i)/factorial(i)*exp(-lambda * dt))
    }
  }
  return(nll)
}

# optimizing algorithm for finding MLE parameters
#                            beta1,        beta2,       beta3,      sigma, lambda,gamma, x∗
initial_beta_jump_bfsg <- c(−0.00447839, 0.03581118, −0.00069272, 0.55070983, 1, 2, 1.8)

result_EM_1_jump_bfsg <- optim(par = initial_beta_jump_bfsg, 
                               fn = EM_nll_jump, 
                               X = df_pro[,2], 
                               dt = 0.05,
                               control=list(maxit=1000),
                               method='L-BFGS-B',
                               lower=c(-5,-5,-5,0.0001,0.0001,-5,-5),
                               upper=c(5,5,5,5,5,5,5))


initial_beta_jump <- result_EM_1_jump_bfsg$par
(result_EM_1_jump <- optim(par = initial_beta_jump_bfsg, 
                           fn = EM_nll_jump, 
                           X = df_pro[,2], 
                           dt = 0.05,
                           control=list(maxit=1000)))

result_jump <- result_EM_1_jump$par

# defining T, dt X_0 and iterations for bootstrap analysis
T <- length(df_pro[,2])*0.05
dt <- 0.05
x0 <- 0
iterations <- 100

# 1D EM-simulation of data
simulate_EM_jump <- function(T, dt, b1, b2, b4, sigma, lambda, alpha, x_val, x0) {
  t <- seq(0, T, by = dt)
  n <- length(t) - 1
  dW <- rnorm(n) * sqrt(dt)  # Correct Wiener increment
  dN <- rpois(n, lambda * dt)  # Poisson increment
  X <- numeric(n + 1)
  X[1] <- x0
  jump <- 0
  for (i in 2:(n + 1)) {
    drift_term <- F_drift(X[i-1], b1, b2, b4)
    jump_term <- F_jump(X[i-1], x_val, alpha)
    jump <- jump + jump_term
    X[i] <- X[i-1] + drift_term * dt + sigma * dW[i-1] + dN[i-1] * jump_term
  }
  return(data.frame(t = T - t, X = X))  # Adjust time vector for correct plotting
}

# plot of simulated data with MLE parameters and subsampling
a <- simulate_EM_jump(T, dt/10, result_jump[1], result_jump[2],result_jump[3], 
                      result_jump[4], result_jump[5], result_jump[6], result_jump[7],x0)

plot(a[seq(1, nrow(a), 10), ]$t, a[seq(1, nrow(a), 10), ]$X,type='l', col='red',
     xlab = "Time before present (kyrs)",ylab = expression(deltaˆ18 ~ O ~ (permil)))

#plot of observed data 
plot(df_pro[,1], df_pro[,2], type='l')

# Euler-Maruyama negative log likelihood for diffusion model
EM_nll <- function(theta, X, dt) {
  N <- length(X) - 1
  beta1<-theta[1]
  beta2<-theta[2]
  beta4<-theta[3]
  sigma<-theta[4]
  
  # Initialize negative log likelihood
  nll <- 0
  
  for (n in 1:N) {
    # Drift function
    F <- F_drift(X[n], beta1, beta2, beta4)
    
    # Diffusion matrix
    SigmaSigma <- sigma^2
    
    # Update negative log likelihood
    nll <- nll + 0.5 * log(SigmaSigma * dt) +
      0.5 * (X[n + 1] - X[n] - F * dt) *
      1/(SigmaSigma * dt) * (X[n + 1] - X[n] - F * dt)
  }
  return(nll)
}

# optimizing for MLE parameters in diffusion model
initial_beta <- c(1,1,1,1)
(result_EM_1 <- optim(par = initial_beta, fn = EM_nll,
                      X = df_pro[,2], dt = 0.05))
result <- result_EM_1$par


# simulating data with the MLE parameters
simulate_EM <- function(T, dt, b1, b2, b4, sigma, x0) {
  t <- seq(0, T, by = dt)
  n<-length(t)-1
  dW <- rnorm(n)  
  X <- numeric(n + 1)
  X[1] <- x0
  for (i in 2:(n + 1)) {
    drift_term <- F_drift(X[i-1], b1, b2, b4)
    X[i] <- X[i-1] + drift_term * dt + sigma * dW[i-1] * sqrt(dt)
  }
  return(data.frame(t = T-t, X = X))
}
# Subsampling
step <- simulate_EM(800, dt/10, result[1], result[2], 
                    result[3],result[4], x0)

# taking every 10'th value since we produces 10 times as many data points per subsampling
result_eul_11d <- step[seq(1, nrow(step), 10), ]

# plot of simulated data versus observed data
plot(df_pro[,1],df_pro[,2], type='l', col='black')
lines(result_eul_11d, type='l', col='darkblue')

# plot of simulated data formatted nicely
plot(result_eul_11d, type='l', col='darkblue', xlab = "Time before present (kyrs)",ylab = expression(deltaˆ18 ~ O ~ (permil)))

# Density plots of both models versus observed data
plot(density(df_pro[,2]), col='darkred', main="", xlab=expression(deltaˆ18 ~ O ~ (permil)))
lines(density(a[seq(1, nrow(a), 10), ]$X), main="")

plot(density(result_eul_11d[1:length(df_pro[,2]),2]),main="", xlab=expression(deltaˆ18 ~ O ~ (permil)))
lines(density(df_pro[,2]), col='darkred')


# qqplot of simulated data versus observed data for both models
qqplot(a[seq(1, nrow(a), 10), ]$X,df_pro[,2], xlab = "Simulated data", ylab = "observed data")

abline(lm(sort(df_pro[,2]) ~ sort(a[seq(0, nrow(a), 10), ]$X)), col = "darkred")

qqplot(result_eul_11d[1:length(df_pro[,2]),2],df_pro[,2], 
       xlab = "Simulated data", ylab = "observed data")

abline(lm(sort(df_pro[,2]) ~ sort(result_eul_11d[1:length(df_pro[,2]),2])), col = "darkred")


# Bootstrap analysis of the MLE parameters 
boot_params <- function(iter, init_params, mult){
  boot_res <- matrix(0, iter, length(init_params)+1)
  boot_res[1,] <- c(init_params,0)
  
  for (i in 2:iter){
    data <- simulate_EM(800, dt/mult, boot_res[1,1], boot_res[1,2], 
                        boot_res[1,3],boot_res[1,4], x0)
    
    optim_data <- data[seq(1, nrow(data), 10), ]
    
    res <- optim(par = boot_res[1,1:4], fn = EM_nll,
                 X = optim_data$X, dt = 0.05)
    
    boot_res[i,] <- c(res$par, res$convergence)
  }
  return(boot_res)
}

Boot <- boot_params(100, result, 10)

boot_params_jump <- function(iter, init_params, mult){
  boot_res <- matrix(0, iter, length(init_params)+1)
  boot_res[1,] <- c(init_params,0)
  
  for (i in 2:iter){
    data <- simulate_EM_jump(800, dt/mult, boot_res[1,1], boot_res[1,2], 
                             boot_res[1,3],boot_res[1,4],
                             boot_res[1,5],boot_res[1,6],
                             boot_res[1,7], x0)
    
    optim_data <- data[seq(1, nrow(data), 10), ]
    
    res <- optim(par = boot_res[1,1:7], fn = EM_nll_jump,
                 X = optim_data$X, dt = 0.05)
    
    boot_res[i,] <- c(res$par, res$convergence)
  }
  return(boot_res)
}
Boot_jump <- boot_params_jump(100, result_jump, 10)
Boot_jump

quantile(Boot_jump[,2], c(0.025, 0.975))
quant_jump <- matrix(0,2,7)
for (i in 1:7){
  quant_jump[,i] <- quantile(Boot_jump[,i], c(0.025, 0.975))
}

quantile(Boot[,2], c(0.025, 0.975))
quant <- matrix(0,2,4)
for (i in 1:4){
  quant[,i] <- quantile(Boot[,i], c(0.025, 0.975))
}

par(mfrow=c(2,4))

plot(density(Boot_jump[,1]), main='beta1')
plot(density(Boot_jump[,2]), main='beta2')
plot(density(Boot_jump[,3]), main='beta3')
plot(density(Boot_jump[,4]), main='sigma')
plot(density(Boot_jump[,5]), main='lambda')
plot(density(Boot_jump[,6]), main='gamma')
plot(density(Boot_jump[,7]), main='x*')

# transition Densities calculated for each datapoint but not visualized in qq plot since i could not figure out how to # compare to convolution of Poisson- and normal-distribution
trans_dens_jump <- function(X,dt,theta){
  n <- length(X)-1
  beta1<-theta[1]
  beta2<-theta[2]
  beta4<-theta[3]
  sigma<-theta[4]
  lambda <- theta[5]
  alpha <- theta[6]
  x_val <- theta[7]
  
  
  jump <- F_jump(X, x_val, alpha)
  drift <- F_drift(X, beta1,beta2,beta4)
  
  mean_val <- X + unlist(lapply(drift, function(df) df * dt))
  
  jump_mean <- unlist(lapply(jump, function(df) df*c(0:10)))
  
  var_val <- sigma^2*dt
  
  density_vals <- matrix(0,nrow = length(X), ncol = 11)
  for (i in 0:10){
    density_vals[,i+1] <- pnorm(q = X, mean = mean_val + jump_mean[seq(1+i, length(jump_mean), 11)], 
                                sd = sqrt(var_val)) * ppois(q = i, lambda * dt)
    
  }
  return(rowSums(density_vals))
}

# calculating transition densities
dens <- trans_dens_jump(df_pro[,2], dt, result_jump)

# Plotting density of transition densities
plot(density(dens))


# calculating the accuracy lost by truncating the likelihood of more then 10 jumps in the MLE for the jump diffusion model with lambda = 59.99457
Pois_trunc <- 1-sum(dpois(0:10,59.99457*0.05))
print(Pois_trunc)
