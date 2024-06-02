# Filepath <- 
# setwd(Filepath)

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

plot(seq(-10,10,0.01),F_drift(seq(-10,10,0.01), 0.03,  0.08, -0.06), type='l')

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

# Euler-Maruyama negative log likelihood
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

#                            β1,             β2,         β3,          σ,      λ    γ,   x∗
initial_beta_jump_bfsg <- c(-0.00447839, 0.03581118, -0.00069272, 0.55070983, 0.1, 0.1, 0.1)

result_EM_1_jump_bfsg <- optim(par = initial_beta_jump_bfsg, 
                          fn = EM_nll_jump, 
                          X = df_pro[,2], 
                          dt = 0.05,
                          control=list(maxit=1000),
                          method='L-BFGS-B',
                          lower=c(-5,-5,-5,0.0001,0.0001,-5,-5),
                          upper=c(5,5,5,5,5,5,5))



?optim
initial_beta_jump <- result_EM_1_jump_bfsg$par
(result_EM_1_jump <- optim(par = initial_beta_jump_bfsg, 
                          fn = EM_nll_jump, 
                          X = df_pro[,2], 
                          dt = 0.05,
                          control=list(maxit=1000)))

result_jump <- result_EM_1_jump$par
result_jump <- c(-0.01073449,  0.04937238,  0.02828193, 0.63363157, 
                 0.04917792,  0.18292174,  0.03550532)

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

simulate_Milstein <- function(T, dt, b1, b2, b4, sigma, lambda, alpha, x_val, x0) {
  t <- seq(0, T, by = dt)
  n <- length(t) - 1
  dW <- rnorm(n) * sqrt(dt)  # Correct Wiener increment
  dN <- rpois(n, lambda * dt)  # Poisson increment
  X <- numeric(n + 1)
  X[1] <- x0
  
  # Define the diffusion function
  diffusion_term <- sigma
  

  for (i in 2:(n + 1)) {
    drift_term <- F_drift(X[i-1], b1, b2, b4)
    jump_term <- F_jump(X[i-1],x_val, alpha)
    
    # Milstein scheme with jump term
    X[i] <- X[i-1] + drift_term * dt + diffusion_term * dW[i-1] +
      0.5 * diffusion_term * (dW[i-1]^2 - dt) + dN[i-1] * jump_term
  }
  
  return(data.frame(t = T - t, X = X))  # Adjust time vector for correct plotting
}

a <- simulate_EM_jump(T, dt/10, result_jump[1], result_jump[2],result_jump[3], 
                 result_jump[4], result_jump[5], result_jump[6], result_jump[7],x0)

plot(a[seq(1, nrow(a), 10), ]$t, a[seq(1, nrow(a), 10), ]$X,type='l', col='red',
     xlab = "Time before present (kyrs)",ylab = expression(deltaˆ18 ~ O ~ (permil)))


b <- simulate_Milstein(T, dt, result_jump[1], result_jump[2],result_jump[3], 
                  result_jump[4], result_jump[5], result_jump[6], result_jump[7]
                  ,x0)
plot(b$t, b$X,type='l', col='red')


plot(df_pro[,1], df_pro[,2], type='l')


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

initial_beta <- c(1,1,1,1)
(result_EM_1 <- optim(par = initial_beta, fn = EM_nll,
                     X = df_pro[,2], dt = 0.05))
result <- result_EM_1$par

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

step <- simulate_EM(800, dt/10, result[1], result[2], 
                              result[3],result[4], x0)

result_eul_11d <- step[seq(1, nrow(step), 10), ]

plot(df_pro[,1],df_pro[,2], type='l', col='black')
lines(result_eul_11d, type='l', col='darkblue')

plot(result_eul_11d, type='l', col='darkblue',
     xlab = "Time before present (kyrs)",ylab = expression(deltaˆ18 ~ O ~ (permil)))

plot(density(df_pro[,2]), col='darkred', main="", xlab=expression(deltaˆ18 ~ O ~ (permil)))
lines(density(a[seq(1, nrow(a), 10), ]$X), main="")

plot(density(result_eul_11d[1:length(df_pro[,2]),2]),main="", xlab=expression(deltaˆ18 ~ O ~ (permil)))
lines(density(df_pro[,2]), col='darkred')



qqplot(a[seq(1, nrow(a), 10), ]$X,df_pro[,2], xlab = "Simulated data", ylab = "observed data")

abline(lm(sort(df_pro[,2]) ~ sort(a[seq(0, nrow(a), 10), ]$X)), col = "darkred")


qqplot(result_eul_11d[1:length(df_pro[,2]),2],df_pro[,2], 
       xlab = "Simulated data", ylab = "observed data")

abline(lm(sort(df_pro[,2]) ~ sort(result_eul_11d[1:length(df_pro[,2]),2])), col = "darkred")



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



trans_dens_jump <- function(X,dt,theta){
  n <- length(X)-1
  beta1<-theta[1]
  beta2<-theta[2]
  beta4<-theta[3]
  sigma<-theta[4]
  lambda <- theta[5]
  alpha <- theta[6]
  x_val <- theta[7]
  
  drifts <- F_drift(X, beta1,beta2,beta4)
  dif <- sigma^2*dt
  jumps <- F_jump(X, x_val, alpha)
  
  trans <- numeric(n)
  for (i in 2:(n+1)){
    jumps_track <- 0
    for (j in 1:10){
      jumps_track <- jumps_track + ((1/(sqrt(2*pi*dif)))*
                                      exp(-(1/2)*(X[i]-X[i-1]-drifts[i-1]*dt-jumps[i-1]*j)^2/(dif)) * 
                                      (((lambda * dt)^j/factorial(j))*exp(-lambda*dt)))
    }
    trans[i-1] <- jumps_track
  }
  return(trans)
}

result_jump

plot(seq(-20,20,0.01), F_drift(seq(-20,20,0.01), result_jump[1],result_jump[2],result_jump[3]))

a <- simulate_EM_jump(T, dt/10, 0.003, 0.05,-0.005, 
                      0.1, 1, 2, 1,x0)

plot(a[seq(1, nrow(a), 10), ]$t, a[seq(1, nrow(a), 10), ]$X,type='l', col='red',
     xlab = "Time before present (kyrs)",ylab = expression(deltaˆ18 ~ O ~ (permil)))


trans_dens <- trans_dens_jump(df_pro[,2], dt, result_jump)
plot(density(trans_dens))

Boot_jump <- boot_params_jump(100, result_jump, 10)
Boot_jump

quantile(Boot_jump[,2], c(0.025, 0.975))
quant_jump <- matrix(0,2,7)
for (i in 1:7){
  quant_jump[,i] <- quantile(Boot_jump[,i], c(0.025, 0.975))
}

quant_jump

plot(density(Boot_jump[,7]))

quantile(Boot[,2], c(0.025, 0.975))
quant <- matrix(0,2,4)
for (i in 1:4){
  quant[,i] <- quantile(Boot[,i], c(0.025, 0.975))
}

quant_jump
names(data.frame(quant))

par(mfrow=c(2,4))


plot(density(Boot_jump[,1]), main='beta1')
plot(density(Boot_jump[,2]), main='beta2')
plot(density(Boot_jump[,3]), main='beta3')
plot(density(Boot_jump[,4]), main='sigma')
plot(density(Boot_jump[,5]), main='lambda')
plot(density(Boot_jump[,6]), main='gamma')
plot(density(Boot_jump[,7]), main='x*')


#Plot af boot data

ggplot(data = data.frame(quant), aes(x = forcats::fct_inorder(Variable))) +
  geom_point( aes(x = forcats::fct_inorder(Variable), y = Optimal), shape = 16, color = "green") +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1, ) +
  geom_point( aes(x = forcats::fct_inorder(Variable), y = Mean), shape = 16, color = "black") +
  labs(x = "Variable", y = "Estimated value") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.margin = unit(1, "lines")) +  # Adjust panel margin to prevent clipping
  theme(axis.ticks.x = element_blank())  # Remove x-axis ticks



a <- c()
for (i in 1:1000){
  a[i] <- sum(rpois(15237,60*dt))
}

1-sum(dpois(0:10,60*0.05))
plot(density(a))
