# Filepath <- 
# setwd(Filepath)

rm(list = ls())

library(yuima)
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
  xlab("Time before present (yrs)") +
  ylab(expression(deltaˆ18 ~ O ~ (permil)))

ggplot(data = df_pro, aes(x = age, y = d18O)) +
  geom_line() +
  xlab("Time before present (yrs)") +
  ylab(expression(deltaˆ18 ~ O ~ (permil)) + 
         theme_minimal())


# Define the SDE parameters and functions
f_drift <- function(x, params) {
  drift <- params[2]*x - params[3]
  return(drift)
}

f_diff <- function(params) {
  diff <- params[4]
  return(diff)
}

f_jump <- function(x, x_val, params) {
  jump <- exp((x_val-x) / params[5])
  return(jump)
}


# Parameters
dt <- 0.05
x_val <- 1.2
T_start <- df_pro[length(df_pro[,1]), 1]
T_end <- df_pro[1, 1]
x0 <- 0
num_steps <- (T_end - T_start) / 0.05
ylim <- c(-6, 6)
parameters <- c(0.5,-0.5,1,0.55)
num_paths <- 10


# Euler-Maruyama negative log likelihood
EM_nll <- function(par, X, dt) {
  N <- length(X) - 1
  
  # Initialize negative log likelihood
  nll <- 0
  
  for (n in 1:N) {
    # Drift function
    drift <- f_drift(X[n], par)
    
    # Diffusion matrix
    SigmaSigma <- f_diff(par)^2

    # Update negative log likelihood
    nll <- nll + 0.5 * log(2*pi*SigmaSigma * dt) +
      0.5 * (X[n + 1] - X[n] - drift * dt) *
      1/(SigmaSigma * dt) * (X[n + 1] - X[n] - drift * dt)
  }
  print(1)
  return(nll)
}

# y-axis limits for plotting
ylim <- c(-6, 6)

# Euler maruyama scheme for optimising values, (0=convergence)
result_EM_1 <- optim(par = parameters, fn = EM_nll,
                     X = df_pro[,2], dt = 0.05)
result_EM_1

result <- result_EM_1$par

T_start <- df_pro[length(df_pro[,1]),1]
T_end <- df_pro[1,1]
x0=0

simulate_EM <- function(T_start, T_end, dt, params, x0) {
  t <- seq(T_start, T_end, by = dt)
  n<-length(t)-1
  dW <- rnorm(n)
  dB <- rpois(n, params[6]*dt)
  X <- numeric(n + 1)
  X[1] <- x0
  dZ <- rnorm(n,params[8],params[7])
  a <- 0
  for (i in 3:(n + 1)) {
    drift_term <- f_drift(X[i-1],params)
    X[i] <- X[i-1] + drift_term * dt + f_diff(params) * dW[i-1] * sqrt(dt)
       + f_jump(X[i], x_val, params) * dB[i] * dZ[i]
  }
  return(data.frame(t = t, X = X))
}
x <- seq(0,50,10)

sim_res <- simulate_EM(T_start, T_end, dt, c(0.901964564, -0.042203700, -0.004701163,  0.550725472, 1,1,1,1), x0)
plot(sim_res$t, sim_res$X, type = 'l', col = 'red', ylim=ylim)

for (i in x){
  sim_res <- simulate_EM(T_start, T_end, dt, c(0.901964564, -0.042203700, -0.004701163,  0.550725472,1 ,1,i,1), x0)
  lines(sim_res$t, sim_res$X, type = 'l', col = 'black', ylim=ylim)
}
sim_res <- simulate_EM(T_start, T_end, dt, c(0.901964564, -0.042203700, -0.004701163,  0.550725472, 1,1,1,1), x0)

plot(sim_res$t, sim_res$X, type = 'l', col = 'red', ylim=ylim)

lines(df_pro$age, df_pro$d18O, type = 'l', col=1, ylim=ylim)


plot(df_pro$age, df_pro$d18O, type = 'l', col=1, ylim=ylim)
lines(sim_res$t, sim_res$X, type = 'l', col = 'red', ylim=ylim)



EM_nll <- function(par, data, dt) {
  N <- length(data) - 1
  
  # Initialize negative log likelihood
  nll <- 0
  
  for (n in 1:N) {
    # Drift function
    drift <- f_drift(data[n], par)
    
    # Diffusion matrix
    Sigma <- f_diff(par)
    
    jump <- f_jump(data[n], x_val, par)
    
    lambda <- par[6]
    
    beta <- par[7]
    
    alpha <- par[8]
    
    pois <- 0
    for (i in 1:7) {  # Adjust the range based on your model
      # Probability of observing i jumps in the interval (dt)
      a <- ((lambda * dt)^(i/jump) / factorial(i/jump) * exp(-lambda * dt)) * 
        (1/(Sigma^2*dt + i*beta^2))*(exp((-1/2)*((data[n+1]-data[n]-drift*dt-i*alpha)^2/(Sigma^2*dt + i*beta^2)))/sqrt(2*pi))
      # Update the Poisson process likelihood
      pois <- pois + a
    }
    nll <- nll - log(pois)
  }
  print(nll)
  return(nll)
}

factorial(1/1000)
parameters = c(1,1,1,1,1,1,1,1)

result_EM_1 <- optim(par = parameters, fn = EM_nll,
                     data = df_pro[,2], dt = 0.05, method = "L-BFGS-B")

result_EM_1

simulate_EM <- function(T_start, T_end, dt, params, x0) {
  t <- seq(T_start, T_end, by = dt)
  n<-length(t)-1
  dW <- rnorm(n)
  dB <- rpois(n, dt)
  dZ <- rnorm(n,params[8],params[7])
  X <- numeric(n + 1)
  X[1] <- x0
  
  a <- 0
  for (i in 2:(n + 1)) {
    drift_term <- f_drift(X[i-1],params)
    X[i] <- X[i-1] + drift_term * dt + f_diff(params) * dW[i-1] * sqrt(dt)
    + f_jump(X[i], x_val, params) * dB[i] * dZ[i]
  }
  return(data.frame(t = t, X = X))
}


a <- simulate_EM(T_start,T_end,dt,result_EM_1$par,x0)
plot(a,type='l')

f_jump(1.5,1.3,c(1,1,1,1,0.2))


tmp <- setYuima(data=setData(zoo(df[,2],order.by=df$age), delta=1/20))

yuima_model <- setModel(drift = "theta2*x-theta3",
                        diffusion = "sigma",
                        jump.coeff = "exp((1.3-x)/beta)",
                        measure = list(intensity = "lambda",
                                       df = list("dnorm(z, alpha, sigma_norm)")),
                                       measure.type="CP",
                                       solve.variable = "x")

str(yuima_model)

yuima <- setYuima(data = tmp@data, model=yuima_model)
str(yuima)

lower <- list(theta2=0.1, theta3=0.1,sigma = 0.1, beta = 0.1, lambda=0.1, alpha=0.1, sigma_norm=0.1)
upper <- list(theta2=10, theta3=10, sigma = 2, beta = 10, lambda=10, alpha=10, sigma_norm=10)
start <- list(theta2=0.7, theta3=1.5, sigma = 0.55, beta = 1, lambda=1, alpha=1, sigma_norm=1)

dt <- 0.05

out <- qmle(yuima, start=start, threshold=sqrt(20), upper=upper, lower=lower, method="L-BFGS-B")

warnings()

out
yuima@model@parameter



Terminal <- 10 
samp <- setSampling(T=Terminal,n=1000)
mod4 <- setPoisson(intensity="beta*(1+sin(lambda*t))", df=list("dconst(z,1)"))
set.seed(123) 
lambda <- 3 
beta <- 5 
y4 <- simulate(mod4, true.par=list(lambda=lambda,beta=beta),sampling=samp) 
par(mfrow=c(2,1)) 
par(mar=c(3,3,1,1)) 
plot(y4) 
f <- function(t) beta*(1+sin(lambda*t)) 
curve(f, 0, Terminal, col="red")




X <- matrix(df_pro[,2])
M <- 1
N <- length(X)

dx <- matrix(0, N, M) # create empty matrix for difference of X states
for (j in 1:M) {
  dx[, j] <- c(X[, j]) # find difference of X states
}

f <- matrix(0, N, 10) # create empty matrix for transition density
estimate <- matrix(0, M, 5)
dx
for (v in 1:M) {
  dif <- dx[, v]
  likelihood <- function(theta, dif, dt) { # create a function for likelihood
    Q1 <- theta[1] # symbolize alpha parameter
    Q2 <- theta[2] # symbolize sigma parameter
    Q3 <- theta[3] # symbolize mu parameter
    Q4 <- theta[4] # symbolize sigma_j parameter
    Q5 <- theta[5] # symbolize lambda parameter
    
    for (ii in 1:N) {
      for (j in 1:10) {
        f[ii, j] <- (exp(-Q5 * dt) * (Q5 * dt)^(j - 1) / factorial(j - 1)) *
          (1 / sqrt(2 * pi * (Q2^2 * dt + (j - 1) * Q4^2))) *
          exp(-((dif[ii] - ((Q1) * dt + (j - 1) * Q3))^2) / (2 * (Q2^2 * dt + (j - 1) * Q4^2)))
      }
    }
    R <- rowSums(f)
    LL <- -sum(log(R)) # find -log likelihood
    print(LL)
    return(LL)
  }
  
  ## Minimize -log likelihood function ##
  estimation <- optim(c(0.5, 2, 0.5, 0.5, 2), likelihood, gr = NULL, dif = dif, dt, method = "L-BFGS-B",
                      lower = c(-Inf, 0, -Inf, 0, -Inf), upper = c(Inf, Inf, Inf, Inf, Inf), hessian = TRUE)
  print(estimation)
  options(scipen = 999)
  estimate[v, ] <- estimation$par # assign each parameter estimation set to estimate matrix
}

parameter_estimations <- colSums(estimate) / M # give the estimated parameter value for each parameter

parameter_estimations
estimation



simulate_EM <- function(T_start, T_end, dt, params, x0) {
  t <- seq(T_start, T_end, by = dt)
  n<-length(t)-1
  dW <- rnorm(n)
  dB <- rpois(n+1, params[6]*dt)
  dZ <- rnorm(n)
  X <- numeric(n + 1)
  X[1] <- x0
  
  for (i in 2:(n + 1)) {
    drift_term <- f_drift(X[i-1],params)
    X[i] <- X[i-1] + drift_term * dt + f_diff(params) * dW[i-1] * sqrt(dt)
    + f_jump(X[i], x_val, params) * (dB[i] * params[7] + sqrt(dB[i]) * params[8] * dZ[i])
  }
  return(data.frame(t = t, X = X, B = dB))
}

a <- simulate_EM(T_start,T_end,dt,c(0.901964564, -0.042203700, -0.004701163,  0.550725472,1, 0.4, 0.6, 1),x0)
plot(a$t, a$X,type='l')







