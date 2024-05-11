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
  xlab("Time before present (yrs)") +
  ylab(expression(deltaˆ18 ~ O ~ (permil)))

ggplot(data = df_pro, aes(x = age, y = d18O)) +
  geom_line() +
  xlab("Time before present (yrs)") +
  ylab(expression(deltaˆ18 ~ O ~ (permil)) + 
         theme_minimal())


# Define the SDE parameters and functions
f_drift <- function(x, params) {
  drift <- -params[1]*x^3 + params[2]*x - params[3]
  return(drift)
}

f_diff <- function(params) {
  diff <- params[4]
  return(diff)
}

f_jump <- function(x, x_val, params) {
  jump <- exp((x - x_val) / params[5])
  return(jump)
}

# Define the Euler-Maruyama approximation
euler_maruyama_step <- function(x, dt, params, dW, dN, x_val) {
  return(x + f_drift(x, params) * dt + f_diff(params) * dW + f_jump(x, x_val, params) * dN)
}

# Generate simulated data
simulate_data <- function(T, dt, params, x0, x_val) {
  timesteps <- T / dt
  t <- seq(0, T, by = dt)
  dW <- rnorm(length(t), mean = 0, sd = sqrt(dt))
  dN <- rpois(length(t), lambda = dt)
  x <- numeric(length(t))
  x[1] <- x0
  for (i in 1:(length(t) - 1)) {
    x[i + 1] <- euler_maruyama_step(x[i], dt, params, dW[i], dN[i], x_val)
  }
  return(list(t = t, x = x))
}

# Define the log-likelihood function for parameter estimation
log_likelihood <- function(params, data) {
  T <- data$T
  dt <- data$dt
  x0 <- data$x0
  x_val <- data$x_val
  observed_data <- data$observed_data
  
  timesteps <- T / dt
  x_simulated <- numeric(length(observed_data))
  x_simulated[1] <- x0
  log_likelihood_value <- 0.0
  for (i in 1:(length(observed_data) - 1)) {
    dW <- rnorm(1, mean = 0, sd = sqrt(dt))
    dN <- rpois(1, lambda = dt)
    x_simulated[i + 1] <- euler_maruyama_step(x_simulated[i], dt, params, dW, dN, x_val)
    log_likelihood_value <- log_likelihood_value - 0.5 * log(2 * pi * dt) - 0.5 * ((observed_data[i + 1] - x_simulated[i + 1]) ^ 2) / dt
    print(log_likelihood_value)
  }
  return(-log_likelihood_value)
}

# Maximum likelihood parameter estimation
estimate_parameters <- function(observed_data, T, dt, x0, x_val, initial_guess) {
  result <- optim(initial_guess, log_likelihood, data = list(T = T, dt = dt, x0 = x0, x_val = x_val, observed_data = observed_data))
  return(result$par)
}

# Example usage
set.seed(42)  # For reproducibility

# Parameters
dt <- 0.05
x_val <- 1.2
T_start <- df_pro[length(df_pro[,1]), 1]
T_end <- df_pro[1, 1]
x0 <- 0
num_steps <- (T_end - T_start) / 0.05
ylim <- c(-6, 6)
parameters <- c(2,2,2,0.55,2,2,2,2)
num_paths <- 10

# Generate simulated data
observed_data <- df_pro[,2]

# Estimate parameters from observed data
estimated_parameters <- estimate_parameters(observed_data, T_end - T_start, dt, x0, x_val, parameters)
print("Estimated parameters:", estimated_parameters)






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
    
    lambda <- par[6]
    
    pois <- 0
    for (i in 1:7){
      a <- (lambda * dt)^(i/f_jump(X[n+1]-X[n], x_val, par))/factorial(i/f_jump(X[n+1]-X[n], x_val, par))*exp(-lambda * dt)
      pois <- pois + a
    }
    # Update negative log likelihood
    nll <- nll + 0.5 * log(2*pi*SigmaSigma * dt) +
      0.5 * (X[n + 1] - X[n] - drift * dt - f_jump(X[n],x_val,par)*lambda) *
      1/(SigmaSigma * dt) * (X[n + 1] - X[n] - drift * dt- f_jump(X[n],x_val,par)*lambda) - log(pois)
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
  dB <- rpois(n, dt)
  X <- numeric(n + 1)
  X[1] <- x0
  
  a <- 0
  for (i in 3:(n + 1)) {
    drift_term <- f_drift(X[i-1],params)
    X[i] <- X[i-1] + drift_term * dt + f_diff(params) * dW[i-1] * sqrt(dt)
       + f_jump(X[i], x_val, params) * dB[i]
  }
  return(data.frame(t = t, X = X))
}

sim_res <- simulate_EM(T_start, T_end, dt, result, x0)

sim_res

plot(sim_res$t, sim_res$X, type = 'l', col = 'red')
lines(df_pro$age, df_pro$d18O, type = 'l', col=1, ylim=ylim)


plot(df_pro$age, df_pro$d18O, type = 'l', col=1, ylim=ylim)
lines(sim_res$t, sim_res$X, type = 'l', col = 'red', ylim=ylim)





EM_nll <- function(par, X, dt) {
  N <- length(X) - 1
  
  # Initialize negative log likelihood
  nll <- 0
  
  for (n in 1:N) {
    # Drift function
    drift <- f_drift(X[n], par)
    
    # Diffusion matrix
    Sigma <- f_diff(par)
    
    lambda <- par[6]
    
    beta <- par[7]
    
    alpha <- par[8]
    
    pois <- 0
    for (i in 1:7){
      a <- (lambda * dt)^(i/f_jump(X[n+1]-X[n], x_val, par))/factorial(i/f_jump(X[n+1]-X[n], x_val, par))*exp(-lambda * dt)
      b <- (1/(sqrt(2*pi)*(Sigma*dt+i*beta)))*exp((-1/2)*((X[n+1]-X[n]-drift*dt-alpha*i)^2/(Sigma*dt+beta*i)^2))
      pois <- pois + a*b
    }
    # Update negative log likelihood
    nll <- nll + log(pois)
  }
  print(1)
  return(nll)
}

EM_nll(parameters,df_pro[,2],dt)

result_EM_1 <- optim(par = parameters, fn = EM_nll,
                     X = df_pro[,2], dt = 0.05,method="L-BFGS-B")
?optim
result_EM_1


simulate_EM <- function(T_start, T_end, dt, params, x0) {
  t <- seq(T_start, T_end, by = dt)
  n<-length(t)-1
  dW <- rnorm(n)
  dB <- rpois(n, dt)
  dZ <- rnorm(n,params[7],params[8])
  X <- numeric(n + 1)
  X[1] <- x0
  
  a <- 0
  for (i in 3:(n + 1)) {
    drift_term <- f_drift(X[i-1],params)
    X[i] <- X[i-1] + drift_term * dt + f_diff(params) * dW[i-1] * sqrt(dt)
    + f_jump(X[i], x_val, params) * dB[i] * dZ[i]
  }
  return(data.frame(t = t, X = X))
}


a <- simulate_EM(T_start,T_end,dt,result_EM_1$par,x0)
plot(a,type='l')

calculate_transition_density <- function(x0, T_start, T_end, dt, params) {

  n_steps <- round((T_end - T_start) / dt)
  
  density <- numeric(n_steps + 1)
  density[1] <- x0
  
  for (i in 1:n_steps) {
    density[i + 1] <- rnorm(1, mean = density[i] + f_drift(density[i], params) *
                              dt, f_diff(params)*sqrt(dt)) 
  }
  return(density)
}

transition_densities <- calculate_transition_density(x0,T_start,T_end,dt,result)

qqnorm(transition_densities)
qqline(transition_densities)

