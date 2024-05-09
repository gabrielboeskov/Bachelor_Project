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
df_pro <- df
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


# Drift function
f_drift <- function(x,params){
  drift <- -params[1]*x^3 + params[2]*x - params[3]
  return(drift)
}

# Diffusion function
f_diff <- function(params){
  diff <- params[4]
  return(diff)
}

# Jump process function
f_jump <- function(x,x_val,params){
  jump <- exp((x-x_val)/params[5])
  return(jump)
}

dt <- 0.05
x_val <- 1.3
T_start <- df_pro[length(df_pro[,1]),1]
T_end <- df_pro[1,1]
x0 <- 0
num_steps <- (T_end-T_start)/0.05
ylim <- c(-6, 6)
parameters <- c(0.03048339, -0.391566545, -0.126585122,  0.5506000193,0,0)
num_paths <- 10

# Step 1: Generate Sample Paths
generate_sample_paths <- function(par, dt, num_paths, num_steps) {
  paths <- matrix(0, nrow = num_paths, ncol = num_steps + 2)
  for (i in 1:num_paths) {
    # Initialize path
    paths[i, 1] <- x0  # Initial condition, starting at 0

    for (j in 1:num_steps) {

      # Generate Brownian increment
      dW <- sqrt(dt) * rnorm(1, mean = 0, sd = 1)
      
      # Generate Poisson increment
      N <- rpois(1, lambda = dt)
      
      # Update path using Euler-Maruyama method
      paths[i, j + 1] <- paths[i, j] + f_drift(paths[i, j], par) * dt + f_diff(par) * dW + f_jump(paths[i, j], x_val, par) * N
    }
  }
  return(paths)
}

# Step 2: Estimate Density
estimate_density <- function(sample_paths, bins = 50) {
  density_estimates <- list()
  for (i in 1:ncol(sample_paths)) {
    density_estimates[[i]] <- density(sample_paths[, i], n = bins)
  }
  return(density_estimates)
}

# Step 3: Evaluate Likelihood
evaluate_likelihood <- function(observed_data, sample_paths) {
  # Assuming observed_data is a vector
  likelihoods <- sapply(2:ncol(sample_paths), function(i) {
    dnorm(observed_data[i], mean = mean(sample_paths[, i]), sd = sd(sample_paths[, i]), log = TRUE)
  })
  return(-likelihoods)
}

likelihood_function <- function(parameters) {
  sample_paths <- generate_sample_paths(parameters, dt = dt, num_paths = num_paths, num_steps = num_steps)
  likelihoods <- evaluate_likelihood(observed_data, sample_paths)
  print(1)
  return(sum(likelihoods))
}


#generate samples
sample_paths <- generate_sample_paths(parameters, dt, 10, num_steps)
plot(sample_paths[1,],type='l')
for (i in 2:num_paths){
  lines(sample_paths[i,],type='l')
}

# Estimate density
density_estimates <- estimate_density(sample_paths)

plot(density_estimates[[1]],type='l')
for (i in 2:num_paths){
  lines(density_estimates[[i]],type='l')
}

# Evaluate likelihood
observed_data <- df_pro[,2]  # Assuming observed data is the first sample path
likelihoods <- evaluate_likelihood(observed_data, sample_paths)
likelihoods

# Optimization
result <- optim(par = parameters, fn = likelihood_function)
print(result$par)

result

sample_paths <- generate_sample_paths(result$par, dt, 10, num_steps)

plot(df_pro[,2],type='l')
plot(sample_paths[1,],type='l')
for (i in 2:num_paths){
  lines(sample_paths[i,],type='l')
}

?optim
install.packages("GenSA")
library(GenSA)

result_gensa <- GenSA(par = parameters, fn = likelihood_function, lower = c(-5,-5,-5,0,-5), upper = c(5,5,5,5,5))

print(result_gensa$par)
































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
    '
    lambda <- par[6]
    
    pois <- 0
    for (i in 1:5){
      a <- (lambda * dt)^(i/f_jump(X[n+1]-X[n], x_val, par))/factorial(i/f_jump(X[n+1]-X[n], x_val, par))*exp(-lambda * dt)
      pois <- pois + a
    }'
    # Update negative log likelihood
    nll <- nll + 0.5 * log(2*pi*SigmaSigma * dt) +
      0.5 * (X[n + 1] - X[n] - drift * dt) *
      1/(SigmaSigma * dt) * (X[n + 1] - X[n] - drift * dt) #* log(pois)
  }
  print(1)
  return(nll)
}


EM_nll(initial_vals, df_pro[,2], dt)

# y-axis limits for plotting
ylim <- c(-6, 6)

initial_vals <- c(2,2,0.2,0.2,0.2,0.2)

# initial parameter study
"
result_loop <- matrix(0, nrow=2, ncol=20)

for (i in range(0.00001,0.001,0.00001)){
  print(i)
  initial_vals <- c(i,i)
  result_EM_1 <- optim(par = initial_vals, fn = EM_nll,
                       X = df_pro[,2], dt = 0.05)
  result_loop[,i+11] <- result_EM_1$par
}
"

# Euler murayama scheme for optimising values, (0=convergence)
result_EM_1 <- optim(par = initial_vals, fn = EM_nll,
                     X = df_pro[,2], dt = 0.05)

result_EM_1
result <- c(result_EM_1$par,mle$par[5:7])
result <- mle$par
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
  dZ <- rnorm(n,mean = params[6], sd = params[7])
  a <- 0
  for (i in 3:(n + 1)) {
    drift_term <- f_drift(X[i-1],params)
    X[i] <- X[i-1] + drift_term * dt + f_diff(params) * dW[i-1] * sqrt(dt)
       + f_jump(X[i], params) * dB[i] * dZ[i]
    a <- a + f_jump(X[i], params) * dB[i] * dZ[i]
  }
  print(a)
  return(data.frame(t = t, X = X))
}

sim_res <- simulate_EM(T_start, T_end, dt, result, x0)

sim_res

plot(sim_res$t, sim_res$X, type = 'l', col = 'red')
lines(df_pro$age, df_pro$d18O, type = 'l', col=1, ylim=ylim)


plot(df_pro$age, df_pro$d18O, type = 'l', col=1, ylim=ylim)
lines(sim_res$t, sim_res$X, type = 'l', col = 'red', ylim=ylim)


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

