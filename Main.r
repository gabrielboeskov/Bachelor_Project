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
ggplot(data = df_pro, aes(x = age, y = d18O)) +
  geom_line() +
  xlab("Time before present (yrs)") +
  ylab(expression(deltaˆ18 ~ O ~ (permil)))
ggplot(data = df, aes(x = age, y = d18O)) +
  geom_line() +
  xlab("Time before present (yrs)") +
  ylab(expression(deltaˆ18 ~ O ~ (permil)))
ggplot(data = df_in, aes(x = age, y = d18O)) +
  geom_line() +
  xlab("Time before present (yrs)") +
  ylab(expression(deltaˆ18 ~ O ~ (permil)))


#Drift function
f_drift <- function(x,params){
  drift <- (params[1]+params[2]*x-params[3]*x^3)
  return(drift*x)
}

f_diff <- function(params){
  diff <- params[4]
  return(diff)
}

dt <- 0.05
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
    nll <- nll + 0.5 * log(SigmaSigma * dt) +
      0.5 * (X[n + 1] - X[n] - drift * dt) *
      1/(SigmaSigma * dt) * (X[n + 1] - X[n] - drift * dt)
  }
  return(nll)
}

initial_vals <- c(5,5,5,5)

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

result <- result_EM_1$par

result

T_start <- df_pro[length(df_pro[,1]),1]
T_end <- df_pro[1,1]
x0=0

simulate_EM <- function(T_start, T_end, dt, params, x0) {
  t <- seq(T_start, T_end, by = dt)
  n<-length(t)-1
  dW <- rnorm(n)  
  X <- numeric(n + 1)
  X[1] <- x0

  
  for (i in 3:(n + 1)) {
    drift_term <- f_drift(X[i-1],params)
    X[i] <- X[i-1] + drift_term * dt + f_diff(params) * dW[i-1] *sqrt(dt)
  }
  return(data.frame(t = t, X = X))
}

sim_res <- simulate_EM(T_start, T_end, dt, result, x0)

plot(sim_res$t, sim_res$X, type = 'l', col = 1)

lines(df_pro$age, df_pro$d18O, type = 'l', col=3)

plot(df_pro$age, df_pro$d18O, type = 'l', col=1)


calculate_transition_density <- function(x0, T_start, T_end, dt, params) {

  n_steps <- round((T_end - T_start) / dt)
  
  density <- numeric(n_steps + 1)
  density[1] <- x0
  
  for (i in 1:n_steps) {
    density[i + 1] <- rnorm(1, mean = density[i] + f_drift(density[i], params) * dt, f_diff(params)*sqrt(dt)) 
  }
  return(density)
}

plot(density(df_pro[,2]))

transition_densities <- calculate_transition_density(x0,T_start,T_end,dt,result)
plot(density(transition_densities))
plot(transition_densities,type='l')

lines(transition_densities, type='l', col=3)

qqnorm(transition_densities)


