# Load required libraries 
library(truncnorm)

# Negative log-likelihood function for the truncated normal distribution
neg_log_lik <- function(params, data, a, b) {
  mu <- params[1]
  sigma <- params[2]
  if (sigma <= 0) return(Inf)
  -sum(log(dtruncnorm(data, a = a, b = b, mean = mu, sd = sigma)))
}

# R0 estimates from Journal of Travel Medicine article: https://academic.oup.com/jtm/article/27/7/taaa156/5901887
data_vector <- c(6, 4.8, 5.2, 7.1, 11, 4.1, 2.38, 3.59, 3.23, 4.21, 1.35)

# Specify truncation bounds
a <- 0
b <- 16

# Initial estimates for mu and sigma based on JTM data
init_mu <- mean(data_vector)
init_sigma <- sd(data_vector)
init_params <- c(init_mu, init_sigma)

# Use optim to minimize the negative log-likelihood
fit <- optim(
  par = init_params,
  fn = neg_log_lik,
  data = data_vector,
  a = a,
  b = b,
  method = "L-BFGS-B",
  lower = c(-Inf, 1e-6)
)

# Extract the fitted parameters
fitted_mu <- fit$par[1]
fitted_sigma <- fit$par[2]

hist(rtruncnorm(10000, a = a, b = b, mean = fitted_mu, sd = fitted_sigma), breaks = 40)
