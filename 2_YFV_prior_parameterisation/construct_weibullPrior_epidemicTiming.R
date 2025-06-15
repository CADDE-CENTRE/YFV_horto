fit_weibull_3quantiles <- function(q025_target = 14, 
                                   q50_target = 29, 
                                   q975_target = 36) {
  
  # Objective: sum of squared differences from target quantiles
  obj_fun <- function(log_par) {
    # Exponentiate to ensure positivity
    shape <- exp(log_par[1])
    scale <- exp(log_par[2])
    
    # Calculate quantiles
    q025  <- qweibull(0.025, shape = shape, scale = scale)
    q50   <- qweibull(0.50,  shape = shape, scale = scale)
    q975  <- qweibull(0.975, shape = shape, scale = scale)
    
    # Sum of squared errors
    sse <- (q025 - q025_target)^2 + 
      (q50  - q50_target)^2  + 
      (q975 - q975_target)^2
    
    return(sse)
  }
  
  init_guess <- c(log(1.5), log(20))  # a rough guess
  fit <- optim(par = init_guess, fn = obj_fun, method = "BFGS")
  
  shape_hat <- exp(fit$par[1])
  scale_hat <- exp(fit$par[2])
  
  # Inspect final quantiles
  q025_hat  <- qweibull(0.025, shape = shape_hat, scale = scale_hat)
  q50_hat   <- qweibull(0.50,  shape = shape_hat, scale = scale_hat)
  q975_hat  <- qweibull(0.975, shape = shape_hat, scale = scale_hat)
  
  return(list(
    shape    = shape_hat,
    scale    = scale_hat,
    sse      = fit$value,
    q025_fit = q025_hat,
    q50_fit  = q50_hat,
    q975_fit = q975_hat,
    converged = fit$convergence == 0
  ))
}

start_date_weibull_fitting <- as.Date("2017-10-09")

## Later Importation Assumptions
lower_late <- as.Date("2017-10-30") - 5
median_late <- as.Date("2017-11-19") - 5
upper_late <-  as.Date("2017-11-30") - 5

weibull_fit_late <- fit_weibull_3quantiles(q025_target = as.numeric(lower_late - start_date_weibull_fitting),
                                            q50_target = as.numeric(median_late - start_date_weibull_fitting),  
                                            q975_target = as.numeric(upper_late - start_date_weibull_fitting))
weibull_shape_late <- weibull_fit_late$shape
weibull_scale_late <- weibull_fit_late$scale

## Earlier Importation Assumptions
lower_early <- as.Date("2017-10-18") - 5
median_early <- as.Date("2017-11-14") - 5
upper_early <-  as.Date("2017-11-27") - 5

weibull_fit_early <- fit_weibull_3quantiles(q025_target = as.numeric(lower_early - start_date_weibull_fitting),
                                           q50_target = as.numeric(median_early - start_date_weibull_fitting),  
                                           q975_target = as.numeric(upper_early - start_date_weibull_fitting))
weibull_shape_early <- weibull_fit_early$shape
weibull_scale_early <- weibull_fit_early$scale

weibull_shape_vector <- c(weibull_shape_early, weibull_shape_late)
weibull_scale_vector <- c(weibull_scale_early, weibull_scale_late)


x <- rweibull(10000, shape = weibull_shape_early, scale = weibull_scale_early)
hist(x, breaks = 20)
