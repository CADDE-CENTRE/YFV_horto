# load required libraries
library(individual); library(dplyr); library(tidyverse); library(doParallel); library(tictoc); 
library(parallel); library(profvis); library(truncnorm)

## Sourcing functions
source("functions/IBM_model.R")
source("functions/particle_filter.R")
source("2_YFV_prior_parameterisation/construct_weibullPrior_epidemicTiming.R")
source("2_YFV_prior_parameterisation/construct_truncNormalPrior_R0.R")

# Loading in fitted parameters
latent_period_fit <- readRDS("outputs/exposure_infectiousDist_stanFit.rds")
latent_period_gamma_shape <- mean(rstan::extract(latent_period_fit, "a")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
latent_period_gamma_rate <- mean(rstan::extract(latent_period_fit, "b")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
infectious_period_fit <- readRDS("outputs/infectious_deathDist_stanFit.rds")
infectious_period_gamma_shape <- mean(rstan::extract(infectious_period_fit, "a")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
infectious_period_gamma_rate <- mean(rstan::extract(infectious_period_fit, "b")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
EIP_gamma_fit <- readRDS("outputs/EIP_adultMice_gammaParams_25degrees.rds")
EIP_gamma_shape <- EIP_gamma_fit$gamma_a
EIP_gamma_rate <- EIP_gamma_fit$gamma_b
death_observation_fit <- readRDS("outputs/deathObservation_distMix_stanFit.rds")

death_observation_gamma_shape <- mean(rstan::extract(death_observation_fit, "a_gamma")[[1]]) 
death_observation_gamma_rate <- mean(rstan::extract(death_observation_fit, "b_gamma")[[1]])
death_observation_exp_rate <- mean(rstan::extract(death_observation_fit, "a_exp")[[1]])
death_observation_prob <- mean(rstan::extract(death_observation_fit, "prob")[[1]])

generation_time <- round(
  (latent_period_gamma_shape / latent_period_gamma_rate) + 
  (infectious_period_gamma_shape / infectious_period_gamma_rate) +
  (EIP_gamma_shape / EIP_gamma_rate) +
  ((1 - death_observation_prob) * (death_observation_gamma_shape / death_observation_gamma_rate) +
     death_observation_prob * (1 / death_observation_exp_rate)))
    
# Loading in and processing Horto/PEAL data for model fitting
horto_df <- readRDS("data/processed_HortoData.rds") %>%
  filter(!is.na(zone_peal)) %>%
  filter(final_yfv_result != "negative")
epi_curve <- incidence::incidence(horto_df$date_collection)
plot(epi_curve)

# Generating incidence data and cutting off first 4 infections
start_date <- as.Date("2017-10-23")
horto_df_fitting <- horto_df %>%
  filter(date_collection > start_date) %>%
  group_by(date_collection) %>%
  summarise(count = n()) %>%
  complete(date_collection = seq.Date(start_date, 
                                      as.Date("2018-01-08"), 
                                      by = "days"),
           fill = list(count = 0))
plot(horto_df_fitting$date_collection, horto_df_fitting$count)
horto_df_fitting$time <- 1:nrow(horto_df_fitting)

## Particle filtering

## Model parameters
N <- 86 - 3 - 3 # (3 negative monkeys - assume the rest not in database killed by yellow fever) and other 3 are predeceased from where we're starting
N_obs <- sum(horto_df_fitting$count)
death_obs_prop <- N_obs / N
dt <- 0.2
initial_infections <- 1
gamma <- 1 / (infectious_period_gamma_shape / infectious_period_gamma_rate)
importations <- 4 / (33 / N_obs) # from the genomic data - 4 clusters containing NHP sequences, from 33/61 NHPs sequenced
importation_last_date <- max(horto_df_fitting$date_collection) - generation_time # upper bound assumed to be 1 generation time before the final monkey death

## Parameters for initial particle filtering to identify parameter regime of highest likelihood
R0_scan <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
start_date_scan <- start_date + seq(0, 36, 3)
transmission_type_scan <- "density_dependent" # c("frequency_dependent", "density_dependent")
exponential_noise_scan <- 1/1e-1
likelihood_distribution_scan <- "negative_binomial" # c("poisson", "negative_binomial") # "negative_binomial"
negative_binomial_size <- 1
R0_prior_function <- function(R0_value) { return(log(dtruncnorm(R0_value, a = 1, b = 16, mean = fitted_mu, sd = fitted_sigma))) }
mrca_prior_scan <- c("early", "late")

iterations <- 10
particles <- 200
cores <- 10

loglikelihood_matrix <- array(data = NA, dim = c(iterations, length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan)))
epilikelihood_matrix <- array(data = NA, dim = c(iterations, length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan)))
importlikelihood_matrix <- array(data = NA, dim = c(iterations, length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan)))
startdatelikelihood_matrix <- array(data = NA, dim = c(iterations, length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan)))
importations_matrix <- array(data = NA, dim = c(iterations, length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan)))
final_size_matrix <- array(data = NA, dim = c(iterations, length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan)))
output_matrix <- array(data = NA, dim = c(iterations, length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan), length(horto_df_fitting$count)))
output_matrix_inc <- array(data = NA, dim = c(iterations, length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan), length(horto_df_fitting$count)))

overall_seed <- 10
set.seed(overall_seed)
simulation_seeds <- array(data = rnbinom(n = iterations * length(R0_scan) * length(start_date_scan) * length(transmission_type_scan) * length(exponential_noise_scan) * length(likelihood_distribution_scan) * length(mrca_prior_scan), 
                                         mu = 10^6, size = 1), 
                          dim = c(length(R0_scan), length(start_date_scan), length(transmission_type_scan), length(exponential_noise_scan), length(likelihood_distribution_scan), length(mrca_prior_scan), iterations))
length(R0_scan) * length(start_date_scan) * length(transmission_type_scan) * length(exponential_noise_scan) * length(likelihood_distribution_scan) * length(mrca_prior_scan)
fresh_run <- TRUE
if (fresh_run) {
  
  ## Looping through R0
  for (i in 1:length(R0_scan)) {
    
    ## Looping through the start dates
    for (j in 1:length(start_date_scan)) {
      
      for (k in 1:length(transmission_type_scan)) {
        
        for (l in 1:length(exponential_noise_scan)) {
          
          for (m in 1:length(likelihood_distribution_scan)) {
            
            for (n in 1:length(mrca_prior_scan)) {
              
              new_sys_time <- Sys.time()
              
              # Selecting the start date and filtering the Horto data to start then
              start_date <- start_date_scan[j]
              data <- horto_df_fitting %>%
                filter(date_collection >= start_date) %>%
                rename(daily_incidence = count)
              steps <- nrow(data) / dt
              days <- nrow(data)
              
              # Calculating the importation rate - note that given the epidemic goes to extinction, we really
              # need to calculate it up to 1 full generation time before the last death (i.e. the time when that last monkey was infected 
              # in our model) - final date in importation_date_range is the actual time when all the susceptibles are depleted
              importation_date_range <- data %>%
                filter(date_collection <= importation_last_date)
              importation_rate <- importations / nrow(importation_date_range) 
              
              # Defining the misc list that supports running the particle filter
              misc <- list(seed = simulation_seeds[i, j, k, l, m, n, ], 
                           steps = steps, 
                           gamma = gamma,
                           particles = particles,
                           dt = dt, 
                           N = N, 
                           start_date = start_date,
                           importation_rate = importation_rate,
                           empirical_importations = importations,
                           transmission_type = transmission_type_scan[k], 
                           exponential_noise_rate = exponential_noise_scan[l],
                           likelihood = c("epidemiological", "importations", "start_date"),
                           initial_infections = initial_infections, 
                           death_obs_prop = death_obs_prop, 
                           initial_run = TRUE, 
                           overall_run_length = steps,
                           latent_period_gamma_shape = latent_period_gamma_shape, 
                           EIP_gamma_shape = EIP_gamma_shape,
                           EIP_gamma_rate = EIP_gamma_rate, 
                           latent_period_gamma_rate = latent_period_gamma_rate,
                           infectious_period_gamma_shape = infectious_period_gamma_shape, 
                           infectious_period_gamma_rate = infectious_period_gamma_rate,
                           death_observation_mixture_gamma_shape = death_observation_gamma_shape, 
                           death_observation_mixture_gamma_rate = death_observation_gamma_rate,
                           death_observation_mixture_exponential_rate = death_observation_exp_rate,
                           death_observation_mixture_prob = death_observation_prob,
                           prior_function = R0_prior_function,
                           weibull_shape = weibull_shape_vector[n],
                           weibull_scale = weibull_scale_vector[n],
                           start_date_weibull_fitting = start_date_weibull_fitting,
                           negative_binomial_size = negative_binomial_size,
                           likelihood_distribution = likelihood_distribution_scan[m])
              
              # Setting up the cluster to run everything in parallel
              cl <- makeCluster(cores)
              clusterExport(cl, varlist = c("r_loglike", "weight_particles_poisson", "weight_particles_negative_binomial", "data", "misc", "run_simulation2", "fitted_mu", "fitted_sigma"))
              clusterEvalQ(cl, {
                library(individual)
                library(truncnorm)
              })
              
              # Running the loglikelihood function in parallel
              R0_temp <- c("R0" = R0_scan[i])
              clusterExport(cl, varlist = c("R0_temp"))
              result_parallel <- parLapply(cl, 1:iterations, function(i) {
                misc_new <- misc
                misc_new$seed <- misc$seed[i]
                temp <- r_loglike(R0_temp, data, misc_new)
                return(temp)
              })
              parallel::stopCluster(cl)
              
              # Storing the output
              padding_zeroes <- rep(0, as.numeric(start_date_scan[j] - start_date_scan[1]))
              for (o in 1:iterations) {
                output_matrix[o, i, j, k, l, m, n, ] <- c(padding_zeroes, result_parallel[[o]]$deaths_trajectory)
                output_matrix_inc[o, i, j, k, l, m, n, ] <- c(padding_zeroes, result_parallel[[o]]$inc_trajectory)
                final_size_matrix[o, i, j, k, l, m, n] <- sum(result_parallel[[o]]$deaths_trajectory)
                loglikelihood_matrix[o, i, j, k, l, m, n] <- result_parallel[[o]]$loglikelihood
                epilikelihood_matrix[o, i, j, k, l, m, n] <- result_parallel[[o]]$likelihood_components$epi
                importlikelihood_matrix[o, i, j, k, l, m, n] <- result_parallel[[o]]$likelihood_components$importations
                startdatelikelihood_matrix[o, i, j, k, l, m, n] <- result_parallel[[o]]$likelihood_components$start
                importations_matrix[o, i, j, k, l, m, n] <- result_parallel[[o]]$importations
              }
              
              print(c(paste0("i = ", i, ", j = ", j, ", k = ", k, ", l = " , l, ", m = ", m, ", n = ", n), paste0("Took ", round((Sys.time() - new_sys_time), 2), " minutes")))

            }
          }
        }
      }
    }
  }
  
  saveRDS(list(output_deaths = output_matrix,
               output_incidence = output_matrix_inc,
               final_size = final_size_matrix, 
               loglike = loglikelihood_matrix, epilikelihood_matrix = epilikelihood_matrix,
               importlikelihood_matrix = importlikelihood_matrix, startdatelikelihood_matrix = startdatelikelihood_matrix, 
               importations = importations_matrix),
          "outputs/parameterScan_hortoEstimation_YesImportations_YesStartDate.rds")
  
} else {
  
  temp <- readRDS("outputs/parameterScan_hortoEstimation_YesImportations_YesStartDate.rds")
  loglikelihood_matrix <- temp$loglike
  epilikelihood_matrix <- temp$epilikelihood_matrix
  importlikelihood_matrix <- temp$importlikelihood_matrix
  startdatelikelihood_matrix <- temp$startdatelikelihood_matrix
  importations <- temp$importations
  output_matrix <- temp$output
}

## Overall loglikelihood
loglik_avg <- apply(loglikelihood_matrix, c(2, 3, 4, 5, 6, 7), mean)
dimnames(loglik_avg) <- list(
  R0 = R0_scan,                   # use your R0_scan vector here
  start_date = paste0("s", as.Date(start_date_scan)),
  transmission_type = transmission_type_scan,
  exponential_noise = exponential_noise_scan,
  likelihood_dist = likelihood_distribution_scan,
  mrca = mrca_prior_scan
)
df_long <- as.data.frame.table(loglik_avg, responseName = "loglikelihood") %>%
  mutate(start_date = as.Date(gsub("s", "", start_date))) %>%
  filter(exponential_noise == "10")
head(df_long)
scales <- c(-90, -70)
overall_likelihood_plot <- ggplot(subset(df_long, likelihood_dist == "negative_binomial" & mrca == "late"), 
                                  aes(x = start_date, y = factor(R0), fill = loglikelihood)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", oob = scales::squish, limits = scales) + 
  labs(x = "Start Date",
       y = expression(R[0]),
       fill = "Avg.\nLoglike") +
  # facet_grid(. ~ mrca) +
  scale_x_date(expand = c(0, 0)) +  # Remove whitespace on the x-axis
  scale_y_discrete(expand = c(0, 0)) +  # Remove whitespace on the y-axis
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

## Sampling parameter combinations
df_long <- data.frame(R0 = R0_scan, loglik_avg) %>%
  pivot_longer(cols = -R0, names_to = "StartDateRaw", values_to = "Value") %>%
  mutate(
    StartDate = as.Date(substr(StartDateRaw, 2, 11), "%Y.%m.%d"), #2, 11 when two distibutions being assessed
    Distribution = case_when(str_detect(StartDateRaw, "poisson") ~ "poisson", 
                             str_detect(StartDateRaw, "negative_binomial") ~ "negative_binomial", 
                             TRUE ~ NA_character_),
    mrca = case_when(str_ends(StartDateRaw, "early") ~ "early", 
                     str_ends(StartDateRaw, "late") ~ "late", 
                     TRUE ~ NA_character_)) %>% # FINISH THIS
  mutate(LogLikelihood_adj = Value - max(Value), 
         Likelihood = exp(LogLikelihood_adj), 
         Probability = Likelihood / sum(Likelihood)) %>%
  select(R0, StartDate, mrca, Distribution, Value, LogLikelihood_adj, Likelihood, Probability)

set.seed(123)
samples <- 10000
mrca <- "early"
sampled_indices <- sample((1:nrow(df_long))[df_long$mrca == mrca], 
                          size = 10000,
                          replace = TRUE,
                          prob = (df_long$Probability)[df_long$mrca == mrca])
sampled_data <- df_long[sampled_indices, c("R0", "StartDate", "Value")] 

## Marginal for R0
avg_R0_values <- sampled_data %>%
  group_by(R0) %>%
  summarise(AvgValue = mean(Value))
sampled_data_R0 <- sampled_data %>%
  left_join(avg_R0_values, by = "R0")
sampled_data_R0$R0 <- as.factor(sampled_data_R0$R0)
R0_marginal_plot <- ggplot(sampled_data_R0, aes(x = R0, fill = AvgValue)) +
  geom_histogram(color = "black", stat = "count") +
  scale_fill_distiller(palette = "RdBu") + # limits = scales, oob = scales::squish) + 
  labs(x = "Marginal R0 Distribution",
       y = "Frequency",
       fill = "Avg.\nLoglike") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

likelihood_plot_legend <- cowplot::get_legend(overall_likelihood_plot)
likelihood_R0_marginal <- cowplot::plot_grid(overall_likelihood_plot + theme(legend.position = "none"),
                                             R0_marginal_plot + theme(legend.position = "none"),
                                             likelihood_plot_legend,
                                             nrow = 1, labels = c("a", "b"), rel_widths = c(1, 1, 0.2))
## Plotting the inferred deaths trajectories
samples <- 10000
sampled_output_matrix <- matrix(nrow = samples, ncol = length(horto_df_fitting$count))
sampled_output_inc_matrix <- matrix(nrow = samples, ncol = length(horto_df_fitting$count))
sampled_Reff_matrix <- matrix(nrow = samples, ncol = length(horto_df_fitting$count))
for (i in 1:samples) {
  R0 <- unlist(sampled_data[i, "R0"])
  implied_beta_sim <- R0 * gamma / N
  R0_index <- which(rownames(loglik_avg) == R0) # paste0("R0=", R0))
  start_date <- sampled_data[i, "StartDate"]
  start_date <- start_date$StartDate
  start_date_index <- which(substr(colnames(loglik_avg), 2, 11) == start_date)
  k <- sample(1:iterations, 1)
  sampled_output_matrix[i, ] <- output_matrix[k, R0_index, start_date_index, 1, 1, which(likelihood_distribution_scan == "negative_binomial"), which(mrca_prior_scan == mrca), ]
  sampled_output_inc_matrix[i, ] <- output_matrix_inc[k, R0_index, start_date_index, 1, 1, which(likelihood_distribution_scan == "negative_binomial"), which(mrca_prior_scan == mrca), ]
  max_num_obs_deaths <- max(cumsum(output_matrix[k, R0_index, start_date_index, 1, 1, which(likelihood_distribution_scan == "negative_binomial"), which(mrca_prior_scan == mrca), ]))
  mod_factor <- N / max_num_obs_deaths
  N_over_time <- N - (cumsum(output_matrix[k, R0_index, start_date_index, 1, 1, which(likelihood_distribution_scan == "negative_binomial"), which(mrca_prior_scan == mrca), ]) * mod_factor) ## THIS IS A TEMPORARY FIX # need to output the inferred deaths, not just the observed ones 
  sampled_Reff_matrix[i, ] <- N_over_time * implied_beta_sim / gamma
}

## Deaths time-series
lower <- apply(sampled_output_matrix, 2, quantile, 0.025)
upper <- apply(sampled_output_matrix, 2, quantile, 0.975)
min <- apply(sampled_output_matrix, 2, min)
max <- apply(sampled_output_matrix, 2, max)
mean <- apply(sampled_output_matrix, 2, mean)
median <- apply(sampled_output_matrix, 2, median)
output_df <- data.frame(time = horto_df_fitting$date_collection, 
                        observed = horto_df_fitting$count, 
                        mean = mean, 
                        median = median,
                        lower = lower, 
                        upper = upper,
                        min = min,
                        max = max)

lower_Reff <- apply(sampled_Reff_matrix, 2, quantile, 0.025)
upper_Reff <- apply(sampled_Reff_matrix, 2, quantile, 0.975)
min_Reff <- apply(sampled_Reff_matrix, 2, min)
max_Reff <- apply(sampled_Reff_matrix, 2, max)
mean_Reff <- apply(sampled_Reff_matrix, 2, mean)
output_df2 <- data.frame(time = horto_df_fitting$date_collection, 
                         observed = horto_df_fitting$count, 
                         mean = mean_Reff, 
                         lower = lower_Reff, 
                         upper = upper_Reff,
                         min = min_Reff,
                         max = max_Reff)

Reff_deaths_plot <- ggplot() +
  geom_ribbon(data = output_df2, aes(x = time, ymin = lower, ymax = upper), fill = "#BF6900", alpha = 0.2) +
  geom_line(data = output_df2, aes(x = time, y = mean), color = "#BF6900", size = 0.75) +
  geom_point(data = output_df2, aes(x = time, y = observed), color = "black", size = 2) +
  geom_ribbon(data = output_df, aes(x = time, ymin = min, ymax = max), fill = "#145C9E", alpha = 0.2) +
  geom_line(data = output_df, aes(x = time, y = mean), color = "#145C9E", size = 0.75) +
  labs(x = "", y = "Daily Reported\nNHP Deaths") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14)) +
  scale_x_date(date_breaks = "1 week",
               limits = c(as.Date("2017-11-06"), as.Date("2018-01-08"))) +
  theme_bw() +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

scale  <- 10 / 12
alt_Reff_deaths_plot <- ggplot() +
  geom_ribbon(data = output_df2, aes(x = time, ymin = lower  * scale, ymax = upper  * scale),
              fill = "#BF6900", alpha = 0.2) +
  geom_line(data  = output_df2, aes(x = time, y = mean  * scale), colour = "#BF6900", size = 0.75) +
  geom_point(data = output_df2, aes(x = time, y = observed), colour = "black", size = 2) +
  geom_ribbon(data = output_df, aes(x = time, ymin = min, ymax = max), fill = "#145C9E", alpha = 0.2) +
  geom_line(data = output_df, aes(x = time, y = mean), colour = "#145C9E", size = 0.75) +
  scale_y_continuous(name = "Daily Reported\nNHP Deaths", 
                     limits = c(0, 10), 
                     breaks = seq(0, 10, 2), 
                     sec.axis = sec_axis(~ . / scale, 
                                         name   = "Reff", 
                                         breaks = seq(0, 14, 2))) +
  scale_x_date(date_breaks = "1 week", limits = c(as.Date("2017-11-06"), as.Date("2018-01-08"))) +
  labs(x = NULL) +
  theme_bw() +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


overall_R0_modelling_plot <- cowplot::plot_grid(likelihood_R0_marginal, alt_Reff_deaths_plot, 
                                                nrow = 2, labels = c("", "c"))

ggsave(plot = overall_R0_modelling_plot,
       filename = "3_horto_YFV_R0_estimation/figures/overall_R0_modelling_plot.pdf",
       height = 7,
       width = 8)
