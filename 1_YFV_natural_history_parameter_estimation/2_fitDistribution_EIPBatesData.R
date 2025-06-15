## Loading required libraries
library(tidyverse); library(rstan); library(loo)

## Recreating data from "The Development of the Virus of Yellow Fever in Haemagogus Mosquitoes" by Bates and Roca-Garcia 1946
degrees25_df_adults <- tibble(days = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20),
                              total= c(27,30,24,21,28,28,28,30,29,30,23, 24, 28, 30, 12, 17 ),    
                              died = c(17,26,0, 3, 0, 2, 2, 0, 17,12, 5, 16, 11, 20, 10, 11),
                              mod =  c(0, 0, 0, 3, 0, 2, 2, 0, 17,12, 5, 16, 11, 20, 10, 11))
degrees25_df_adults$perc <- degrees25_df_adults$mod/degrees25_df_adults$total

## Adults at 25 degrees
model_gamma2 <- stan_model("2_YFV_natural_history_parameter_estimation_real/models/EIP_gamma_model2.stan")
data_stan_adults25 <- list(N = length(degrees25_df_adults$days),
                           day = degrees25_df_adults$days,
                           infected = degrees25_df_adults$total,
                           died = degrees25_df_adults$mod,
                           a_1 = 0.1,
                           a_2 = 10,
                           b_1 = 0.1,
                           b_2 = 10,
                           min_p_death_prior_mean = 0.02,
                           min_p_death_prior_sd = 0.02,
                           max_p_death_prior_mean = 0.8 - 0,
                           max_p_death_prior_sd = 0.1)
fit_adults25 <- sampling(model_gamma2, data = data_stan_adults25, iter = 2000, chains = 1)
summary(fit_adults25)
# pairs(fit_adults25)

posterior_samples <- rstan::extract(fit_adults25, "died_rep")
died_rep_array <- posterior_samples$died_rep
df_posterior_died_summary <- data.frame(
  day = degrees25_df_adults$days,
  observed = degrees25_df_adults$mod,
  predicted_mean = apply(died_rep_array, 2, mean),
  predicted_lower = apply(died_rep_array, 2, quantile, probs = 0.025),
  predicted_upper = apply(died_rep_array, 2, quantile, probs = 0.975))
died_posterior_plot <- ggplot(df_posterior_died_summary, aes(x = day)) +
  geom_point(aes(y = observed), color = "red", size = 2) +
  geom_line(aes(y = predicted_mean), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = predicted_lower, ymax = predicted_upper), 
              alpha = 0.2, fill = "blue") +
  labs(title = "Observed vs. Posterior Predictive Death Counts",
       y = "Deaths",
       x = "Days") +
  theme_minimal()

posterior_samples <- rstan::extract(fit_adults25, "mortality")
mortality_array <- posterior_samples$mortality
df_posterior_mortality_summary <- data.frame(
  day = degrees25_df_adults$days,
  observed = degrees25_df_adults$perc,
  predicted_mean = apply(mortality_array, 2, mean),
  predicted_lower = apply(mortality_array, 2, quantile, probs = 0.025),
  predicted_upper = apply(mortality_array, 2, quantile, probs = 0.975))
mortality_posterior_plot <- ggplot(df_posterior_mortality_summary, aes(x = day)) +
  geom_point(aes(y = observed), color = "black", size = 2) +
  geom_line(aes(y = predicted_mean), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = predicted_lower, ymax = predicted_upper), 
              alpha = 0.2, fill = "blue") +
  labs(title = "Observed vs. Posterior Predictive Death Proportion",
       y = "Proportion Dying",x = "Days") +
  theme_bw()

cowplot::plot_grid(mortality_posterior_plot, 
                   died_posterior_plot,
                   labels = c("A", "B"))

modelled_a <- mean(rstan::extract(fit_adults25, "a")[[1]])
modelled_b <- mean(rstan::extract(fit_adults25, "b")[[1]])
modelled_a / modelled_b

saveRDS(object = data.frame(gamma_a = modelled_a, gamma_b = modelled_b),
        file = "outputs/EIP_adultMice_gammaParams_25degrees.rds")

a <- ggplot(data = df_posterior_mortality_summary) +
  geom_point(x = 0, y = 0.63, shape = 21, color = "black", fill = "grey", size = 4, alpha = 0.1) +
  geom_point(x = 1, y = 0.87, shape = 21, color = "black", fill = "grey", size = 4, alpha = 0.1) + 
  geom_line(aes(x = day, y = predicted_mean), color = "#DB3069", size = 1) +
  geom_ribbon(aes(x = day, ymin = predicted_lower, ymax = predicted_upper), 
              alpha = 0.2, fill = "#DB3069") +
  geom_point(aes(x = day, y = observed), shape = 21, color = "black", fill = "#CDC776", size = 4) +
  lims(y = c(0, 1)) +
  labs(x = "Days Since Mosquito Infected", y = "Proportion of Mosquito-Bitten\nMice Dying") +
  theme_bw(base_family = "sans") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

df_gamma <- lapply(seq_along(1:1000), function(i) {
  x_vals <- seq(0, 20, by = 0.5)
  data.frame(i = i, x = x_vals,
    pdfval = dgamma(x_vals, shape = rstan::extract(fit_adults25, "a")[[1]][i], 
                            rate = rstan::extract(fit_adults25, "b")[[1]][i]))}) %>%
  bind_rows()

b <- ggplot(df_gamma, aes(x = x, y = pdfval, group = i)) +
  geom_line(alpha = 0.05, col = "#DB3069", size = 0.2) +
  labs(x = "Days Since Mosquito Infected", y = "Proportion of Mosquitoes\nBecoming Infectious") +
  theme_bw(base_family = "sans") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

overall <- cowplot::plot_grid(a, b, labels = c("a", "b"), nrow = 2)
ggsave(filename = "1_YFV_natural_history_parameter_estimation/figures/SI_EIP_figure.pdf",
       plot = overall,
       width = 5, height = 7)

overall <- cowplot::plot_grid(a, b, labels = c("a", "b"), nrow = 1)
ggsave(filename = "1_YFV_natural_history_parameter_estimation/figures/SI_EIP_figure_alt.pdf",
       plot = overall,
       width = 10, height = 3.5)


