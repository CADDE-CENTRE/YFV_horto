## Loading required libraries
library(tidyverse); library(rstan)

# Recreating table 1 from "The susceptibility of Howler monkeys to yellow fever virus" Laemmert & Klumm, 1950.

## Interpret (e.g.) 6th day after infection as being 5 days post-infection; assume monkey 4 dies on day 10
days_death_post_infection <- c(5, 8, 2, 10, 6, 4, 3, 4, 3, 9)
## Interpret - "trace" doesn't count; if detectable on 1st day, then 0.5; count 6th day after infection as being 5 days post-infection
days_virus_post_infection <- c(2, 2, 1, 2, 2, 1, 1, 1, 0.5, 2)

df <- tibble(id = seq_along(days_virus_post_infection),
             days_virus_post_infection = days_virus_post_infection,
             days_death_post_infection = days_death_post_infection,
             days_death_post_virus = days_death_post_infection - days_virus_post_infection)  %>% 
  mutate(death_truncated = if_else(is.na(days_death_post_infection), 1, 0)) %>% 
  mutate(days_death_post_infection = if_else(is.na(days_death_post_infection), -99, days_death_post_infection))

## Fitting the time between infection and death
model <- stan_model("2_YFV_natural_history_parameter_estimation_real/models/laemmert_monkey_death.stan")
data_stan <- list(N = nrow(df),
                  days = df$days_death_post_infection,
                  truncated = df$death_truncated,
                  days_truncated = 10,
                  a_1 = 1,
                  a_2 = 0.5,
                  b_1 = 1,
                  b_2 = 2)
fit <- sampling(model, data=data_stan, iter=2000, chains=4)
hist(rstan::extract(fit, "days_simulated")[[1]], breaks = 50)
summary(fit)
saveRDS(fit, "outputs/infection_deathDist_stanFit.rds")

df_gamma <- lapply(seq_along(1:1000), function(i) {
  x_vals <- seq(0, 15, by = 0.5)
  data.frame(i = i, x = x_vals,
             pdfval = dgamma(x_vals, shape = rstan::extract(fit, "a")[[1]][i], 
                             rate = rstan::extract(fit, "b")[[1]][i]))}) %>%
  bind_rows()
df_real <- tibble(days_death_post_infection = df$days_death_post_infection[!is.na(days_death_post_infection)], 
                  iteration = 1:10,
                  type = "actual")  %>%
  count(days_death_post_infection) %>% # Count occurrences of each day
  mutate(normalized_count = n / sum(n)) # Normalize by total counts

a <- ggplot() +
  geom_bar(data = df_real,
           aes(x = days_death_post_infection, y = normalized_count, fill = "empirical"), 
           stat = "identity", alpha = 1, width = 0.7, fill = "#948D9B") +
  geom_line(data = df_gamma,
            aes(x = x, y = pdfval, group = i), alpha = 0.05, 
            size = 1, adjust = 1, col = "#80C498") +
  scale_x_continuous(breaks = seq(0, 15, 1)) +
  scale_y_continuous(name = "Density / Normalized Count") +
  xlab("Days from infection to death") +
  ylab("Density (Simulated) or Normalized Count (Empirical)") +
  theme_bw(base_family = "sans") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")

## Fitting the time between infection and death
data_stan <- list(N = nrow(df),
                  days = df$days_virus_post_infection,
                  truncated = df$death_truncated,
                  days_truncated = 10,
                  a_1 = 1,
                  a_2 = 0.5,
                  b_1 = 1,
                  b_2 = 2)
fit <- sampling(model, data=data_stan, iter=2000, chains=4)
hist(rstan::extract(fit, "days_simulated")[[1]], breaks = 50)
summary(fit)
saveRDS(fit, "outputs/exposure_infectiousDist_stanFit.rds")

df_gamma <- lapply(seq_along(1:1000), function(i) {
  x_vals <- seq(0, 15, by = 0.5)
  data.frame(i = i, x = x_vals,
             pdfval = dgamma(x_vals, shape = rstan::extract(fit, "a")[[1]][i], 
                             rate = rstan::extract(fit, "b")[[1]][i]))}) %>%
  bind_rows()
df_real <- tibble(days_virus_post_infection = df$days_virus_post_infection[!is.na(df$days_virus_post_infection)], 
                  iteration = 1:10,
                  type = "actual") %>%
  count(days_virus_post_infection) %>% # Count occurrences of each day
  mutate(normalized_count = n / sum(n)) # Normalize by total counts

b <- ggplot() +
  geom_bar(data = df_real,
           aes(x = days_virus_post_infection, y = normalized_count, fill = "empirical"), 
           stat = "identity", alpha = 1, width = 0.7, fill = "#948D9B") +
  geom_line(data = df_gamma,
            aes(x = x, y = pdfval, group = i), alpha = 0.05, 
            size = 1, adjust = 1, col = "#D6BD94") +
  scale_x_continuous(breaks = seq(0, 15, 1)) +
  scale_y_continuous(name = "Density / Normalized Count") +
  xlab("Days from infection to detectable virus") +
  ylab("Density (Simulated) or Normalized Count (Empirical)") +
  theme_bw(base_family = "sans") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")

## Fitting the time between infectiousness (virus detectable) and death
data_stan <- list(N = nrow(df),
                  days = df$days_death_post_virus,
                  truncated = df$death_truncated,
                  days_truncated = 10,
                  a_1 = 1,
                  a_2 = 0.5,
                  b_1 = 1,
                  b_2 = 2)
fit <- sampling(model, data=data_stan, iter=2000, chains=4)
hist(rstan::extract(fit, "days_simulated")[[1]], breaks = 50)
summary(fit)
saveRDS(fit, "outputs/infectious_deathDist_stanFit.rds")

df_gamma <- lapply(seq_along(1:1000), function(i) {
  x_vals <- seq(0, 15, by = 0.5)
  data.frame(i = i, x = x_vals,
             pdfval = dgamma(x_vals, shape = rstan::extract(fit, "a")[[1]][i], 
                             rate = rstan::extract(fit, "b")[[1]][i]))}) %>%
  bind_rows()

df_real <- tibble(days_death_post_virus = df$days_death_post_virus, 
                  iteration = 1:10,
                  type = "actual") %>%
  count(days_death_post_virus) %>% # Count occurrences of each day
  mutate(normalized_count = n / sum(n)) # Normalize by total counts

c <- ggplot() +
  geom_bar(data = df_real,
           aes(x = days_death_post_virus, y = normalized_count, fill = "empirical"), 
           stat = "identity", alpha = 1, width = 0.7, fill = "#948D9B") +
  geom_line(data = df_gamma,
            aes(x = x, y = pdfval, group = i), alpha = 0.05, 
            size = 1, adjust = 1, col = "#D67663") +
  scale_x_continuous(breaks = seq(0, 15, 1)) +
  scale_y_continuous(name = "Density / Normalized Count") +
  xlab("Days from detectable virus to death") +
  ylab("Density (Simulated) or Normalized Count (Empirical)") +
  theme_bw(base_family = "sans") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")

first_part <- cowplot::plot_grid(a + theme(legend.position = "none"), 
                                 b + theme(legend.position = "none"), 
                                 c, 
                                 labels = c("a", "b", "c"),
                                 nrow = 1)

ggsave(filename = "2_YFV_natural_history_parameter_estimation_real/figures/SI_NHP_NatHist.pdf",
       plot = first_part,
       width = 12, height = 4)


