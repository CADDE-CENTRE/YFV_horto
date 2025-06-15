data {
  int N;
  int day[N];
  int infected[N];
  int died[N];
  
  // priors for the gamma distribution
  real a_1;
  real a_2;
  real b_1;
  real b_2;
  real min_p_death_prior_mean;
  real max_p_death_prior_mean;
  real min_p_death_prior_sd;
  real max_p_death_prior_sd;
}

parameters {
  real<lower=0> a; // gamma dist parameter
  real<lower=0> b; // gamma dist parameter
  real<lower=0,upper=1> p_death_min; // baseline mortality parameter
  real<lower=p_death_min,upper=1> p_death_max; // baseline mortality parameter
}

model {
  p_death_min ~ normal(min_p_death_prior_mean, min_p_death_prior_sd);
  p_death_max ~ normal(max_p_death_prior_mean, max_p_death_prior_sd);
  a ~ uniform(a_1, a_2);
  b ~ uniform(b_1, b_2);
  for (i in 1:N) {
    real temp_day = day[i];
    real mortality_p = p_death_min + (p_death_max - p_death_min) * gamma_cdf(temp_day, a, b);
    died[i] ~ binomial(infected[i], mortality_p);
  }
}

generated quantities {
  real days_simulated = gamma_rng(a, b);
}
