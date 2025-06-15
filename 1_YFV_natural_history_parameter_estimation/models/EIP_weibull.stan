data {
  int N;
  int day[N];
  int infected[N];
  int died[N];
  
  // priors for the weibull distribution
  real a_1;
  real a_2;
  real b_1;
  real b_2;
  real min_p_death_prior;
  real max_p_death_prior;
}

parameters {
  real<lower=a_1,upper=a_2> a; // weibull dist parameter
  real<lower=b_1,upper=b_2> b; // weibull dist parameter
  real<lower=0.01,upper=1> p_death_min; // baseline mortality parameter
  real<lower=0.01,upper=1> p_death_max; // baseline mortality parameter
}

model {
  p_death_min ~ normal(p_death_prior, 0.02);
  a ~ uniform(a_1, a_2);
  b ~ uniform(b_1, b_2);
  for (i in 1:N) {
    real temp_day = day[i];
    if (temp_day == 0) {
      temp_day = temp_day + 0.01;
    }
    real mortality_p = p_death_min + (1 - p_death_min) * weibull_cdf(temp_day, a, b);
    died[i] ~ binomial(infected[i], mortality_p);
  }
}

generated quantities {
  real days_simulated = weibull_rng(a, b);
}
