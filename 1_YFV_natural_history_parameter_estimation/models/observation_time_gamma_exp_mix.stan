data {
  int N;
  real days[N];
  
  // priors
  real a_1;
  real a_2;
  real b_1;
  real b_2;
}

parameters {
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;
  real<lower=0> a_exp;
  real<lower=0,upper=1> prob;
}

model {
  a_gamma ~ gamma(a_1, a_2);
  b_gamma ~ gamma(b_1, b_2);
  a_exp ~ uniform(0.1, 10);
  prob ~ uniform(0, 1);
  for (i in 1:N) {
    target += log_mix(prob,
                      exponential_lpdf(days[i] | a_exp),
                      gamma_lpdf(days[i] | a_gamma, b_gamma));
  }
}

generated quantities {
  vector[N] days_simulated;
  real is_min;
  real is_max;
  
  // Generate mixture samples
  for (i in 1:N) {
    if (bernoulli_rng(prob) == 1)
      days_simulated[i] = exponential_rng(a_exp);
    else
      days_simulated[i] = gamma_rng(a_gamma, b_gamma);
  }
  
  // Checks: compare the min and max of simulated days to observed days
  is_min = min(days_simulated) <= min(days);
  is_max = max(days_simulated) >= max(days);
}
