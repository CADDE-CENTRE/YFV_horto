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
  real<lower=0> a;
  real<lower=0> b;
}

model {
  days ~ gamma(a, b);
  a ~ gamma(a_1, a_2);
  b ~ gamma(b_1, b_2);
}

generated quantities {
  real is_min;
  real is_max;
  real days_simulated = gamma_rng(a, b);
  {
    vector[N] days_sim;
    for(i in 1:N)
      days_sim[i] = gamma_rng(a, b);
    is_min = min(days_sim) <= min(days);
    is_max = max(days_sim) >= max(days);
  }
}
