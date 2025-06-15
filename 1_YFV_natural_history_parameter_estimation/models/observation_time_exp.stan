data {
  int N;
  real days[N];
  
  // priors
  real a_1;
}

parameters {
  real<lower=0> a;
}

model {
  days ~ exponential(a);
  a ~ uniform(0.1, 10);
}

generated quantities {
  real is_min;
  real is_max;
  real days_simulated = exponential_rng(a);
  {
    vector[N] days_sim;
    for(i in 1:N)
      days_sim[i] = exponential_rng(a);
    is_min = min(days_sim) <= min(days);
    is_max = max(days_sim) >= max(days);
  }
}
