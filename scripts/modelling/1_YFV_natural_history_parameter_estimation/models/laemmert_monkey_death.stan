data {
  int N;
  real days[N];
  int truncated[N];
  real days_truncated;
  
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
  for(i in 1:N) {
    if(truncated[i] == 0)
      days[i] ~ gamma(a, b);
    else
      target += gamma_lccdf(days_truncated| a, b);
  }
  a ~ gamma(a_1, a_2);
  b ~ gamma(b_1, b_2);
}

generated quantities {
  
  // no is_max because max is truncated
  real is_min;
  real is_truncated;
  real days_simulated = gamma_rng(a, b);
  {
    vector[N] days_sim;
    for(i in 1:N)
      days_sim[i] = gamma_rng(a, b);
    is_min = min(days_sim) <= min(days);
    is_truncated = max(days_sim) >= days_truncated;
  }
}
