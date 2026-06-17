functions {
  real pk2(real x, real c1, real c2, real mu1, real mu2, real sigma1, real sigma2) {
    return c1 * exp(-square((x - mu1) / sigma1)) + c2 * exp(-square((x - mu2) / sigma2));
  }  
}
data {
  int<lower=0> N;
  int<lower=0> A;
  vector<lower=0>[N] B;
  vector<lower=0>[N] E;
  array[N] int<lower=1, upper=A> age;
}
parameters {
  real log_c1;
  real log_c2;
  real<lower=0> mu1;
  real<lower=0> d_mu2;
  real log_sigma1;
  real log_sigma2;
}
transformed parameters {
  vector[A] fertility_schedule;
  for(i in 1:A) {
    fertility_schedule[i] = pk2(i, exp(log_c1), exp(log_c2), mu1, mu1 + d_mu2, exp(log_sigma1), exp(log_sigma2));
  }
}
model {
  log_c1 ~ normal(-3.0, 1.0);
  log_c2 ~ normal(-2, 1.0);
  
  mu1 ~ normal(20, 2);
  d_mu2 ~ normal(20, 10);
  
  log_sigma1 ~ normal(1, 0.5);
  log_sigma2 ~ normal(1, 0.5);
  
  // Poisson likelihood
  vector[N] log_rate = log(E) + log(fertility_schedule[age]);
  target += sum(B .* log_rate - exp(log_rate));
}

