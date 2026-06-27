functions {
  real pk2(real x, real c1, real c2, real mu1, real mu2, real sigma1, real sigma2) {
    return c1 * exp(-square((x - mu1) / sigma1)) 
      + c2 * exp(-square((x - mu2) / sigma2));
  }  
}
data {
  int<lower=0, upper=1> likelihood_on;
  int<lower=0> N;
  int<lower=0> A;
  vector<lower=0>[N] B;
  vector<lower=0>[N] E;
  array[N] int<lower=1, upper=A> age;
  vector[A] ages;
  
  // For each schedule parameter k in {c1, c2, mu1, mu2, sigma1, sigma2}:
  //    P_k:                number of coefficients (>= 1)
  //    X_k:                N x P_k design matrix
  //    n_penalties_k:      number of penalties (>= 0)
  //    penalty_sizes_k:    size of each penalty block
  //    penalty_starts_X_k: starting column in X_k for each penalty's coefficients
  //    penalty_starts_S_k: starting row/col in S_packed_k for each penalty
  //    S_packed_k:         block-diagonal stack of penalty matrices
  //    n_terms_k:          number of logical terms (intercept + parametric + smooth)
  //    term_starts_k:      starting column in X_k for each term
  //    term_sizes_k:       number of columns each term spans
  //    term_scales_k:      per-term proir SD (0 means "skip direct prior, penalty only")
  //    term_centers_k:     per-term prior mean
  //    tau_scales_k:       per-penalty tau prior scales
  
  int<lower=1> P_c1;
  matrix[N, P_c1] X_c1;
  int<lower=0> n_penalties_c1;
  array[n_penalties_c1] int<lower=1> penalty_sizes_c1;
  array[n_penalties_c1] int<lower=1> penalty_starts_X_c1;
  array[n_penalties_c1] int<lower=1> penalty_starts_S_c1;
  matrix[n_penalties_c1 > 0 ? sum(penalty_sizes_c1) : 0,
         n_penalties_c1 > 0 ? sum(penalty_sizes_c1) : 0] S_packed_c1;
  int<lower=0> n_terms_c1;
  array[n_terms_c1] int<lower=1> term_starts_c1;
  array[n_terms_c1] int<lower=1> term_sizes_c1;
  vector<lower=0>[n_terms_c1] term_scales_c1;
  vector[n_terms_c1] term_centers_c1;
  vector<lower=0>[n_penalties_c1] tau_scales_c1;
         
  int<lower=1> P_c2;
  matrix[N, P_c2] X_c2;
  int<lower=0> n_penalties_c2;
  array[n_penalties_c2] int<lower=1> penalty_sizes_c2;
  array[n_penalties_c2] int<lower=1> penalty_starts_X_c2;
  array[n_penalties_c2] int<lower=1> penalty_starts_S_c2;
  matrix[n_penalties_c2 > 0 ? sum(penalty_sizes_c2) : 0,
         n_penalties_c2 > 0 ? sum(penalty_sizes_c2) : 0] S_packed_c2;
  int<lower=0> n_terms_c2;
  array[n_terms_c2] int<lower=1> term_starts_c2;
  array[n_terms_c2] int<lower=1> term_sizes_c2;
  vector<lower=0>[n_terms_c2] term_scales_c2;
  vector[n_terms_c2] term_centers_c2;
  vector<lower=0>[n_penalties_c2] tau_scales_c2;
         
  int<lower=1> P_mu1;
  matrix[N, P_mu1] X_mu1;
  int<lower=0> n_penalties_mu1;
  array[n_penalties_mu1] int<lower=1> penalty_sizes_mu1;
  array[n_penalties_mu1] int<lower=1> penalty_starts_X_mu1;
  array[n_penalties_mu1] int<lower=1> penalty_starts_S_mu1;
  matrix[n_penalties_mu1 > 0 ? sum(penalty_sizes_mu1) : 0,
         n_penalties_mu1 > 0 ? sum(penalty_sizes_mu1) : 0] S_packed_mu1;
  int<lower=0> n_terms_mu1;
  array[n_terms_mu1] int<lower=1> term_starts_mu1;
  array[n_terms_mu1] int<lower=1> term_sizes_mu1;
  vector<lower=0>[n_terms_mu1] term_scales_mu1;
  vector[n_terms_mu1] term_centers_mu1;
  vector<lower=0>[n_penalties_mu1] tau_scales_mu1;
         
  int<lower=1> P_mu2;
  matrix[N, P_mu2] X_mu2;
  int<lower=0> n_penalties_mu2;
  array[n_penalties_mu2] int<lower=1> penalty_sizes_mu2;
  array[n_penalties_mu2] int<lower=1> penalty_starts_X_mu2;
  array[n_penalties_mu2] int<lower=1> penalty_starts_S_mu2;
  matrix[n_penalties_mu2 > 0 ? sum(penalty_sizes_mu2) : 0,
         n_penalties_mu2 > 0 ? sum(penalty_sizes_mu2) : 0] S_packed_mu2;
  int<lower=0> n_terms_mu2;
  array[n_terms_mu2] int<lower=1> term_starts_mu2;
  array[n_terms_mu2] int<lower=1> term_sizes_mu2;
  vector<lower=0>[n_terms_mu2] term_scales_mu2;
  vector[n_terms_mu2] term_centers_mu2;      
  vector<lower=0>[n_penalties_mu2] tau_scales_mu2;
  
  int<lower=1> P_sigma1;
  matrix[N, P_sigma1] X_sigma1;
  int<lower=0> n_penalties_sigma1;
  array[n_penalties_sigma1] int<lower=1> penalty_sizes_sigma1;
  array[n_penalties_sigma1] int<lower=1> penalty_starts_X_sigma1;
  array[n_penalties_sigma1] int<lower=1> penalty_starts_S_sigma1;
  matrix[n_penalties_sigma1 > 0 ? sum(penalty_sizes_sigma1) : 0,
         n_penalties_sigma1 > 0 ? sum(penalty_sizes_sigma1) : 0] S_packed_sigma1;
  int<lower=0> n_terms_sigma1;
  array[n_terms_sigma1] int<lower=1> term_starts_sigma1;
  array[n_terms_sigma1] int<lower=1> term_sizes_sigma1;
  vector<lower=0>[n_terms_sigma1] term_scales_sigma1;
  vector[n_terms_sigma1] term_centers_sigma1;             
  vector<lower=0>[n_penalties_sigma1] tau_scales_sigma1;
         
  int<lower=1> P_sigma2;
  matrix[N, P_sigma2] X_sigma2;
  int<lower=0> n_penalties_sigma2;
  array[n_penalties_sigma2] int<lower=1> penalty_sizes_sigma2;
  array[n_penalties_sigma2] int<lower=1> penalty_starts_X_sigma2;
  array[n_penalties_sigma2] int<lower=1> penalty_starts_S_sigma2;
  matrix[n_penalties_sigma2 > 0 ? sum(penalty_sizes_sigma2) : 0,
         n_penalties_sigma2 > 0 ? sum(penalty_sizes_sigma2) : 0] S_packed_sigma2;
  int<lower=0> n_terms_sigma2;
  array[n_terms_sigma2] int<lower=1> term_starts_sigma2;
  array[n_terms_sigma2] int<lower=1> term_sizes_sigma2;
  vector<lower=0>[n_terms_sigma2] term_scales_sigma2;
  vector[n_terms_sigma2] term_centers_sigma2;      
  vector<lower=0>[n_penalties_sigma2] tau_scales_sigma2;
}
parameters {
  vector[P_c1]     beta_c1;
  vector[P_c2]     beta_c2;
  vector[P_mu1]    beta_mu1;
  vector[P_mu2]    beta_mu2;
  vector[P_sigma1] beta_sigma1;
  vector[P_sigma2] beta_sigma2;
  
  // Smoothness hyperparameters: one per penalty
  vector<lower=0>[n_penalties_c1]     tau_c1;
  vector<lower=0>[n_penalties_c2]     tau_c2;
  vector<lower=0>[n_penalties_mu1]    tau_mu1;
  vector<lower=0>[n_penalties_mu2]    tau_mu2;
  vector<lower=0>[n_penalties_sigma1] tau_sigma1;
  vector<lower=0>[n_penalties_sigma2] tau_sigma2;
}
transformed parameters {
  // Cell-level schedule parameters, on the working scale
  vector[N] eta_c1     = X_c1     * beta_c1;
  vector[N] eta_c2     = X_c2     * beta_c2;
  vector[N] eta_mu1    = X_mu1    * beta_mu1;
  vector[N] eta_d_mu2  = X_mu2    * beta_mu2;
  vector[N] eta_sigma1 = X_sigma1 * beta_sigma1;
  vector[N] eta_sigma2 = X_sigma2 * beta_sigma2;
  
  vector[N] log_rate;
  vector[N] schedule;
  for(i in 1:N) {
    schedule[i] = pk2(
      ages[age[i]],
      exp(eta_c1[i]),
      exp(eta_c2[i]),
      eta_mu1[i],
      eta_mu1[i] + exp(eta_d_mu2[i]),
      exp(eta_sigma1[i]),
      exp(eta_sigma2[i])
    );
    log_rate[i] = log(schedule[i]) + log(E[i]);
  }
}
model {
  
  // Per-term priors
  // For each term: if scale > 0, apply a normal prior to its coefficient block.
  // If scale == 0, the term gets no direct prior (its smooth's penalty contribution below
  // handles it).
  for(t in 1:n_terms_c1) {
    if(term_scales_c1[t] > 0) {
      segment(beta_c1, term_starts_c1[t], term_sizes_c1[t]) ~ normal(term_centers_c1[t], term_scales_c1[t]);
    }
  }
  for(t in 1:n_terms_c2) {
    if(term_scales_c2[t] > 0) {
      segment(beta_c2, term_starts_c2[t], term_sizes_c2[t]) ~ normal(term_centers_c2[t], term_scales_c2[t]);
    }
  }
  for(t in 1:n_terms_mu1) {
    if(term_scales_mu1[t] > 0) {
      segment(beta_mu1, term_starts_mu1[t], term_sizes_mu1[t]) ~ normal(term_centers_mu1[t], term_scales_mu1[t]);
    }
  }
  for(t in 1:n_terms_mu2) {
    if(term_scales_mu2[t] > 0) {
      segment(beta_mu2, term_starts_mu2[t], term_sizes_mu2[t]) ~ normal(term_centers_mu2[t], term_scales_mu2[t]);
    }
  }
  for(t in 1:n_terms_sigma1) {
    if(term_scales_sigma1[t] > 0) {
      segment(beta_sigma1, term_starts_sigma1[t], term_sizes_sigma1[t]) ~ normal(term_centers_sigma1[t], term_scales_sigma1[t]);
    }
  }
  for(t in 1:n_terms_sigma2) {
    if(term_scales_sigma2[t] > 0) {
      segment(beta_sigma2, term_starts_sigma2[t], term_sizes_sigma2[t]) ~ normal(term_centers_sigma2[t], term_scales_sigma2[t]);
    }
  }
  
  
  // Penalty contributions
  for(k in 1:n_penalties_c1) {
    int sz = penalty_sizes_c1[k];
    int sX = penalty_starts_X_c1[k];
    int sS = penalty_starts_S_c1[k];
    target += -0.5 / square(tau_c1[k]) * quad_form(block(S_packed_c1, sS, sS, sz, sz), segment(beta_c1, sX, sz));
  }
  
  for(k in 1:n_penalties_c2) {
    int sz = penalty_sizes_c2[k];
    int sX = penalty_starts_X_c2[k];
    int sS = penalty_starts_S_c2[k];
    target += -0.5 / square(tau_c2[k]) * quad_form(block(S_packed_c2, sS, sS, sz, sz), segment(beta_c2, sX, sz));
  }
  
  for(k in 1:n_penalties_mu1) {
    int sz = penalty_sizes_mu1[k];
    int sX = penalty_starts_X_mu1[k];
    int sS = penalty_starts_S_mu1[k];
    target += -0.5 / square(tau_mu1[k]) * quad_form(block(S_packed_mu1, sS, sS, sz, sz), segment(beta_mu1, sX, sz));
  }
  
  for(k in 1:n_penalties_mu2) {
    int sz = penalty_sizes_mu2[k];
    int sX = penalty_starts_X_mu2[k];
    int sS = penalty_starts_S_mu2[k];
    target += -0.5 / square(tau_mu2[k]) * quad_form(block(S_packed_mu2, sS, sS, sz, sz), segment(beta_mu2, sX, sz));
  }
  
  for(k in 1:n_penalties_sigma1) {
    int sz = penalty_sizes_sigma1[k];
    int sX = penalty_starts_X_sigma1[k];
    int sS = penalty_starts_S_sigma1[k];
    target += -0.5 / square(tau_sigma1[k]) * quad_form(block(S_packed_sigma1, sS, sS, sz, sz), segment(beta_sigma1, sX, sz));
  }
  
  for(k in 1:n_penalties_sigma2) {
    int sz = penalty_sizes_sigma2[k];
    int sX = penalty_starts_X_sigma2[k];
    int sS = penalty_starts_S_sigma2[k];
    target += -0.5 / square(tau_sigma2[k]) * quad_form(block(S_packed_sigma2, sS, sS, sz, sz), segment(beta_sigma2, sX, sz));
  }
  
  // Hyperpriors
  for(k in 1:n_penalties_c1)     tau_c1[k]     ~ student_t(3, 0, tau_scales_c1[k]);
  for(k in 1:n_penalties_c2)     tau_c2[k]     ~ student_t(3, 0, tau_scales_c2[k]);
  for(k in 1:n_penalties_mu1)    tau_mu1[k]    ~ student_t(3, 0, tau_scales_mu1[k]);
  for(k in 1:n_penalties_mu2)    tau_mu2[k]    ~ student_t(3, 0, tau_scales_mu2[k]);
  for(k in 1:n_penalties_sigma1) tau_sigma1[k] ~ student_t(3, 0, tau_scales_sigma1[k]);
  for(k in 1:n_penalties_sigma2) tau_sigma2[k] ~ student_t(3, 0, tau_scales_sigma2[k]);
  
  if(likelihood_on == 1) {
    target += sum(B .* log_rate - exp(log_rate));
  }
}
