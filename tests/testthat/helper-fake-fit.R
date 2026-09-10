# Shared test fixtures. Loaded automatically by testthat before tests run.

#' Small synthetic dataset for tests (age x group counts of births/exposure)
make_synthetic_data <- function(n_age = 12, groups = c("A", "B")) {
  ages <- seq(15, 49, length.out = n_age)
  df <- expand.grid(age = ages, group = groups, stringsAsFactors = FALSE)
  df$group <- factor(df$group, levels = groups)
  df$exposure <- 1000
  set.seed(42)
  df$births <- rpois(nrow(df), lambda = 50)
  df
}

#' Hand-build a pk_fit-classed object with synthetic draws, bypassing Stan.
#'
#' Runs the real process_effects()/build_stan_data() pipeline (pure R), but
#' fills in beta_draws/tau_draws with small synthetic matrices instead of
#' calling cmdstanr, so S3 methods can be exercised quickly and
#' deterministically.
make_fake_pk_fit <- function(
  effects = pk_effects(),
  data = make_synthetic_data(),
  n_draws = 25,
  with_prior = FALSE
) {
  set.seed(123)
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")
  processed <- process_effects(effects, data)

  make_draws <- function() {
    beta <- stats::setNames(
      lapply(param_names, function(p) {
        n_coef <- processed[[p]]$n_coef
        matrix(rnorm(n_draws * n_coef, sd = 0.2), nrow = n_draws, ncol = n_coef)
      }),
      param_names
    )
    tau <- stats::setNames(
      lapply(param_names, function(p) {
        n_pen <- length(processed[[p]]$penalties)
        matrix(abs(rnorm(n_draws * n_pen, sd = 0.2)), nrow = n_draws, ncol = n_pen)
      }),
      param_names
    )
    list(beta = beta, tau = tau)
  }

  posterior <- make_draws()
  prior <- if (with_prior) make_draws() else list(beta = NULL, tau = NULL)

  fake_cmdstan_fit <- list(
    metadata = function() {
      list(num_chains = 1L, iter_sampling = n_draws, iter_warmup = n_draws)
    },
    time = function() list(total = 1.0)
  )

  structure(
    list(
      fit = fake_cmdstan_fit,
      prior_fit = if (with_prior) list() else NULL,
      processed = processed,
      beta_draws = posterior$beta,
      tau_draws = posterior$tau,
      prior_beta_draws = prior$beta,
      prior_tau_draws = prior$tau,
      sample_prior = if (with_prior) "yes" else "no",
      data = data,
      effects = effects,
      age_col = "age",
      births_col = "births",
      exposure_col = "exposure",
      diagnostics = list(max_rhat = 1.0, min_ess = 1000, n_divergent = 0)
    ),
    class = "pk_fit"
  )
}
