#' Fit the Peristera-Kostaki Model
#'
#' @param data data frame
#' @param births column name of number of births
#' @param exposure column name of exposure
#' @param age column name of age
#' @param effects model specification, specify via \code{pk_effects()}
#'
#' @export
pk_fit <- function(
  data,
  births,
  exposure,
  age,
  effects = pk_effects(),
  sample_prior = c("no", "yes", "only"),
  ...
) {
  sample_prior <- match.arg(sample_prior)

  checkmate::assert_data_frame(data)
  checkmate::assert_string(births)
  checkmate::assert_string(exposure)
  checkmate::assert_string(age)
  checkmate::assert_choice(births, choices = names(data))
  checkmate::assert_choice(exposure, choices = names(data))
  checkmate::assert_choice(age, choices = names(data))
  checkmate::assert_class(effects, "pk_effects")

  processed <- process_effects(effects, data)
  stan_data <- build_stan_data(processed, data, births, exposure, age)

  if (Sys.getenv("BAYESFERTILITY_DEV") == "TRUE") {
    model <- cmdstanr::cmdstan_model("src/stan/pk2.stan")
  } else {
    model <- instantiate::stan_package_model(
      name = "pk2",
      package = "bayesfertility"
    )
  }

  make_init_fun <- function(stan_data) {
    expand_centers <- function(p) {
      centers <- numeric(p$n_coef)
      for (i in seq_along(p$terms_info)) {
        centers[
          p$terms_info[[i]]$col_idx
        ] <- p$prior_centers[i]
      }
      centers
    }

    function() {
      list(
        beta_c1 = array(expand_centers(processed$c1), dim = stan_data$P_c1),
        beta_c2 = array(expand_centers(processed$c2), dim = stan_data$P_c2),
        beta_mu1 = array(expand_centers(processed$mu1), dim = stan_data$P_mu1),
        beta_mu2 = array(expand_centers(processed$mu2), dim = stan_data$P_mu2),
        beta_sigma1 = array(
          expand_centers(processed$sigma1),
          dim = stan_data$P_sigma1
        ),
        beta_sigma2 = array(
          expand_centers(processed$sigma2),
          dim = stan_data$P_sigma2
        ),
        tau_c1 = rep(0.3, stan_data$n_penalties_c1),
        tau_c2 = rep(0.3, stan_data$n_penalties_c2),
        tau_mu1 = rep(0.3, stan_data$n_penalties_mu1),
        tau_mu2 = rep(0.3, stan_data$n_penalties_mu2),
        tau_sigma1 = rep(0.3, stan_data$n_penalties_sigma1),
        tau_sigma2 = rep(0.3, stan_data$n_penalties_sigma2)
      )
    }
  }

  stan_data$likelihood_on <- 1L
  fit <- NULL
  beta_draws <- NULL
  tau_draws <- NULL
  if (sample_prior != "only") {
    fit <- model$sample(stan_data, init = make_init_fun(stan_data), ...)

    draws <- extract_all_draws(fit, processed)
    beta_draws <- draws$beta
    tau_draws <- draws$tau
  }

  prior_fit <- NULL
  prior_beta_draws <- NULL
  prior_tau_draws <- NULL
  if (sample_prior == "yes" || sample_prior == "only") {
    stan_data_prior <- stan_data
    stan_data_prior$likelihood_on <- 0L
    prior_fit <- model$sample(
      stan_data_prior,
      init = make_init_fun(stan_data_prior),
      ...
    )

    prior_draws <- extract_all_draws(prior_fit, processed)
    prior_beta_draws <- prior_draws$beta
    prior_tau_draws <- prior_draws$tau
  }

  structure(
    list(
      fit = fit,
      prior_fit = prior_fit,
      processed = processed,
      beta_draws = beta_draws,
      tau_draws = tau_draws,
      prior_beta_draws = prior_beta_draws,
      prior_tau_draws = prior_tau_draws,
      sample_prior = sample_prior,
      data = data,
      effects = effects,
      age_col = age,
      births_col = births,
      exposure_col = exposure,
      diagnostics = if (!is.null(fit)) {
        pk_quick_diagnostics(fit)
      } else {
        NULL
      }
    ),
    class = "pk_fit"
  )
}

#' Extract posterior draws of beta_<param> from a cmdstanr fit
#' @noRd
extract_beta_draws <- function(fit, param_name, n_coef) {
  var_names <- paste0("beta_", param_name, "[", seq_len(n_coef), "]")
  draws <- fit$draws(variables = var_names, format = "draws_matrix")
  unclass(draws)
}

#' Extract posterior draws of tau_<param> from a cmdstanr fit
#' @noRd
extract_tau_draws <- function(fit, param_name, n_penalties) {
  if (n_penalties == 0) {
    n_draws <- nrow(fit$draws(variables = "lp__", format = "draws_matrix"))
    return(matrix(numeric(0), nrow = n_draws, ncol = 0))
  }
  var_names <- paste0("tau_", param_name, "[", seq_len(n_penalties), "]")
  draws <- fit$draws(variables = var_names, format = "draws_matrix")
  unclass(draws)
}

extract_all_draws <- function(fit, processed) {
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")
  beta <- setNames(
    lapply(param_names, function(p) {
      extract_beta_draws(fit, p, processed[[p]]$n_coef)
    }),
    param_names
  )
  tau <- setNames(
    lapply(param_names, function(p) {
      extract_tau_draws(fit, p, length(processed[[p]]$penalties))
    }),
    param_names
  )
  list(beta = beta, tau = tau)
}

pk_quick_diagnostics <- function(cmdstanr_fit) {
  summ <- cmdstanr_fit$summary(variables = NULL)

  list(
    max_rhat = max(summ$rhat, na.rm = TRUE),
    min_ess = min(summ$ess_bulk, na.rm = TRUE),
    n_divergent = sum(cmdstanr_fit$diagnostic_summary()$num_divergent)
  )
}
