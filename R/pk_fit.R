#' Fit the Peristera-Kostaki Model
#'
#' @description
#' The Peristera-Kostaki (PK) model is a parametric model of
#' age-specific fertility rates. The model has a two-component parametric schedule
#' in which each "hump" is governed by an
#' amplitude, a location (age), and a width:
#' \deqn{ASFR(age) = c_1 \exp(-((age - \mu_1)/\sigma_1)^2) +
#'   c_2 \exp(-((age - \mu_2)/\sigma_2)^2)}
#' Observed births are modeled as Poisson counts with rate
#' `exposure * ASFR(age)`. Each of the six schedule parameters (`c1`, `c2`,
#' `mu1`, `mu2`, `sigma1`, `sigma2`) can depend on covariates via a formula
#' supplied through \code{pk_effects()}, including smooth/random-effect terms,
#' which are fit with partial pooling via Stan.
#'
#' @details
#' Model fitting is delegated to a pre-compiled Stan model (`pk2.stan`) via
#' \pkg{cmdstanr}. When `sample_prior` requests it, the model is fit twice:
#' once with the likelihood enabled (the posterior fit) and once with the
#' likelihood switched off (the prior fit, useful for prior predictive checks
#' and prior/posterior comparisons via \code{pk_prior_posterior()}).
#'
#' @param data A data frame containing the columns named by `births`,
#'   `exposure`, and `age`, plus any covariates referenced in `effects`.
#' @param births String giving the name of the column in `data` holding
#'   observed birth counts (the Poisson response).
#' @param exposure String giving the name of the column in `data` holding
#'   exposure (person-time at risk) for each row.
#' @param age String giving the name of the column in `data` holding age.
#' @param effects Model specification, created via \code{pk_effects()},
#'   giving a formula for each of the six PK schedule parameters (`c1`, `c2`,
#'   `mu1`, `mu2`, `sigma1`, `sigma2`) plus any prior customizations. Defaults
#'   to intercept-only formulas for all six parameters.
#' @param sample_prior One of `"no"` (default; only fit the posterior),
#'   `"yes"` (fit both the posterior and, separately, the prior with the
#'   likelihood disabled), or `"only"` (fit only the prior).
#' @param ... Additional arguments passed on to cmdstanr's
#'   `CmdStanModel$sample()` method (e.g. `chains`, `iter_warmup`,
#'   `iter_sampling`, `parallel_chains`, `seed`), used for both the posterior
#'   and, when requested, the prior fit.
#'
#' @return An object of class `pk_fit`, a list containing:
#'   \describe{
#'     \item{fit}{The posterior `CmdStanMCMC` fit object, or `NULL` if
#'       `sample_prior = "only"`.}
#'     \item{prior_fit}{The prior `CmdStanMCMC` fit object (likelihood
#'       disabled), or `NULL` if `sample_prior = "no"`.}
#'     \item{processed}{Per-schedule-parameter design matrices, penalty
#'       matrices, and resolved priors produced by `process_effects()`.}
#'     \item{beta_draws, tau_draws}{Named lists (one entry per schedule
#'       parameter) of posterior draws for the coefficients and smoothing
#'       hyperparameters, or `NULL` if no posterior fit was run.}
#'     \item{prior_beta_draws, prior_tau_draws}{The corresponding prior
#'       draws, or `NULL` if no prior fit was run.}
#'     \item{sample_prior}{The `sample_prior` argument as supplied.}
#'     \item{data}{The original `data` argument, used as the default
#'       `newdata` in downstream methods such as `predict.pk_fit()`.}
#'     \item{effects}{The `pk_effects` object as supplied.}
#'     \item{age_col, births_col, exposure_col}{The `age`, `births`, and
#'       `exposure` column-name strings as supplied.}
#'     \item{diagnostics}{A list with `max_rhat`, `min_ess`, and
#'       `n_divergent` summarizing convergence of the posterior fit, or
#'       `NULL` if no posterior fit was run.}
#'   }
#'
#' @examples
#' \dontrun{
#' # `data` is typically pre-aggregated: one row per age (and, if `effects`
#' # includes covariates, per age/covariate combination), giving weighted
#' # birth counts and exposure, e.g.:
#' #   data |>
#' #     dplyr::group_by(age, race) |>
#' #     dplyr::summarize(
#' #       births = sum(weight * fert),
#' #       exposure = sum(weight)
#' #     )
#'
#' # Fit a single national schedule (intercept-only for every parameter)
#' national_fit <- pk_fit(
#'   data,
#'   "births",
#'   "exposure",
#'   "age",
#'   effects = pk_effects(),
#'   sample_prior = "yes",
#'   chains = 4,
#'   parallel_chains = 4,
#'   iter_warmup = 500,
#'   iter_sampling = 500,
#'   adapt_delta = 0.999,
#'   max_treedepth = 14
#' )
#'
#' # Let every schedule parameter vary by a covariate such as race
#' race_fit <- pk_fit(
#'   data,
#'   "births",
#'   "exposure",
#'   "age",
#'   effects = pk_effects(
#'     c1 = ~1 + race,
#'     c2 = ~1 + race,
#'     mu1 = ~1 + race,
#'     mu2 = ~1 + race,
#'     sigma1 = ~1 + race,
#'     sigma2 = ~1 + race
#'   ),
#'   sample_prior = "yes",
#'   parallel_chains = 4
#' )
#'
#' # Downstream: predicted schedules/parameters and total fertility rate
#' predict(national_fit, type = "schedule_curve")
#' posterior_tfr(national_fit)
#' }
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

    z_size <- function(n_pen, sizes) {
      if (n_pen == 0) 0 else sum(sizes)
    }

    function() {
      list(
        beta_c1_raw = array(
          expand_centers(processed$c1),
          dim = stan_data$P_c1
        ),
        beta_c2_raw = array(
          expand_centers(processed$c2),
          dim = stan_data$P_c2
        ),
        beta_mu1_raw = array(
          expand_centers(processed$mu1),
          dim = stan_data$P_mu1
        ),
        beta_mu2_raw = array(
          expand_centers(processed$mu2),
          dim = stan_data$P_mu2
        ),
        beta_sigma1_raw = array(
          expand_centers(processed$sigma1),
          dim = stan_data$P_sigma1
        ),
        beta_sigma2_raw = array(
          expand_centers(processed$sigma2),
          dim = stan_data$P_sigma2
        ),
        tau_c1 = rep(0.3, stan_data$n_penalties_c1),
        tau_c2 = rep(0.3, stan_data$n_penalties_c2),
        tau_mu1 = rep(0.3, stan_data$n_penalties_mu1),
        tau_mu2 = rep(0.3, stan_data$n_penalties_mu2),
        tau_sigma1 = rep(0.3, stan_data$n_penalties_sigma1),
        tau_sigma2 = rep(0.3, stan_data$n_penalties_sigma2),

        z_c1 = rep(
          0,
          z_size(stan_data$n_penalties_c1, stan_data$penalty_sizes_c1)
        ),
        z_c2 = rep(
          0,
          z_size(stan_data$n_penalties_c2, stan_data$penalty_sizes_c2)
        ),
        z_mu1 = rep(
          0,
          z_size(stan_data$n_penalties_mu1, stan_data$penalty_sizes_mu1)
        ),
        z_mu2 = rep(
          0,
          z_size(stan_data$n_penalties_mu2, stan_data$penalty_sizes_mu2)
        ),
        z_sigma1 = rep(
          0,
          z_size(stan_data$n_penalties_sigma1, stan_data$penalty_sizes_sigma1)
        ),
        z_sigma2 = rep(
          0,
          z_size(stan_data$n_penalties_sigma2, stan_data$penalty_sizes_sigma2)
        )
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
#'
#' @param fit a `CmdStanMCMC` fit object
#' @param param_name one of "c1", "c2", "mu1", "mu2", "sigma1", "sigma2"
#' @param n_coef number of coefficients for this schedule parameter
#'
#' @return an `n_draws` x `n_coef` matrix of draws
#' @noRd
extract_beta_draws <- function(fit, param_name, n_coef) {
  var_names <- paste0("beta_", param_name, "[", seq_len(n_coef), "]")
  draws <- fit$draws(variables = var_names, format = "draws_matrix")
  unclass(draws)
}

#' Extract posterior draws of tau_<param> from a cmdstanr fit
#'
#' @param fit a `CmdStanMCMC` fit object
#' @param param_name one of "c1", "c2", "mu1", "mu2", "sigma1", "sigma2"
#' @param n_penalties number of penalty terms (smoothness hyperparameters)
#'   for this schedule parameter
#'
#' @return an `n_draws` x `n_penalties` matrix of draws, or an `n_draws` x 0
#'   matrix if `n_penalties` is 0
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

#' Extract beta and tau draws for all six schedule parameters
#'
#' @param fit a `CmdStanMCMC` fit object
#' @param processed output of `process_effects()`, used to determine
#'   `n_coef`/number of penalties for each schedule parameter
#'
#' @return a list with elements `beta` and `tau`, each a named list (one
#'   entry per schedule parameter) of draws matrices from
#'   `extract_beta_draws()`/`extract_tau_draws()`
#' @noRd
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

#' Compute quick convergence diagnostics for a cmdstanr fit
#'
#' @param cmdstanr_fit a `CmdStanMCMC` fit object
#'
#' @return a list with elements `max_rhat`, `min_ess` (minimum bulk ESS), and
#'   `n_divergent` (total divergent transitions across chains)
#' @noRd
pk_quick_diagnostics <- function(cmdstanr_fit) {
  summ <- cmdstanr_fit$summary(variables = NULL)

  list(
    max_rhat = max(summ$rhat, na.rm = TRUE),
    min_ess = min(summ$ess_bulk, na.rm = TRUE),
    n_divergent = sum(cmdstanr_fit$diagnostic_summary()$num_divergent)
  )
}
