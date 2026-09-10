#' Predict fertility schedule quantities from a fitted PK model
#'
#' @description
#' Computes model-implied quantities from a fitted \code{pk_fit} object at a
#' set of covariate values, using posterior (or prior) draws. Four kinds of
#' quantity can be requested via `type`:
#' \itemize{
#'   \item `"asfr"`: the age-specific fertility rate at specific (row, age)
#'     combinations in `newdata`.
#'   \item `"linear_predictor"`: the working-scale linear predictor `eta` for
#'     each of the six schedule parameters (`c1`, `c2`, `mu1`, `mu2`,
#'     `sigma1`, `sigma2`), one value per row of `newdata`.
#'   \item `"schedule_params"`: the same six schedule parameters transformed
#'     back to their natural scale (see \code{pk_effects()} for the
#'     working-scale transforms), one value per row of `newdata`.
#'   \item `"schedule_curve"`: the full fertility schedule evaluated over
#'     `age_grid`, for each distinct covariate combination ("cell") implied
#'     by `newdata`.
#' }
#'
#' @param object A fitted `pk_fit` object, as returned by \code{pk_fit()}.
#' @param newdata A data frame giving the covariate values to predict at
#'   (and, for `type = "asfr"`, ages as well, unless `age` is supplied
#'   separately). Defaults to the training data (`object$data`) if `NULL`.
#' @param type One of `"asfr"`, `"schedule_params"`, `"linear_predictor"`, or
#'   `"schedule_curve"` (see Description).
#' @param source Whether to predict from `"posterior"` (default) or `"prior"`
#'   draws; `"prior"` requires the fit to have been created with
#'   `sample_prior = "yes"` or `"only"`.
#' @param age Only used when `type = "asfr"`. One of `NULL` (use the age
#'   column named in `object$age_col`, looked up in `newdata`), a string
#'   naming a column of `newdata` to use as age, or a numeric vector of ages
#'   with length `nrow(newdata)`.
#' @param age_grid Only used when `type = "schedule_curve"`. A numeric vector
#'   of ages at which to evaluate the schedule; defaults to `seq(15, 49, by = 1)`.
#' @param summarize If `TRUE` (default), collapse draws into a point estimate
#'   (mean and median) and a `conf.level` credible interval. If `FALSE`,
#'   return the raw draws.
#' @param conf.level Credible interval width used when `summarize = TRUE`.
#' @param ndraws Optional number of draws to randomly subsample before
#'   summarizing or returning; `NULL` (default) uses all available draws.
#' @param ... Unused; present for S3 method consistency.
#'
#' @return The shape depends on `type` and `summarize`:
#'   \describe{
#'     \item{`type = "asfr"`}{If `summarize`, a tibble with one row per row of
#'       `newdata` and columns `.fitted`, `.median`, `.lower`, `.upper`. If
#'       not, an `n_draws` x `nrow(newdata)` matrix of ASFR draws.}
#'     \item{`type = "linear_predictor"` or `"schedule_params"`}{If
#'       `summarize`, a tibble with one row per row of `newdata` and, for each
#'       of the six schedule parameters, columns named
#'       `<parameter>.fitted`/`.median`/`.lower`/`.upper`. If not, a named
#'       list (one element per schedule parameter) of `n_draws` x
#'       `nrow(newdata)` matrices.}
#'     \item{`type = "schedule_curve"`}{If `summarize`, a long-format tibble
#'       with one row per (cell, age) combination, including the covariate
#'       columns from `newdata`, `age`, and `.fitted`/`.median`/`.lower`/
#'       `.upper`. If not, an `n_draws` x `n_cells` x `n_ages` array (with
#'       `age_grid` values as the third dimension's names).}
#'   }
#'
#' @examples
#' \dontrun{
#' predict(fit, type = "asfr")
#' predict(fit, newdata = pred_levels, type = "schedule_params")
#' predict(fit, type = "schedule_curve", conf.level = 0.95)
#' }
#'
#' @export
predict.pk_fit <- function(
  object,
  newdata = NULL,
  type = c("asfr", "schedule_params", "linear_predictor", "schedule_curve"),
  source = c("posterior", "prior"),
  age = NULL,
  age_grid = NULL,
  summarize = TRUE,
  conf.level = 0.95,
  ndraws = NULL,
  ...
) {
  source <- match.arg(source)
  checkmate::assert_class(object, "pk_fit")
  checkmate::assert_data_frame(newdata, null.ok = TRUE)
  checkmate::assert_flag(summarize)
  checkmate::assert_number(conf.level, lower = 0, upper = 1)
  checkmate::assert_count(ndraws, positive = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(
    age_grid,
    null.ok = TRUE,
    any.missing = FALSE,
    finite = TRUE,
    min.len = 2
  )
  type <- match.arg(type)

  if (is.null(newdata)) {
    newdata <- object$data
  }

  if (type == "asfr") {
    age_values <- resolve_age(age, newdata, object)
    draws <- compute_asfr_draws(object, newdata, age_values, source)
    draws <- maybe_subsample_draws(draws, ndraws)

    if (summarize) {
      summarize_draws_matrix(draws, conf.level = conf.level)
    } else {
      draws
    }
  } else if (type == "linear_predictor") {
    eta_draws <- compute_eta_draws(object, newdata, source)
    eta_draws <- maybe_subsample_draws_list(eta_draws, ndraws)

    if (summarize) {
      summarize_draws_list(eta_draws, conf.level = conf.level)
    } else {
      eta_draws
    }
  } else if (type == "schedule_params") {
    eta_draws <- compute_eta_draws(object, newdata, source)
    param_draws <- eta_to_params(eta_draws)
    param_draws <- maybe_subsample_draws_list(param_draws, ndraws)

    if (summarize) {
      summarize_draws_list(param_draws, conf.level = conf.level)
    } else {
      param_draws
    }
  } else if (type == "schedule_curve") {
    if (is.null(age_grid)) {
      age_grid <- seq(15, 49, by = 1)
    }

    eta_draws <- compute_eta_draws(object, newdata, source)
    eta_draws <- maybe_subsample_draws_list(eta_draws, ndraws)
    param_draws <- eta_to_params(eta_draws)

    curve_draws <- evaluate_schedule_grid(param_draws, age_grid)

    if (summarize) {
      return(summarize_curve_draws(
        curve_draws,
        newdata,
        age_grid,
        conf.level = conf.level
      ))
    } else {
      dimnames(curve_draws) <- list(
        draw = NULL,
        cell = NULL,
        age = as.character(age_grid)
      )
      return(curve_draws)
    }
  }
}

#' Compute working-scale linear predictors for all 6 schedule parameters
#'
#' @param fit a `pk_fit` object
#' @param newdata data frame of covariate values to predict at
#' @param source "posterior" (default) or "prior"; which beta draws to use
#'
#' @return a named list with elements c1, c2, mu1, mu2, sigma1, sigma2. Each
#'   is an n_draws x n_rows matrix of draws on the working scale
#'
#' @noRd
compute_eta_draws <- function(fit, newdata, source = c("posterior", "prior")) {
  source <- match.arg(source)
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")

  beta_source <- switch(
    source,
    posterior = fit$beta_draws,
    prior = fit$prior_beta_draws
  )

  if (is.null(beta_source)) {
    cli::cli_abort(c(
      "No {source} draws available in this {.cls pk_fit} object.",
      "i" = if (source == "prior") {
        "Refit with {.code sample_prior = \"yes\"} to enable prior comparison."
      } else {
        "This fit was created with {.code sample_prior = \"only\"}. No posterior draws are available."
      }
    ))
  }

  # Build prediction matrices
  X_pred <- lapply(param_names, function(p) {
    build_prediction_matrix(
      processed_param = fit$processed[[p]],
      formula = fit$effects[[p]],
      newdata = newdata,
      training_data = fit$data
    )
  })

  names(X_pred) <- param_names

  eta <- lapply(param_names, function(p) {
    beta_source[[p]] %*% t(X_pred[[p]])
  })
  names(eta) <- param_names
  eta
}

#' Convert eta from working scale to parameter scale
#'
#' @param eta output of `compute_eta_draws()`: a named list with elements
#'   c1, c2, mu1, mu2, sigma1, sigma2, each an n_draws x n_rows matrix
#'
#' @return a named list of the same shape as `eta`, with each schedule
#'   parameter transformed to its natural scale (exponentiated for c1/c2/
#'   sigma1/sigma2; mu1 unchanged; mu2 as `mu1 + exp(eta$mu2)`)
#' @noRd
eta_to_params <- function(eta) {
  mu1 <- eta$mu1
  list(
    c1 = exp(eta$c1),
    c2 = exp(eta$c2),
    mu1 = mu1,
    mu2 = mu1 + exp(eta$mu2),
    sigma1 = exp(eta$sigma1),
    sigma2 = exp(eta$sigma2)
  )
}

#' Build a prediction-time design matrix for one schedule parameter
#'
#' @param processed_param one element of fit$processed
#' @param formula original formula for this parameter
#' @param newdata data frame at which to predict
#' @param training_data the original training data
#'
#' @return a matrix with nrow(newdata) rows and processed_params$n_coef columns
#'
#' @noRd
build_prediction_matrix <- function(
  processed_param,
  formula,
  newdata,
  training_data
) {
  # Check for whether the formula only has parametric components (no smooths)
  if (is.null(processed_param$pregam)) {
    xlevels <- stats::.getXlevels(stats::terms(formula), training_data)
    return(stats::model.matrix(formula, data = newdata, xlev = xlevels))
  }

  pregam <- processed_param$pregam

  if (!is.null(pregam$pterms)) {
    xlevels_p <- pregam$xlevels %||%
      stats::.getXlevels(pregam$pterms, training_data)

    newdata[["..dummy_y.."]] <- rep(0, nrow(newdata))
    X_param <- stats::model.matrix(
      pregam$pterms,
      data = newdata,
      xlev = xlevels_p
    )
  } else {
    X_param <- matrix(numeric(0), nrow = nrow(newdata), ncol = 0)
  }

  smooth_blocks <- lapply(pregam$smooth, function(sm) {
    mgcv::PredictMat(sm, data = newdata)
  })

  X_pred <- do.call(cbind, c(list(X_param), smooth_blocks))

  if (ncol(X_pred) != processed_param$n_coef) {
    stop(sprintf(
      "Prediction matrix has %d columns but training matrix has %d. Check that newdata has the same factor levels as training data.",
      ncol(X_pred),
      processed_param$n_coef
    ))
  }

  X_pred
}

#' Compute ASFR posterior draws at a set of (cell, age) combinations
#'
#' @param fit a `pk_fit` object
#' @param newdata data frame of covariate values, one row per cell
#' @param age_values numeric vector of ages, length `nrow(newdata)`, giving
#'   the age to evaluate ASFR at for each row of `newdata`
#' @param source "posterior" or "prior"; which draws to use
#'
#' @return n_draws x n_rows matrix of ASFR draws
#' @noRd
compute_asfr_draws <- function(fit, newdata, age_values, source) {
  eta <- compute_eta_draws(fit, newdata, source)
  params <- eta_to_params(eta)

  age_mat <- matrix(
    age_values,
    nrow = nrow(params$c1),
    ncol = length(age_values),
    byrow = TRUE
  )

  asfr <- params$c1 *
    exp(-((age_mat - params$mu1) / params$sigma1)^2) +
    params$c2 * exp(-((age_mat - params$mu2) / params$sigma2)^2)

  asfr
}

#' Evaluate Peristera-Kostaki Model 2 schedule on a grid of ages for each cell
#'
#' @param params output of eta_to_params() (list of n_draws x n_cells matrices)
#' @param age_grid numeric vector of ages
#' @return n_draws x n_cells x n_ages array
#' @noRd
evaluate_schedule_grid <- function(params, age_grid) {
  n_draws <- nrow(params$c1)
  n_cells <- ncol(params$c1)
  n_ages <- length(age_grid)

  result <- array(0, dim = c(n_draws, n_cells, n_ages))

  for (a_idx in seq_along(age_grid)) {
    age <- age_grid[a_idx]
    result[,, a_idx] <- params$c1 *
      exp(-((age - params$mu1) / params$sigma1)^2) +
      params$c2 * exp(-((age - params$mu2) / params$sigma2)^2)
  }

  result
}

#' Resolve the age values to use for prediction
#'
#' @param age NULL (use `object$age_col` looked up in `newdata`), a column
#'   name in `newdata`, or a numeric vector of length `nrow(newdata)`
#' @param newdata data frame to resolve age from
#' @param object a `pk_fit` object
#'
#' @return a numeric vector of length `nrow(newdata)`
#' @noRd
resolve_age <- function(age, newdata, object) {
  # Resolve age values
  if (is.null(age)) {
    age_col <- object$age_col
    checkmate::assert_choice(
      age_col,
      choices = names(newdata),
      .var.name = "age column from fit"
    )
    return(newdata[[age_col]])
  } else if (is.character(age)) {
    checkmate::assert_string(age)
    checkmate::assert_choice(age, choices = names(newdata))
    return(newdata[[age]])
  } else if (is.numeric(age)) {
    checkmate::assert_numeric(
      age,
      len = nrow(newdata),
      any.missing = FALSE,
      finite = TRUE
    )
    return(age)
  }

  cli::cli_abort(c(
    "{.arg age} must be {.code NULL}, a column name, or a numeric vector.",
    "x" = "You supploied an object of class {.cls {class(age)}}."
  ))
}
