#' @importFrom rstantools posterior_predict
#' @export
rstantools::posterior_predict

#' Draw from the posterior predictive distribution of births
#'
#' For each posterior draw and each (cell, age) combination in `newdata`,
#' simulates an integer count from Poisson(exposure * ASFR) where ASFR is
#' drawn from the posterior (or prior, if requested).
#'
#' @param object a pk_fit object
#' @param newdata data frame at which to predict; defaults to training data
#' @param age column name in newdata, numeric vector of ages, or NULL to use the
#'   age column from the original fit
#' @param exposure column name in newdata, numeric vector, or NULL to use the
#'   exposure column from the original fit
#' @param source character; "posterior" (default) or "prior"
#' @param ndraws optional integer; subsample draws for speed
#' @param ... unused
#'
#' @return an n_draws x n_obs integer matrix of posterior predictive
#'   draws of births
#'
#' @export
posterior_predict.pk_fit <- function(
  object,
  newdata = NULL,
  age = NULL,
  exposure = NULL,
  source = c("posterior", "prior"),
  ndraws = NULL,
  ...
) {
  checkmate::assert_class(object, "pk_fit")
  checkmate::assert_data_frame(newdata, null.ok = TRUE)
  checkmate::assert_count(ndraws, positive = TRUE, null.ok = TRUE)
  source <- match.arg(source)

  if (is.null(newdata)) {
    newdata <- object$data
  }

  exposure_values <- resolve_exposure(exposure, newdata, object)

  rate_draws <- predict(
    object,
    newdata = newdata,
    type = "asfr",
    age = age,
    source = source,
    summarize = FALSE,
    ndraws = ndraws
  )

  n_draws <- nrow(rate_draws)
  n_obs <- ncol(rate_draws)
  expected_mat <- rate_draws *
    matrix(exposure_values, nrow = n_draws, ncol = n_obs, byrow = TRUE)

  yrep <- matrix(
    stats::rpois(n_draws * n_obs, lambda = as.vector(expected_mat)),
    nrow = n_draws,
    ncol = n_obs
  )

  yrep
}

#' @noRd
resolve_exposure <- function(exposure, newdata, object) {
  if (is.null(exposure)) {
    exposure_col <- object$exposure_col
    checkmate::assert_choice(
      exposure_col,
      choices = names(newdata),
      .var.name = "exposure column from fit"
    )
    return(newdata[[exposure_col]])
  }

  if (is.character(exposure)) {
    checkmate::assert_string(exposure)
    checkmate::assert_choice(exposure, choices = names(newdata))
    return(newdata[[exposure]])
  } else if (is.numeric(exposure)) {
    checkmate::assert_numeric(
      exposure,
      len = nrow(newdata),
      lower = 0,
      any.missing = FALSE,
      finite = TRUE
    )
    return(exposure)
  }

  cli::cli_abort(c(
    "{.arg exposure} must be {.code NULL}, a column name, or a numeric vector.",
    "x" = "You supplied an object of class {.cls {class(exposure)}}."
  ))
}
