#' Calculate the posterior total fertility rate (TFR)
#'
#' @description
#' Computes the total fertility rate implied by a fitted PK model: the
#' age-specific fertility rate (ASFR) schedule is evaluated on a grid over
#' `age_range` spaced by `age_step` and integrated (via the trapezoidal rule)
#' to give TFR, separately for each covariate combination ("cell") in
#' `newdata`.
#'
#' @param object A fitted `pk_fit` object, as returned by \code{pk_fit()}.
#' @param newdata A data frame giving the covariate values ("cells") to
#'   compute TFR for. Defaults to the distinct covariate combinations implied
#'   by the model's effects formulas in the training data (see the internal
#'   `unique_cells()` helper) if `NULL`.
#' @param age_range Length-2 numeric vector giving the lower and upper age
#'   bounds to integrate the schedule over. Defaults to `c(15, 49)`.
#' @param age_step Spacing of the age grid used for the trapezoidal
#'   integration. Defaults to `0.5`.
#' @param source Whether to compute TFR from `"posterior"` (default) or
#'   `"prior"` draws; `"prior"` requires the fit to have been created with
#'   `sample_prior = "yes"` or `"only"`.
#' @param summarize If `TRUE` (default), collapse draws into a point estimate
#'   (mean and median) and a `conf.level` credible interval. If `FALSE`,
#'   return the raw draws.
#' @param conf.level Credible interval width used when `summarize = TRUE`.
#' @param ndraws Optional number of draws to randomly subsample before
#'   computing TFR; `NULL` (default) uses all available draws.
#' @param ... Currently unused; accepted for interface consistency with
#'   \code{predict.pk_fit()}.
#'
#' @return If `summarize`, a tibble with one row per row of `newdata`,
#'   including its covariate columns plus `.fitted`, `.median`, `.lower`, and
#'   `.upper` giving the TFR estimate and its credible interval. If not, an
#'   `n_draws` x `nrow(newdata)` matrix of TFR draws.
#'
#' @examples
#' \dontrun{
#' posterior_tfr(fit)
#' posterior_tfr(fit, newdata = pred_levels, conf.level = 0.95)
#' }
#'
#' @export
posterior_tfr <- function(
  object,
  newdata = NULL,
  age_range = c(15, 49),
  age_step = 0.5,
  source = c("posterior", "prior"),
  summarize = TRUE,
  conf.level = 0.95,
  ndraws = NULL,
  ...
) {
  source <- match.arg(source)
  checkmate::assert_class(object, "pk_fit")
  checkmate::assert_data_frame(newdata, null.ok = TRUE)
  checkmate::assert_numeric(age_range, len = 2, finite = TRUE, sorted = TRUE)
  checkmate::assert_number(age_step, lower = 0, finite = TRUE)
  checkmate::assert_flag(summarize)
  checkmate::assert_number(conf.level, lower = 0, upper = 1)
  checkmate::assert_count(ndraws, positive = TRUE, null.ok = TRUE)

  if (is.null(newdata)) {
    newdata <- unique_cells(object)
  }

  age_grid <- seq(age_range[1], age_range[2], by = age_step)

  curve_draws <- predict(
    object,
    newdata = newdata,
    type = "schedule_curve",
    age_grid = age_grid,
    source = source,
    summarize = FALSE,
    ndraws = ndraws
  )

  # Trapezoidal integration over age
  n_ages <- dim(curve_draws)[3]
  endpoints <- curve_draws[,, 1] + curve_draws[,, n_ages]
  totals <- apply(curve_draws, c(1, 2), sum)
  tfr_draws <- age_step * (totals - 0.5 * endpoints)

  if (summarize) {
    summary_cols <- summarize_draws_matrix(tfr_draws, conf.level = conf.level)
    dplyr::bind_cols(tibble::as_tibble(newdata), summary_cols)
  } else {
    tfr_draws
  }
}
