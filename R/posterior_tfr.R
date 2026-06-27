#' Calculate posterior TFR
#' @param object object
#' @export
posterior_tfr <- function(
  object,
  newdata = NULL,
  age_range = c(15, 49),
  age_step = 0.5,
  source = c("prior", "posterior"),
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
