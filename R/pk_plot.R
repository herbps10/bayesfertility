#' Plot a fitted PK model
#'
#' @description
#' Plots the fitted fertility schedule (age-specific fertility rate vs. age,
#' with a credible ribbon) implied by a \code{pk_fit} object, via
#' \code{predict.pk_fit(type = "schedule_curve")}. If the model's effects
#' formulas reference any covariates, the plot is automatically faceted over
#' the distinct covariate combinations found in the training data. For a
#' comparison against observed data, see \code{pp_check(type = "ribbon")}
#' instead.
#'
#' @param x A fitted `pk_fit` object, as returned by \code{pk_fit()}.
#' @param newdata A data frame giving the covariate values to facet over.
#'   Defaults to the training data (`x$data`) if `NULL`; only the columns
#'   referenced in `x$effects` are used, deduplicated to distinct
#'   combinations.
#' @param age_grid Numeric vector of ages at which to evaluate the schedule;
#'   defaults to `seq(15, 49, by = 1)` (see `predict.pk_fit()`).
#' @param source Whether to plot the `"posterior"` (default) or `"prior"`
#'   fitted schedule; `"prior"` requires the fit to have been created with
#'   `sample_prior = "yes"` or `"only"`.
#' @param conf.level Credible interval width for the ribbon.
#' @param ... Additional arguments, currently unused.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot(fit)
#' plot(fit, source = "prior", conf.level = 0.8)
#' }
#'
#' @import ggplot2
#' @export
plot.pk_fit <- function(
  x,
  newdata = NULL,
  age_grid = NULL,
  source = c("posterior", "prior"),
  conf.level = 0.95,
  ...
) {
  source <- match.arg(source)
  checkmate::assert_class(x, "pk_fit")
  checkmate::assert_data_frame(newdata, null.ok = TRUE)
  checkmate::assert_number(conf.level, lower = 0, upper = 1)

  if (is.null(newdata)) {
    newdata <- x$data
  }

  cell_vars <- schedule_cell_vars(x)
  cell_data <- if (length(cell_vars) > 0) {
    dplyr::distinct(newdata[, cell_vars, drop = FALSE])
  } else {
    data.frame(row.names = 1)
  }

  curve_df <- predict(
    x,
    newdata = cell_data,
    type = "schedule_curve",
    age_grid = age_grid,
    summarize = TRUE,
    conf.level = conf.level,
    source = source
  )

  plot_schedule_curve(
    curve_df,
    facet_vars = cell_vars,
    title = "Fitted fertility schedule"
  )
}

#' Build a ggplot of a fitted schedule curve, optionally with an observed overlay
#'
#' @param curve_df summarized schedule-curve tibble from
#'   `predict.pk_fit(type = "schedule_curve")` (columns: facet covariates,
#'   age, .fitted, .median, .lower, .upper)
#' @param obs_df optional tibble of observed points/intervals to overlay
#'   (columns: age, asfr, .lower, .upper, and any facet covariates); `NULL`
#'   to omit the observed layer
#' @param facet_vars character vector of column names to facet by;
#'   `character(0)` for a single unfaceted panel
#' @param title plot title
#'
#' @return a ggplot object
#' @noRd
plot_schedule_curve <- function(
  curve_df,
  obs_df = NULL,
  facet_vars = character(0),
  title = NULL
) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = curve_df,
      ggplot2::aes(x = age, ymin = .lower, ymax = .upper),
      fill = "#3F7EBA",
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = age, y = .fitted),
      color = "#3F7EBA"
    )

  if (!is.null(obs_df)) {
    p <- p +
      ggplot2::geom_errorbar(
        data = obs_df,
        ggplot2::aes(x = age, ymin = .lower, ymax = .upper),
        width = 0,
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        data = obs_df,
        ggplot2::aes(x = age, y = asfr),
        size = 0.8
      )
  }

  if (length(facet_vars) > 0) {
    p <- p +
      ggplot2::facet_wrap(
        stats::as.formula(paste("~", paste(facet_vars, collapse = " + ")))
      )
  }

  p +
    ggplot2::labs(
      title = title,
      x = "Age",
      y = "Age-specific fertility rate"
    ) +
    ggplot2::theme_minimal()
}
