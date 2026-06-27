#' @importFrom bayesplot pp_check
#' @export
bayesplot::pp_check

#' Posterior predictive checks for a pk_fit object
#'
#' @param object a pk_fit object
#' @param type which check to perform: "dens_overlay" (default),
#'   "hist", "intervals", "scatter_avg", "stat", "stat_grouped", "ribbon"
#' @param newdata data frame; defaults to training data
#' @param group optional column name for grouped checks
#' @param stat for type = "stat" / "stat_grouped": a function or character
#'   name of a function (e.g., "mean", "sd", "max")
#' @param ndraws number of replicate datasets to overlay (default 50 for
#'   density-style plots, NULL = all for other types)
#' @param yrep precomputed posterior predictive matrix (optional;
#'   avoids recomputation)
#' @param source character; "posterior" (default) or "prior"
#' @param ... additional arguments passed to the underlying bayesplot function
#'
#' @return a ggplot object
#'
#' @export
pp_check.pk_fit <- function(
  object,
  type = c(
    "dens_overlay",
    "hist",
    "intervals",
    "scatter_avg",
    "stat",
    "stat_grouped",
    "ribbon"
  ),
  newdata = NULL,
  group = NULL,
  stat = "mean",
  ndraws = NULL,
  yrep = NULL,
  source = c("posterior", "prior"),
  ...
) {
  checkmate::assert_class(object, "pk_fit")
  type <- match.arg(type)
  source <- match.arg(source)

  if (is.null(newdata)) {
    newdata <- object$data
  }

  if (is.null(yrep)) {
    yrep <- posterior_predict(
      object,
      newdata = newdata,
      source = source,
      ndraws = ndraws
    )
  }

  y_obs <- newdata[[object$births_col]]

  if (type %in% c("dens_overlay", "hist") && is.null(ndraws)) {
    n_show <- min(50, nrow(yrep))
    yrep_show <- yrep[sample.int(nrow(yrep), n_show), , drop = FALSE]
  } else {
    yrep_show <- yrep
  }

  switch(
    type,
    dens_overlay = bayesplot::ppc_dens_overlay(y_obs, yrep_show, ...),
    hist = bayesplot::ppc_hist(
      y_obs,
      yrep_show[1:min(8, nrow(yrep_show)), ],
      ...
    ),
    intervals = bayesplot::ppc_intervals(y_obs, yrep, ...),
    scatter_avg = bayesplot::ppc_scatter_avg(y_obs, yrep, ...),
    stat = pp_check_stat(y_obs, yrep, stat = stat, ...),
    stat_grouped = pp_check_stat_grouped(y_obs, yrep, group, stat, ...),
    ribbon = pp_check_ribbon(object, newdata, yrep, source, group, ...)
  )
}

#' @noRd
resolve_stat <- function(stat) {
  if (is.function(stat)) {
    return(stat)
  }
  if (is.character(stat)) {
    checkmate::assert_string(stat)
    fn <- match.fun(stat)
    return(fn)
  }

  cli::cli_abort("{.arg stat} must be a function or a character name of one.")
}

#' @noRd
pp_check_stat <- function(y_obs, yrep, stat, ...) {
  stat_fn <- resolve_stat(stat)
  bayesplot::ppc_stat(y_obs, yrep, stat = stat_fn, ...)
}

#' @noRd
pp_check_stat_grouped <- function(y_obs, yrep, newdata, group, stat, ...) {
  if (is.null(group)) {
    cli::cli_abort(
      "{.arg group} must be supplied for {.code type = \"stat_grouped\"}."
    )
  }
  checkmate::assert_choice(group, choices = names(newdata))

  stat_fn <- resolve_stat(stat)
  group_values <- newdata[[group]]
  bayesplot::ppc_stat_grouped(
    y_obs,
    yrep,
    group = group_values,
    stat = stat_fn,
    ...
  )
}

#' @noRd
pp_check_ribbon <- function(
  object,
  newdata,
  yrep,
  source = c("prior", "posterior"),
  group = NULL,
  conf.level = 0.95,
  age_grid = NULL,
  ...
) {
  source <- match.arg(source)
  if (is.null(age_grid)) {
    age_grid <- seq(15, 49, by = 0.5)
  }
  alpha <- (1 - conf.level) / 2

  exposure_col <- object$exposure_col
  age_col <- object$age_col
  births_col <- object$births_col

  exposures <- newdata[[exposure_col]]
  yrep_rates <- sweep(yrep, 2, exposures, FUN = "/")

  obs_summary <- tibble::tibble(
    age = newdata[[age_col]],
    asfr = newdata[[births_col]] / exposures,
    .lower = matrixStats::colQuantiles(yrep_rates, probs = alpha),
    .upper = matrixStats::colQuantiles(yrep_rates, probs = 1 - alpha)
  )

  cell_vars <- unique(unlist(lapply(object$effects[1:6], all.vars)))
  cell_vars <- setdiff(cell_vars, age_col)
  cell_data <- if (length(cell_vars) > 0) {
    dplyr::distinct(newdata[, cell_vars, drop = FALSE])
  } else {
    data.frame(row.names = 1)
  }

  curve_df <- predict(
    object,
    newdata = cell_data,
    type = "schedule_curve",
    age_grid = age_grid,
    summarize = TRUE,
    conf.level = conf.level,
    source = source
  )

  if (!is.null(group)) {
    checkmate::assert_choice(group, choices = names(newdata))
    obs_summary[[group]] <- newdata[[group]]
  }

  p <- ggplot2::ggplot()

  if (!is.null(group)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = curve_df,
        ggplot2::aes(
          x = age,
          ymin = .lower,
          ymax = .upper
          #fill = .data[[group]]
        ),
        alpha = 0.2
      ) +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(
          x = age,
          y = .fitted,
          #color = .data[[group]]
        )
      ) +
      ggplot2::geom_errorbar(
        data = obs_summary,
        ggplot2::aes(
          x = age,
          ymin = .lower,
          ymax = .upper,
          #color = .data[[group]]
        ),
        width = 0,
        alpha = 0.6
      ) +
      ggplot2::geom_point(
        data = obs_summary,
        ggplot2::aes(
          x = age,
          y = asfr,
          #color = .data[[group]]
        ),
        size = 0.8
      ) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", group)))
  } else {
    p <- p +
      ggplot2::geom_ribbon(
        data = curve_df,
        ggplot2::aes(x = age, ymin = .lower, ymax = .upper),
        alpha = 0.2,
        fill = "#3F7EBA"
      ) +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = age, y = .fitted),
        color = "#3F7EBA"
      ) +
      ggplot2::geom_errorbar(
        data = obs_summary,
        ggplot2::aes(x = age, ymin = .lower, ymax = .upper),
        width = 0,
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        data = obs_summary,
        ggplot2::aes(x = age, y = asfr)
      )
  }

  p +
    ggplot2::labs(
      title = "Posterior predictive ribbon: ASFRs vs. observed",
      x = "Age",
      y = "Age-specific fertility rate"
    ) +
    ggplot2::theme_minimal()
}
