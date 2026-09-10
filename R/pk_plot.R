#' Plot a fitted PK model
#'
#' @description
#' Intended to visualize a fitted \code{pk_fit} object (e.g. its fertility
#' schedule curve and/or posterior diagnostics). Not yet implemented --- this
#' method is currently a placeholder with an empty body and returns `NULL`.
#' Use \code{predict.pk_fit()} with `type = "schedule_curve"` together with
#' \pkg{ggplot2} directly, or \code{pk_prior_posterior()}, in the meantime.
#'
#' @param x A fitted `pk_fit` object, as returned by \code{pk_fit()}.
#'
#' @return `NULL`, invisibly.
#'
#' @import ggplot2
#' @export
plot.pk_fit <- function(x) {}
