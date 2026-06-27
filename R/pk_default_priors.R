#' Default per-parameter prior scales and centers
#'
#' @noRd
pk_default_prior_scales <- function(param_name) {
  list(
    "(Intercept)" = switch(
      param_name,
      c1 = 0.5,
      c2 = 0.5,
      mu1 = 0.7,
      mu2 = 0.2,
      sigma1 = 0.15,
      sigma2 = 0.2
    ),
    parametric = switch(
      param_name,
      c1 = 1.0,
      c2 = 1.0,
      mu1 = 0.5,
      mu2 = 0.3,
      sigma1 = 0.15,
      sigma2 = 0.3
    ),
    smooth = NULL
  )
}

pk_default_prior_centers <- function(param_name) {
  list(
    "(Intercept)" = switch(
      param_name,
      c1 = -3.0,
      c2 = -2.0,
      mu1 = 20.0,
      mu2 = log(8),
      sigma1 = log(3.5),
      sigma2 = log(7.5)
    ),
    parametric = 0,
    smooth = 0
  )
}

pk_default_tau_scales <- function(param_name) {
  switch(
    param_name,
    c1 = 0.5,
    c2 = 0.5,
    mu1 = 0.2,
    mu2 = 0.3,
    sigma1 = 0.2,
    sigma2 = 0.2
  )
}
