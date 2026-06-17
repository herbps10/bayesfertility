pk_effects <- function(
  c1 = ~1,
  c2 = ~1,
  mu1 = ~1,
  mu2 = ~1,
  sigma1 = ~1,
  sigma2 = ~1
) {
  list(
    c1 = c1,
    c2 = c2,
    mu1 = mu1,
    mu2 = mu2,
    sigma1 = sigma1,
    sigma2 = sigma2
  )
}

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
  ...
) {
  if (Sys.getenv("BAYESFERTILITY_DEV") == "TRUE") {
    model <- cmdstanr::cmdstan_model("src/stan/pk2.stan")
  } else {
    model <- instantiate::stan_package_model(
      name = "pk2",
      package = "bayesfertility"
    )
  }

  A <- max(data[[age]])

  stan_data <- list(
    A = A,
    N = nrow(data),
    B = data[[births]],
    E = data[[exposure]],
    age = data[[age]]
  )

  fit <- model$sample(stan_data, ...)

  fit
}
