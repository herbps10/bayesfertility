#' Build the Stan data list from the processed effects
#'
#' @param processed output of process_effects()
#' @param data the original data frame
#' @param births,exposure,age column names
#'
#' @noRd
build_stan_data <- function(processed, data, births, exposure, age) {
  param_names <- c("c1", "c2", "mu1", "mu2", "sigma1", "sigma2")
  N <- nrow(data)

  ages <- sort(unique(data[[age]]))

  stan_data <- list(
    A = length(ages),
    N = N,
    B = data[[births]],
    E = data[[exposure]],
    ages = ages,
    age = match(data[[age]], ages)
  )

  for (p in param_names) {
    pp <- processed[[p]]

    # Design matrix
    stan_data[[paste0("P_", p)]] <- pp$n_coef
    stan_data[[paste0("X_", p)]] <- pp$X

    # Per-term bookkeeping (parametric + smooth terms)
    terms_info <- pp$terms_info
    n_terms <- length(terms_info)

    term_starts <- vapply(terms_info, function(t) t$col_idx[1], integer(1))
    term_sizes <- vapply(terms_info, function(t) length(t$col_idx), integer(1))

    term_scales <- ifelse(is.na(pp$prior_scales), 0, pp$prior_scales)
    term_centers <- pp$prior_centers

    stan_data[[paste0("n_terms_", p)]] <- n_terms
    stan_data[[paste0("term_starts_", p)]] <- term_starts
    stan_data[[paste0("term_sizes_", p)]] <- term_sizes
    stan_data[[paste0("term_scales_", p)]] <- term_scales
    stan_data[[paste0("term_centers_", p)]] <- term_centers

    # Penalty bookkeeping
    n_pen <- length(pp$penalties)
    stan_data[[paste0("n_penalties_", p)]] <- n_pen

    pen_is_re <- vapply(
      pp$penalties,
      function(e) as.integer(e$is_re),
      integer(1)
    )
    stan_data[[paste0("penalty_is_re_", p)]] <- pen_is_re

    if (n_pen > 0) {
      pen_sizes <- vapply(
        pp$penalties,
        function(e) length(e$col_idx),
        integer(1)
      )
      pen_starts_X <- vapply(
        pp$penalties,
        function(e) e$col_idx[1],
        integer(1)
      )
      pen_starts_S <- c(1L, head(cumsum(pen_sizes), -1) + 1L)
      pen_smooth_idx <- vapply(
        pp$penalties,
        function(e) e$smooth_idx,
        integer(1)
      )

      stan_data[[paste0("penalty_sizes_", p)]] <- pen_sizes
      stan_data[[paste0("penalty_starts_X_", p)]] <- pen_starts_X
      stan_data[[paste0("penalty_starts_S_", p)]] <- pen_starts_S
      stan_data[[paste0("penalty_smooth_", p)]] <- pen_smooth_idx
      stan_data[[paste0("n_smooth_groups_", p)]] <- length(unique(
        pen_smooth_idx
      ))
      stan_data[[paste0("tau_scales_", p)]] <- pp$tau_scales

      # Pack penalty matrices into a block-diagonal matrix
      total <- sum(pen_sizes)
      S_packed <- matrix(0, nrow = total, ncol = total)
      offset <- 0
      for (j in seq_len(n_pen)) {
        sz <- pen_sizes[j]
        S_packed[
          (offset + 1):(offset + sz),
          (offset + 1):(offset + sz)
        ] <- pp$penalties[[j]]$S
        offset <- offset + sz
      }
      stan_data[[paste0("S_packed_", p)]] <- S_packed
    } else {
      stan_data[[paste0("penalty_sizes_", p)]] <- integer(0)
      stan_data[[paste0("penalty_starts_X_", p)]] <- integer(0)
      stan_data[[paste0("penalty_starts_S_", p)]] <- integer(0)
      stan_data[[paste0("penalty_smooth_", p)]] <- integer(0)
      stan_data[[paste0("n_smooth_groups_", p)]] <- 0L
      stan_data[[paste0("tau_scales_", p)]] <- numeric(0)
      stan_data[[paste0("S_packed_", p)]] <- matrix(0, 0, 0)
    }
  }

  stan_data
}
