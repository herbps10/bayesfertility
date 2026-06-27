maybe_subsample_draws <- function(draws_mat, ndraws) {
  if (is.null(ndraws) || ndraws >= nrow(draws_mat)) {
    return(draws_mat)
  }
  idx <- sample.int(nrow(draws_mat), ndraws)
  draws_mat[idx, , drop = FALSE]
}

maybe_subsample_draws_list <- function(draws_list, ndraws) {
  if (is.null(ndraws)) {
    return(draws_list)
  }
  n_draws <- nrow(draws_list[[1]])
  if (ndraws >= n_draws) {
    return(draws_list)
  }

  idx <- sample.int(n_draws, ndraws)
  lapply(draws_list, function(mat) mat[idx, , drop = FALSE])
}

summarize_draws_matrix <- function(draws_mat, conf.level = 0.95) {
  alpha <- (1 - conf.level) / 2
  tibble::tibble(
    .fitted = matrixStats::colMeans2(draws_mat),
    .median = matrixStats::colMedians(draws_mat),
    .lower = matrixStats::colQuantiles(draws_mat, probs = alpha),
    .upper = matrixStats::colQuantiles(draws_mat, probs = 1 - alpha)
  )
}

#' @importFrom matrixStats colMeans2 colMedians colQuantiles
#' @importFrom dplyr bind_cols
#' @noRd
summarize_draws_list <- function(draws_list, conf.level = 0.95) {
  alpha <- (1 - conf.level) / 2
  summaries <- lapply(names(draws_list), function(nm) {
    mat <- draws_list[[nm]]
    tibble::tibble(
      "{nm}.fitted" := matrixStats::colMeans2(mat),
      "{nm}.median" := matrixStats::colMedians(mat),
      "{nm}.lower" := matrixStats::colQuantiles(mat, probs = alpha),
      "{nm}.upper" := matrixStats::colQuantiles(mat, probs = 1 - alpha)
    )
  })

  dplyr::bind_cols(summaries)
}

#' Summarise schedule_curve draws into a long-format tibble
#'
#' @param curve_draws n_draws x n_cells x n_ages array
#' @param newdata the cell-level data frame
#' @param age_grid numeric vector of ages
#' @param conf.level credible interval width
#'
#' @return long-format tibble with one row per (cell, age) combination
#' @noRd
summarize_curve_draws <- function(
  curve_draws,
  newdata,
  age_grid,
  conf.level = 0.95
) {
  alpha <- (1 - conf.level) / 2
  n_cells <- dim(curve_draws)[2]
  n_ages <- dim(curve_draws)[3]

  # Compute summaries per (cell, age): collapse over draws (dimension 1)
  mean_mat <- apply(curve_draws, c(2, 3), mean)
  median_mat <- apply(curve_draws, c(2, 3), stats::median)
  lower_mat <- apply(curve_draws, c(2, 3), stats::quantile, probs = alpha)
  upper_mat <- apply(curve_draws, c(2, 3), stats::quantile, probs = 1 - alpha)

  # Expand newdata to one per (cell, age) combination
  expanded <- newdata[rep(seq_len(n_cells), each = n_ages), , drop = FALSE]
  expanded$age <- rep(age_grid, times = n_cells)

  expanded$.fitted <- as.vector(t(mean_mat))
  expanded$.median <- as.vector(t(median_mat))
  expanded$.lower <- as.vector(t(lower_mat))
  expanded$.upper <- as.vector(t(upper_mat))

  tibble::as_tibble(expanded)
}


#' Deduplicate training data based on variables that appear in any of the effects formulas
#' @importFrom dplyr distinct
#' @noRd
unique_cells <- function(object) {
  cell_vars <- unique(unlist(lapply(object$effects, all.vars)))
  cell_vars <- setdiff(cell_vars, object$age_col)
  dplyr::distinct(object$data[, cell_vars, drop = FALSE])
}
