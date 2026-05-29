#' @keywords internal
# simple functions used in mHMM
dif_matrix <- function(rows, cols, data = NA) {
  return(matrix(data, ncol = cols, nrow = rows))
}

#' Obtain stationary distribution of a transition probability matrix
#'
#' @param gamma
#'
#' @returns Vector with marginal probabilities.
#'
#' @export
stat_dist <- function(gamma) {
  m <- nrow(gamma)
  delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  return(delta)
}

#' Return relative transition probabilities
#'
#' @param gamma transition probability matrix
#'
#' @returns Relative transition probability matrix
#'
#' @export
gamma_relative <- function(gamma) {
  m <- nrow(gamma)
  gamma_rel <- gamma / (rowSums(gamma) - diag(gamma))
  diag(gamma_rel) <- 0
  return(gamma_rel)
}
