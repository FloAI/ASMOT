#' @title Normalize Vector
#' @keywords internal
normalize_vec <- function(x) {
  if (sum(x) == 0) return(x)
  return(x / sum(x))
}

#' @title Auto-Tune Penalty Parameter (rho)
#' @export
estimate_rho <- function(M, quantile_cutoff = 0.5) {
  return(stats::quantile(M, probs = quantile_cutoff))
}

#' @title Correlation to Distance
#' @description Converts correlation matrix to a metric distance (1 - abs(cor)).
#' @keywords internal
cor_to_dist <- function(mat) {
  C <- cor(mat)
  C[is.na(C)] <- 0 
  return(1 - abs(C))
}
