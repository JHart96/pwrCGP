#' Estimate point of diminishing returns using diminishing returns/elbow method
#'
#' This function estimates the point at which small increases in correlation
#' come at increasingly high costs of sampling effort. If correlation is
#' suitably high, this can be used to determine how large window sizes should be
#' when subsetting data. However a low required sampling effort will not fix
#' poor correlation, and this should be considered when carrying out analyses.
#'
#' @param soc_diff Positive real value describing the social differentiation as
#'   estimated using `net_cor`.
#' @param rho_max Real value between 0 and 1 describing the maximum correlation
#'   of practical relevance.
#' @param I_max Maximum feasible sampling effort (interactions seen per dyad
#'   under even sampling). Leave it as is unless higher sampling is possible.
#'   Setting this too high is not a major problem, setting it too low might lead
#'   to unreliable results.
#'
#' @return List containing the optimal sampling effort and the correlation at
#'   that sampling effort.
#' @export
#'
#' @examples
pwr_elbow <- function(soc_diff, rho_max, I_max=1000) {
  I <- seq(0, I_max, 0.01) # Maximum sampling effort.

  rho_ <- (soc_diff * sqrt(I))/(sqrt(1 + I * soc_diff^2))

  rho_argmax <- which.min(abs(rho_ - rho_max))
  max_rho_ <- rho_[rho_argmax]
  max_I <- I[rho_argmax]

  theta <- atan2(max_rho_, max_I)
  co = cos(theta)
  si = sin(theta)
  rotation_matrix = matrix(c(co, -si, si, co), nrow=2, ncol=2)

  rotated_data <- cbind(rho_, I) %*% rotation_matrix

  rho_prime <- rho_[which.max(rotated_data[, 1])] # Optimal rho, by diminishing returns principle.

  I_prime <- I[which.max(rotated_data[, 1])] # Corresponding optimal sampling effort.

  list(sampling_effort=I_prime, correlation=rho_prime)
}
