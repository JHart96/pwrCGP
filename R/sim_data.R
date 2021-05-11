#' Simulate networks based on Gamma-Poisson model.
#'
#' Simulates networks based on the Gamma-Poisson model and outputs corresponding
#' adjacency matrices.
#'
#' @param n Integer greater than zero describing number of nodes in the network.
#' @param mean_sampling Real number greater than zero describing the mean amount
#'   of time spent sampling each dyad.
#' @param S Real number greater than zero describing the desired social
#'   differentiation of the network.
#' @param mu Real number greater than zero describing the mean number of
#'   interactions observed per dyad.
#' @param directed TRUE/FALSE indicating if the network should be directed.
#'
#' @return A list of three `n x n` matrices. `X` is the interaction observation
#'   matrix, `D` is the sampling time, and `A` is a matrix of the true
#'   interaction rates.
#' @export
#'
#' @examples
simulate_data_gp <- function (n, mean_sampling, S, mu, directed=FALSE) {
  # Convert social differentiation and mean interaction rate to Gamma parameters a and b.
  a <- 1/S**2
  b = a/mu

  # Variable sampling effort.
  d <- rpois(n^2, mean_sampling) # Individual sampling times.
  d[d == 0] = 1 # Every dyad must be sampled for a non-zero amount of time.

  # Generate association rates.
  alpha <- rgamma(n^2, a, b) # Actual rate.

  # Generate matrices.
  A <- matrix(alpha, n, n)
  X <- matrix(rpois(n^2, alpha*d), n, n)
  D <- matrix(d, n, n)

  # If directed only remove the diagonals, if undirected remove the
  # lower triangle and replace with upper triangle to get symmetry.
  if (directed) {
    A <- A * upper.tri(A) + A * lower.tri(A)
    X <- X * upper.tri(X) + X * lower.tri(X)
    D <- D * upper.tri(D) + D * lower.tri(D)
  } else {
    A <- A * upper.tri(A); A <- A + t(A)
    X <- X * upper.tri(X); X <- X + t(X)
    D <- D * upper.tri(D); D <- D + t(D)
  }

  list(X=X, D=D, A=A)
}
