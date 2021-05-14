#' Estimate network correlation and social differentiation for interaction
#' networks.
#'
#' Estimates properties of an estimated network with `n` nodes including network
#' correlation and social differentiation.
#'
#' By default this function plots a QQ diagnostic plot to indicate how well the
#' raw data fit the Gamma-Poisson model assumed by this method. If most of the
#' points are close to the diagonal line, a good fit can be assumed. Sometimes
#' the tails of the distribution may curve away from the line, but provided this
#' is limited to the tails it is unlikely to significantly affect the results.
#'
#' Note: This function will only work with interaction rate networks.
#'
#' @param X n x n integer-valued matrix with zero-valued diagonals where each
#'   entry (i,j) gives the number of observed interactions between i and j.
#' @param D n x n real-valued positive matrix with zero-valued diagonals where
#'   each entry (i, j) gives the amount of time spent sampling i and j.
#' @param directed TRUE/FALSE indicating if the network is directed or not.
#' @param ci Real value between 0 and 1 giving the width of the confidence
#'   interval.
#' @param show_diagnostic TRUE/FALSE indicating if to show the diagnostic QQ
#'   plot.
#'
#' @return Summary table of network properties including network correlation and
#'   social differentiation.
#' @export
#'
#' @examples
#' # Simulate observation data.
#' sim_data <- simulate_data_gp(20, 10, 0.25, 0.5)
#' X <- sim_data$X
#' D <- sim_data$D
#'
#' # Estimate network sampling properties.
#' net_cor(X, D)
net_cor <- function(X, D, directed=FALSE, ci=0.95, show_diagnostic=TRUE) {
  # Likelihood function for negative binomial model.
  lk_nbinom <- function(par, x, d) {
    par <- exp(par)
    r <- par[1]
    p <- par[2]/(par[2] + d)
    sum(dnbinom(x, r, p, log=TRUE))
  }

  # Convert matrices to list using upper triangle.
  x <- X[upper.tri(X)]
  d <- D[upper.tri(D)]

  # If network is directed, include lower triangle too.
  if (directed) {
    x <- c(x, X[lower.tri(X)])
    d <- c(d, D[lower.tri(D)])
  }

  # Estimate parameters using optim.
  optim_obj <- optim(c(0, 0), function(par) -lk_nbinom(par, x, d), hessian=TRUE)

  # Generate samples of parameters using quadratic approximation.
  parameters <- exp(MASS::mvrnorm(1e5, optim_obj$par, solve(optim_obj$hessian)))
  a_ <- parameters[,1]
  b_ <- parameters[,2]

  # Diagnostic QQ plot
  if (show_diagnostic) {
    q_ideal <- as.vector(quantile(rnbinom(length(d) * 1000, median(a_), median(b_)/(median(b_) + rep(d, 1000))), probs=seq(0, 1, 0.01)))
    q_true <- as.vector(quantile(x, probs=seq(0, 1, 0.01)))

    plot(q_ideal, q_true, xlim=c(0, max(q_true)), main="Gamma-Poisson QQ Plot", xlab="Theoretical Quantiles", ylab="Observed Quantiles")
    abline(coef=c(0, 1))
  }

  # Compute social differentiation and interaction rate.
  S <- 1/sqrt(a_)
  mu <- a_/b_

  # Compute sampling effort and network correlation.
  I <- mu * .hmean(d)
  rho <- (S * sqrt(I))/sqrt(1 + S^2 * I)

  # Calculate CI bounds.
  ci_lower <- 0.5 * (1 - ci)
  ci_upper <- 0.5 * (1 - ci) + ci

  # Standard error and quantiles of sampling effort and network correlation.
  S_se <- sd(S)
  rho_se <- sd(rho)
  mu_se <- sd(mu)
  S_q <- quantile(S, probs=c(ci_lower, 0.5, ci_upper))
  rho_q <- quantile(rho, probs=c(ci_lower, 0.5, ci_upper))
  mu_q <- quantile(mu, probs=c(ci_lower, 0.5, ci_upper))

  summary <- matrix(nrow=6, ncol=4)
  rownames(summary) <- c("Observed Social Differentiation", "Mean Interaction Rate", "Sampling Effort", "Est. Interaction Rate", "Est. Social Differentiation", "Est. Correlation")
  colnames(summary) <- c("Estimate", "SE", "Lower CI", "Upper CI")
  summary[1, 1] <- sd(x)/mean(x)
  summary[2, 1] <- mean(x/d)
  summary[3, 1] <- mean(x/d) * .hmean(d)
  summary[4, ] <- c(mu_q[2], mu_se, mu_q[1], mu_q[3])
  summary[5, ] <- c(S_q[2], S_se, S_q[1], S_q[3])
  summary[6, ] <- c(rho_q[2], rho_se, rho_q[1], rho_q[3])

  summary <- signif(summary, 3)
  summary
}

# Harmonic mean.
.hmean <- function(x) {
  length(x)/sum(1/x)
}
