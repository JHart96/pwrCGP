#' Power analysis for nodal regression.
#'
#' This function computes the power of nodal regression where a node-level
#' network metric (resonse, Y) is regressed against a nodal covariate
#' (predictor, X). Networks are ususally imperfectly sampled, which can affect
#' the reliability of power estimates. This function takes this into account
#' using social differentiation `soc.diff` and sampling effort `samp.eff`, as
#' calculated by `net.cor`.
#'
#' @param nodes Positive integer number of nodes.
#' @param effect Real value between 0 and 1 describing the effect size (correlation coefficient) between response and predictor in the true underlying network, before considering sampling.
#' @param soc_diff Positive real value descriing the social differentiation estimated using `net_cor`.
#' @param int_rate Positive real value describing interaction rate estimated using `net_cor`.
#' @param samp_times Real value describing mean sampling time per dyad or squares matrix of sampling times for each pair of nodes.
#' @param sig_level Real valued significance level of power analysis. Defaults to 0.05.
#' @param metric_fn Metric function used to compute node-level response. Choose from presets or provide a function that accepts igraph graphs as input and outputs a vector of the correct size.
#' @param directed TRUE/FALSE indicating if the network is directed.
#' @param num_iters Number of iterations used to calculate power.
#'
#' @return A list of three `n x n` matrices. `X` is the interaction observation
#'   matrix, `D` is the sampling time, and `A` is a matrix of the true
#'   interaction rates.
#' @export
#'
#' @examples
#' #' # Simulate observation data.
#' sim_data <- simulate_data_gp(20, 10, 2, 0.5)
#' X <- sim_data$X
#' D <- sim_data$D
#'
#' # Estimate social differentiation and interaction rate.
#' net_cor(X, D) # Pretend this gives us exactly soc_diff = 2.0 and int_rate=0.5.
#'
#' # Run node regression power analysis.
#' pwr_nodereg(20, 0.5, 2, 0.5, D)
#'
pwr_nodereg <- function (nodes, effect, soc_diff, int_rate, samp_times, sig_level=0.05, metric_fn=c("strength", "eigenvector", "closeness", "betweenness"), directed=FALSE, num_iters=1000) {
  if (length(metric_fn) > 1) {
    metric_fn = "strength"
  }
  if (metric_fn == "strength") {
    metric_fn <- igraph::strength
  } else if (metric_fn == "eigenvector") {
    metric_fn <- function(g) igraph::eigen_centrality(g)$vector
  } else if (metric_fn == "closeness") {
    metric_fn <- function(g) igraph::closeness(g)
  } else if (metric_fn == "betweenness") {
    metric_fn <- function(g) igraph::betweenness(g)
  }
  if (is.null(metric_fn)) {
    metric_fn <- igraph::strength
  }

  if (length(samp_times) == 1) {
    sampling_times <- matrix(rpois(nodes^2, samp_times), nodes, nodes)
    if (directed) {
      sampling_times <- sampling_times * upper.tri(sampling_times) + sampling_times * lower.tri(sampling_times)
      mode = "directed"
    } else{
      sampling_times <- sampling_times * upper.tri(sampling_times)
      sampling_times <- sampling_times + t(sampling_times)
      mode = "undirected"
    }
  } else {
    sampling_times <- samp_times
  }

  a <- 1/soc_diff^2
  b <- 1/(soc_diff^2 * int_rate)

  pvals <- rep(0, 1)

  for (i in 1:num_iters) {
    alpha <- matrix(rgamma(nodes^2, a, b), nodes, nodes) # Actual rate.

    if (directed) {
      alpha <- alpha * upper.tri(alpha) + alpha * lower.tri(alpha)
      mode = "directed"
    } else{
      alpha <- alpha * upper.tri(alpha)
      alpha <- alpha + t(alpha)
      mode = "undirected"
    }

    net <- igraph::graph_from_adjacency_matrix(alpha, mode=mode, weighted=TRUE)
    metric <- metric_fn(net)

    r <- effect
    effect_size <- r * sqrt(1/(1 - r))/(sd(metric) * sqrt(r + 1))

    trait <- 1 + effect_size * metric + rnorm(nodes, sd=1)

    alpha <- as.vector(alpha)
    sampling_times <- as.vector(sampling_times)
    alpha_hat <- rpois(nodes^2, alpha*sampling_times)/sampling_times # Sampled rate.
    alpha_hat[is.nan(alpha_hat)] <- 0

    # Build network
    net_est <- igraph::graph_from_adjacency_matrix(matrix(alpha_hat, nodes, nodes), weighted=TRUE)
    metric_est <- metric_fn(net_est)

    pvals[i] <- summary(lm(trait ~ metric_est))$coefficients[2, 4]
  }
  power <- mean(pvals < 0.05)
  list(nodes=nodes, effect=effect, soc_diff=soc_diff, int_rate=int_rate, power=power)
}
