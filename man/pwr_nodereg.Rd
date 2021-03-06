% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pwr_reg.R
\name{pwr_nodereg}
\alias{pwr_nodereg}
\title{Power analysis for nodal regression.}
\usage{
pwr_nodereg(
  nodes,
  effect,
  soc_diff,
  int_rate,
  samp_times,
  sig_level = 0.05,
  metric_fn = c("strength", "eigenvector", "closeness", "betweenness"),
  directed = FALSE,
  num_iters = 1000
)
}
\arguments{
\item{nodes}{Positive integer number of nodes.}

\item{effect}{Real value between 0 and 1 describing the effect size (correlation coefficient) between response and predictor in the true underlying network, before considering sampling.}

\item{soc_diff}{Positive real value descriing the social differentiation estimated using `net_cor`.}

\item{int_rate}{Positive real value describing interaction rate estimated using `net_cor`.}

\item{samp_times}{Real value describing mean sampling time per dyad or squares matrix of sampling times for each pair of nodes.}

\item{sig_level}{Real valued significance level of power analysis. Defaults to 0.05.}

\item{metric_fn}{Metric function used to compute node-level response. Choose from presets or provide a function that accepts igraph graphs as input and outputs a vector of the correct size.}

\item{directed}{TRUE/FALSE indicating if the network is directed.}

\item{num_iters}{Number of iterations used to calculate power.}
}
\value{
A list of three `n x n` matrices. `X` is the interaction observation
  matrix, `D` is the sampling time, and `A` is a matrix of the true
  interaction rates.
}
\description{
This function computes the power of nodal regression where a node-level
network metric (resonse, Y) is regressed against a nodal covariate
(predictor, X). Networks are ususally imperfectly sampled, which can affect
the reliability of power estimates. This function takes this into account
using social differentiation `soc.diff` and sampling effort `samp.eff`, as
calculated by `net.cor`.
}
\examples{
#' # Simulate observation data.
sim_data <- simulate_data_gp(20, 10, 2, 0.5)
X <- sim_data$X
D <- sim_data$D

# Estimate social differentiation and interaction rate.
net_cor(X, D) # Pretend this gives us exactly soc_diff = 2.0 and int_rate=0.5.

# Run node regression power analysis.
pwr_nodereg(20, 0.5, 2, 0.5, D)

}
