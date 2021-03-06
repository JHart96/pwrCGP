% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pwr_elbow.R
\name{pwr_elbow}
\alias{pwr_elbow}
\title{Estimate point of diminishing returns using diminishing returns/elbow method}
\usage{
pwr_elbow(soc_diff, rho_max, I_max = 1000)
}
\arguments{
\item{soc_diff}{Positive real value describing the social differentiation as
estimated using `net_cor`.}

\item{rho_max}{Real value between 0 and 1 describing the maximum correlation
of practical relevance.}

\item{I_max}{Maximum feasible sampling effort (interactions seen per dyad
under even sampling). Leave it as is unless higher sampling is possible.
Setting this too high is not a major problem, setting it too low might lead
to unreliable results.}
}
\value{
List containing the optimal sampling effort and the correlation at
  that sampling effort.
}
\description{
This function estimates the point at which small increases in correlation
come at increasingly high costs of sampling effort. If correlation is
suitably high, this can be used to determine how large window sizes should be
when subsetting data. However a low required sampling effort will not fix
poor correlation, and this should be considered when carrying out analyses.
}
