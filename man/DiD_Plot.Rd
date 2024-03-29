% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DiD_Plot.R
\name{DiD_Plot}
\alias{DiD_Plot}
\title{DiD_Plot}
\usage{
DiD_Plot(data, n_iterations = 1000, alpha = 0.05)
}
\arguments{
\item{data}{The matched data. Should be a SparsePanelMatch object.}

\item{n_iterations}{Integer representing number of iterations to use for bootstrapping (defaults to 1000).}

\item{alpha}{The alpha level to use when calculating the confidence intervals (defaults to 0.05)}
}
\value{
A ggplot object
}
\description{
Creates a classic two by two Difference-in-Difference plot of the means (with percentile bootstrapped confidence intervals) for the treated and untreated groups, pre and post treatment.
}
\examples{
MatchedData <- Sparse_PanelMatch(data = CMP, time = "date", unit = "party", treatment = "wasingov", outcome = "sdper103", treatment_lags = 3, outcome_leads = 0, time_window_in_months = 60, match_missing = TRUE, covs = c("pervote", "lag_sd_rile"), qoi = "att", refinement_method = "mahalanobis", size_match = 5, use_diagonal_covmat = TRUE)
DiD_Plot(MatchedData)
}
