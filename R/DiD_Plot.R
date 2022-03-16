#' DiD_Plot
#' 
#' Creates a classic two by two Difference-in-Difference plot of the means (with percentile bootstrapped confidence intervals) for the treated and untreated groups, pre and post treatment.
#' 
#' @param data The matched data. Should be a SparsePanelMatch object.
#' @param n_iterations Integer representing number of iterations to use for bootstrapping (defaults to 1000).
#' @param alpha The alpha level to use when calculating the confidence intervals (defaults to 0.05)
#' @return A ggplot object
#' @examples
#' MatchedData <- Sparse_PanelMatch(data = CMP, time = "date", unit = "party", treatment = "wasingov", outcome = "sdper103", treatment_lags = 3, outcome_leads = 0, time_window_in_months = 60, match_missing = TRUE, covs = c("pervote", "lag_sd_rile"), qoi = "att", refinement_method = "mahalanobis", size_match = 5, use_diagonal_covmat = TRUE)
#' DiD_Plot(MatchedData)


DiD_Plot <- function(data, n_iterations = 1000, alpha = 0.05) {
  if(class(data) != "SparsePanelMatch"){stop("Data is not SparsePanelMatch object")}
  
  output <- data
  
  # Calculate qois
  df1 <- data.table::setDT(output$summary)
  plotdata <- data.frame(Outcome = c(weighted.mean(df1$lag_outcome[df1$treatment == 1], w = df1$weight[df1$treatment == 1]),
                                     weighted.mean(df1$lag_outcome[df1$treatment == 0], w = df1$weight[df1$treatment == 0]),
                                     weighted.mean(df1$outcome[df1$treatment == 1], w = df1$weight[df1$treatment == 1]),
                                     weighted.mean(df1$outcome[df1$treatment == 0], w = df1$weight[df1$treatment == 0])))
  
  boots <- matrix(NA, nrow = n_iterations, ncol = 4)
  for (k in 1:n_iterations) {
    clusters <- unique(df1[["unit"]])
    units <- sample(clusters, size = length(clusters), replace = T)
    df.bs <- lapply(units, function(x) which(df1[, "unit"] == x))
    d.sub1 <- df1[unlist(df.bs), ]
    boots[k,] <- c(weighted.mean(d.sub1$lag_outcome[d.sub1$treatment == 1], w = d.sub1$weight[d.sub1$treatment == 1]),
                   weighted.mean(d.sub1$lag_outcome[d.sub1$treatment == 0], w = d.sub1$weight[d.sub1$treatment == 0]),
                   weighted.mean(d.sub1$outcome[d.sub1$treatment == 1], w = d.sub1$weight[d.sub1$treatment == 1]),
                   weighted.mean(d.sub1$outcome[d.sub1$treatment == 0], w = d.sub1$weight[d.sub1$treatment == 0]))
  }
  
  plotdata$bootstrap_low <- apply(boots, 2, quantile, probs = alpha/2)
  plotdata$bootstrap_high <- apply(boots, 2, quantile, probs = (1 - (alpha/2)))
  plotdata$Time <- c('Pre','Pre','Post','Post')
  plotdata$Treated <- c('Treated','Untreated','Treated','Untreated')
  
  ggplot2::ggplot(plotdata, aes(x=Time, y=Outcome, ymin = bootstrap_low, ymax = bootstrap_high, colour = Treated)) + 
    geom_pointrange(position = position_dodge(width=0.05)) + ylab('Outcome') + 
    scale_x_discrete(limits = rev) + xlab(NULL) +
    geom_line(aes(group = Treated), position = position_dodge(width=0.05)) +
    labs(title = "Difference-in-Difference plot with bootstrapped standard errors") +
    theme_bw() + 
    theme(
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background= element_blank(),
      legend.box.background = element_blank(),
      legend.title=element_blank(),
      legend.position = "bottom",
      strip.background = element_blank())
}
