#' Sparse_PanelEstimate 
#' 
#' Estimate DiD estimator for Sparse_PanelMatch object with bootstrapped standard errors
#'
#' @param data A SparsePanelMatch object.
#' @param n_iterations The number of iterations to use for bootstrapping (defaults to 1000).
#' @param alpha The alpha level to use when calculating the confidence intervals (defaults to 0.05).
#' @return A SparsePanelEstimate object.
#' @examples
#' MatchedData <- Sparse_PanelMatch(data = CMP, time = "date", unit = "party", treatment = "wasingov", outcome = "sdper103", treatment_lags = 3, outcome_leads = 0, time_window_in_months = 60, match_missing = TRUE, covs = c("pervote", "lag_sd_rile"), qoi = "att", refinement_method = "mahalanobis", size_match = 5, use_diagonal_covmat = TRUE)
#' Sparse_PanelEstimate(data = MatchedData, n_iterations = 1000, alpha = 0.05)
Sparse_PanelEstimate <- function(data, n_iterations = 1000, alpha = 0.05) {
  
  if(class(data) != "SparsePanelMatch"){stop("Data is not SparsePanelMatch object")}
  
  output <- data
  
  # Calculate qois
  df1 <- data.table::as.data.table(output$summary)
  df1$ideological_change <- (df1$outcome - df1$lag_outcome)
  df1 <- df1[treatment==0, weight := (weight*-1)]
  Estimate <- tibble(coefs = sum(df1$ideological_change*df1$weight)/length(unique(df1$group)))
  Estimate$lead <- 't+0'
  
  if (output$outcome_leads > 0){
    Estimate <- bind_rows(Estimate, tibble(coefs = sapply(1:output$outcome_leads, function(x) sum(sapply(sapply(1:output$outcome_leads, function (x) paste0('lead_outcome_',x)), function(x) (df1[[x]] - df1$lag_outcome))[,x]*df1$weight)/length(unique(df1$group))),
                                           lead = sapply(1:output$outcome_leads, function (x) paste0('t+',x))))
    Estimate <- Estimate[,c(2,1)]
  }
  cat("Coefficients calculated. Beginning bootstrapping\n")                         
  # Bootstrap SD
  boots <- matrix(NA, nrow = n_iterations, ncol = (output$outcome_leads+1))
  for (k in 1:n_iterations) {
    clusters <- unique(df1[['unit']])
    units <- sample(clusters, size = length(clusters), replace = T)
    df.bs <- lapply(units, function(x) which(df1[,'unit'] == x)) # creates index of where in main dataset they match each element in 'units'
    d.sub1 <- df1[unlist(df.bs),] # take those indexes to create new dataset
    boots[k,1] <- sum((d.sub1$outcome - d.sub1$lag_outcome)*d.sub1$weight)/length(unique(d.sub1$group))
    if (output$outcome_leads > 0){
      boots[k,2:ncol(boots)] <- sapply(1:output$outcome_leads, function(x) sum(sapply(sapply(1:output$outcome_leads, function (x) paste0('lead_outcome_',x)), function(x) (d.sub1[[x]] - d.sub1$lag_outcome))[,x]*d.sub1$weight)/length(unique(d.sub1$group)))
    }
  }
  
  Estimate$bootstrap_coefs <- apply(boots, 2, function (x) mean(x))
  Estimate$bootstrap_sd <- apply(boots, 2, function (x) sd(x))
  Estimate$bootstrap_low <-  apply(boots, 2, quantile, probs = alpha/2) # percentile confidence interval
  Estimate$bootstrap_high <-  apply(boots, 2, quantile, probs = (1-(alpha/2))) # percentile confidence interval
  Estimate$bootstrap_coefs_BC <- 2*Estimate$coefs - colMeans(boots) # bias corrected point estimate (Efron & Tibshirani 1993 p138)
  Estimate$bootstrap_low_BC <- apply((2*matrix(nrow = n_iterations, ncol = length(Estimate$coefs), Estimate$coefs, byrow = TRUE) - boots), 2, quantile,
                                     probs = alpha/2)
  Estimate$bootstrap_high_BC <- apply((2*matrix(nrow = n_iterations, ncol = length(Estimate$coefs), Estimate$coefs, byrow = TRUE) - boots), 2, quantile,
                                      probs = (1-(alpha/2)))# bias corrected confidence interval
  output <- list(summary = Estimate,
                 qoi = output$qoi,
                 covs = output$covs,
                 treatment_lags = output$treatment_lags,
                 outcome_leads = output$outcome_leads,
                 refinement_method = output$refinement_method,
                 time_window_in_months = output$time_window_in_months,
                 outcome = output$outcome,
                 n_iterations = n_iterations,
                 alpha = alpha)
  class(output) <- 'SparsePanelEstimate'
  return(output)
}


summary.SparsePanelEstimate <- function(object) {
  cat(" Matched DiD estimate of",toupper(object$qoi),'with refinement by',object$refinement_method,'\n With bootstrapped standard errors (using',object$n_iterations,"iterations, alpha =",object$alpha,")\n Bias corrected confidence intervals also available.\n\n")
  print(object$summary %>% select(lead, coefs, bootstrap_sd, bootstrap_low, bootstrap_high))
}


plot.SparsePanelEstimate <- function(object) {
  plotdata <- object$summary
  ylim <- c(min(min(plotdata$bootstrap_low),min(plotdata$bootstrap_low_BC)), max(max(plotdata$bootstrap_high),max(plotdata$bootstrap_high_BC)))
    graphics::plot(x = 1:(nrow(plotdata)),y = plotdata$coefs, frame = TRUE, pch = 16, cex = 1.5,
                   xaxt = "n", ylab = paste0("Estimated ",toupper(object$qoi)), xlab = "Time", ylim = ylim,
                   main = "Estimated Effects of Treatment Over Time")
    graphics::axis(side = 1, at = 1:nrow(plotdata), labels = plotdata$lead)
    graphics::segments(1:(nrow(plotdata)), plotdata$bootstrap_low, 1:(nrow(plotdata)), plotdata$bootstrap_high)
    graphics::abline(h = 0, lty = "dashed")
}
