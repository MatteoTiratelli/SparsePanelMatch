#' Sparse_PanelMatch
#' 
#' Find matching units across irregular panel data:
#' 1. Match each treated observation (unit-time) to untreated observations occurring within a user-defined time window 
#' 2. Limit that matched set to control observations with exactly the same treatment history over the last n observations
#' 3. Within that matched set, refine further or weight using Propensity Scores, Covariate Balancing Propensity Scores or Mahalanobis distance
#' 
#' @param data The name of the dataset
#' @param time String of characters representing the name of the time variable. Should be in Date format.
#' @param unit String of characters representing the name of the unit variable. Should be a numeric variable.
#' @param treatment String of characters representing the name of the treatment variable. Should be numeric variable where 0 = untreated and 1 = treated.
#' @param outcome String of characters representing the name of the outcome variable. Should be a numeric variable.
#' @param treatment_lags Integer representing the number of treatment lags to use for the exact matching on treatment history.
#' @param outcome_leads Integer representing the number of outcome leads for which to calculate DiD estimands (each compared to Yt-1).
#' @param time_window_in_months Integer representing the size of the time window in months within which to search for control observations.
#' @param match_missing Whether or not to match missing values within treatment history (defaults to TRUE).
#' @param covs List of covariates to use in matching refinement procedures (should be a list of strings of characters).
#' @param qoi The wuantity of interest: either 'ATT' (Average Treatment effect on Treated) or 'ATC' (Average Treatment effect on Controls).
#' @param refinement_method The refinement method used to improve covariate balancing after the exact matching procedure. Can take one of: "none"; "CBPS.weight" (adjusts weights by Covariate Balancing Propensity Scores); "CBPS.match" (reduces the size of each control set to size_match by lowest CBPS); "ps.weight" and "ps.match" do the same for standard propensity scores; "mahalanobis" reduces set size to size_match by Mahalanobis distance. For more details see Imai, Kim & Wang (2018).
#' @param size_match Integer representing the maximum number of control observations within each refined set. For use with 'CBPS.match', 'ps.match' and 'mahalanobis'. Has no impact otherwise.
#' @param use_diagonal_covmat Whether or not mahalanobis distance calculation should use a regular covariance matrix (FALSE), or a diagonal matrix with only covariate variances (TRUE). Imai, Kim & Wang (2018) report that, "in many cases, setting this to TRUE can lead to better covariate balance, especially when there is high correlation between variables". Default is FALSE.
#' @return A Sparse_PanelMatch object
#' @examples
#' MatchedData <- Sparse_PanelMatch(data = CMP, time = "date", unit = "party", treatment = "wasingov", outcome = "sdper103", treatment_lags = 3, outcome_leads = 0, time_window_in_months = 60, match_missing = TRUE, covs = c("pervote", "lag_sd_rile"), qoi = "att", refinement_method = "mahalanobis", size_match = 5, use_diagonal_covmat = TRUE)


Sparse_PanelMatch <- function(data, time, unit, treatment, outcome,
                              treatment_lags, outcome_leads, time_window_in_months,
                              match_missing = TRUE,
                              covs, qoi = c("att", "atc"),
                              refinement_method = c("none","CBPS.weight", "CBPS.match", "ps.weight", "ps.match", "mahalanobis"),
                              size_match,
                              use_diagonal_covmat = FALSE) {
  
  ### Prepare dataset
  # Rename vars
  df1 <- data.table::setDT(data)
  df1 <- setnames(df1, treatment, "treatment")
  df1 <- setnames(df1, outcome, "outcome")
  df1 <- setnames(df1, time, "time")
  df1 <- setnames(df1, unit, "unit")
  df1 <- setnames(df1, covs, sapply(1:length(covs), function (x) paste0("control", x)))                                   
  
  if(class(df1$unit) != "numeric"){stop("Unit variable is not numeric. Please convert")}
  if(class(df1$treatment) != "numeric"){stop("Treatment variable is not numeric. Please convert")}
  if(class(df1$time) != "Date"){stop("Time variable is not in Date format. Please convert")}
  if(class(df1$outcome) != "numeric"){stop("Outcome variable is not numeric. Please convert")}
  if(sum(is.na(df1$treatment)) > 0){stop("Treatment variable contains missing values.")}
  if(sum(is.na(df1$unit)) > 0){stop("Unit variable contains missing values.")} 
  if(sum(is.na(df1$time)) > 0){stop("Time variable contains missing values.")}                                   
  
  
  # create outcome lag
  df1 <- df1[order(time), "lag_outcome" := shift(outcome, 1) , unit]
  
  # create treatment lags
  df1 <- df1[order(time), sapply(1:treatment_lags, function (x) paste0("lag_treatment_", x)) := shift(treatment, 1:treatment_lags) , unit]
  
  # create outcome leads
  if (outcome_leads > 0) {
    df1 <- df1[order(time), sapply(1:outcome_leads, function (x) paste0("lead_outcome_", x)) := shift(outcome, 1:outcome_leads, type=c("lead")) , unit]
  }
  
  # Deal with missing treatment history values
  lagindex <- grep('lag_treatment_', colnames(df1))                                 
  if (match_missing == FALSE) {
    treatmentslist <- sapply(1:treatment_lags, function (x) paste0("lag_treatment_", x))                               
    df1 %>% drop_na(all_of(treatmentslist)) -> df1
  }
  if (match_missing == TRUE) {
    df1 %>%
      mutate(across(starts_with("lag_treatment_"), ~replace_na(., 99))) -> df1
  }
  
  df1 <- as_tibble(df1)
  
  df1 %>% drop_na(outcome, lag_outcome) -> df1
  
  
  ### Exact matching on treatment history  
                             
  if(qoi == "att"){
    # For each unit, find the dates when Treatment = 1, but was 0 at previous observation
    units <- unique(df1$unit[df1$treatment == 1 & df1$lag_treatment_1 == 0])
  }
  
  if(qoi == "atc"){
    # For each unit, find the dates when Treatment = 0, but was 1 at previous observation
    units <- unique(df1$unit[df1$treatment == 0 & df1$lag_treatment_1 == 1])
  }
  
  find_controls <- function (x, y) { # For each treatment unit, find matching control unites
    
    # create list of treatment history for given treated observation
    list1 <- as.vector(df1[, lagindex][df1$unit == y & df1$time == x,])
    
    # Subset by matching treatment history
    index <- which(apply(df1[, lagindex], 1, function(x) all(x == list1)))
    temp <- df1[index,]
    
    # Refine subset by finding untreated observations in time window
    tw <- time_window_in_months/2
    x %m-% months(tw) -> start
    x %m+% months(tw) -> end
    temp <- temp[temp$time %in% seq.Date(start, end, by = "month") & temp$treatment == 0 & temp$unit != y,]
    
    if(nrow(temp) > 0) { # If there are matching control units, collate and weigh them
      
      control <- temp
      control$group <- paste0(y,' ',x)
      control$treatment <- 0
      control$weight <- 1/nrow(temp)
      
      treated <- df1[df1$time == x & df1$unit == y,]
      treated$group <- paste0(y,' ',x)
      treated$treatment <- 1
      treated$weight <- 1
      
      set <- bind_rows(treated, control)
      return(set)
    }
  }                             
  
  find_exact_matches <- function (z) { # Identify treatment units
    
    if(qoi == "att"){ # For each unit, find the dates when Treatment = 1, but was 0 at previous observation
      listofdates <- df1$time[df1$unit == z & df1$treatment == 1 & df1$lag_treatment_1 == 0]
    }
    
    if(qoi == "atc"){ # For each unit, find the dates when Treatment = 0, but was 1 at previous observation
      listofdates <- df1$time[df1$unit == z & df1$treatment == 0 & df1$lag_treatment_1 == 1]
    }
    
    lapply(listofdates, find_controls, y = z) # Find matching controls
  }
  
  sets <- lapply(units, find_exact_matches)                         
  output <- bind_rows(sets)                         
  cat("Exact matching complete. Starting refinement with", refinement_method,'\n')                         
  
                         
  ### Refinement
  
  if(refinement_method %in% c("CBPS.weight", "CBPS.match", "ps.weight", "ps.match")) {
    
    controlslist <- sapply(1:length(covs), function (x) paste0("control", x))                                 
    output %>% drop_na(all_of(controlslist)) -> output
    
    # CBPS and PS matching
    if(refinement_method == "CBPS.weight" | refinement_method == "CBPS.match") {
      quiet <- capture.output(fit0 <- CBPS::CBPS(reformulate(response = 'treatment', termlabels = sapply(1:length(covs), function (x) paste0("control", x))),
                          family = binomial(link = "logit"), data = output))
    }
    
    if(refinement_method == "ps.weight" | refinement_method == "ps.match") {
      fit0 <- glm(reformulate(response = 'treatment', termlabels = sapply(1:length(covs), function (x) paste0("control", x))),
                  family = binomial(link = "logit"), data = output)
    }
    
    # Calculate initial inverse propensity score weights (Hirano et al. 2003) [same method for CBPS]
    inverse_PS_weighting <- function (x, B) {
      xx <- cbind(1, as.matrix(x[, sapply(1:length(covs), function (y) paste0("control", y))]))
      x[, (ncol(x) + 1)] <- as.vector(1 - 1/(1+exp(xx %*% fit0$coefficients)))
      names(x)[ncol(x)] <- "ps"
      return(x)
    }
    
    sets <- split(output, f = output$group)
    sets_with_ps <- lapply(sets, inverse_PS_weighting, B = fit0$coefficients)
    
    
    # Calculate final normalised weights, or match by size_match
    
    if(refinement_method == "CBPS.weight" | refinement_method == "ps.weight") {
      adjust_weights <- function(set) {
        set <- arrange(set, treatment)
        control.ps.set <- set[set$treatment == 0,]
        if(nrow(control.ps.set) == 1) {
          set$weight <- c(1,1)}
        vec.ratio <- control.ps.set$ps / (1 - control.ps.set$ps) #just for clarity
        if(sum(vec.ratio) == 0) {
          set$weight <- rep(1 / nrow(control.ps.set), nrow(control.ps.set))
        }
        if(sum(vec.ratio) > 0 & nrow(control.ps.set) > 1) {
          set$weight <- c((vec.ratio)/sum(vec.ratio), 1)
        }
        return(set)
      }
      sets <- lapply(sets_with_ps, adjust_weights)
      output <- bind_rows(sets)
    }
    
    if(refinement_method == "CBPS.match" | refinement_method == "ps.match"){
      restrict_sets <- function(set, size_match) {
        if((nrow(set)-1) > size_match){
          treated.ps <- set$ps[set$treatment == 1]
          set$dist <- abs(treated.ps - set$ps)
          dist.to.beat <- max(utils::head(x = sort(set$dist), n = size_match+1))
          set <- set[set$dist <= dist.to.beat,]
          set$weight[set$treatment == 0] <- (1/size_match)
        }
        return(set)
      }
      sets <- lapply(sets_with_ps, restrict_sets, size_match = size_match)
      output <- bind_rows(sets)
    }
  }
  
  if(refinement_method == "mahalanobis") {
    
    controlslist <- sapply(1:length(covs), function (x) paste0("control", x))                                 
    output %>% drop_na(all_of(controlslist)) -> output
    sets <- split(output, f = output$group)
    
    # For each matched set, use unit ids to find ind all other observations of those units
    build_maha_sets <- function(set){
      if((nrow(set)-1) > size_match){
        listofunits <- unique(set$unit[set$treatment == 0])
        expandedset <- df1[df1$unit %in% listofunits,]
        expandedset <- add_row(expandedset, cbind(set[set$treatment == 1,][, sapply(1:length(covs), function (y) paste0("control", y))], unit = 999999999))
        return(expandedset)
      }
    }
    
    maha_sets <- lapply(sets, build_maha_sets)
    maha_sets <- maha_sets[lengths(maha_sets) > 0]
    
    # For each of the new Maha sets, calculate mahalanobis distance for each observation compared to the treated observation from original matched set, then take average for each unit
    maha_calculations <- function(set, size.match = size_match, use.diagonal.covmat = use_diagonal_covmat) {
      if(nrow(set) > (size.match+1)) {
        
        center.data <- set[set$unit == 999999999,][, sapply(1:length(covs), function (y) paste0("control", y))]
        set <- set[set$unit != 999999999,]
        cov.data <- set[, sapply(1:length(covs), function (y) paste0("control", y))]
        
        if(use.diagonal.covmat == TRUE) {
          cov.matrix <- diag(apply(cov.data, 2, var), ncol(cov.data), ncol(cov.data))
        }
        if (use.diagonal.covmat == FALSE) {
          cov.matrix <- cov(cov.data)
        }
        
        set$maha <- tryCatch({
          mahalanobis(x = as.matrix(cov.data), center = as.matrix(center.data), cov = as.matrix(cov.matrix))
        }, warning = function(w) {
          
        }, error = function(e) {
          cov.matrix <- cov(cov.data)
          mahalanobis(x = as.matrix(cov.data), center = as.matrix(center.data), cov = as.matrix(cov.matrix), inverted = TRUE)
        })
        set %>% group_by(unit) %>% summarise(maha = mean(maha)) -> set
        return(set)
      }
    }
    
    maha_sets_distance <- lapply(maha_sets, maha_calculations)
    
    # Merge to add mean maha distance to original sets
    listofdfs <- names(maha_sets_distance)
    for (i in 1:length(listofdfs)) {
      x <- maha_sets_distance[[listofdfs[[i]]]]
      y <- sets[[listofdfs[[i]]]]
      sets[[listofdfs[[i]]]] <- merge(x,y, by = 'unit', all = TRUE)
    }
    
    # Restrict size of each matched set to size_match
    restrict_sets_maha <- function(set, size.match = size_match) {
      treated_maha <- set[set$treatment == 1,]
      control_maha <- set[set$treatment == 0,]
      control_maha <- arrange(control_maha, maha)[1:size.match,]
      control_maha$weight <- (1/size.match)
      set <- rbind(treated_maha, control_maha)
      set$maha <- NULL
      return(set)
    }
    
    for (i in 1:length(listofdfs)) {
      sets[[listofdfs[[i]]]] <- restrict_sets_maha(sets[[listofdfs[[i]]]])
    }
    
    output <- bind_rows(sets)
  }
  
  # Generate output
  output <- list(summary = as_tibble(output),
                 qoi = qoi,
                 covs = covs,
                 treatment_lags = treatment_lags,
                 outcome_leads = outcome_leads,
                 time_window_in_months = time_window_in_months,
                 outcome = outcome,
                 refinement_method = refinement_method,
                 size_match = size_match)
  class(output) <- 'SparsePanelMatch'
  return(output)
}



summary.SparsePanelMatch <- function(object) {
  cat(" Matched DiD for Time-Series Cross-Sectional Data (Imai, Kim & Wang 2018)\n Method adapted by matching treated to untreated observations within",object$time_window_in_months,
      "month window\n Exact matching using treatment history over",
      object$treatment_lags,"periods\n Matches refined using",object$refinement_method,
      "with covariates:",paste(object$covs, collapse = ', '),'\n ',length(unique(object$summary$group)),'matched sets. Overall n = ',nrow(matches$summary),'\n\n')
  object$summary
}
