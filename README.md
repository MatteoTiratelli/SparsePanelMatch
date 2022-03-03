# SparsePanelMatch
Applies Imai, Kim &amp; Wang (2018) ["Matching Methods for Causal Inference with Time-Series Cross-Sectional Data"](https://imai.fas.harvard.edu/research/tscs.html) to sparse panel data where observations are irregular and there are many missing values.

**Generalised procedure:**
1. Match each treated observation (unit-time) with control observations occurring within a user-defined time window (this is the only significant difference to the original method, where it is assumed that the panel data is well-ordered and regular, meaning that each observation is matched with every other observation at that year/month/date)
2. Within that matched set, find control observations with exactly the same treatment history over the last n period (e.g. over the last 3 election cycles)
3. Within that matched set, refine further using either Propensity Scores, Covariate Balancing Propensity Scores or Mahalanobis distance
4. Calculate a Difference-in-Difference estimator across the matched set (with bootstrapped standard errors)

**Installation:**
``` r
devtools::install_github("https://github.com/MatteoTiratelli/SparsePanelMatch")
```

## Motivating example
[Imai, Kim &amp; Wang (2018)](https://imai.fas.harvard.edu/research/tscs.html) have developed a matching procedure to facilitate causal inference with time-series cross-sectional data. Although their approach is more general, their software package ([PanelMatch](https://github.com/insongkim/PanelMatch)) only works with regular panel data where there are repeated observations of each unit at identical and equally spaced moments in time. A canonical example might be repeated observations of countries each year.
There are times, however, when observations are irregular by design. For example, elections across countries rarely coincide, particularly when we look at more granular measures of time. Data of this sort might look like the graph below, with irregular observations and much missing data.

![Graph showing irregular panel data](https://github.com/MatteoTiratelli/matteotiratelli.github.io/raw/master/Files/Irregular.png)

## Matching procedure

**Exact matching:**
1. Match each treated observation (unit-time) with control observations occurring within a user-defined time window (this is the only significant difference to the original method, where it is assumed that the panel data is well-ordered and regular, meaning that each observation is matched with every other observation at that year/month/date).
2. Within that matched set, find control observations with exactly the same treatment history over the last n period (e.g. over the last 3 election cycles).

![Matching procedure](https://github.com/MatteoTiratelli/matteotiratelli.github.io/raw/master/Files/matching.png)

**Refinement:**

Within that matched set, refine further using either Propensity Scores, Covariate Balancing Propensity Scores or Mahalanobis distance. This can improve covariate balancing.
``` r
  Sparse_PanelMatch(data = cmp, time = "date", unit = "party", 
  treatment = "wasingov", outcome = "sdper103", 
  treatment_lags = 3, outcome_leads = 2, 
  time_window_in_months = 60, match_missing = TRUE, 
  covs = c("pervote", "lag_sd_rile"), qoi = "att", 
  refinement_method = "CBPS.weight", size_match = 5, 
  use_diagonal_covmat = TRUE)
```
## Estimation procedure
In order to calculate the Average Treatment effect on Treated/Control, we calculate a Difference-in-Difference estimator. Within each matched set, we compare the difference in the treated unit with the (weighted) mean difference in control units. These individual Difference-in-Difference scores are then averaged across all of the available matched sets. Standard errors are calculated by block bootstrapping (resampling across units).

``` r
# sum((Yt - Yt-1) - mean(Y't - Y't-1)) / number_matched_sets
plot(Sparse_PanelEstimate(data = matches, n_iterations = 1000, alpha = 0.05))
```

![Plot of effects over time](https://github.com/MatteoTiratelli/matteotiratelli.github.io/raw/master/Files/plot_zoom_png.png)
