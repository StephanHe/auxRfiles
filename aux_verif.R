## Calculation of PIT values of ensemble forecasts
## Inputs:
## - NxR matrix of ensemble forecasts
## - vector of observations of size N
## Outputs:
## - vector of size N of PIT values

pit <- function(ens, obs) {
  storeseed <- .Random.seed
  on.exit({ .Random.seed <<- storeseed })
  set.seed(pi)
  N <- dim(ens)[1]
  R <- dim(ens)[2] 
  obsranks <- apply(cbind(obs, ens), 1, function(x) if(anyNA(x)) return(NA) else rank(x, ties.method = "random")[1])
  ## lower <- pmax(obsranks - 1.5,0) / R
  ## upper <- pmin(obsranks - 0.5,R) / R
  lower <- (obsranks - 1) / (R+1)
  upper <- (obsranks) / (R+1)
  runif(N, lower, upper)
}


## Calculation of widths of prediction intervals from ensemble forecasts
## --> used to evaluate sharpness
## Inputs:
## - NxR matrix of ensemble forecasts
## - for centered two-sided prediction intervals: the upper and lower quantile to be considered
## - for one-sided prediction intervals (applies to bounded variables like precipitation or wind speed): the upper
##   quantile to be considered (the lower is zero)
## Outputs:
## - vector of size N of prediction interval widths


calcpredwidth <- function(ens, CDFVals, quantiles)  if(anyNA(ens)) return(NA) else diff(approx(CDFVals, sort(ens), xout = quantiles)$y) ## auxiliary function for predwidth

predwidth <- function(ens, quantiles, center = TRUE)  {
  R <- dim(ens)[2]
  CDFVals <- {1:R}/{R+1}
  if(center) {
      stopifnot(length(quantiles) == 2)
      apply(ens, 1, calcpredwidth, CDFVals = CDFVals, quantiles = quantiles)
  } else {
    stopifnot(length(quantiles) == 1)
    apply(ens, 1, function(x) if(anyNA(x)) return(NA) else quantile(x, probs = quantiles))  
    }
}
