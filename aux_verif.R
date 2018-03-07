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

predwidth <- function(ens, obs, quantiles = c(0.125,0.875), center = TRUE)  {
  ## obs needed only for compatibility with veriApply function {easyVerification}
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

## predwidth(ens, c(0.125, 0.875), TRUE)
## predwidth(ens, c(0.125, 0.875), FALSE)
## predwidth(ens, quantiles = c(0.875), TRUE)
## predwidth(ens, quantiles = c(0.875), FALSE)


## Diebold Mariano test wrapper function
## --> applies univariate Diebold-Mariano tests to arrays of arbitrary size
## --> last dimension of the array refers to the realizations

 dm.test.arr <- function(forc, ref) {
    tdim <- dim(forc)
    ndim <- length(tdim)

    forc1 <- apply(forc, 1:{ndim-1}, function(x) list(c(x)))
    ref1 <- apply(ref, 1:{ndim-1}, function(x) list(c(x)))
    
    forc2 <- lapply(forc1, function(x) x[[1]])
    ref2 <- lapply(ref1, function(x) x[[1]])
    
    f.dm.test <- function(x,y, alternative, h, power) {
      if(all(is.na(x)) | all(is.na(y))) return(NA)
      dm.test(x,y, alternative, h, power)$p.val
    }


    testres <- mapply(function(x,y) f.dm.test(x, y, "two.sided", h = 1, power = 1),
                    x = forc2, y = ref2)
    
    return(array(testres, dim=tdim[-ndim]))
    
  }

