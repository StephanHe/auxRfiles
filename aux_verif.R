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
      outp <- try(dm.test(x,y, alternative, h, power)$p.val, silent=TRUE)
      if(class(outp) == "try-error") return(NA)
      return(outp)
    }

    testres <- mapply(function(x,y) f.dm.test(x, y, "two.sided", h = 1, power = 1),
                    x = forc2, y = ref2)
    
    return(array(testres, dim=tdim[-ndim]))
    
  }


## Function that calculates bootstrap ci for the mean score difference when including
## different number of lags (1; 1,2; 1,2,3; 1,...,n), where 1 are the youngest, and n are the oldest runs
## of a lagged ensemble
## ens: matrix of hindcasts
## obs: reference vectors
## nlags: number of different initializations/lags to be considered (is 4 for UK Met Office GloSea5)
## ref.ind: reference indicices for calculation of the observed climatology
## prob: vector of quantile levels 


f_veriflagged <- function(ens, obs, nlags, scorefunct, ref.ind, prob = NULL) {

  out <- list()
  
  if(deparse(substitute(scorefunct)) == "EnsCrps") {
    linds <- sapply(seq(1, nmemb, by = nmemb/nlags), function(x) 1:{x+nmemb/nlags-1}, simplify = FALSE)
    score <- sapply(linds, function(ii) scorefunct(ens[,ii], obs))
 } else if(deparse(substitute(scorefunct)) == "EnsRps") {
   obs_prob <-  f_probforc(matrix(obs, ncol = 1), ref.ind = ref.ind, prob = prob)
   linds <- sapply(seq(1, nmemb, by = nmemb/nlags), function(x) 1:{x+nmemb/nlags-1}, simplify = FALSE)
   lforcprobs <- sapply(linds, function(ii) f_probforc(ens[,ii], ref.ind = ref.ind, prob = prob),
                        simplify = FALSE)
   score <- sapply(lforcprobs, function(x) scorefunct(x, obs_prob, format = "members"))   
   } else { stop("no valid scoring function defined") }


  out[[1]] <- apply(score, 2, mean)
  
  bs <- function(data, ind) {
    mean(data[ind,2] - data[ind,1])
  }

  
  pairlist <- combn(1:nlags, 2, simplify = FALSE)
   
  out[[2]] <- t(sapply(pairlist, function(ii) {
    bsres <- boot(data = score[,ii], statistic = bs, R = 1000)
    ci <- boot.ci(bsres, conf=0.95, type="basic", index = 1)
    c(ii, bsres$t0, ci$basic[4:5])
  }, simplify=TRUE))

  colnames(out[[2]]) <- c("lag_01", "lag_02", "meandiff", "ci95_low", "ci95_high")
  names(out) <- c("meanscores", "BS conf ints")
  out
}  


## auxliary function needed in f_veriflagged
f_probforc <- function(x, ref.ind, prob) {
  if(length(dim(x)) == 1) { 
    convert2prob(x, prob = prob, ref.ind = ref.ind)
  } else { ## for ensemble forecasts
    nens <- max(apply(x, 1, function(x) sum(!is.na(x))))
    count2prob(convert2prob(x[,1:nens], prob = prob, ref.ind = ref.ind), type = 4)
  }
}
