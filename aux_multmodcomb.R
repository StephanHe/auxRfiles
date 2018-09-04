########################################################################################
#
# Auxiliary R functions needed for multi-model combination
#
#
########################################################################################

library(R.utils)
library(nor1mix)

f_mat2list <- function(x) lapply(seq_len(ncol(x)), function(jj) x[,jj])


f_multmod2 <- function(mu, sigma, probs, weights) {
  if(anyNA(mu)) return(rep(NA, length(probs)))

  
  ## handle cases where the standard root method gets stuck
  out <- try(withTimeout(qnorMix(probs, norMix(mu = mu, sigma = sigma, w = weights),
                                               tol = .Machine$double.eps^0.19), timeout = 1,
                         onTimeout="error"), silent = TRUE)
  if(class(out) == "try-error") {
    out <- try(qnorMix(probs, norMix(mu = mu, sigma = sigma, w = weights),
                                     tol = .Machine$double.eps^0.19, method = "eachRoot"),
               silent = TRUE)
    if(class(out) == "try-error") { 
      out <- rep(NA, length(probs))
    }
  }
  return(out)
}



f_multmod2_sst <- function(mu, sigma, probs, weights, nn) {
  ## variant of multmod2 based on resampling that can  handle sst without getting stuck due to
  ## numerical issues
  if(anyNA(mu)) return(rep(NA, length(probs)))

  ## nn <- 1e4
  
  tsample <- NULL
  for(mm in 1:length(mu)) {       
    tsample <- c(tsample, rnorm(round(nn * weights[mm], 0), mean = mu[mm], sd = sigma[mm]))
  }
  quantile(tsample, probs = probs)
}



f_multmod <- function(input, weights, rsampling = FALSE, nresamp = 1e4) {
  ## input: list of input ensembles matrices
  ## weights: matrix of weights needed for the multimodel combination


  nyears <- dim(input[[1]])[1]
  minsize <- min(sapply(input, function(x) sum(apply(is.na(x), 2, sum) < nyears)))
  
  if(!minsize) return(NA) ## returns NA if minsize is 0, i.e. a ensemble w/o any values
  
  nmemb <- minsize * length(input) ## ensures that output ensemble is of the same size as the simple multimodel
  
  probs <- {0:(nmemb - 1) + 0.5} / nmemb
  probslist <- replicate(nyears, probs, simplify=FALSE)
  
  MU1 <- matrix(unlist(sapply(input, function(x) apply(x, 1, mean, na.rm = T),
                              simplify = FALSE)),  ncol = nyears, byrow = TRUE)
  SIGMA1 <- matrix(unlist(sapply(input, function(x) apply(x, 1, sd, na.rm = T),
                                 simplify = FALSE)), ncol = nyears, byrow = TRUE)


  SIGMA1[SIGMA1 < 1e-6] <- 1e-6 ## needed in order to avoid zero variances (and as a consequence hanging up f_multmod2)
  if(rsampling) {   
    return(t(mapply(f_multmod2_sst, f_mat2list(MU1),
                    f_mat2list(SIGMA1), probslist, f_mat2list(weights),
                    MoreArgs=list(nn=nresamp))))
  
  } else {
    return(t(mapply(f_multmod2, f_mat2list(MU1), f_mat2list(SIGMA1), probslist, f_mat2list(weights))))
  }
  
}



## ensemble mean difference:
f_md <- function(x) {

  if(all(is.na(x))) return(NA) 

  Nna <- sum(is.na(x))
  
  sum(abs(outer(x, x, FUN="-")), na.rm = T) / (length(x)-Nna)^2
}
