########################################################################################
#
# Auxiliary R functions needed for multi-model combination
#
#
########################################################################################

library(nor1mix)

f_mat2list <- function(x) lapply(seq_len(ncol(x)), function(jj) x[,jj])

f_multmod2 <- function(mu, sigma, probs, weights) {
  if(anyNA(mu)) return(rep(NA, length(probs)))
  qnorMix(probs, norMix(mu = mu, sigma = sigma, w = weights))
}


f_multmod <- function(input, weights) {
  ## input: list of input ensembles matrices
  ## weights: matrix of weights needed for the multimodel combination


  nyears <- dim(input[[1]])[1]
  minsize <- min(sapply(input, function(x) sum(apply(is.na(x), 2, sum) < nyears)))
  nmemb <- minsize * length(input) ## ensures that output ensemble is of the same size as the simple multimodel
  
  probs <- {0:(nmemb - 1) + 0.5} / nmemb
  probslist <- replicate(nyears, probs, simplify=FALSE)
  
  MU1 <- matrix(unlist(sapply(input, function(x) apply(x, 1, mean, na.rm = T),
                              simplify = FALSE)),  ncol = nyears, byrow = TRUE)
  SIGMA1 <- matrix(unlist(sapply(input, function(x) apply(x, 1, sd, na.rm = T),
                                 simplify = FALSE)), ncol = nyears, byrow = TRUE)

  return(t(mapply(f_multmod2, f_mat2list(MU1), f_mat2list(SIGMA1), probslist, f_mat2list(weights))))
  
}



## ensemble mean difference:
f_md <- function(x) {

  if(all(is.na(x))) return(NA) 

  Nna <- sum(is.na(x))
  
  sum(abs(outer(x, x, FUN="-")), na.rm = T) / (length(x)-Nna)^2
}
