
## function that defines contrasts: 
f_contrasts <- function(nbin, type = NULL) {

##   if(nbin %% 2 == 0) stop("function only implemented for odd number of bins")
  if(!{type %in% c("linear", "ends", "v_shape", "u_shape")}) stop("please select a contrast type: 'linear', 'ends', 'v_shape', 'u_shape'")

  if(nbin %% 2 == 1) {
    ## odd number of bins
  
    if(type == "linear") {
      h <- (nbin - 1)/2
      a <- -h
      b <- 1
      contr0 <- seq(a, a + (nbin - 1) * b, by = b)
      const_normalize <- sqrt(sum(contr0^2))
      return(contr0 / const_normalize)
    }

    if(type == "ends") {
      h <- (nbin - 1)/2
      a <- (2*h - 1)
      b <- 2
      contr0 <- c(a, rep(-b, nbin - 2), a)
      const_normalize <- sqrt(sum(contr0^2))
      return(contr0 / const_normalize)    
    }

    if(type == "v_shape") {
      h <- (nbin - 1)/2
      a <- h^2
      b <- 2*h+1
      contr0 <- a - c(0, {1:{h-1}}*b, h*b, {{h-1}:1}*b, 0)
      const_normalize <- sqrt(sum(contr0^2))
      return(contr0 / const_normalize)       
    }
    
    if(type == "u_shape") {
      h <- (nbin - 1)/2
      b <- h * (h+1) / 3
      contr0 <- c(seq(h,1,by = -1)^2, 0, seq(1,h,by = 1)^2) - b
      const_normalize <- sqrt(sum(contr0^2))
      return(contr0 / const_normalize)
    }

  } else {
    ## even number of bins

    if(type == "linear") {
      h <- nbin/2
      a <- -(2 * h - 1)
      b <- 2
      contr0 <- seq(a, a + (nbin - 1) * b, by = b)
      const_normalize <- sqrt(sum(contr0^2))
      return(contr0 / const_normalize)
    }

    if(type == "ends") {
      h <- nbin/2
      a <- h - 1
      b <- 1
      contr0 <- c(a, rep(-b, nbin - 2), a)
      const_normalize <- sqrt(sum(contr0^2))
      return(contr0 / const_normalize)    
    }
    if(type == "v_shape") {
      h <- nbin/2
      a <- h-1
      b <- 2
      contr0 <- a - c({0:(h-1)}*b, {(h-1):0}*b)
      const_normalize <- sqrt(sum(contr0^2))
      return(contr0 / const_normalize)       
    }
    if(type == "u_shape") {
      h <- nbin/2
      b <- (4 * h^2 - 1)/3
      contr0 <- c(seq(2 * h - 1, 1, -2)^2, seq(1, 2 * h - 1, 2)^2) - b
      const_normalize <- sqrt(sum(contr0^2))
      return(contr0 / const_normalize)       
    }
}
  
 
} ## end of f_contrasts



## function that performs Chi squared tests for flat rank histograms:
f_chi2rank <- function(ens = NULL, obs = NULL, tranks = NULL, contrasttypes) {
  ## needs either pairs of ens and obs or tranks as input
  ## ens: matrix of ensemble forecasts
  ## obs: vector of corresponding observations
  ## tranks: vector of ranks of the observations

  ## set extra random seed for function in order to ensure reproducability
  old <- try(.Random.seed, silent = TRUE)
  if(class(old) == "try-error") {runif(1)
                                 old <- .Random.seed
                               }
  on.exit( { .Random.seed <<- old }) ## ensures that seed applies only within the function
  set.seed(0)
  
  stopifnot(length(contrasttypes) == 2)

  if(is.null(tranks)) {
    ## map ens to nbin number of bins:
    tranks <- apply(cbind(obs, ens), 1, function(x) rank(x, ties.method = "random")[1])
  }

  probs <- {0:(mm - 1) + 0.5} / mm
  cdf_intervals <- cbind(c(0, probs), c(probs, 1))

  cdfvals_at_obs <- sapply(1:nn, function(ii) runif(1, min = cdf_intervals[tranks[ii], 1], max = cdf_intervals[tranks[ii], 2]))

  ## set bin limits: 
  bin_limits <- t(replicate(nn, {{1:nbin}/nbin}[-nbin]))

  ## map cdfvals to the bins:
  binvals_at_obs<- apply(cbind(cdfvals_at_obs, bin_limits), 1, function(x) rank(x, ties.method = "min")[1])


  ## calculate x_i's
  obsperbin <-  sapply(1:nbin, function(x) sum(binvals_at_obs == x))
  
  expect_obsperbin <- nn / nbin

  x_vec <- (obsperbin - expect_obsperbin)/sqrt(expect_obsperbin)

  pvals <- rep(NA, 4)
  names(pvals) <- c("full", "1st contrast", "2nd contrast", "residual")
  
  T_full <- sum(x_vec^2)
  pvals[1] <- 1 - pchisq(T_full, df = nbin - 1)
  
  ## sum((obsperbin - expect_obsperbin)^2/expect_obsperbin) ## check if this is equal to T_full --> seems to be OK

  ## 1st contrast:
  (T_1 <- sum(f_contrasts(nbin, type = contrasttypes[1]) * x_vec)^2)
  pvals[2] <- 1 - pchisq(T_1, df = 1)


  ## 2nd contrast:
  (T_2 <- sum(f_contrasts(nbin, type = contrasttypes[2]) * x_vec)^2)
  pvals[3] <- 1 - pchisq(T_2, df = 1)

  
  T_resid <- T_full - T_1 - T_2
  pvals[4] <-  1 - pchisq(T_resid, df = nbin - 3)

  return(pvals)
  
} ## end of f_chi2rank
