
PsiNorm <- function(zval, optionalzcdf = NULL){

  # (c) 2004 Stephen E. Chick, all rights reserved
  # This file is for input to matlab, and does calculations
  # to support the 'selecting the best system' paper in the
  # bayesian environment.
  # Edited 2014 to include improvement in computation value from Mills ratio
  #
  # Converted to R code by Michael Dymock 2024
  #
  # Code is 'as is' and no guarantees for correctness.
  #
  # INPUTS:
    # zval : test statistic values (vector)
  # optionalzcdf : optional argument, if passed, it is assumed to be a
  # precomputed version of normcdf(zval). Typically it will not be passed.
  #
  # OUTPUTS:
    # rval : vector of outputs
  #

  ZVALLIM <- 10
  zpdf <- dnorm(zval)
  zbig <- zval > ZVALLIM
  zsmall <- zval <= ZVALLIM
  rbig <- zbig*(zpdf/(zval^2 + 1)) # uses Mills ratio for improving stability when z is big

  if(is.null(optionalzcdf)){
    rsmall <- zsmall*(zpdf - zval*pnorm(-zval))
  } else {
    rsmall <- zsmall*(zpdf - zval*(1 - optionalzcdf))
  }
  rval <- rbig + rsmall

  return(rval)
}

# following naive version is mathematically correct but computationally unstable for large |zval|
# rval <- dnorm(zval) - zval*pnorm(- zval)
# rval <- max(rval,0)