
TerminalRewardFunction <- function(valuevec, discountfactor, predvarsample = NULL, nochangeallowed = FALSE){

  # (c) 2014 Chick, Forster, Pertile
  # This file is for input to R, and does calculations
  # to support the 'selecting the best system' paper in the
  # bayesian environment.
  # Edited 2014 to include improvement in computation value from Mills ratio
  # Converted to R code by Michael Dymock 2024
  #
  # Code is 'as is' and no guarantees for correctness.
  #
  # INPUTS:
    # valuevec : vector of means of posterior distribution
    # discountfactor: discountfactor to apply to reward after the patients are seen
    # predvarsample: a scalar with the predictive variance of the posterior mean
    # nochangeallowed: optional, SHOULD DEFAULT TO FALSE
  #
  # OUTPUTS:
    # rval : vector of outputs

  if(is.null(predvarsample)) stop("TerminalRewardFunction: not enough arguments")
  if(nochangeallowed) predvarsample <- 0 # if no change is allowed, then it is equivalent to saying that there is no variance in posterior decision, or pred variance is 0

  if(predvarsample > 0){
    zvec <- -valuevec/sqrt(predvarsample)
    rval <- discountfactor^(1 - nochangeallowed)*sqrt(predvarsample)*PsiNorm(zvec)
  } else if(predvarsample == 0){ # this line indicates that if no change is allowed, the decision can be implemented without waiting for the delay of tau patients
    rval <- (discountfactor^(1 - nochangeallowed))*pmax(valuevec, 0)
  } else {
    rval <- 0*valuevec
    warning("TerminalRewardFunction: called with negative predvarsample")
  }

  return(rval)
}