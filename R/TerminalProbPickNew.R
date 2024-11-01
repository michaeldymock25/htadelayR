
TerminalProbPickNew <- function(valuevec, predvarsample, nochangeallowed = FALSE){

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
  # predvarsample: a scalar with the predictive variance of the posterior mean
  # nochangeallowed: optional, SHOULD DEFAULT TO FALSE
  #
  # OUTPUTS:
    # rval : vector of outputs

  # if no change is allowed, then it is equivalent to saying that there is no variance in posterior decision, or pred variance is 0
  if(nochangeallowed) predvarsample <- 0
  if(predvarsample > 0){
    rval <- pnorm(valuevec/sqrt(predvarsample))
  } else if(predvarsample == 0){  # this line indicates that if no change is allowed, the decision can be implemented without waiting for the delay of tau patients
    rval <- as.numeric(valuevec >= 0)
  } else {
    stop("TerminalProbPickNew: called with negative predvarsample")
  }

  return(rval)
}