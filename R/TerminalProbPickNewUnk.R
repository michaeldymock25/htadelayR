
TerminalProbPickNewUnk <- function(valuevec, predvarsample, nochangeallowed = FALSE, dof = NULL){

  # (c) 2014 Chick, Forster, Pertile
  # Converted to R code by Michael Dymock 2024
  # This file is for input to R, and does calculations
  # to support the 'selecting the best system' paper in the
  # bayesian environment.
  # Call if variance is unknown
  #
  # Code is 'as is' and no guarantees for correctness.
  #
  # INPUTS:
    # valuevec : vector of means of posterior distribution
  # predvarsample: a scalar with the predictive variance of the posterior mean
  # nochangeallowed: optional, SHOULD DEFAULT TO FALSE
  # dof: degrees of freedom for relevant t distribution
  #
  # OUTPUTS:
    # rval : vector of outputs
  #

  if(is.null(dof)){
    rval <- 0
    warning("TerminalProbPickNew: not enough arguments")
  }

  # if no change is allowed, then it is equivalent to saying that there is no variance in posterior decision, or pred variance is 0
  if(nochangeallowed) predvarsample <- 0
  if(predvarsample > 0){
    zvec <- valuevec/sqrt(predvarsample)
    rval <- pt(zvec, dof)
  } else if(predvarsample == 0){  # this line indicates that if no change is allowed, the decision can be implemented without waiting for the delay of tau patients
    rval <- as.numeric(valuevec >= 0)
  } else {
    rval <- 0
    stop("TerminalProbPickNewUnk: called with negative predvarsample")
  }

  return(rval)
}