
PsiNormUV <- function(zval, dof, optionaltcdf = NULL){

  # (c) 2004 Stephen E. Chick, all rights reserved
  # This file is for input to matlab, and does calculations
  # to support the 'selecting the best system' paper in the
  # bayesian environment.
  #
  # Converted to R code by Michael Dymock 2024
  #
  # Assumes unknown sampling variance in loss function.
  #
  # Code is 'as is' and no guarantees for correctness.
  #
  # INPUTS:
    # zval : test statistic values (vector)
  # dof : degrees of freedom.
  # optionalzcdf : optional argument, if passed, it is assumed to be a
  # precomputed version of tcdf(zval). Typically it will not be passed.
  #
  # OUTPUTS:
    # rval : vector of outputs
  #\frac{\nu + s^2}{\nu-1} \phi_\nu(s) - s (1- \Phi_\nu(s) )
  #
  # For BIG values of zval, one might try to use formula
  # based on Mills Ratio and Soms (1976) as cited in Soms (Jasa, 1980, vol 75, number 370, page 438, equation 2.1), or as in the 2001 Machine
  # Learning article about student distribution mill's ratio approximations.
#
# NOTE: THIS APPROX IS TO BE VALIDATED.

  tpdfvals <- dt(zval, dof)

  if(is.null(optionaltcdf)){
    rval <- (dof + zval^2)*tpdfvals/(dof - 1) - zval*pt(-zval, dof)
  } else {
    rval <- tpdfvals - zval*(1 - optionaltcdf)
  }

  return(rval)
}