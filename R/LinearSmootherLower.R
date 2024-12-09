
LinearSmootherLower <- function(tvec, valvec, tighten = FALSE) {

  # Converted to R code by Michael Dymock 2024
  #LinearSmootherLower: smoothens the stopping boundary from a grid based
  #setting to a linearly smoothed, by taking the convex hull below the points
  # INPUTS:
    #   tvec: times at which the vector should be smoothed
  #   valvec: vector containing the discretized stopping values
  #   tighten: if false, then don't try to make time vector sparser. if true,
  #      then try to 'sparsen' the time vector
  # OUTPUTS:
  #   tout, valout: the smoothed time values and observation values
  #
  # behavior is undefined if valvec is not nonincreasing, and if tvec is not nondecreasing.
  #
  # 1) Error checking and setup

  tlen <- length(tvec)
  vlen <- length(valvec)

  if(tlen != vlen){
    warning("LinearSmootherLower: time and value vectors are not of the same dimension")
    return(list(tout = tvec, valout = valvec))
  }

  # 2) Find the vector of 'jumps', and the times they jump.

  tout <- tvec[tlen]
  valout <- valvec[tlen]

  for(i in (tlen-1):1){
    if(valvec[i] != valvec[i+1]){
      tout <- c(tvec[i], tout)
      valout <- c(valvec[i], valout)
    }
  }
  # if the first time step is not at a jump, then need to do something a bit special: here, we extrapolate linearly from preceding steps
  if(tout[1] != tvec[1]){
    i <- 1
    if(length(tout) > 2){
      ratio <- (valout[2] - valout[1])/(tout[2] - tout[1])
      newval <- valout[2] - ratio*(tout[2] - tvec[1])
      tout <- c(tvec[i], tout)
      valout <- c(newval, valout)
    } else {
      warning("LinearSmootherLower: potentially odd smoothing at t0")
      tout <- c(tvec[i], tout)
      valout <- c(valvec[i], valout)
    }
  }
  # SECOND: Go back and find the linear interpolations if smoothening is needed
  # cleanup: doublecheck if last entries are the same
  toutlen <- length(tout)
  if(toutlen > 2 & tout[toutlen-1] == tout[toutlen]){
    tout <- tout[1:(toutlen-1)]
    valout[toutlen-1] <- min(valout[(toutlen-1):toutlen])
    valout <- valout[1:(toutlen-1)]
  }

  if(!tighten){
    valout <- approx(tout, valout, tvec, method = "linear")$y
    tout <- tvec
  }

  return(list(tout = tout, valout = valout))
}