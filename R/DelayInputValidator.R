
#' @title DelayInputValidator
#' @description Checks the basic and advanced data structures for validity
#' @param basic A list containing basic input parameters created by the DelayInputConstructor function
#' @param advanced A list containing advanced input parameters (e.g., computational settings) created by the DelayInputConstructor function
#' @return A list containing validated basic and advanced data structures and an error list containing messages if any errors were found
#' @rdname DelayInputValidator
#' @export
DelayInputValidator <- function(basic, advanced){

  routinename <- "DelayInputValidator"
  errorlist <- c()
  MINBIGX <- 5  # mimimum number of dmu per sigma for numerical analysis
  rval <- TRUE

  defaults <- DelayInputConstructor()  # find the default values used in the construtor function.
  defaultbasic <- defaults$basic
  defaultadvanced <- defaults$advanced

  if(basic$sigma <= 0){
    basic$sigma <- defaultbasic$sigma # standard deviation of sample differences (between people in two trial arms)
    errorlist <- c(errorlist, sprintf("#s: sampled std dev was negative: sigma reset to #f", routinename, basic$sigma))
    rval <- FALSE
  }

  if(basic$c < 0){
    basic$c <- defaultbasic$c # cost per sample
    errorlist <- c(errorlist, sprintf("#s: negative sampling cost, c reset to #f", routinename, basic$c))
    rval <- FALSE
  }

  if(basic$theta > 1 || basic$theta <= 0){
    tmpval <- basic$theta
    basic$theta <- defaultbasic$theta # discount factor
    errorlist <- c(errorlist, sprintf("#s: invalid discount factor (#f), theta reset to #f", routinename, tmpval, basic$theta))
    rval <- FALSE
  }

  if(basic$PPatients <= 0){
    basic$PPatients <- defaultbasic$PPatients  # number of patients covered post-trial
    errorlist <- c(errorlist, sprintf("#s: PPatients should be positive, reset to #f", routinename, basic$PPatients))
    rval <- FALSE
  }

  if(basic$ICost < 0){
    # basic$ICost <- defaultbasic$ICost  # fixed cost if adoption of new technology is made
    errorlist <- c(errorlist, sprintf("#s: normally ICost is not negative - tread with care", routinename))
  }

  if(basic$tau < 0){
    basic$tau <- defaultbasic$tau # delay (in number of patients treated) until data is observed
    errorlist <- c(errorlist, sprintf("#s: tau must be nonnegative, and was reset to #f", routinename, basic$tau))
    rval <- FALSE
  }

  if(basic$TMax <= basic$tau){
    basic$TMax <- 2*basic$tau # maximum number of samples (patient pairs) to be allowed during trial
    errorlist <- c(errorlist, sprintf("#s: TMax too small relative to tau, and was reset to #f", routinename, basic$TMax))
    rval <- FALSE
  }

  # basic$.online <- TRUE    # default: online learning, meaning results of patients in trial are counted in expected reward

  if(basic$t0 <= 1.0){
    basic$t0 <- defaultbasic$t0 # effective number of samples in prior distribution for unknown mean
    errorlist <- c(errorlist, sprintf("#s: t0 should exceed 1.0, reset to #f", routinename, basic$t0))
    rval <- FALSE
  }

  if(basic$mu0 > basic$mumax | basic$mu0 < basic$mumin){
    basic$mu0 <- (basic$mumax + basic$mumin)/2
    errorlist <- c(errorlist, sprintf("#s: warning, potentially invalid value of mu0, reset to #f", routinename, basic$mu0))
  }

  if(advanced$UnkVariance){ # if variance is unknown, then validate scale parameter for unknown variance
    if(advanced$UnkVarianceShape == -1.0){ # if we don't use the three parameter form for the unknown variance,
      if(basic$t0 < 4){ # xi_0 <- (t0-1)/2, and var unknown mean has 2xi_0 = (t0-1) dof, and want var to exist
        rval <- FALSE
        errorlist <- c(errorlist, sprintf("#s: error, basic$t0 is #f but should be at least 4 if variance is unknown", routinename, basic$t0))
      }
    } else if(advanced$UnkVarianceShape < 2){
      rval <- FALSE
      errorlist <- c(errorlist, sprintf("#s: error, advanced$UnkVarianceShape is #f but should be at least 2 if variance is unknown",
                                        routinename, advanced$UnkVarianceShape))
    }
    # NOTE: if variance is unknown,
    if(advanced$DOPDE){
      errorlist <- c(errorlist, sprintf("#s: warning, if sampling variance is UNKnown, normally one should set advanced$DOPDE to false",
                                        routinename))
    }
  } else {
    if(!advanced$DOPDE){
      errorlist <- c(errorlist, sprintf("#s: warning, if sampling variance is Known, normally one should set advanced$DOPDE to true",
                                        routinename))
    }
  }

  ## below needs fixing at some point
  # if(!is.null(advanced$DistributionType)){
  #   if(class(advanced$Distribution) != advanced$DistributionType){
  #     errorlist <- c(errorlist, sprintf("#s: warning, advanced$DistributionType is #s but instance advanced$Distribution is of type #s",
  #                                       routinename, advanced$DistributionType, class(advanced$Distribution)))
  #   }
  # }

  if(advanced$verbose){
    newval <- 20  # in fact, to get a reasonable contour plot this should be at least 20 or 40
    if(advanced$NumGridsForContours <= newval){
      errorlist <- c(errorlist, sprintf("#s: NumGridsForContours might be too small, suggest at least #f", routinename, newval))
    }
  }

  # for optimal stage 1 sampling budget, set to 0 or less to check at all integers, use positive
  # number, such as 50, to specify that one should check at s = tau/50, 2 tau/50, 3 tau/50, ..., tau
  if(advanced$StageOneChecks > 0){
    if(ceiling(advanced$StageOneChecks) != floor(advanced$StageOneChecks)){
      advanced$StageOneChecks <- -1
      errorlist <- c(errorlist, sprintf("#s: StageOneChecks should be negative or a positive integer, was reset to #f",
                                        routinename, advanced$StageOneChecks))
    }
  }

  #    advanced$fixedP <- TRUE # set to true if expected reward on stopping is for P * expected reward per patient, false if patients not tested due to early stopping can also benefit from better alterantive
  #    advanced$nochangeallowed <- FALSE   # default to false, so that one waits til all data in (tau samples more), before selecting. If true, then choice must be made before tau outstanding samples are observed
  #    advanced$verbose <- TRUE # print out summary information during the run
  #    advanced$RegretPenalty <- 0   # 0 for no regret penalty, positive for regret penalty
  #    advanced$smoothed <- TRUE   # smoothen the stopping boundaries if true, otherwise do not smoothen them.
  if(advanced$MAXPU <= 0.1 | advanced$MAXPU > 0.49){
    advanced$MAXPU <- defaultadvanced$MAXPU
    errorlist <- c(errorlist, sprintf("#s: MAXPU was reset to #f", routinename, advanced$MAXPU))
  }

  if(is.na(advanced$RegretPenalty)){
    advanced$RegretPenalty <- defaultadvanced$RegretPenalty
    errorlist <- c(errorlist, sprintf("#s: NA no longer supported for RegretPenalty: reset to #f", routinename, advanced$RegretPenalty))
  }

  if(advanced$RegretPenalty < 0){
    advanced$RegretPenalty <- max(0, defaultadvanced$RegretPenalty)
    errorlist <- c(errorlist, sprintf("#s: RegretPenalty needs to be nonnegative: reset to #f", routinename, advanced$RegretPenalty))
  }

  # Monte Carlo data
  #    advanced$simNumReps <- 400 # 0 for no simulations, >2 for running sample paths of the trials, and simulation estimates of power curves
  #    advanced$CRN <- TRUE       # true to use CRN for noise across Bayes and frequentist estimations, across MAT/experiments
  #    advanced$CRNAcrossExperiment <- 0 # seed to use if CRN is desired across the MAT experiment
  #    advanced$CRNAcrossBayesMu <- TRUE # true to use CRN acruss mu values for a given MAT in TestDelayIterate experiments
  #    advanced$keepAllOutput <- FALSE # true to keep all output from all replications, false if only means and std of various items are to be kept
  #    advanced$NumPointsQuadrature <- 81 # SHOULD BE ODD INTEGER: number of points for quadrature for computing expected regret once tau samples come in
  if(advanced$simNumReps < 0){
    advanced$simNumReps <- 0
    errorlist <- c(errorlist, sprintf("#s: simNumReps must be at least 0 (0 to run no reps)", routinename))
  }

  if(advanced$NumPointsQuadrature < 1) {
    advanced$NumPointsQuadrature <- defaultadvanced$NumPointsQuadrature
    errorlist <- c(errorlist, sprintf("#s: NumPointsQuadrature must be at least 1 (reset to #f)", routinename, advanced$NumPointsQuadrature))
  }

  # advanced$z <- 0 # test statsitic for frequentist tests of whether new treatment exceeds the old.
  # should default to 0 if only means are to be used, set to 1.96 for example, if stronger evidence of new better than old is needed
  # next few lines figure out how big dt and dw should be, to have pu = pd = MAXPU when we get back to start of time horizon.
  # For this, see page 2 of notes from Steve on 22 Mar 09 from Arlotto Chick
  # Gans project on leraning curves.
  MINTzero <- basic$t0
  if(advanced$UnkVariance){
    if(advanced$UnkVarianceShape == -1){
      tmpxi <- (basic$t0)/2 # use if there is a noninformative prior for mu, sigma^2 of sampling distribution
      # tmpxi <- (basic$t0 - 1)/2  # use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    } else {
      tmpxi <- (advanced$UnkVarianceShape - 1)/2 # use if the prior is otherwise specified
    }
    fudgefactor <- sqrt((2*tmpxi)/(2*tmpxi - 2))
  } else {
    fudgefactor <- 1
  }

  if(advanced$MinGridPerStdev >= 1){
    bigX <- max(advanced$MinGridPerStdev, sqrt(1.01*2*advanced$MAXPU*MINTzero)) # how many dmu per sigma do we need?
    advanced$dmu <- basic$sigma/bigX
    advanced$dt <- (2*advanced$MAXPU*MINTzero^2/bigX^2/fudgefactor)/(1 - (2*advanced$MAXPU*MINTzero/bigX^2/fudgefactor))
  } else if(advanced$MinGridPerStdev < 0){ # try to set dt to the negative of MinGridPerStdev
    advanced$dt <- -advanced$MinGridPerStdev
    bigX <- sqrt(2*advanced$MAXPU*MINTzero^2/advanced$dt) # how many dmu per sigma do we need?
    # advanced$dmu <- basic$sigma/bigX
    advanced$dmu <- fudgefactor*basic$sigma*sqrt(advanced$dt/(2*advanced$MAXPU*MINTzero*(MINTzero + advanced$dt)))
  } else { # try dt = 1
    advanced$dt <- 1
    # bigX <- sqrt(2*advanced$MAXPU*MINTzero^2/advanced$dt)   # how many dmu per sigma do we need?
    # advanced$dmu <- basic$sigma/bigX
    advanced$dmu <- fudgefactor*sqrt(basic$sigma^2*advanced$dt/(2*advanced$MAXPU*MINTzero*(MINTzero + advanced$dt)))
    bigX <- basic$sigma/advanced$dmu
  }

  if(bigX < MINBIGX){
    bigX <- MINBIGX
    advanced$dmu <- basic$sigma/bigX
    advanced$dt <- 2*advanced$MAXPU*MINTzero^2/bigX^2/fudgefactor

    errorlist <- c(errorlist, sprintf("#s: warning, grid parameters might be inaccurate: dmu reset to #f, dt reset to #f",
                                      routinename, advanced$dmu, advanced$dt))
  }

  # These are some new parameters which might not have been set up in by
  # DelayInputConstructor. Indeed, they are best set by DelayInputValidator
  # rather than by the constructor, as they may depend in a complicated way
  # on the other parameters. Be sure to call UtilMakeTvecMuvec in conjuction
  # with these issues.
  # Places where UnkVariance is handled differently than known variance
  # UtilMakeTvecMuvec()
  # DelayCurvesRecur()
  # DelaySimComputer()

  if(advanced$UnkVariance){
  # FIX: Can set dt to 1 for the unknown variance case, but then the PCS and expected number of samples
  # might not be computed correctely. Would need fixing / numerical stability checks in DelayCurvesRecur
  # advanced$dt <- 1   # SEC: Note, this is due to KG* formulation rather than PDE for unknown variance case.
  }
  tlen <- 1 + max(1, ceiling((basic$TMax - basic$tau)/advanced$dt))

  # find grid of mu values over full hi-low range, and with grid point at mu0
  # this grid allows us to return a value function estimation at precisely
  # the value of mu0 without having to interpolate
  advanced$mushift <- basic$mu0%%advanced$dmu
  # advanced$mushift <- (basic$ICost/basic$PPatients)%%advanced$dmu
  basic$mumax <- advanced$dmu*ceiling((advanced$mushift + basic$mumax)/advanced$dmu)
  basic$mumin <- advanced$dmu*floor((-advanced$mushift + basic$mumin)/advanced$dmu)

  # The following code reconstructs the time grid and mu grid in full. Better
  # to call it where it is needed: UtilMakeTvecMuvec
  # advanced$tvec <- (basic$t0 + basic$tau) + advanced$dt*(0:(advanced$tlen-1))
  # handle integer number of grid points, cover range of mu needed
  # advanced$muvec <- advanced$mushift + seq(basic$mumin, basic$mumax, by = advanced$dmu)

  if(!(advanced$simNumReps == 0 | advanced$simNumReps > 2)){
    advanced$simNumReps <- defaultadvanced$simNumReps  # effective number of samples in prior distribution for unknown mean
    errorlist <- sprintf("%s\n%s: simNumReps should be 0 or greater than 2, reset to %f\n", errorlist, routinename, advanced$simNumReps)
    rval <- 0
  }

  if("simFreqDeltaVec" %in% names(advanced)){
    if(advanced$simNumReps == 0){
      advanced$simNumReps <- defaultadvanced$simNumReps  # effective number of samples in prior distribution for unknown mean
      errorlist <- sprintf("%s\n%s: simNumReps should be greater than 2 if simFreqDeltaVec non-empty, simNumReps reset to %f\n",
                           errorlist, routinename, advanced$simPowerSteps)
      rval <- 0
    }
  } else {
    if(!("numinsimFreqDeltaVec" %in% names(advanced))){
      numinsimFreqDeltaVec <- 200  # can change the number of points in freqdeltavec
    } else {
      numinsimFreqDeltaVec <- advanced$numinsimFreqDeltaVec
    }
    advanced$simFreqDeltaVec <- (basic$sigma/sqrt(basic$t0))*(11/3)* seq(-ceiling(numinsimFreqDeltaVec/2),
                                                                         ceiling(numinsimFreqDeltaVec/2))/ceiling(numinsimFreqDeltaVec/2)
  }

  if(length(errorlist) > 0){
    warning(paste(errorlist, collapse = "\n"))
  }

  if(rval){
    return(list(basic = basic, advanced = advanced, rval = rval, errorlist = errorlist))
  } else {
    return(NULL)
  }
}