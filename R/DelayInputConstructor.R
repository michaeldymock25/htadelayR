
#' @title DelayInputConstructor
#' @description Creates the initial basic and advanced data structures which are used as inputs for the model
#' @param basicarray A list containing any basic input parameters
#' @param advancedarray A list containing any advanced input parameters (e.g., computational settings)
#' @return A list containing the initial basic and advanced data structures
#' @rdname DelayInputConstructor
#' @export
DelayInputConstructor <- function(basicarray = list(), advancedarray = list()){

  basic <- list(c = 0,             # marginal cost per sampled pair
                theta = 0.999,     # discount factor per sampled pair
                PPatients = 10000, # number of future patients affected
                ICost = 0,         # fixed investment cost
                tau = 40,          # delay of outcomes
                TMax = 1000,       # maximum number of patient pairs
                numpairs = 700,    # number of patient pairs to be tested for simulation
                online = TRUE,     # default: online learning, meaning results of patients in trial are counted in expected reward
                # If there is a sampling distribution which is normally distributed, and
                # which has unknown mean and known variance, and the PDE calculations can
                # rely on knowing the variance, then use the following:
                mu0 = 0,           # mean of prior distribution for unknown mean reward
                t0 = 5,            # effective number of samples in prior distribution for unknown mean reward
                sigma = 100)       # standard deviation of sample differences (between people in two trial arms)

  # careful: UnkVariance is true when the unknown variance is to be
  # sampled prior to each sample path for the trial.
  # UnkVarBound is true when the unknown variance prior is to be used for
  # the plugin estimator for the inference process.
  # these can be set independently. Careful!

  advanced <- list(UnkVariance = FALSE,      # false if variance is known, true if unknown (default to known variance)
                   UnkVarianceShape = -1,    # shape parameter for unknown mean for use with plug-in estimator when variance is assumed unknown
                   UnkVarBound = FALSE,      # false if variance is known, true if unknown (default to known variance)
                   # Set UnkVarianceShape = -1 to have a 3-parameter conjugate prior
                   #   for the unknown mean and variance, as in Chick and Inoue (Oper Res 2001), eq. (2) and (3)
                   # Set UnkVarianceShape = shape parameter of inverted gamma distribution for
                   #   unknown variance, as in \xi_{i,0} of Chick and Frazier (Man Sci 2012), eq. (20).
                   #   In this case, \chi_{i,0} is set so that basic$sigma = \chi_{i,0} / (\xi_{i,0} - 1), meaning that basic$sigma gives the
                   #   a priori mean value of the unknown variance
                   DistributionType = NULL,  # set to empty string unless object used for sampling distribution
                   Distribution = NULL,
                   # If there is a sampling distribution which is normally distributed, and
                   # which has unknown mean and known variance, and the PDE calculations can
                   # NOT rely on knowing the variance, then use the following:
                     # basic$mu0 = 0          # mean of prior distribution for uknown mean reward
                     # basic$t0 = 5           # effective number of samples in prior distribution for unknown mean reward
                     # basic$sigma = 100       # standard deviation of sample differences (between people in two trial arms)
                     # advanced$UnkVariance = TRUE # identify if the stopping boundary will be adjusted to account for sampling variance or not
                     # advanced$DistributionType = dnorm # identify the type of sampling distribution object
                     # create an instance of the object with the right hyperparameters
                     # advanced$Distribution = advanced$DistributionType(basic$mu0, basic$t0, basic$sigma)

                   # If there is a sampling distribution which is normally distributed, and
                   #% which has unknown mean and UNknown variance, and the PDE calculations can
                   # NOT rely on knowing the variance, then use the following:
                     # basic$mu0 = 0         # mean of prior distribution for uknown mean reward
                     # basic$t0 = 5          # effective number of samples in prior distribution for unknown mean reward
                     # basic$sigma = 100     # standard deviation of sample differences (between people in two trial arms)
                     # xi0 = 20              # shape parameter for unknown variance
                     # advanced$UnkVariance = TRUE # identify if the stopping boundary will be adjusted to account for sampling variance or not
                     # advanced$DistributionType = dnorm #identify the type of sampling distribution object
                     # create an instance of the object with the right hyperparameters
                     # advanced$Distribution = advanced.DistributionType(basic$mu0, basic$t0, basic$sigma, xi0)
                   mushift = 0,
                   MinGridPerStdev = 30,     # minimum number of delta-mus per standard deviation
                   dt = NA,                  # delta t for the recursion: set to NA if it should be computed from dmu
                   fixedP = TRUE,            # set to true if expected reward on stopping is for P * expected reward per patient, false if
                                             # patients not tested due to early stopping can also benefit from better alternative
                   nochangeallowed = FALSE,  # default to false, so that one waits til all data in (tau samples more), before selecting. If true,
                                             # then choice must be made before tau outstanding samples are observed
                   verbose = TRUE,           # print out summary information during the run
                   RegretPenalty = 0,        # positive for regret, or 0 for standard calculation (no regret penalty)
                   smoothed = TRUE,          # smoothen the stopping boundaries if true, otherwise do not smoothen them
                   MAXPU = 0.475,            # maximum probability of going up (or down) in the trinomial tree which is to be constructed
                   StageOneChecks = -1,      # for optimal stage 1 sampling budget, set to 0 or less to check at all integers, use positive
                   DOPDE = TRUE,             # use the PDE mechanism to compute stage II if true (the default), or a KG-type otherwise
                   # KGNoPCS = FALSE # don't compute the PCS and ENumSamps, ONLY active if DOPDE is false (so that KG* computation is active)
                   # number, such as 50, to specify that one should check at s = tau/50, 2 tau/50, 3 tau/50, ..., tau
                   simNumReps = 200 ,        # 0 for no simulations, > 2 for running sample paths of the trials and power curves
                   CRN = TRUE,               # true to use CRN for noise across Bayes and frequentist estimations, across MAT
                   CRNAcrossExperiment = 37, # seed for use with common random numbers (if CRN is true)
                   CRNAcrossBayesMu = TRUE,  # true to use CRN across mu in TestDelayIterate experiments too
                   keepAllOutput = TRUE,     # true to keep all output from all replications, false if only for means and std of various items
                   DoRegretIntegral = FALSE, # TRUE if regret should be done to high accuracy (integral) or FALSE for quadrature approximation
                   NumPointsQuadrature = 80, # If DoRegretIntegral is false, then INTEGER number of points for quadrature approximation for regret
                   z = 0)                    # used as a test statistic for frequentist routines

  # Now that basic parameters are set up, check to see if values have been
  # passed by end user. Note that the parameters further down which are
  # computed especially by parameters above need to be handled specially.
  # Same for SOME parameters which are set by the validator routine, e.g.
  # allow end user to specify simFreqDeltaVec,but not tlen, PlotCheck,
  # mushift as those are selected for internal consistency purposes.

  modified <- DelayInputModifier(basic = basic, advanced = advanced, basicarray = basicarray, advancedarray = advancedarray)
  rval <- modified$rval
  basic <- modified$basic
  advanced <- modified$advanced

  # COMPUTED PARAMETERS: If more are put here, please fix the code to insure
  # that they don't overwrite parameter values which may have been entered by
  # the end user.

  if(!rval){
    return(list(basic = list(), advanced = list()))
  } else {
    if(advanced$UnkVariance){
      if(advanced$UnkVarianceShape == -1){
        tmpxi <- basic$t0/2         # use if there is a noninformative prior for mu, sigma^2 of sampling distribution
        # tmpxi <- (basic.t0 - 1)/2 # use if there is a noninformative prior for mu, sigma^2 of sampling distribution
      } else {
        tmpxi <- (advanced$UnkVarianceShape - 1)/2 # use if the prior is otherwise specified
      }
      fudgefactor <- sqrt((2*tmpxi)/(2*tmpxi - 2))
    } else {
      fudgefactor <- 1
    }
    # the increment used for the grid for the posterior mean
    if("dmu" %in% names(advanced)) advanced$dmu <- fudgefactor*basic$sigma/advanced$MinGridPerStdev
    # the upper value to be included in the grid for the posterior mean
    if("mumax" %in% names(basic)) basic$mumax <- 15*basic$sigma/sqrt(basic$t0)
    # the lower value to be included in the grid for the posterior mean
    if("mumin" %in% names(basic)) basic$mumin <- -basic$mumax
    return(list(basic = basic, advanced = advanced))
  }
}