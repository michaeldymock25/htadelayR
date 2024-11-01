
#' @title DelayCurvesRecur
#' @description Implements the Stage II recursion
#' @import stats
#' @param basic A list containing basic input parameters validated by the DelayInputValidator function
#' @param advanced A list containing advanced input parameters (e.g., computational settings) validated by the DelayInputValidator function
#' @return A list containing estimated model details
#' @rdname DelayCurvesRecur
#' @export
DelayCurvesRecur <- function(basic, advanced){

  #======= BLOCK 1 ========
    # Validate the input parameters, resetting some values if needed, and put
  # the values into the local space as needed. Those adjusted values can be
  # returned to the original calling routine

  validate <- DelayInputValidator(basic, advanced)
  if(!validate$rval | length(validate$messages) > 0) warning(validate$messages)

  # take values from parameter structures and set local variables for use.
  online <- basic$online
  sigma <- basic$sigma
  tau <- basic$tau
  #mu0 <- basic$mu0
  t0 <- basic$t0
  c <- basic$c
  P <- basic$PPatients
  I <- basic$ICost
  theta <- basic$theta
  #TMax <- basic$TMax
  tvec_muvec <- UtilMakeTvecMuvec(basic, advanced) # tvec should run from t0+tau to approximately t0+Tmax in small increments
  tvec <- tvec_muvec$tvec
  muvec <- tvec_muvec$muvec
  tlen <- length(tvec)
  dt <- advanced$dt
  dmu <- advanced$dmu
  DOPDE <- advanced$DOPDE

  if(theta == 1){ # put in discounted costs for added samples
    thetadtfactor <- dt
  } else {
    thetadtfactor <- -(1 - theta^dt)/log(theta)
  }

  # NOTE FOR SPECIAL CASE OF UKNOWN VARIANCE:
  # Using UnkVariance to modify pu and pd (probability in trinomial tree) helps with getting right variance in unknown mean below, but the
  # code below does not account for the learning about the unknown variance
  # with samples. When running code with pu and pd influencing mean and not
  # variance, the upper and lower boundaries have odd behaviors (they curve back in toward I/P, and even cross each other in some cases),
  # and the VOI is underestimated - implying that Stage I sampling, when correctly
  # computed, has more value relative to the estimate from Stage II sampling.
  # This incoherency means that it is probably better to plug in an estimator
  # of the variance for the upper and lower boundaries when handling the case
  # of Unknown Variance. This is best avoided in BOTH Phase II (DelayCurvesRecur.R) and
  # Phase I (DelayStageOne.R) routines.
  # Thus, we issue a warning message here to indicate that even though code
  # is here for pu and pd to handle the mean, we do not use it in
  # computations. It is left here in order to enable further development on
  # the calculations for unknown variance at a later date.

  if(advanced$UnkVariance) {
    if(DOPDE){
      warning("DelayCurvesRecur: PDE approach with unknown variance may underestimate value of stage II sampling")
    } else {
      warning("DelayCurvesRecur: KG* approach to unknown variance may underestimate value of stage II sampling")
    }
  }

  #====== BLOCK 2 ========
    # allocate parameters, vectors et for local computation

  musize <- length(muvec)
  flipmu <- rev(muvec)
  rawboundtop <- rep(0, length(tvec)) # use this to store the upper 'boundary' of the stopping region, from grid
  rawboundbot <- rep(0, length(tvec)) # use this to store the lower 'boundary' of the stopping region, from grid
  dtdmu2ratio <- dt/dmu^2             # ratio which is used frequently

  # INITIALIZE TERMINAL CONDITIONS

  discountfactor <- theta^tau         # amount of discounting from a given time until the data for that time arrives
  i <- tlen
  curt <- tvec[i]
  numpeopletreated <- P
  predvarsample <- (numpeopletreated^2*sigma^2*tau)/((curt - tau)*curt) # note: at time curt, the effective number of samples is curt-tau,
                                                                        # because tvec starts at time t0+tau, and we have tau samples
                                                                        # not yet arrived

  if(advanced$UnkVariance){                                   # don't define enddof if variance of sampling distribution is known
    if(advanced$UnkVarianceShape == -1){
      enddof <- curt - tau                                    #P?? use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    } else {
      enddof <- curt + 2*advanced$UnkVarianceShape - tau - t0 #P?? use if the prior is otherwise specified
    }
  }

  # This 'RegretPenalty' code, if changed, needs to be changed in several
  # places: one place in DelayStageOne, two places in DelayCurvesRecur

  if(advanced$UnkVariance){ # terminal regret and expected reward with perfect info are different if sampling variance is unknown
    Cinit <- TerminalRewardFunctionUnk(muvec*numpeopletreated - I, discountfactor, predvarsample, advanced$nochangeallowed, enddof)
    postvar <- (numpeopletreated^2*sigma^2)/(enddof + tau*(1 - advanced$nochangeallowed)) # posterior variance after all samples arrive which
                                                                                          # can be seen by time of decision
    # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
    # information, given the information state at the time of stopping, assuming one can continue
    reward <- TerminalRegretUnk(muvec*numpeopletreated - I, discountfactor^(1 - advanced$nochangeallowed), predvarsample,
                                postvar, advanced, enddof)
    RewardTmaxIII <- reward$ereward
    pcsstopnow <- reward$pcs
    # Get expected reward with perfect information at time t=Tmax, assuming no discounting
    RewardPITmaxIII <- TerminalRewardFunctionUnk(muvec*numpeopletreated - I, 1.0, (numpeopletreated^2*sigma^2)/(curt - tau), FALSE, enddof)
  } else {
    Cinit <- TerminalRewardFunction(muvec*numpeopletreated - I, discountfactor, predvarsample, advanced$nochangeallowed)
    postvar <- (numpeopletreated^2*sigma^2)/(curt + tau*(1 - advanced$nochangeallowed)) # posterior variance after all samples arrive which
                                                                                        # can be seen by time of decision
    # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
    # information, given the information state at the time of stopping, assuming one can continue
    reward <- TerminalRegret(muvec*numpeopletreated - I, discountfactor^(1 - advanced$nochangeallowed), predvarsample, postvar, advanced)
    RewardTmaxIII <- reward$ereward
    pcsstopnow <- reward$pcs
    # Get expected reward with perfect information at time t=Tmax, assuming no discounting
    RewardPITmaxIII <- TerminalRewardFunction(muvec*numpeopletreated - I, 1.0, (numpeopletreated^2*sigma^2)/(curt - tau), FALSE)
  }

  pcsstopprev <- pcsstopnow
  #           mat$RewardPITmaxIII # Expected reward with perfect information at time t=Tmax
  #           mat$RewardTmaxIII   # Expected reward, assuming tau samples pending at time t=Tmax
  # compute expected regret of a decision, assuming stopping. This would be the expected value of perfect information, given the information
  # state at the time of stopping, discounted to tau time steps later if needed, for commensurability with Cinit.
  RegretIII <- RewardPITmaxIII - RewardTmaxIII

  Cinit <- Cinit - abs(advanced$RegretPenalty)*RegretIII # set initial estimate of reward to go at terminal time tHoriz
  # stoprewardvec <- Cinit
  Cin <- Cinit    # set initial estimate of reward to go at terminal time tHoriz
  Cintmp <- Cinit # done just to size Cintmp correctly - used in intermediate calculations

  ## HERE BACKWARD INDUCTION STARTS

  # Issue: Difference between accept option to decision in tau time steps or not to allow waiting for tau time
  # steps / versus the stopping boundary issue which says wait for tau time steps before selecting the alternative
  # find first element with positive value for expected reward at end of horizon
  # tstval <- max(Cinit > 0)
  # tstindx <- which.max(Cinit > 0)
  # if(tstval > 0){
  #    rawboundbot[i] <- muvec[tstindx]     # use this to store the lower 'boundary' of the stopping region, from grid
  #    rawboundtop[i] <- muvec[tstindx]     # use this to store the upper 'boundary' of the stopping region, from grid
  # } else {
  rawboundtop[i] <- I/numpeopletreated
  rawboundbot[i] <- I/numpeopletreated
  #
  # the following three are other interesting values of interest to a
  # Bayesian. These computations are valid everywhere in the stopping region,
  # but they will require tweaks below to account for their dynamics in the
  # continuation region. Keep these all as column vectors. The as.numeric() before a
  # logical expression converts the 0s and 1s to double precision rather than logical values.

  Pupperin <- as.numeric(muvec >= rawboundtop[i]) # if on or above upper boundary, prob of stopping on/above top is 1, otherwise it is 0,
                                                  # at end of time horizon
  Pupperin[muvec == rawboundtop[i]] <- 0.5        # call it 50/50 at equality
  # find probability of picking new alternative
  if(advanced$UnkVariance){                       # handle case of unknown variance and known variance separately
    Ppicknewin <- TerminalProbPickNewUnk(muvec*numpeopletreated - I, predvarsample, advanced$nochangeallowed, enddof)
  } else {
    Ppicknewin <- TerminalProbPickNew(muvec*numpeopletreated - I, predvarsample, advanced$nochangeallowed)
  }

  # initialize the expected number of samples until stopping to 0: at end of horizon no more samples can be taken
  ENumSampsin <- rep(0, length(muvec))

  # start to build up the matrices with return values. Time values in columns, different mu values in rows
  B0mat <- matrix(Cinit, nrow = length(Cinit), ncol = 1)
  Puppermat <- matrix(Pupperin, nrow = length(Pupperin), ncol = 1)
  Pnewmat <- matrix(Ppicknewin, nrow = length(Ppicknewin), ncol = 1)
  PCSmat <- matrix(pcsstopnow, nrow = length(pcsstopnow), ncol = 1)
  maxindx <- 1
  minindx <- length(Cintmp)
  nonmonotoniclower <- 0
  nonmonotonicupper <- 0
  RewardPIMatII <- matrix(RewardPITmaxIII, nrow = length(RewardPITmaxIII), ncol = 1) # Expected reward with perfect information at time t=Tmax
  RewardMatII <- matrix(RewardTmaxIII, nrow = length(RewardTmaxIII), ncol = 1)    # Expected reward, assuming tau samples pending at time t=Tmax
  tmat <- tvec[i]

  for(i in (tlen - 1):1){
    curt <- tvec[i]  # set current time index starting from t0+tau up to t0+TMAX
    numpeopletreated <- P + (1 - advanced$fixedP)*(tvec[tlen] - curt)
    predvarsample <- (numpeopletreated^2*sigma^2*tau)/((curt - tau)*curt)  # note: at time curt, the effective number of samples is curt-t0,
                                                                           # because tvec starts at time t0, and we have tau samples
                                                                           # not yet arrived
    if(advanced$UnkVariance){   # don't define enddof if variance of sampling distribution is known
      if(advanced$UnkVarianceShape == -1){
        enddof <- curt - tau  # use if there is a noninformative prior for mu, sigma^2 of sampling distribution
      } else {
        enddof <- curt + 2*advanced$UnkVarianceShape - tau - t0 # use if the prior is otherwise specified
      }
      fudgefactor <- sqrt(enddof/(enddof - 2))
      postvar <- (numpeopletreated^2*sigma^2)/(enddof + tau*(1 - advanced$nochangeallowed)) # posterior variance after all samples arrive which
                                                                                            # can be seen by time of decision
    } else {
      fudgefactor <- 1
      postvar <- (numpeopletreated^2*sigma^2)/(curt + tau*(1 - advanced$nochangeallowed)) # posterior variance after all samples arrive which
                                                                                          # can be seen by time of decision
    }

    # This 'RegretPenalty' code, if changed, needs to be changed in several
    # places: one place in DelayStageOne, two places in DelayCurvesRecur

    if(advanced$UnkVariance){ # don't define enddof if variance of sampling distribution is known
      ## WHAT FOLLOWS IS THE FUNDAMENTAL STEP HERE
      stoprewardvec <- TerminalRewardFunctionUnk(muvec*numpeopletreated - I, discountfactor, predvarsample, advanced$nochangeallowed, enddof)
      # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
      # information, given the information state at the time of stopping, assuming one can continue
      reward <- TerminalRegretUnk(muvec*numpeopletreated - I, discountfactor^(1 - advanced$nochangeallowed),
                                  predvarsample, postvar, advanced, enddof)
      Rewardmax <- reward$ereward
      pcsstopnow <- reward$pcs
      if(abs(advanced$RegretPenalty) != 0){
        RewardPI <- TerminalRewardFunctionUnk(muvec*numpeopletreated - I, 1.0, (numpeopletreated^2*sigma^2)/(curt - tau), FALSE, enddof)
        Regretcontout <- RewardPI - Rewardmax
        stoprewardvec <- stoprewardvec - abs(advanced$RegretPenalty)*Regretcontout # set initial estimate of reward to go at terminal time tHoriz
      }
    } else {
      stoprewardvec <- TerminalRewardFunction(muvec*numpeopletreated - I, discountfactor, predvarsample, advanced$nochangeallowed)
      # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
      # information, given the information state at the time of stopping, assuming one can continue
      reward <- TerminalRegret(muvec*numpeopletreated - I, discountfactor^(1 - advanced$nochangeallowed), predvarsample, postvar, advanced)
      Rewardmax <- reward$ereward
      pcsstopnow <- reward$pcs
      if(abs(advanced$RegretPenalty) != 0){
        RewardPI <- TerminalRewardFunction(muvec*numpeopletreated - I, 1.0, (numpeopletreated^2*sigma^2)/(curt - tau), FALSE)
        Regretcontout <- RewardPI - Rewardmax
        stoprewardvec <- stoprewardvec - abs(advanced$RegretPenalty)*Regretcontout # set initial estimate of reward to go at terminal time tHoriz
      }
    }

    pu <- fudgefactor*dtdmu2ratio*sigma^2/((curt - tau)*((curt - tau) + dt))/2 # probability of going up depends on t now
    pd <- pu                                                                   # probability of going down
    psame <- 1 - pu - pd                                                       # probability of going straight less a factor from the discounting

    if(DOPDE){
      # Use free boundary iteration to compute Cintmp, the vector containing value of continuing
      # at present, this code is debugged for the case of known variance.
      # for the case of unknown variance, the results are not computing
      # well. THus, we have set DOPDE false above so that the KG* naive
      # one-step look-ahead approximation is used.
      if(psame < 0.001) warning("DelayCurvesRecur: oops, psame smaller than expected")
      Cintmp[2:(musize - 1)] <- pu*Cin[3:musize] + psame*Cin[2:(musize - 1)] + pd*Cin[1:(musize - 2)]
      Cintmp <- theta^dt*Cintmp + thetadtfactor*(muvec*online - c) # and add in sampling cost and the online sampling reward if appropriate
    } else { # Use KG* type approach to compute Cintmp, the vector containing value of continuing
      numrepsleft <- basic$t0 + basic$TMax - curt
      # FIX: KGset hardcode below could probably be replaced with 2 or 3
      # paramters in the advanced structure, but
      #        KGset <- unique(c(0, 2^(0:min(2, max(1, floor(log2(numrepsleft))))/3), numrepsleft)) # check powers of 2 up to number of remaining samples, and all remaining samples
      #        KGset <- unique(c(dt, 2^seq(-0.5, max(1.5, ceiling(log2(min(128, numrepsleft)))), by = 0.5))) # check powers of 2 up to number of remaining samples, and all remaining samples
      if(i > 2){
        KGset <- unique(c(dt, 2^seq(-1, max(1.5, ceiling(log(min(1024, numrepsleft), base = 2))), by = -0.25)))
      } else { # if i is small, we are near beginning of time horizon, check a lot more 'look aheads' for this.
        KGset <- unique(c(dt, 2^seq(-1, max(1.5, ceiling(log(min(1024, numrepsleft), base = 2))), by = -0.125)))
      }
      for(j in seq_along(KGset)){ # see what happens if we take an extra KGset(j) samples before stopping to observe the remaining samples
        # note: at time curt, the effective number of samples is curt-t0, because tvec starts at time t0, and we have tau samples not yet arrived
        predtst <- (numpeopletreated^2*sigma^2*(tau + KGset[j]))/((curt - tau)*(curt + KGset[j]))
        if(advanced$UnkVariance){
          # posterior variance after all samples arrive which can be seen by time of decision
          posttst <- (numpeopletreated^2*sigma^2)/(enddof + KGset[j] + tau*(1 - advanced$nochangeallowed))
          # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
          # information, given the information state at the time of stopping, assuming one can continue
          OneShotReward <- TerminalRewardFunctionUnk(muvec*numpeopletreated - I, theta^KGset[j]*discountfactor^(1 - advanced$nochangeallowed),
                                                     predtst, advanced$nochangeallowed, enddof)
        } else {
          # posterior variance after all samples arrive which can be seen by time of decision
          posttst <- (numpeopletreated^2*sigma^2)/(KGset[j] + curt + tau*(1 - advanced$nochangeallowed))
          # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
          # information, given the information state at the time of stopping, assuming one can continue
          OneShotReward <- TerminalRewardFunction(muvec*numpeopletreated - I, theta^KGset[j]*discountfactor^(1 - advanced$nochangeallowed),
                                                  predtst, advanced$nochangeallowed)
        }

        if(theta == 1){ # put in discounted costs for added samples
          OneShotReward <- OneShotReward + KGset[j]*(muvec*online - c)
        } else { # -(1-basic$theta^dt)/log(basic$theta) is integral of theta^t dt from 0 to KGset[j]
          OneShotReward <- OneShotReward - ((1 - theta^KGset[j])/log(theta))*(muvec*online - c)
        }

        if(j == 1){
          Cintmp <- OneShotReward
        } else {
          Cintmp <- pmax(Cintmp, OneShotReward)
        }
      }
    }

    Cintmp[musize] <- stoprewardvec[musize] # now, the discounting and sampling cost might given incorrect values for the extremes of mu,
    Cintmp[1] <- stoprewardvec[1]           # so we again must reset their values to indure that we remain in the stopping set
    # find stopping boundary
    Cout <- pmax(Cintmp, stoprewardvec)     # implement maximizer for bellman equation
    stopped <- (Cintmp <= stoprewardvec)    # stopped[i] is true if it is optimal to stop when mean is muvec[i]
    oldmaxindx <- maxindx
    oldminindx <- minindx
    minval <- min(stopped)
    minindx <- which.min(stopped)

    if(minval > 0){ # if there is not a nonempty stopping region
      minindx <- which.min(muvec* numpeopletreated - I < 0)
      maxindx <- minindx # length(muvec)+1-minindx # Thanks Paolo! (SEC touched lines 203 and 205 on 7 oct 2014)
      rawboundbot[i] <- muvec[minindx] # I/numpeopletreated    # use this to store the lower boundary of the stopping region, from grid
      rawboundtop[i] <- muvec[maxindx] # flipmu(maxindx)I/numpeopletreated use this to store the upper boundary of the stopping region, from grid
    } else {  # there is a nontrivial stopping region
      rawboundbot[i] <- muvec[minindx] # grid point in bottom boundary is in contin set
      flipstop <- rev(stopped)
      maxval <- min(flipstop)
      maxindx <- which.min(flipstop)   # this finds first point in the continuation set
      if(maxval > 0){                  # if there not a value in the continuation set
        # assume stop bound is as high as possible this condition should never happen, as there is at least one 1 in stopvec
        rawboundtop[i] <- muvec[length(muvec)]
      } else {
        rawboundtop[i] <- flipmu[maxindx] # get first grid point in continuation set
      }
      maxindx <- length(muvec) + 1 - maxindx
    }
    # some test code to see if enforcing monotonicity of continuation set results in an overcoming of some numerical stability issues
    if(minindx > oldminindx){
      nonmonotoniclower <- nonmonotoniclower + 1
      minindx <- min(minindx, oldminindx)
      rawboundbot[i] <- muvec[minindx]
    }
    if(maxindx < oldmaxindx){
      nonmonotonicupper <- nonmonotonicupper + 1
      maxindx <- max(maxindx, oldmaxindx)
      rawboundtop[i] <- muvec[maxindx]
    }

    # compute the prob(accept new technology) values, and prob(stop on or above upper boundary)
    # if above upper boundary, prob of stopping above top is 1, otherwise it is 0, at end of time horizon
    # find probability of picking new alternative
    Pupperout <- as.numeric(muvec >= rawboundtop[i])
    if(advanced$UnkVariance){ # handle case of unknown variance and known variance separately
      if(advanced$UnkVarianceShape == -1){
        enddof <- curt # use if there is a noninformative prior for mu, sigma^2 of sampling distribution
      } else {
        enddof <- curt + 2*advanced$UnkVarianceShape - t0 # use if the prior is otherwise specified
      }
      Ppicknewout <- TerminalProbPickNewUnk(muvec*numpeopletreated - I, predvarsample, advanced$nochangeallowed, enddof)
    } else {
      Ppicknewout <- TerminalProbPickNew(muvec*numpeopletreated - I, predvarsample, advanced$nochangeallowed)
    }

    # by default, expected num samps more to do is 0, reset below for case of contin region
    # for values in continuation region, compute expectation via
    # conditioning on what happens at time t+dt
    ENumSampsout <- numeric(length(muvec))
    Pupperout[minindx:maxindx] <- DelayContinExpectation(minindx, maxindx, pu, pd, psame, Pupperin)
    Ppicknewout[minindx:maxindx] <- DelayContinExpectation(minindx, maxindx, pu, pd, psame, Ppicknewin)

    if(advanced$UnkVariance) {
      ENumSampsout[minindx:maxindx] <- dt + DelayContinExpectation(minindx, maxindx, pu, pd, psame, ENumSampsin)
    } else {
      ENumSampsout[minindx:maxindx] <- dt + DelayContinExpectation(minindx, maxindx, pu, pd, psame, ENumSampsin)
    }

    # now, update the variables for the recursion
    Pupperin <- Pupperout
    Ppicknewin <- Ppicknewout
    ENumSampsin <- ENumSampsout
    pcsout <- pcsstopnow
    pcsout[minindx:maxindx] <- DelayContinExpectation(minindx, maxindx, pu, pd, psame, pcsstopprev)

    if(i == tlen - 1){
      # This stores the matrices for the output values, and sets a new
      # time vector to manage them all. FIX: maybe some speed optimization can be done here by preallocating these things?
      B0mat <- cbind(Cout, B0mat)
      Puppermat <- cbind(Pupperin, Puppermat)
      Pnewmat <- cbind(Ppicknewin, Pnewmat)
      PCSmat <- cbind(pcsstopnow, PCSmat)

      if(abs(advanced$RegretPenalty) == 0){ # this was computed earlier for the case of regretpenalty ~= 0
        if(advanced$UnkVariance){ # handle case of unknown variance and known variance separately
          RewardPI <- TerminalRewardFunctionUnk(muvec*numpeopletreated - I, 1.0, (numpeopletreated^2*sigma^2)/(curt - tau), FALSE, enddof)
        } else {
          RewardPI <- TerminalRewardFunction(muvec*numpeopletreated - I, 1.0, (numpeopletreated^2*sigma^2)/(curt - tau), FALSE)
        }
      }

      RewardPIMatII <- cbind(RewardPI, RewardPIMatII) # Expected reward with perfect information at time t=Tmax
      RewardMatII <- cbind(Rewardmax, RewardMatII)    # Expected reward, assuming tau samples pending at time t=Tmax
      tmat <- c(tvec[i], tmat)
    }

    pcsstopprev <- pcsout
    Cin <- Cout
  }

  ## COMPUTE BOUNDARY (with smoothing if needed)

  if(!advanced$smoothed){  # no smoothing desired
    mintvec <- tvec
    bndupper <- rawboundtop
    bndlower <- rawboundbot
  } else { # do smoothing, using the sparser boundary as a basis
    tighten <- TRUE # false : keep all time values, true: try to reduce number of time values
    if(advanced$UnkVariance){
      result <- LinearSmoothifier(tvec, -rawboundtop, tighten)
      mintvec1 <- result$tout
      bndupper <- -result$valout
    } else {
      result <- LinearSmootherLower(tvec, rawboundtop, tighten)
      mintvec1 <- result$tout
      bndupper <- result$valout
    }
    result <- LinearSmootherLower(tvec, rawboundbot, tighten)
    mintvec2 <- result$tout
    bndlower <- result$valout

    mintvec <- sort(unique(c(mintvec1, mintvec2)))
    bndupper <- approx(mintvec1, bndupper, mintvec, rule = 2)$y
    bndlower <- approx(mintvec2, bndlower, mintvec, rule = 2)$y
    gridest <- max(1, 1 + round((min(bndlower) - basic$mumin)/advanced$dmu))
    ENumSampsin[gridest] <- max(ENumSampsin[gridest], 0.01)
    gridest <- min(length(ENumSampsin), 1 + round((max(bndupper) - basic$mumin)/advanced$dmu))
    ENumSampsin[gridest] <- max(ENumSampsin[gridest], 0.01)
  }

  # Finish setting output values
  # B0vec <- Cin
  # Build a structure which has all the extra outputs, if desired
  mat <- list()
  mat$B0mat <- B0mat                # value function
  mat$Puppermat <- Puppermat        # Probability of stopping due to hitting/exceeding upper boundary
  mat$Pnewmat <- Pnewmat            # Probability that if one stopped at a given (t, mu), that new technology is adopted
  mat$PCSmat <- PCSmat              # Probability that if one stopped at a given (t, mu), that new technology is adopted
  mat$ENumSamps <- ENumSampsin
  mat$Bayespcsvec <- pcsout
  mat$tmat <- tmat                  # should be a row vector, gives time indices for the above matrices for contour plots

  mat$muvec <- as.vector(muvec)     # should be a column vector of the mu values
  mat$tvec <- mintvec               # tvec is for the boundaries: may be different from tmat, as tvec tries to compress
  mat$bndupper <- bndupper          # info requirements to get at curves of boundary, whereas tmat is equally spaced, in order
  mat$bndlower <- bndlower          # to help contour plot
  mat$B0vec <- Cin                  # first column of B0mat (might trim this as it is extra unneeded info
  # Get expected reward with perfect information at time t=Tmax, assuming no discounting
  mat$RewardPITmaxIII <- as.vector(RewardPITmaxIII)
  # the 'min' corrects numerical stability issues with RewardTmaxIII
  mat$RewardTmaxIII <- pmin(as.vector(RewardTmaxIII), as.vector(RewardPITmaxIII))
  mat$RewardPIMatII <- RewardPIMatII  # Expected reward with perfect information at time t in stage II
  mat$RewardMatII <- RewardMatII      # Expected reward, assuming tau samples pending at time t in stage II

  retval <- (max(bndupper) < max(muvec)) | (min(bndlower) > min(muvec))

  if(nonmonotoniclower + nonmonotonicupper > 0) warning("forced boundaries to retain monotonicity properties")

  return(list(retval = retval, mat = mat))
}