
#' @title DelayStageOne
#' @import stats
#' @description Implements the Stage I recursion
#' @param basic A list containing basic input parameters validated by the DelayInputValidator function
#' @param advanced A list containing advanced input parameters (e.g., computational settings) cvalidated by the DelayInputValidator function
#' @param mat A list containing estimated model details created by the DelayCurvesRecur function
#' @return A list containing a modified copy of the mat input list
#' @rdname DelayStageOne
#' @export
DelayStageOne <- function(basic, advanced, mat){

  ######################################################
    # STEP 1: Initialize variables
  ######################################################

  muvec <- mat$muvec
  sigma <- basic$sigma
  t0 <- basic$t0
  theta <- basic$theta
  tau <- basic$tau
  TMAX <- basic$TMax
  #PPatients <- basic$PPatients
  ICost <- basic$ICost
  c <- basic$c
  B0vec <- mat$B0vec
  advanced2 <- advanced
  advanced2$DoRegretIntegral <- TRUE # force integral for regret calculations for one stage procedure

  ######################################################
    # Step 2: compute the value of having perfect information at time 0 and making the best decision
  ######################################################

  numpatients <- basic$PPatients + (1 - advanced$fixedP)*TMAX
  predvarsample <- numpatients^2*sigma^2/t0

  #NOTE FOR SPECIAL CASE OF UKNOWN VARIANCE:
  # Using UnkVariance to modify pu and pd (probability in trinomial tree) helps with getting right variance in unknown mean below, but the
  # code below does not account for the learning about the unknown variance
  # with samples. When running code with pu and pd influencing mean and not
  # variance, the upper and lower boundaries have odd behaviors (they curve back in toward I/P, and even cross each other in some cases),
  # and the VOI is underestimated - implying that Stage I sampling, when correctly
  # computed, has more value relative to the estimate from Stage II sampling.
  # This incoherency means that it is probably better to plug in an estimator
  # of the variance for the upper and lower boundaries when handling the case
  # of Unknown Variance. This is best avoided in BOTH Phase II (DelayCurvesRecur.m) and
  # Phase I (DelayStageOne.R) routines.
  # Thus, we issue a warning message here to indicate that even though code
  # is here for pu and pd to handle the mean, we do not use it in
  # computations. It is left here in order to enable further development on
  # the calculations for unknown variance at a later date.

  # if(advanced$UnkVariance){
  #    warning("DelayCurvesRecur: Unknown variance not implemented in Phase I analysis, fudge factor applied later in monte carlo simulations")
  #    advanced$UnkVariance <- FALSE
  # }

  if(advanced$UnkVariance){ # don't define enddof if variance of sampling distribution is known
    if(advanced$UnkVarianceShape == -1){
      dof <- t0  # use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    } else {
      dof <- 2*advanced$UnkVarianceShape # use if the prior is otherwise specified
    }
    RewardPIT0I <- TerminalRewardFunctionUnk(muvec*numpatients - ICost, 1.0, predvarsample, FALSE, dof)
  } else {
    RewardPIT0I <- TerminalRewardFunction(muvec*numpatients - ICost, 1.0, predvarsample, FALSE)
  }

  ######################################################
    # Step 3: For cases of delay or no delay in samples, check for best one-stage allocation.
  ######################################################

  if(tau > 0){  # IF tau > 0, then we have the opprtunity to have a delay
    if(advanced$StageOneChecks <= 0){ # SEC: Updated 19 nov 2014 to allow for user to specify number of checks for optimal stage 1 sampling size
      deltas <- 1 # Give 0 or a negative number to check for s=1, 2, 3, ..., \tau.
                  # Give positive integer to specify the number of equally spaced values of time to check
    } else {
      deltas <- tau/advanced$StageOneChecks # check at times t0 + [ tau, 2 tau, ... tau*tau] / tau
    }
    svec <- seq(t0 + deltas, t0 + tau, by = deltas) # if sampling to occur in stage 1, must take at least 1 sample
    svec[length(svec)] <- t0 + tau
    PPatientvec <- rep(basic$PPatients, length(svec)) + (1 - advanced$fixedP)*(TMAX - svec) # number of patients to be treated following adoption
    discountvector <- theta^(svec - t0 + tau) # decision will be made tau units of time after last of s observations made
    predvarvec <- (PPatientvec^2*sigma^2)*(svec - t0)/(t0*svec) # predictive variance of posterior mean to be seen after svec samples
    PPMax <- basic$PPatients + (1 - advanced$fixedP)*TMAX # total number of patients in contract plus max in trial who could be switched
    G0base <- matrix(0, nrow = length(muvec), ncol = length(svec)) # preallocate size of matrix for G0

    for(i in seq_along(svec)){ # for each time value s for sampling, compute terminal expected reward for each
      # This 'RegretPenalty' code, if changed, needs to be changed in several
      # places: one place in DelayStageOne, two places in DelayCurvesRecur
      postvar <- (PPatientvec[i]^2*sigma^2)/svec[i] # posterior variance after all samples arrive which can be seen by time of decision
      if(advanced$UnkVariance){
        G0base[,i] <- TerminalRewardFunctionUnk(muvec*PPatientvec[i] - ICost, discountvector[i], predvarvec[i], FALSE, dof)
        # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
        # information, given the information state at the time of stopping, assuming one can continue
        RewardI <- TerminalRegretUnk(muvec*PPatientvec[i] - ICost, discountvector[i]^(1 - advanced$nochangeallowed),
                                     predvarvec[i], postvar, advanced2, dof)$ereward
      } else {
        G0base[,i] <- TerminalRewardFunction(muvec*PPatientvec[i] - ICost, discountvector[i], predvarvec[i], FALSE)
        # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
        # information, given the information state at the time of stopping, assuming one can continue
        RewardI <- TerminalRegret(muvec*PPatientvec[i] - ICost, discountvector[i]^(1 - advanced$nochangeallowed),
                                  predvarvec[i], postvar, advanced2)$ereward
      }
      G0base[,i] <- G0base[,i] - abs(advanced$RegretPenalty)*(RewardPIT0I - RewardI) # set initial estimate of reward to go at terminal time tHoriz
    }

    if(theta < 1){
      B0hat <- -c*(1 - theta^tau)/(1 - theta) + theta^tau*B0vec
      G0initcost <- matrix(-c*(1 - theta^(svec - t0))/(1 - theta), nrow = length(muvec), ncol = length(svec), byrow = TRUE)
      if(basic$online) {
        B0hat <- B0hat + muvec*(1 - theta^tau)/(1 - theta)
        G0initcost <- G0initcost + t(muvec)%*%(1 - theta^(svec - t0)/(1 - theta))
      }
    } else {
      B0hat <- -c*tau + B0vec
      G0initcost <- matrix(-c*(svec - t0), nrow = length(muvec), ncol = length(svec), byrow = TRUE)
      if(basic$online){
        B0hat <- B0hat + muvec*tau
        G0initcost <- G0initcost + t(muvec)%*%(svec - t0)
      }
    }

    G0 <- G0initcost + G0base
    beststage1 <- apply(G0, 1, max)
    bests <- apply(G0, 1, which.max)
    bestsvec <- pmin(svec[bests] - t0, tau) # get optimal number of samples

    # check if taking 0 samples is optimal among policies which do not advance to Stage II.
    if(advanced$nochangeallowed){ # In stage I, no change allowed implies that stopping now before new data arrives requires that one pick the
                                  # best. One can only change if one continues to stage II.
      beststage1 <- pmax(muvec*PPMax - ICost, 0) # check if s=0 is even better than positive s in (0, tau]
      bestsvec[beststage1 == pmax(muvec*PPMax - ICost, 0)] <- 0
    } else {
      # The default is advanced.nochangeallowed=false, the following code, which permits stopping
      # before stage II after a few samples have been taken, and allows one to decide which alternative is best.
      beststage1 <- pmax(beststage1, pmax(muvec*PPMax - ICost, 0)) # check if s=0 is even better than positive s in (0, tau]
      bestsvec[beststage1 == pmax(muvec*PPMax - ICost, 0)] <- 0
    }

    # if optimal to continue into stage II sampling, then set optimal number of samples to tau + 1.
    # Note: formally this is correct when the variance is known. When the variance is unknown we use a kludge: we test to see if the VOI from
    # one stage is maximizes on [0, tau] at tau and consider that to be a decision to go to stage II (the reason is that in the variance
    # unknown case, the VOI in stage I is calculated with student distributions, where as for stage II diffusion we assume known
    # variance (variance is known for diffusion sampled on a continuous interval). Thus, the B0hat is not directly commensurate with the
    # stage I VOI due to the known/unknown variance mismatch.

    if(advanced$UnkVariance & advanced$DOPDE){
      bestsvec[bestsvec > tau*2/3] <- tau + 1 # this is the kludge: go to stage II if it appears VOI is increasing at tau samples
    } else {
      bestsvec[beststage1 < B0hat] <- tau + 1
    }

    B0hat <- pmax(B0hat, beststage1) # again this is right for known variance, an approx (lower bound) for unknown variance.

    if(advanced$verbose){
      bestsrealzero <- (bestsvec == 0)
      bestszero <- (bestsvec <= tau & bestsvec > 0)
      beststau <- (bestsvec > tau)
      cat("Number of bestsrealzero: ", sum(bestsrealzero), "\n")
      cat("Number of bestszero: ", sum(bestszero), "\n")
      cat("Number of beststau: ", sum(beststau), "\n")
    }

    # set output variables, including adding fields to the output structure.
    matout <- mat
    matout$B0hat <- B0hat
    matout$bestsvec <- bestsvec
    matout$optENumSamps <- bestsvec
    # compute mean optimal number of samples over both stages
    matout$optENumSamps[bestsvec == (tau + 1)] <- tau + mat$ENumSamps[bestsvec == (tau + 1)]
  } else { # tau = 0: this means the data from stage II is what matters most - there is no stage I
    matout <- mat
    matout$B0hat <- B0vec # no cost in stage 1
    bestsvec <-  rep(0, length = length(mat$ENumSamps))
    matout$bestsvec <- bestsvec
    matout$bestsvec[mat$ENumSamps > 0] <- tau + 1
    matout$optENumSamps <- mat$ENumSamps
  }

  ######################################################
    # Step 4: Now compute statistics for one stage trial of all TMax patients
  ######################################################

  PPMax <- basic$PPatients
  discountfactor <- theta^(TMAX + (1 - advanced$nochangeallowed)*tau)
  predvarsample <- (PPMax^2*sigma^2*(TMAX - advanced$nochangeallowed*tau))/((t0 + TMAX - advanced$nochangeallowed*tau)*t0)
  # posterior variance after all samples arrive which can be seen by time of decision
  postvar <- (PPMax^2*sigma^2)/(t0 + TMAX - advanced$nochangeallowed*tau)

  if(advanced$UnkVariance){
    if(advanced$UnkVarianceShape == -1){
      dof <- t0 + TMAX - advanced$nochangeallowed*tau # use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    } else {
      dof <- 2*advanced$UnkVarianceShape + TMAX - advanced$nochangeallowed*tau # use if the prior is otherwise specified
    }
    # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
    # information, given the information state at the time of stopping, assuming one can continue
    reward <- TerminalRegretUnk(muvec*PPMax - ICost, discountfactor^(1 - advanced$nochangeallowed), predvarsample, postvar, advanced2, dof)
    OneShotReward <- reward$ereward
    OneshotPCS <- reward$pcs
  } else {
    # compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect
    # information, given the information state at the time of stopping, assuming one can continue
    reward <- TerminalRegret(muvec*PPMax - ICost, discountfactor^(1 - advanced$nochangeallowed), predvarsample, postvar, advanced2)
    OneShotReward <- reward$ereward
    OneshotPCS <- reward$pcs
  }

  OneshotRegret <- RewardPIT0I - OneShotReward
  # set initial estimate of reward to go at terminal time tHoriz
  BayesOneshotExpectedReward <- OneShotReward - abs(advanced$RegretPenalty)*OneshotRegret

  # need to take out sampling costs!
  if(theta < 1){ # subtract off the sampling costs, depending on whether there is discounting or not
    rewarddelta2 <- (basic$online*muvec - c)*(1 - theta^basic$TMax)/(1 - theta)
  } else {
    rewarddelta2 <- (basic$online*muvec - c)*basic$TMax
  }
  BayesOneshotExpectedReward <- BayesOneshotExpectedReward + rewarddelta2

  matout$OneshotRegret <- OneshotRegret
  matout$OneshotPCS <- OneshotPCS
  matout$BayesOneshotExpectedReward <- BayesOneshotExpectedReward
  zfreq <- abs(muvec*PPMax - ICost)/(PPMax*sigma/sqrt(basic$TMax)) # SEC: fixed one shot pcs and eoc, 30 dec 2014.
  if(advanced$UnkVariance){
    if(advanced$UnkVarianceShape == -1){
      dof <- t0
    } else {
      dof <- 2*advanced$UnkVarianceShape
    }
    matout$FreqOneShotPCS <- pt(zfreq, dof)
    matout$FreqOneShotEOC <- PPMax*(sigma/sqrt(basic$TMax))*PsiNormUV(zfreq, dof)
  } else {
    matout$FreqOneShotPCS <- pnorm(zfreq)
    matout$FreqOneShotEOC <- PPMax*(sigma/sqrt(basic$TMax))*PsiNorm(zfreq)
  }

  matout$Threshpoint <- UtilGetABCDVector(bestsvec, tau)
  if(advanced$verbose){
    cat('bestsvec\n')
    print(bestsvec)
    cat('matout$Threshpoint\n')
    print(matout$Threshpoint)
  }
  matout$RewardPIT0I <- RewardPIT0I
  tvectotest <- seq(0, basic$TMax)
  matout <- DelayOptimalBayesOneStage(basic, advanced, tvectotest, matout)

  tstval <- (matout$bestsvec[1] == 0) & (matout$bestsvec[length(matout$bestsvec)] == 0)
  if(!tstval){
    warning('DelayStageOne: range of nonzero optimal allocations may exceed [mumin, mumax], consider making basic$mumax larger and/or basic$mumin smaller')
  }

  return(matout)
}