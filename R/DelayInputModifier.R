
DelayInputModifier <- function(basic, advanced, basicarray = list(), advancedarray = list()){

  # Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
  # (alpha order).
  #
  # (c) 2014, S Chick
  # Created: 14 April 2014
  # Last touched: 14 April 2014
  # Converted to R code by Michael Dymock 2024
  #
  # The first two batches of parameters are most critical in terms of defining the
  # problem structure of the clinical trial.

  rval <- TRUE
  if(length(basicarray) > 0){
    for(nm in names(basicarray)){
      if(nm %in% names(basic) | nm %in% c("mumax", "mumin")){
        basic[[nm]] <- basicarray[[nm]]
      } else {
        warning(sprintf("invalid basic field: %s", nm))
        rval <- FALSE
      }
    }
  }
  if(length(advancedarray) > 0){
    for(nm in names(advancedarray)){
      if(nm %in% names(advanced) | nm %in% c("dmu", "simFreqDeltaVec")){
        advanced[[nm]] <- advancedarray[[nm]]
      } else {
        warning(sprintf("invalid advanced field: %s", nm))
        rval <- FALSE
      }
    }
  }

  return(list(rval = rval, basic = basic, advanced = advanced))
}