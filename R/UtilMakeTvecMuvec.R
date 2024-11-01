
UtilMakeTvecMuvec <- function(basic, advanced){

  # UtilMakeTvecMuvec returns the time grid and mu grid, given the problem structures defined in the parameters basic, advanced.
  # Converted to R code by Michael Dymock 2024

  tlen <- 1 + max(1, ceiling((basic$TMax - basic$tau)/advanced$dt))
  tvec <- (basic$t0 + basic$tau) + advanced$dt*seq(0, tlen - 1)
  # handle integer number of grid points, cover range of mu needed
  muvec <- advanced$mushift + advanced$dmu*(round(basic$mumin/advanced$dmu):round(basic$mumax/advanced$dmu))

  return(list(tvec = tvec, muvec = muvec))
}