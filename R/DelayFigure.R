
#' @title DelayFigure
#' @description Visualises stopping boundaries
#' @import ggplot2
#' @param basic A list containing basic input parameters validated by the DelayInputValidator function
#' @param mat A list containing estimated model details created by the DelayStageOne function
#' @param xlims Limits for the horizontal axis
#' @param xnbreaks Number of breaks for the horizontal axis
#' @param ylims Limits for the vertical axis
#' @param ynbreaks Number of breaks for the vertical axis
#' @param thresholds If TRUE, plot threshold boundaries for a one stage trial
#' @param points If TRUE, plot points for a series of one stage trials
#' @return A list containing estimated model details
#' @rdname DelayFigure
#' @export
DelayFigure <- function(basic, mat, xlims = NULL, xnbreaks = NULL, ylims = NULL, ynbreaks = NULL, thresholds = TRUE, points = TRUE){

  n <- E_INMB <- bound <- NULL

  # Visualises optimal bounds for the value-based sequential design
  # Created by Michael Dymock 2024

  dat_line <- data.frame(n = rep(mat$tvec, 2),
                         E_INMB = c(mat$bndlower, mat$bndupper),
                         bound = rep(c("lower", "upper"), each = length(mat$tvec)))

  if(points){
    dat_points <- data.frame(n = c(mat$bestsvec[mat$Threshpoint[1]:mat$Threshpoint[3]],
                                   mat$bestsvec[mat$Threshpoint[4]:mat$Threshpoint[2]]),
                             E_INMB = c(mat$muvec[mat$Threshpoint[1]:mat$Threshpoint[3]],
                                        mat$muvec[mat$Threshpoint[4]:mat$Threshpoint[2]]))
    fig <- ggplot(dat_points) +
      geom_point(aes(x = n, y = E_INMB), colour = "#88CCEE")
  } else {
    fig <- ggplot(dat_line)
  }

  if(thresholds){
    fig <- fig +
      geom_segment(x = 0, xend = basic$tau + basic$t0, y = mat$muvec[mat$Threshpoint[1]], linetype = "dashed") +
      geom_segment(x = 0, xend = basic$tau + basic$t0, y = mat$muvec[mat$Threshpoint[2]], linetype = "dashed") +
      geom_segment(x = 0, xend = basic$tau + basic$t0, y = mat$muvec[mat$Threshpoint[3]], linetype = "dashed") +
      geom_segment(x = 0, xend = basic$tau + basic$t0, y = mat$muvec[mat$Threshpoint[4]], linetype = "dashed")
  }

  fig <- fig +
    geom_line(data = dat_line, aes(x = n, y = E_INMB, group = bound), colour = "#CC6677") +
    geom_vline(xintercept = basic$tau + basic$t0, linetype = "dashed") +
    geom_vline(xintercept = max(mat$tvec), linetype = "dashed") +
    scale_x_continuous("Pairwise allocations", n.breaks = xnbreaks, limits = xlims) +
    scale_y_continuous("Prior/Posterior mean of E[INMB]", n.breaks = ynbreaks, limits = ylims)

  return(fig)
}