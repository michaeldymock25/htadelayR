# htadelayR

This repository contains an R implementation of the Matlab code provided at https://github.com/sechick/htadelay/ by Stephen E. Chick.

The code is based on the model developed in the paper: 

'A Bayesian Decision-Theoretic Model of Sequential Experimentation with Delayed Response'

by Stephen E. Chick, Martin Forster, and Paolo Pertile (published in JRSSB and available at https://doi.org/10.1111/rssb.12225). 

The functions provided in the code are mostly the same as in the Matlab implementation (visualisations are not produced within DelayCurvesRecur and DelayStageOne functions). Instead, the main visualisation of interest visualises the stopping boundaries and is produced by the DelayFigure function.

The R implementation requires the user to install the htadelayR R package via devtools::install_github("michaeldymock25/htadelayR").

To estimate the model for a given trial design, the DelayInputConstructor function should be called (with optional trial specific inputs), followed by the DelayInputValidator function. The validated output should be passed to the DelayCurvesRecur function followed by the DelayStageOne function (including the intermediate output from DelayCurvesRecur).