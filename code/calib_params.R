## calib_params.R: ACMG-calibration priors. Defines CALIB_PARAM only; sourced by
## calibrate.R and plot_calibration.R. Edit here to recalibrate.
##   alpha = prior P(pathogenic) per class; OP = odds-of-pathogenicity ladder base.
CALIB_PARAM <- list(deletion  = list(alpha = 0.046, OP = 1051),
                    insertion = list(alpha = 0.008, OP = 12397))
