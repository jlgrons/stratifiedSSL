# Updated: 2021-08-30

#' Apparent estimates for brier score (MSE) and misclassification rate (OMR).
#'
#' @param X_labeled Covariate matrix for labeled data set.
#' @param beta_SL Numeric vector of regression coefficients.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of regression coefficients.
#' @param samp_prob Numeric vector of weights.
#' @param resamp_weight Numeric vector of resampling weights.
#' @param threshold Threshold for overall misclassification rate.
#' @return Supervised MSE and OMR.
#'
SupervisedApparentAccuracy <- function(X_labeled, y, beta_SL, samp_prob,
                                           resamp_weight = NULL,
                                           threshold = 0.5){

  if(is.null(resamp_weight)){resamp_weight <- rep(1, length(y))}
  weight <- resamp_weight / samp_prob / mean(resamp_weight/samp_prob)

  MSE_SL <- MeanSquaredError(X_labeled, beta_SL, y, weight = weight)
  OMR_SL <- AbsoluteError(X_labeled, beta_SL, y, weight = weight,
                          threshold = threshold)

  return(list(mse_sl = MSE_SL, omr_sl = OMR_SL))
}


# Note: Can add these to simulation code.
# mse.naive = mean((Yt - lp.t.ssl)^2)
# ae.naive = mean(abs(Yt - lp.t.ssl.ind))

# mse.dr = mean((Yt - lp.t.dr)^2 * weight)
# ae.dr = mean(abs(Yt - lp.t.dr.ind) * weight)
