# Updated: 2021-04-22

#' Apparent estimates for brier score (MSE) and misclassification rate (OMR).
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param basis_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param beta_SSL Numeric vector of regression coefficients.
#' @param beta_imp Numeric vector of regression coefficients for imputation.
#' @param samp_prob Numeric vector of weights.
#' @param ind_est Optional numeric vector to indicate inds for refitting model.
#' @param resamp_weight Numeric vector of resampling weights.
#' @param threshold Threshold for overall misclassification rate.
#' @export
#' @return Semi-supervised MSE and OMR.
#'
SemiSupervisedApparentAccuracy <- function(basis_labeled, basis_unlabeled,
                                           X_labeled, X_unlabeled, y, beta_SSL,
                                           beta_imp, samp_prob, ind_est = NULL,
                                           resamp_weight = NULL,
                                           threshold = 0.5){

  if(is.null(resamp_weight)){resamp_weight <- rep(1, length(y))}
  if(is.null(ind_est)){ind_est <- 1:length(y)}

  weight <- resamp_weight[ind_est] / samp_prob[ind_est] / mean(
    resamp_weight[ind_est] / samp_prob[ind_est])

  MSE_SSL <- SemiSupervisedAccuracy(basis_labeled, basis_unlabeled,
                                    X_labeled, X_unlabeled, y, beta_SSL,
                                    beta_imp, ind_est,
                                    weight = weight,
                                    threshold = NULL)

  OMR_SSL <- SemiSupervisedAccuracy(basis_labeled, basis_unlabeled,
                                    X_labeled, X_unlabeled, y, beta_SSL,
                                    beta_imp, ind_est,
                                    weight = weight,
                                    threshold = threshold)

  return(list(mse_ssl = MSE_SSL, omr_ssl = OMR_SSL))
}
