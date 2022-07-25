# Updated: 2021-04-22

#' Apparent estimates for brier score (MSE) and misclassification rate (OMR).
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param basis_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param num_folds Number of folds for cross-validation.
#' @param lambda Regularization parameter for the imputation model.
#' @export
#' @return Cross-validated residuals.
#'

CrossValResids <- function(basis_labeled, basis_unlabeled, X_labeled,
                           X_unlabeled, y, samp_prob, num_folds, lambda = NULL){

  n_labeled <- nrow(X_labeled)
  p_basis <- ncol(basis_labeled)
  ind_cv <- split(1:n_labeled, sample(rep(1:num_folds,
                                          floor(n_labeled / num_folds))))

  if(is.null(lambda)){
    lambda <- log(p_basis) / (floor((num_folds - 1) * n_labeled/num_folds))^1.5
    }

  # Regression estimates with the kth fold removed.
  cv_ests <- lapply(1:num_folds, function(kk) {
    inds_fold <- as.vector(unlist(ind_cv[-kk]))
    y_fold <- y[inds_fold]
    samp_prob_fold <- samp_prob[inds_fold]
    SemiSupervisedRegression(basis_labeled[inds_fold,], basis_unlabeled,
                             X_labeled[inds_fold,], X_unlabeled, y_fold,
                             samp_prob_fold, lambda = lambda)
  })

  beta_ssl_cv <- sapply(cv_ests, "[[", 1)
  beta_sl_cv <- sapply(cv_ests, "[[", 2)
  gamma_cv <- sapply(cv_ests, "[[", 3)
  beta_dr_cv <- sapply(cv_ests, "[[", 4)

  # Residuals based on the kth fold and the beta with the kth fold removed.
  res_cv <- lapply(1:num_folds, function(kk){
    inds_fold <- as.vector(unlist(ind_cv[kk]))
    beta_ssl_fold <- beta_ssl_cv[,kk]
    beta_sl_fold <- beta_sl_cv[,kk]
    gamma_fold <- gamma_cv[,kk]
    beta_dr_fold <- beta_dr_cv[,kk]
    samp_prob_fold <- samp_prob[inds_fold]
    y_fold <- y[inds_fold]
    pred_prob_ssl <- Expit(cbind(1, X_labeled[inds_fold, ]) %*% beta_ssl_fold)
    pred_prob_sl <- Expit(cbind(1, X_labeled[inds_fold, ]) %*% beta_sl_fold)
    pred_imp <- Expit(cbind(1, basis_labeled[inds_fold, ]) %*% gamma_fold)
    pred_prob_dr <- Expit(cbind(1, X_labeled[inds_fold, ]) %*% beta_dr_fold)
    resids <- list(beta_ssl = (y_fold - pred_prob_ssl) * (1 / samp_prob_fold),
                beta_sl = (y_fold - pred_prob_sl) * (1 / samp_prob_fold),
                gamma = (y_fold - pred_imp) * (1 / samp_prob_fold),
                beta_dr = (y_fold - pred_prob_dr) * (1 / samp_prob_fold))
  })

  inds_cv_ordered <- order(unlist(ind_cv))
  mean_weight <- mean(1/samp_prob);

  # Reorder the residuals and divide by the mean of the weights.
  resids_beta_ssl <- c(unlist(sapply(res_cv,
                                     "[[", 1)))[inds_cv_ordered] / mean_weight
  resids_beta_sl <- c(unlist(sapply(res_cv,
                                    "[[", 2)))[inds_cv_ordered] / mean_weight
  resids_gamma <- c(unlist(sapply(res_cv,
                                     "[[", 3)))[inds_cv_ordered] / mean_weight
  resids_beta_dr <- c(unlist(sapply(res_cv,
                                    "[[", 4)))[inds_cv_ordered] / mean_weight

  return(list(resids_beta_ssl = resids_beta_ssl,
              resids_beta_sl = resids_beta_sl,
              resids_gamma = resids_gamma,
              resids_beta_dr = resids_beta_dr,
              beta_ssl_cv = beta_ssl_cv,
              beta_sl_cv = beta_sl_cv,
              gamma_cv = gamma_cv,
              beta_dr_cv = beta_dr_cv))
}



