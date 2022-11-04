# Updated: 2021-08-30

#' CV estimates for brier score (MSE) and misclassification rate (OMR).
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param basis_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param min_var_weight Numeric vector of minimum variance weights.
#' @param num_folds Scalar indicating number of folds for CV.
#' @param reps Scalar indicating number of repitions for CV.
#' @param threshold Threshold for overall misclassification rate.
#' @param lambda0 Initial lambda for imputation model.
#' @export
#' @return CV semi-supervised and supervised MSE and OMR.
#'
# computes the CV estimates for model evaluation parameters
CrossValAccuracy <- function(basis_labeled, basis_unlabeled,
                              X_labeled, X_unlabeled, y, samp_prob,
                              min_var_weight, num_folds = 3, reps = 10,
                              threshold = 0.5, lambda0 = NULL){

  pp <- ncol(basis_labeled)
  n_labeled <- nrow(X_labeled)

  mse_cv_ssl <- matrix(NA, nrow = num_folds, ncol = reps)
  omr_cv_ssl <- matrix(NA, nrow = num_folds, ncol = reps)

  mse_cv_sl <- matrix(NA, nrow = num_folds, ncol = reps)
  omr_cv_sl <- matrix(NA, nrow = num_folds, ncol = reps)

  mse_cv_naive <- matrix(NA, nrow = num_folds, ncol = reps)
  omr_cv_naive <- matrix(NA, nrow = num_folds, ncol = reps)

  mse_cv_dr <- matrix(NA, nrow = num_folds, ncol = reps)
  omr_cv_dr <- matrix(NA, nrow = num_folds, ncol = reps)


  if(is.null(lambda0)){lambda0 = log(pp) /
    (floor((num_folds-1)*n_labeled/num_folds))^1.5}

  for(j in 1:reps){
    # set.seed(j)
    ind.cv <- split(1:n_labeled,
                    sample(rep(1:num_folds, floor(n_labeled/num_folds))))

    for(i in 1:num_folds){

      inds_val <- as.vector(unlist(ind.cv[i]))
      inds_tr <- as.vector(setdiff(as.vector(unlist(ind.cv)), inds_val))

      wg_val <- samp_prob[inds_val]
      wg_tr <- samp_prob[inds_tr]

      y_val <- y[inds_val]
      X_labeled_val <- X_labeled[inds_val, ]
      basis_val <- basis_labeled[inds_val, ]

      y_tr <- y[inds_tr]
      X_labeled_tr <- X_labeled[inds_tr, ]
      basis_tr <- basis_labeled[inds_tr, ]

      beta_tr <- SemiSupervisedRegression(basis_tr, basis_unlabeled,
                                          X_labeled_tr, X_unlabeled,
                                          y_tr, wg_tr, lambda = lambda0)

      beta_ssl_tr <- beta_tr$beta_SSL
      beta_sl_tr <- beta_tr$beta_SL
      gamma_tr <-  beta_tr$gamma
      beta_mv_tr <- min_var_weight*beta_ssl_tr +
        (1-min_var_weight)*beta_sl_tr
      beta_dr_tr <- beta_tr$beta_DR
      beta_naive_tr <- beta_tr$beta_SL_unweighted

      # Supervised estimates.
      acc_sl_val <- SupervisedApparentAccuracy(X_labeled_val, y_val,
                                               beta_sl_tr, wg_val,
                                               resamp_weight = NULL,
                                               threshold)

      mse_cv_sl[i,j] <- acc_sl_val$mse_sl
      omr_cv_sl[i,j] <- acc_sl_val$omr_sl

      # Semi-supervised estimates.
      acc_ssl_val <- SemiSupervisedApparentAccuracy(basis_labeled,
                                                    basis_unlabeled,
                                                    X_labeled,
                                                    X_unlabeled,
                                                    y,
                                                    beta_mv_tr,
                                                    gamma_tr,
                                                    samp_prob,
                                                    inds_val,
                                                    resamp_weight = NULL,
                                                    threshold)

      mse_cv_ssl[i,j] <- acc_ssl_val$mse_ssl
      omr_cv_ssl[i,j] <- acc_ssl_val$omr_ssl

      # Naive supervised estimates.
      acc_naive_val <- SupervisedApparentAccuracy(X_labeled_val, y_val,
                                               beta_sl_tr,
                                               rep(1, length(y_val)),
                                               resamp_weight = NULL,
                                               threshold = threshold)

      mse_cv_naive[i,j] <- acc_naive_val$mse_sl
      omr_cv_naive[i,j] <- acc_naive_val$omr_sl

      # DR estimates.
      acc_dr_val <- SupervisedApparentAccuracy(X_labeled_val, y_val,
                                               beta_dr_tr,
                                               wg_val,
                                               resamp_weight = NULL,
                                               threshold = threshold)

      mse_cv_dr[i,j] = acc_dr_val$mse_sl
      omr_cv_dr[i,j] = acc_dr_val$omr_sl
    }
  }

  return(list(mse_ssl = mean(mse_cv_ssl, na.rm = T),
              omr_ssl = mean(omr_cv_ssl, na.rm = T),
              mse_sl = mean(mse_cv_sl, na.rm = T),
              omr_sl = mean(omr_cv_sl, na.rm = T),
              mse_dr = mean(mse_cv_dr, na.rm = T),
              omr_dr = mean(omr_cv_dr, na.rm = T),
              mse_naive = mean(mse_cv_naive, na.rm = T),
              omr_naive = mean(omr_cv_naive, na.rm = T),
              mse_ssl_all = mse_cv_ssl,
              omr_ssl_all = omr_cv_ssl))

}
