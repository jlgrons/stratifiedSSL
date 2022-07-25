# Updated: 2021-04-19

#' Component-wise minimum variance semi-supervised regression.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param basis_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param min_var_weight Minimum variance weight for semi-supervised estimate.
#' @param beta_SL Supervised regression coefficient vector.
#' @param beta_MV MinVar Semi-supervised regression coefficient vector.
#' @param beta_DR Density ratio regression coefficient vector.
#' @param resids_beta_SL Residuals from the supervised regression model.
#' @param resids_beta_imp Residuals from the imputation model.
#' @param resids_beta_dr Residuals from the density ratio estimator.
#' @param proj_dr Projection from density ratio estimator.
#' @param inverse_information Inverse information  matrix.
#' @param num_resamples Number of resamples.
#' @param threshold Threshold for over misclassification rate.
#' @export
#' @return Pertrubed estimates.
#' @importFrom stats coef glm quantile rbeta rbinom rlogis rnorm
AccuracyStdErrorEstimation <- function(basis_labeled, basis_unlabeled,
                                       X_labeled, X_unlabeled, y,
                                       samp_prob, min_var_weight,
                                       beta_SL, beta_MV, beta_DR,
                                       resids_beta_SL,
                                       resids_beta_imp, resids_beta_dr,
                                       proj_dr, inverse_information,
                                       num_resamples = 500, threshold = 0.5){
  set.seed(34)
  n_labeled <- length(y)
  resamp_weight <- sapply(1:num_resamples, function(kk) 4*rbeta(n_labeled,
                                                              0.5, 1.5))

  resids_beta_SL_weighted <- ((resamp_weight - 1) * resids_beta_SL)
  resids_beta_imp_weighted <- ((resamp_weight - 1) * resids_beta_imp)
  resids_beta_dr_weighted <- ((resamp_weight - 1) * resids_beta_dr)
  proj_dr_pert <- t(proj_dr) %*% (resamp_weight - 1)

  # Note: double check correctness of this.
  X_labeled_intercept <- cbind(1, X_labeled)
  IF_beta_imp <- inverse_information %*% t(
    X_labeled_intercept) %*% resids_beta_imp_weighted / n_labeled
  IF_beta_SL <- inverse_information %*% t(
    X_labeled_intercept) %*% resids_beta_SL_weighted / n_labeled


  beta_SSL_pert <- beta_MV + diag(1-min_var_weight) %*% IF_beta_SL + diag(
    min_var_weight) %*% IF_beta_imp
  beta_SL_pert <- beta_SL + IF_beta_SL

  T_1.dr.p <- inverse_information %*% t( X_labeled_intercept) %*% resids_beta_dr_weighted / n_labeled
  T_2.dr.p <- inverse_information %*% proj_dr_pert / n_labeled
  beta_DR_pert <- beta_DR + T_1.dr.p - T_2.dr.p

  beta_imp_pert <- sapply(1:num_resamples, function(kk){
    RidgeRegression(basis_labeled, y,
                    weights = resamp_weight[,kk] / samp_prob/mean(
                      resamp_weight[,kk] / samp_prob),
             lambda = log(ncol(basis_labeled)) / (n_labeled^1.5))
  })

  perturbations_sl <- lapply(1:num_resamples, function(kk){
    SupervisedApparentAccuracy(X_labeled, y, beta_SL_pert[,kk], samp_prob,
                               resamp_weight = resamp_weight[ ,kk],
                               threshold = threshold)})

  perturbations_dr <- lapply(1:num_resamples, function(kk){
    SupervisedApparentAccuracy(X_labeled, y, beta_DR_pert[,kk], samp_prob,
                               resamp_weight = resamp_weight[ ,kk],
                               threshold = threshold)})

  perturbations_ssl <- lapply(1:num_resamples, function(kk){
    SemiSupervisedApparentAccuracy(basis_labeled, basis_unlabeled,
                                   X_labeled, X_unlabeled,
                                   y, beta_SSL_pert[,kk],  beta_imp_pert[, kk],
                                   samp_prob,
                                   resamp_weight = resamp_weight[ ,kk],
                                   threshold = threshold)})

  ssl_pert_mse <- sapply(perturbations_ssl, "[[", 1)
  ssl_pert_omr <- sapply(perturbations_ssl, "[[", 2)

  sl_pert_mse <- sapply(perturbations_sl, "[[", 1)
  sl_pert_omr <- sapply(perturbations_sl, "[[", 2)

  dr_pert_mse <- sapply(perturbations_dr, "[[", 1);
  dr_pert_omr <- sapply(perturbations_dr, "[[", 2);

  return(list(ssl_pert_mse = ssl_pert_mse,
              sl_pert_mse = sl_pert_mse,
              dr_pert_mse = dr_pert_mse,
              ssl_pert_omr = ssl_pert_omr,
              sl_pert_omr = sl_pert_omr,
              dr_pert_omr = dr_pert_omr,
              beta_imp_pert = beta_imp_pert,
              beta_SSL_pert = beta_SSL_pert,
              beta_SL_pert= beta_SL_pert))

}
