#' Runs the full analysis
#'
#' @param X_labeled matrixof n*p
#' @param X_unlabeled matrix of N*p
#' @param S_labeled vector of n*1
#' @param S_unlabeled matrix of N*1
#' @param y Numeric outcome vector
#' @param samp_prob vector n*1
#' @param n_lab integer
#' @param num_knots number of knots for basis expansion (default is 3)
#' @param basis_type type of basis to use (default is 'NS_basis')
#' @param num_folds number of folds for cross-validations (default is 3)
#' @param my_threshold threshold to determine classification (default is 0.5)
#' @param reps number of replications for cross-validation (default is 10)
#' @param num_perts number of replications for perturbation (default is 500)
#' @return Perturbed estimates and corresponding residuals
#' @export
#' 

RunStratifiedSSL <- function(X_labeled, X_unlabeled, S_labeled, S_unlabeled, y, samp_prob, n_lab,
                             num_knots = 3, basis_type = 'NS_basis', num_folds = 3, my_threshold = 0.5, reps = 10, num_perts = 500){
  
  
  # Get basis expansion
  if(basis_type == 'NS_basis'){
    
    my_basis <- NaturalSplineBasis(rbind(X_labeled, X_unlabeled),
                                   c(S_labeled, S_unlabeled),
                                   num_knots = num_knots)
    
  }else{
    
    my_basis <- AlternativeBasis(rbind(X_labeled, X_unlabeled),
                                 c(S_labeled, S_unlabeled),
                                 num_knots = num_knots, basis_type = basis_type) 
    
  }
  
  basis_labeled <- my_basis[1:n_lab, ]
  basis_unlabeled <- my_basis[(n_lab+1):nrow(my_basis), ]
  
  # Fit the regression model
  regression_result <- SemiSupervisedRegression(basis_labeled,
                                                basis_unlabeled,
                                                X_labeled,
                                                X_unlabeled,
                                                y,
                                                samp_prob,
                                                lambda = NULL)
  
  # Supervised beta.
  beta_sl <- regression_result$beta_SL
  # Semi-supervised beta.
  beta_ssl <- regression_result$beta_SSL
  # Naive beta with no sampling probability.
  beta_naive <- regression_result$beta_SL_unweighted
  # Gamma from imputation.
  gamma <- regression_result$gamma
  # Beta from density ratio method.
  beta_dr <- regression_result$beta_DR
  # DR projection.
  proj_dr <- regression_result$proj_DR
  
  num_folds <- 3
  cv_residuals <- CrossValResids(basis_labeled, basis_unlabeled, X_labeled,
                                 X_unlabeled, y, samp_prob, num_folds)
  
  # Save the CV residuals.
  resids_sl <- cv_residuals$resids_beta_sl
  resids_ssl <- cv_residuals$resids_beta_ssl
  resids_dr <- cv_residuals$resids_beta_dr
  resids_gamma <- cv_residuals$resids_gamma
  
  # Standard error estimates for supervised and semi-supervised estimates.
  beta_sl_se_obj <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                       beta_sl, resids_sl)
  beta_ssl_se_obj <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                        beta_ssl, resids_gamma)
  beta_dr_se_obj <- StdErrorEstimationDR(X_labeled, X_unlabeled, y,
                                         beta_dr, resids_dr, proj_dr)
  
  # Supervised SE and inverse info.
  beta_sl_se <- beta_sl_se_obj$std_error
  beta_sl_inv_info <- beta_sl_se_obj$inverse_information
  # Semi-supervised SE and inverse info.
  beta_ssl_se <- beta_ssl_se_obj$std_error
  beta_ssl_inv_info <- beta_ssl_se_obj$inverse_information
  # Semi-supervised SE and inverse info.
  beta_dr_se <- beta_dr_se_obj
  
  # Minimum Variance Estimator (here the component-wise optimal estimator).
  beta_minvar <- SemiSupervisedMinVarRegression(beta_ssl, beta_sl,
                                                beta_ssl_se,
                                                beta_sl_se,
                                                resids_gamma,
                                                resids_sl,
                                                X_labeled,
                                                X_unlabeled,
                                                epsilon = NULL)
  
  # Save the results.
  beta_mv <- beta_minvar$beta_SSL_min_var
  mv_weight <- beta_minvar$min_var_weight
  beta_mv_se <- beta_minvar$se_beta_weight
  
  
  acc_sl <- SupervisedApparentAccuracy(X_labeled, y, beta_sl, samp_prob,
                                       resamp_weight = NULL,
                                       threshold = my_threshold)
  
  acc_ssl <- SemiSupervisedApparentAccuracy(basis_labeled, basis_unlabeled,
                                            X_labeled, X_unlabeled,
                                            y, beta_mv, gamma, samp_prob,
                                            resamp_weight = NULL,
                                            threshold = my_threshold)
  
  acc_dr <- SupervisedApparentAccuracy(X_labeled, y, beta_dr, samp_prob,
                                       resamp_weight = NULL,
                                       threshold = my_threshold)
  
  acc_naive <- SupervisedApparentAccuracy(X_labeled, y, beta_sl,
                                          rep(1, length(y)),
                                          resamp_weight = NULL,
                                          threshold = my_threshold)
  
  acc_ap_omr <- c(acc_ssl$omr_ssl, acc_sl$omr_sl, acc_dr$omr_sl, acc_naive$omr_sl)
  acc_ap_mse <- c(acc_ssl$mse_ssl, acc_sl$mse_sl, acc_dr$mse_sl, acc_naive$mse_sl)

  # acc_ap_omr
  # acc_ap_mse
  
  
  acc_cv <- CrossValAccuracy(basis_labeled, basis_unlabeled,
                             X_labeled, X_unlabeled, y, samp_prob,
                             mv_weight[,1],
                             num_folds, reps,
                             my_threshold, lambda0 = NULL)
  
  acc_cv_omr <- c(acc_cv$omr_ssl, acc_cv$omr_sl, acc_cv$omr_dr, acc_cv$omr_naive)
  acc_cv_mse <- c(acc_cv$mse_ssl, acc_cv$mse_sl, acc_cv$mse_dr, acc_cv$mse_naive)
  
  # acc_cv_omr
  # acc_cv_mse
  
  cv_weight <- num_folds / (2 * num_folds - 1)
  acc_en_omr <- cv_weight * acc_ap_omr + ((1-cv_weight) * acc_cv_omr)
  acc_en_mse <- cv_weight * acc_ap_mse + ((1-cv_weight) * acc_cv_mse)
  
  # acc_en_omr
  # acc_en_mse
  
  
  acc_pert <- AccuracyStdErrorEstimation(basis_labeled, basis_unlabeled,
                                         X_labeled, X_unlabeled, y,
                                         samp_prob,
                                         mv_weight[,1],
                                         beta_sl,
                                         beta_mv,
                                         beta_dr,
                                         resids_sl,
                                         resids_gamma,
                                         resids_dr,
                                         proj_dr,
                                         beta_ssl_inv_info,
                                         num_resamples = num_perts,
                                         threshold = my_threshold)
  
  # Save all the residuals.
  resids_all <- cbind(resids_ssl, resids_sl, resids_dr, resids_gamma)
  
   
  # Save all the standard errors.
  se_all <- cbind(beta_ssl_se,  beta_mv_se, beta_sl_se, beta_dr_se)
  
  # Save all the beta coefficients.
  beta_all <- cbind(beta_ssl, beta_mv, beta_sl, beta_dr, beta_naive)
  
  # Save results.
  ssl_pert_mse <- acc_pert$ssl_pert_mse
  sl_pert_mse <- acc_pert$sl_pert_mse
  dr_pert_mse <- acc_pert$dr_pert_mse
  
  ssl_pert_omr <- acc_pert$ssl_pert_omr
  sl_pert_omr <- acc_pert$sl_pert_omr
  dr_pert_omr <- acc_pert$dr_pert_omr
  
  # Save all the results.
  pert_mse_all <- cbind(ssl_pert_mse, sl_pert_mse, dr_pert_mse)
  pert_omr_all <- cbind(ssl_pert_omr, sl_pert_omr, dr_pert_omr)
  
  return(list(resids_sl = resids_sl,
              resids_ssl = resids_ssl,
              resids_dr = resids_dr,
              resids_gamma = resids_gamma,
              resids_all = resids_all,
              beta_mv = beta_mv,
              mv_weight = mv_weight,
              beta_mv_se = beta_mv_se,
              se_all = se_all,
              beta_all = beta_all,
              ssl_pert_mse = ssl_pert_mse,
              sl_pert_mse = sl_pert_mse,
              dr_pert_mse = dr_pert_mse,
              ssl_pert_omr = ssl_pert_omr,
              sl_pert_omr = sl_pert_omr,
              dr_pert_omr = dr_pert_omr,
              pert_mse_all = pert_mse_all,
              pert_omr_all = pert_omr_all
              
  ))
}

