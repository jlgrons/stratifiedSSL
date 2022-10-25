# Updated: 2021-04-19

#' Component-wise minimum variance semi-supervised regression.
#'
#' @param beta_SSL Semi-supervised regression coefficient vector.
#' @param beta_SL Supervised regression coefficient vector.
#' @param beta_SSL_se Semi-supervised regression standard error.
#' @param beta_SL_se Supervised regression standard error.
#' @param resids_imp Residuals from the imputation model.
#' @param resids_beta_SL Residuals from the supervised regression model.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param epsilon Small offset to help with numerical stability.
#' @return List with parameter and SE est as well as estimated min_var_weights.
#'
SemiSupervisedMinVarRegression <- function(beta_SSL, beta_SL,
                                           beta_SSL_se, beta_SL_se,
                                           resids_imp, resids_beta_SL,
                                           X_labeled, X_unlabeled,
                                           epsilon = NULL){

  X_all <- cbind(1, rbind(X_labeled, X_unlabeled))
  X_labeled_int <- cbind(1, X_labeled)

  n_all <- nrow(X_all)
  n_labeled <- nrow(X_labeled_int)
  ones <- c(1, 1)
  p <- ncol(X_labeled)
  min_var_weight <- matrix(NA, p + 1, 2)

  if(is.null(epsilon)){epsilon <-  (n_labeled*(beta_SL_se^2 +
                                                 beta_SSL_se^2))/(2*n_labeled^0.6)}

  # Compute the minimum variance estimator.
  # Note: The residuals have been divided by mean of min_var_weights.
  SSL_pred_prob_deriv <- c(ExpitDerivative(X_all %*% beta_SSL))
  SSL_info_matrix <- solve(t(X_all) %*% (X_all * SSL_pred_prob_deriv))
  SSL_scaled_info_matrix <- SSL_info_matrix * n_all

  imp_IF <- SSL_scaled_info_matrix %*% t(X_labeled_int * c(resids_imp))
  beta_SL_IF <- SSL_scaled_info_matrix %*% t(X_labeled_int * c(resids_beta_SL))

  for (i in 1:(p + 1)){
    comb_IF <- cbind(imp_IF[i,], beta_SL_IF[i,])
    cov <- solve(t(comb_IF)%*%comb_IF  + epsilon[i]*diag(2)) / n_labeled
    denominator <- t(ones) %*% cov %*% ones
    numerator <- t(ones) %*% cov
    min_var_weight[i,] <- numerator / c(denominator)

    # Check with Molei about this.
    if (NA %in% min_var_weight[i,]){
      min_var_weight[i,] <- c(1, 0)
    }
  }

  # Componentwise minimum variance estimator.
  beta_SSL_min_var <- beta_SSL * min_var_weight[,1] +
    beta_SL * min_var_weight[,2]

  # Standard error of minimum variance estimator.
  SSL_mv_pred_prob_deriv <- c(ExpitDerivative(X_all %*% beta_SSL_min_var))
  SSL_mv_info_matrix <- solve(t(X_all) %*% (X_all * SSL_mv_pred_prob_deriv))
  SSL_mv_scaled_info_matrix <- SSL_mv_info_matrix * n_all

  beta_SSL_mv_weight <- diag(min_var_weight[,1])
  beta_SL_mv_weight <-  diag(min_var_weight[,2])

  resids_mv_weight <- beta_SSL_mv_weight %*% t(X_labeled_int*c(resids_imp)) +
    beta_SL_mv_weight %*% t(X_labeled_int*c(resids_beta_SL))
  scaled_resids_mv_weight <- resids_mv_weight / n_labeled

  se_beta_weight <- sqrt(diag(
    SSL_mv_scaled_info_matrix %*%
      (scaled_resids_mv_weight %*% t(scaled_resids_mv_weight)) %*%
      SSL_mv_scaled_info_matrix))

  return(list(beta_SSL_min_var = beta_SSL_min_var,
              min_var_weight = min_var_weight,
              se_beta_weight = se_beta_weight))
}

