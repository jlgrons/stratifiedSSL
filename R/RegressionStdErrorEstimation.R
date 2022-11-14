# Updated: 2021-08-30

#' Standard error estimation based on influence function expansion.
#'
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param beta Regression parameter estimate to compute the score function.
#' @param residual Estimated residual.
#' @return Standard error estimate.
#'
StdErrorEstimation <- function(X_labeled, X_unlabeled, y, beta, residual){

  X_labeled_intercept <- cbind(1, X_labeled)
  X_all_intercept <- cbind(1, rbind(X_labeled, X_unlabeled))

  n_all <- nrow(X_all_intercept)
  n_labeled <- nrow(X_labeled_intercept)

  score <-  t(X_labeled_intercept) %*% (
    X_labeled_intercept * c(residual)^2) / n_labeled

  expit_deriv <- c(ExpitDerivative(X_all_intercept %*% beta))
  # Ask Molei about this.
  expit_deriv[which(is.na(expit_deriv))] <- 0

  information <- t(X_all_intercept) %*% (X_all_intercept * expit_deriv) / n_all
  inverse_information <- solve(information)
  variance_est <-  inverse_information %*% score %*% inverse_information

  return(list(std_error = sqrt(diag(variance_est) / n_labeled),
              inverse_information = inverse_information))
}


#' Standard error estimation based on influence function expansion
#' for DR method.
#'
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param beta_dr DR regression parameter estimate to compute the score function.
#' @param resid Estimated residual.
#' @param proj_dr DR projection.
#' @return Standard error estimate for DR method.
#'
#'
StdErrorEstimationDR <- function(X_labeled, X_unlabeled, y, beta_dr,
                                 resid, proj_dr){

  X_labeled_intercept <- cbind(1, X_labeled)
  X_all <- cbind(1, rbind(X_labeled, X_unlabeled))

  n_all <- nrow(X_all)
  n_labeled <- nrow(X_labeled_intercept)

  U_dr <- diag(as.vector(resid)) %*% X_labeled_intercept
  S <-  t(U_dr - proj_dr) %*% (U_dr - proj_dr) /n_labeled

  ddd <- c(ExpitDerivative(X_all %*% beta_dr))
  ddd[which(is.na(ddd))] = 0

  I <- t(X_all) %*% (X_all * ddd)/n_all

  A <- solve(I)
  vars <-  A %*% S %*% A

  return(sqrt(diag(vars)/n_labeled))

}

