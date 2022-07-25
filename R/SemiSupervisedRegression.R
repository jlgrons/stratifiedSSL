# Updated: 2021-08-30

#' Semi-supervised regression.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param basis_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param lambda Penalization parameter for initial ridge estimator.
#' @export
#' @return Vector containing regression coefficients.
#' @importFrom stats coef glm quantile rbeta rbinom rlogis rnorm
SemiSupervisedRegression <- function(basis_labeled, basis_unlabeled, X_labeled,
                                     X_unlabeled, y, samp_prob, lambda = NULL){

  num_labeled <- length(y)
  p_basis <- ncol(basis_labeled)
  p_X <- ncol(X_labeled)

  # Default lambda if nothing input by the user.
  if(is.null(lambda)){lambda <- log(p_basis)/num_labeled^1.5}

  # Standardized weights.
  weights <- 1/samp_prob/mean(1/samp_prob)

  # Step 1: Basis function regression for imputation.
  gamma <- RidgeRegression(basis_labeled, y, weights = weights,
                              lambda = lambda)

  # Step 2: Semi-supervised regression using the imputations from Step 1.
  basis_all <- rbind(basis_labeled, basis_unlabeled)
  X_all <- rbind(X_labeled, X_unlabeled)
  basis_imps <- Expit(cbind(1, basis_all) %*% gamma);

  beta_SSL <- tryCatch(glm(basis_imps ~ X_all,
                          family = 'binomial')$coeff,
                      error = function(e) rep(NA, p_X+1))

  # Fit supervised estimate for comparison.
  beta_SL <- tryCatch(glm(y ~ X_labeled, family = 'binomial',
                         weights = weights)$coeff,
                     error = function(e) rep(NA, p_X+1))

  # Fit supervised estimate, without weights, for comparison.
  beta_SL_unweighted <- tryCatch(glm(y ~ X_labeled, family = 'binomial')$coeff,
                        error = function(e) rep(NA, p_X+1))

  # Fit the density ratio estimate for comparison.
  beta_DR <- DensityRatioRegression(basis_labeled, basis_unlabeled, X_labeled,
                                    X_unlabeled, y, samp_prob, lambda = NULL)

  return(list(beta_SSL = beta_SSL,
              beta_SL = beta_SL,
              gamma = gamma,
              beta_DR = beta_DR$beta_dr,
              beta_SL_unweighted = beta_SL_unweighted,
              proj_DR = beta_DR$proj_dr,
              lambda = lambda))
}
