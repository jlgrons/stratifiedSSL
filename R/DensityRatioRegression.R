
Newton_DR <- function(phi.t, phi.v, weights,
                      max.iter = 100, tol = 1e-5, initial = rep(0, ncol(phi.t)),
                      lambda0 = NULL){

  n.t <- nrow(phi.t);
  n.v <- nrow(phi.v);
  if(is.null(lambda0)){lambda0 = log(ncol(phi.t))/n.t^1.5};

  error <- Inf
  iter <- 0
  gamma <- initial
  phi.v.mean <- rowMeans(t(phi.v))

  while(iter < max.iter & error > tol){
    gamma_old <- gamma
    z <- as.vector(phi.t %*% gamma)
    w <- exp(z) * weights
    phiT.phi <- crossprod(phi.t, w * phi.t) / n.t
    E.phi <- t(phi.t) %*% w / n.t

    gamma <- gamma + solve(phiT.phi + diag(rep(lambda0, ncol(phiT.phi)))) %*% (phi.v.mean - E.phi)
    error <- sqrt(mean((gamma - gamma_old)^2))
  }
  return(gamma)

}


#' Density ratio regression_
#'
#' @param basis_labeled Basis matrix for labeled data set_
#' @param basis_unlabeled Basis matrix for unlabeled data set_
#' @param X_labeled Covariate matrix for labeled data set_
#' @param X_unlabeled Covariate matrix for unlabeled data set_
#' @param y Numeric outcome vector_
#' @param samp_prob Numeric vector of weights_
#' @param lambda Penalization parameter for initial ridge estimator_
#' @export
#' @return Vector containing regression coefficients_
#'
DensityRatioRegression <- function(basis_labeled, basis_unlabeled, X_labeled,
                                     X_unlabeled, y, samp_prob, lambda = NULL){

  p_X <- ncol(X_labeled)
  num_labeled <- nrow(basis_labeled)
  basis_all <- rbind(basis_labeled, basis_unlabeled)

  # A typo fixed in this line:

  n_all <- length(basis_all[,1])

  # Standardized weights.
  weights <- 1/samp_prob/mean(1/samp_prob)

  basis_phi <- cbind(rep(1, n_all), basis_all)
  phi_labeled <- basis_phi[1:num_labeled, ]
  phi_unlabeled <- basis_phi[-c(1:num_labeled), ]
  dim_basis <- ncol(basis_phi)

  ##### I made a change here on solving for the density ratio:
  gamma.dr <- Newton_DR(phi_labeled, basis_phi, weights, lambda0 = lambda)

  density_ratio <- as.vector(exp(phi_labeled %*% gamma.dr))
  weights_dr <- density_ratio * weights

  #phiT_phi <- t(basis_phi) %*% basis_phi / (num_labeled + n_all)
  #E_phi <- t(phi_labeled) %*% weights / num_labeled

  #theta_ratio_1 <- solve(phiT_phi + diag(rep(0.01, ncol(phiT_phi))))
  #theta_ratio_2 <- (rowMeans(t(phi_unlabeled)) - E_phi)
  #theta_ratio <- theta_ratio_1 %*% theta_ratio_2

  #density_ratio <- exp(phi_labeled %*% theta_ratio)
  #weights_dr <- density_ratio * weights

  beta_dr <- tryCatch(glm(y~X_labeled, family = 'binomial',
                          weights = weights_dr)$coeff,
                      error = function(e) rep(NA, p_X + 1));

  # I made a change on the way to calculate phiT_phi:

  phiT_phi <- crossprod(phi_labeled, weights_dr * phi_labeled) / num_labeled

  # "as_vector" fixed as "as.vector", "g_logit" fixed as "Expit":

  u_dr_1 <- diag(as.vector(y - Expit(cbind(1, X_labeled) %*% beta_dr)))

  u_dr <- t(u_dr_1 %*% cbind(1, X_labeled))
  E_uT_phi <- u_dr %*% diag(as.vector(weights)) %*% phi_labeled / sum(weights)

  # May replace this 0.01 with lambda, should give a value to lambda if it is set as NULL.

  proj_coef_dr <- E_uT_phi %*% solve(phiT_phi + diag(rep(0.01, ncol(phiT_phi))))
  proj_dr <- diag(as.vector(weights)) %*% phi_labeled %*% t(proj_coef_dr)

  # Add the estimated "density_ratio" as one of the returned outcomes:

  return(list(beta_dr = beta_dr, proj_dr = proj_dr, DR_est = density_ratio))
}



