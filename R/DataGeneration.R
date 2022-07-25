
CovariateGen <- function(n, mu, Sigma){

  return(MASS::mvrnorm(n, mu, Sigma))

}


ARoneCovMat <- function(p, rho){

  id_matrix <- diag(p)
  ar_one_matrix <- rho^abs(row(id_matrix)-col(id_matrix))

  return(ar_one_matrix)
}


StratificationVar <- function(covariates, num_strata = 2){

  N <- length(covariates[ ,1])

  strat_var_1 <- I(covariates[ ,1] + rnorm(N) < 0.5)
  strat_var_3 <- I(covariates[ ,3] + rnorm(N) < 0.5)

  if (num_strata == 2){
    strat_var <- ifelse(strat_var_1, 1, 0)
  }

  if (num_strata == 4){
    strat_var <- rep(0, N)
    strat_var[intersect(which(strat_var_1 == 1), which(strat_var_3 == 1))] <- 1
    strat_var[intersect(which(strat_var_1 == 1), which(strat_var_3 == 0))] <- 2
    strat_var[intersect(which(strat_var_1 == 0), which(strat_var_3 == 1))] <- 3
  }

  return(strat_var)
}


#' Generate Data for Simulation Studies
#'
#' @param n_lab Size of labeled data.
#' @param n_unlab Size of unlabeled data.
#' @param p Number of covariates.
#' @param rho Covariance parameter for AR-1 covariance structure.
#' @param signal Values of nonzero regression parameters.
#' @param model_specification Choice of model specification.
#' @param num_strata Number of strata for stratified sampling.
#' @export
#' @return List of relevant data objects.
DataGeneration <- function(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                     model_specification = 'outcome_correct_imputation_correct',
                     num_strata = 2){

  ## model_specification:
  ## 'outcome_correct_imp_correct'  = 'CC'
  ## 'outcome_wrong_imp_correct'  = 'IC'
  ## 'outcome_wrong_imp_wrong'  = 'II'
  ## Supplement:
  ## 'outcome_wrong_imp_correct'  = 'IC1'
  ## 'outcome_wrong_imp_wrong'  = 'II1'
  ## Gaussian mixture

  # Total data size.
  N <- n_lab + n_unlab

  # Regression parameter.
  signal <- c(signal, rep(0, p - length(signal)))

  # Covariance for covariates.
  sigma <- 3*ARoneCovMat(p = p, rho = rho)

  # Covariates.
  covariates <- CovariateGen(N, mu = rep(0, p), sigma)

  # Linear predictor.
  lin_pred <- c(covariates %*% signal)

  # Outcome and stratification variable.
  if (model_specification == 'outcome_correct_imputation_correct'){

    strat_var <- StratificationVar(covariates, num_strata = num_strata)
    Y <- I(lin_pred + rlogis(N) > 0)

  }

  if (model_specification == 'outcome_incorrect_imputation_correct'){

    strat_var <- StratificationVar(covariates, num_strata = num_strata)

    if (num_strata == 2){
      gamma_coef <- c(signal, c(0.5, 0, 0, 0.5),
                      rep(0, 8), - 0.5, rep(0, 4), -0.5)
      basis <- InteractionBasis(covariates, strat_var)
    }

    if (num_strata == 4){
      gamma_coef <- c(signal, c(0.5, 0, 0, 0.5), rep(0, 8),
                      -0.5, rep(0, 4), -0.5, 0, 0)
      basis <- InteractionBasis(covariates, strat_var)
    }

    Y <- I(c(basis %*% gamma_coef) + rlogis(N) > 0)
  }


  if (model_specification == 'outcome_incorrect_imputation_correct_supp'){

    strat_var_1 <- I(covariates[,1] + covariates[,2] + rnorm(N) > 1.5)
    strat_var <- ifelse(strat_var_1, 1, 0)

    gamma_coef <- c(signal, c(0.5, 0, 0, 0.5), rep(0, 8),
                    -0.5, rep(0, 4), -0.5)
    basis <- InteractionBasis(covariates, strat_var)

    Y0 <- basis %*% gamma_coef * strat_var +
      (0.8 * basis %*% gamma_coef - 5) * (1 - strat_var) + rlogis(N)
    Y <- I(Y0 > 1)

  }


  if (model_specification == 'outcome_incorrect_imputation_incorrect'){

    strat_var <- StratificationVar(covariates, num_strata = num_strata)
    incorrect_signal <- c(-2, rep(0,3), -3, -3, 0, 0, rep(0,2))
    Y0 <- lin_pred + covariates[,1]^2 + covariates[,3]^2 +
      rgumbel(N, -2, 0.3)*exp(covariates%*% incorrect_signal)
    Y <- I(Y0 > 0)
  }


  if (model_specification == 'outcome_incorrect_imputation_incorrect_supp'){

    strat_var_1 <- I(covariates[,1] + covariates[,2] + rnorm(N) > 1.5)
    strat_var <- ifelse(strat_var_1, 1, 0)

    incorrect_signal = c(-2, rep(0,3), -3, -3, 0, 0, rep(0,2))
    Y0 <- (lin_pred + covariates[,1]^2 + covariates[,3]^2) * strat_var +
      (0.8 * lin_pred - 5)* (1 - strat_var) +
      rgumbel(N, -2, 0.3)*exp(c(covariates%*% incorrect_signal))
    Y <- I(Y0 > 1)
  }


  if (model_specification == 'gaussian_mixture'){

    mu_diff <- c(0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.1, -0.1, rep(0, p - 8))
    sigma_diff <- ARoneCovMat(p, rho = 0.3) - diag(rep(0.4, p)) + matrix(0.2, p, p)

    Y <- rbinom(N, 1, 0.5)
    covariates <- matrix(0, N, p)

    mu1 <- rep(0, p)
    sigma1 <- ARoneCovMat(p = p, rho = rho)
    covariates[which(Y == 0), ] <- CovariateGen(length(which(Y == 0)),
                                                mu1, sigma1)

    mu2 <- mu1 + mu_diff
    sigma2 <- sigma1 + sigma_diff
    covariates[which(Y == 1), ] <- CovariateGen(length(which(Y == 1)),
                                                mu2, sigma2)

    strat_var <- StratificationVar(covariates, num_strata = num_strata)

    covariates[which(Y == 1), 3] <- covariates[which(Y == 1), 3] +
      0.12 * covariates[which(Y == 1), 3]^3
    covariates[which(Y == 1), 4] <- covariates[which(Y == 1), 4] +
      0.12 * covariates[which(Y == 1), 4]^3

    covariates[which(Y == 1), 7] <- covariates[which(Y == 1), 7] +
      0.12 * covariates[which(Y == 1), 7]^3
    covariates[which(Y == 1), 8] <- covariates[which(Y == 1), 8] +
      0.12 * covariates[which(Y == 1), 8]^3


  }

  # Size of strata
  n_each <- n_lab/num_strata

  if (num_strata == 2){

    # Stratum 1,
    stratum_1 <- which(strat_var == 1)
    ind_1 <- sample(stratum_1, size = n_each)

    # Stratum 2.
    stratum_2 <- which(strat_var == 0)
    ind_2 <- sample(stratum_2, size = n_each)

    # Indices of labeled and unlabeled data.
    ind_lab <- c(ind_1, ind_2)
    ind_unlab <- setdiff(1:N, ind_lab)

    # Full data.
    dat_full <- cbind(Y, covariates)

    # Name labeled and unlabeled data sets.
    dat_lab <- dat_full[ind_lab, ]
    S_lab <- strat_var[ind_lab]
    dat_unlab <- dat_full[ind_unlab, ]
    S_unlab <- strat_var[ind_unlab]
    S_full <- c(S_lab, S_unlab)

    # Random sampling (only used in Section S4).
    dat_random_samp <- dat_full[1:n_lab, ]
    Y_random_samp <- dat_random_samp[,1]
    covariates_random_samp <- dat_random_samp[,-1]
    covariates_unlab_random_samp <- dat_full[-c(1:n_lab), -1]

    # Relevant labeled data.
    Y_lab <- dat_lab[,1]
    covariates_lab <- dat_lab[,-1]

    # Relevant unlabeled data.
    covariates_unlab <- dat_unlab[,-1]

    # Sampling probabilities.
    ind_S1_unlab <- (which(strat_var > 0))
    ind_S2_unlab  <- (which(strat_var <= 0))
    N_1 <- length(ind_S1_unlab)
    N_2 <- length(ind_S2_unlab)
    samp_prob <- c(rep(n_each/N_1, n_each), rep(n_each/N_2, n_each))
  }

  if (num_strata == 4){

    # Stratum 1.
    stratum_1 <- which(strat_var == 0)
    ind_1 <- sample(stratum_1, size = n_each)

    # Stratum 2.
    stratum_2 <- which(strat_var == 1)
    ind_2 <- sample(stratum_2, size = n_each)

    # Stratum 3.
    stratum_3 <- which(strat_var == 2)
    ind_3 <- sample(stratum_3, size = n_each)

    # Stratum 4.
    stratum_4 <- which(strat_var == 3)
    ind_4 <- sample(stratum_4, size = n_each)

    # Indices of labeled and unlabeled data.
    ind_lab <- c(ind_1, ind_2, ind_3, ind_4)
    ind_unlab <- setdiff(1:N, ind_lab)

    # Full data.
    dat_full <- cbind(Y, covariates)

    # Labeled and Unlabeled data sets.
    dat_lab <- dat_full[ind_lab, ]
    S_lab <- strat_var[ind_lab]
    dat_unlab <- dat_full[ind_unlab, ]
    S_unlab <- strat_var[ind_unlab]
    S_full <- c(S_lab, S_unlab)

    # Random sampling (only used in Section S4).
    dat_random_samp <- dat_full[1:n_lab, ]
    Y_random_samp <- dat_random_samp[,1]
    X_random_samp <- dat_random_samp[,-1]
    X_unlab_random_samp <- dat_full[-c(1:n_lab), -1]

    # Relevant labeled data.
    Y_lab <- dat_lab[,1]
    X_lab <- dat_lab[,-1]

    # Relevant unlabeled data.
    X_unlab <- dat_unlab[,-1]

    # Sampling probabilities.
    ind_S1_unlab <- (which(strat_var == 0))
    ind_S2_unlab  <- (which(strat_var == 1))
    ind_S3_unlab <- (which(strat_var == 2))
    ind_S4_unlab  <- (which(strat_var == 3))

    N_1 <- length(ind_S1_unlab)
    N_2 <- length(ind_S2_unlab)
    N_3 <- length(ind_S3_unlab)
    N_4 <- length(ind_S4_unlab)

    samp_prob = c(rep(n_each/N_1, n_each), rep(n_each/N_2, n_each),
                  rep(n_each/N_3, n_each), rep(n_each/N_4, n_each))

  }

  return(list(covariates = covariates, Y = Y, strat_var = S_full,
              S_lab = S_lab, S_unlab = S_unlab, covariates_lab = covariates_lab,
              covariates_unlab = covariates_unlab, Y_lab = Y_lab,
              samp_prob = samp_prob, signal = signal,
              covariates_random_samp = covariates_random_samp,
              Y_random_samp = Y_random_samp,
              S_random_samp = strat_var,
              covariates_unlab_random_samp = covariates_unlab_random_samp))
}
