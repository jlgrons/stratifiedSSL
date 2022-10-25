#' Simulated data
#'
#' A dataset used to replicate results from the publication. We assume correct outcome and correct imputation
#' Generated using the function
#' simulation_data <- DataGeneration(n_lab = 400, n_unlab = 20000, p = 10, rho = 0.4 signal = c(1, 1, 0.5, 0.5), model_specification = 'outcome_correct_imputation_correct',num_strata = 2)
#' \itemize{
#'   \item covariates_lab. labeled set X
#'   \item covariates_unlab. unlabeled set X
#'   \item S_lab. labeled S
#'   \item S_unlab. unlabeled S
#'   \item Y_lab. outcome 
#'   \item samp_prob. probabilities
#' }
#' @docType data
#' @keywords datasets
#' @name simulation_data
#' @usage data(simulation_data)
#' @format A list of parameters and vectors
#' @export
NULL