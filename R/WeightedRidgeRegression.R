# Updated: 2021-08-30


#' @importFrom stats coef glm quantile rbeta rbinom rlogis rnorm
WeightedRidgeRegression <- function(X, y, weights, weights_mom, indx_mom,
                            lambda0 = 1e-04, initial = rep(0, 1 + ncol(X))){

  if(is.null(weights)){weights = rep(1, length(y))}

  # multiple lambdas to make sure convergence happens
  gamma = NewtonGlmnet(X, y, weights, weights_mom, indx_mom, lambda0 =lambda0,
                        initial = initial)

  return(gamma)

}
