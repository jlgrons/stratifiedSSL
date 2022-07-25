
RidgeRegression <- function(X, y, weights = NULL, lambda = 1e-04,
                            family = 'binomial'){

  if(is.null(weights)){weights = rep(1, length(y))}

  # Multiple lambdas to make sure convergence happens.
  gamma = tryCatch((coef(glmnet(X, y, weights = weights, alpha = 0,
                                lambda = seq(lambda, 100*lambda, length.out = 100),
                                family = family))),
                   error = function(e) rep(NA, ncol(X) + 1 ))

  # Take the value at specified value of lambda.
  gamma = gamma[, ncol(gamma)]

  return(gamma)
}
