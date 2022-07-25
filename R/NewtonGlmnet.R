
NewtonGlmnet <- function(X, y, weights, weights_mom, indx_mom,
                         lambda0 = 1e-04, max_iter = 100, tol = 1e-4,
                         initial = rep(0, 1 + ncol(X))){
  error <- Inf
  iter <- 0
  gamma <- initial
  X <- cbind(1, X)
  n <- nrow(X)

  sqloss <- mean(weights * (y - Expit(as.vector(X %*% gamma)))^2)
  indx_mom <- c(1, 1 + indx_mom)

  while(iter < max_iter & error > tol){

    iter <- iter + 1
    gamma_old <- gamma
    sqloss_old <- sqloss

    # Update the minimization:

    z <- as.vector(X %*% gamma)
    y_ <- y - Expit(z) + ExpitDerivative(z) * z
    x_ <- as.vector(ExpitDerivative(z)) * X
    xTx <- crossprod(x_, as.vector(weights) * x_) / n + lambda0 / 2 * diag(rep(1, ncol(x_)))
    xTy <- t(x_) %*% (as.vector(y_) * as.vector(weights)) / n
    C <- crossprod(X[,indx_mom], as.vector(weights_mom) * x_) / n
    b <- as.vector(t(X[,indx_mom]) %*% (as.vector(y_) * as.vector(weights_mom)) / n)
    mat_bind <- cbind(rbind(xTx, C), rbind(t(C), matrix(0, nrow(C), nrow(C))))
    vec_bind <- c(xTy, b)
    solution_bind <- solve(mat_bind) %*% vec_bind
    gamma <- solution_bind[1:ncol(X)]
    sqloss <- mean(weights * (y - Expit(as.vector(X %*% gamma)))^2)

    if (sqloss_old < sqloss){
      gamma <- gamma_old
    }
    error <- sqrt(mean((gamma - gamma_old)^2))
  }

  return(gamma)
}
