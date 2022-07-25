# Updated: 2021-04-19

#' Computes the natural spline basis.
#'
#' @param X Covariate matrix.
#' @param S Vector with stratum id.
#' @param num_knots Number of knots.
#' @export
#' @return Matrix containing natural spline basis.
#' @importFrom stats coef glm quantile rbeta rbinom rlogis rnorm
#'
NaturalSplineBasis <- function(X, S, num_knots){

  X <- as.matrix(X)
  basis.X <- c()

  for(i in 1:ncol(X)){
    X_i <- X[,i]

    # Quantiles to determine appropriate knots.
    knots <- quantile(X_i, seq(0, 1, length = num_knots))

    # Changes quantiles if there aren't enough unique values.
    j <- 0

    while(length(unique(knots)) != num_knots){
      j <- j + 1

      knots <- unique(quantile(X_i, seq(0, 1, length = num_knots  + j)))
    }

    # Compute the natural spline basis.
    d_k <- (TruncatedCubic(X_i, knots[num_knots-1])) / (knots[num_knots] - knots[num_knots-1])
    evals <- sapply(1:(num_knots-2), function(ii){
    d_i <- (TruncatedCubic(X_i, knots[ii])) / (knots[num_knots] - knots[ii]);
        basis.new <- d_i - d_k})

    # Bind original variable and basis.
    basis.X = cbind(basis.X, cbind(X_i, evals))

  }

  # One hot encoding of stratification variable.
  basis.S <- OneHotEncoding(S)

  # Return basis including everything.
  basis <- cbind(basis.X, basis.S)

  return(basis)
}

