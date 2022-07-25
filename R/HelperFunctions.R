# Updated: 2021-04-22
# Note: This now has to many functions and needs to be split up.

# -----------------------------------------------------------------------------

Logit <- function(x) {
  return(log(x / (1 - x)))
}


Expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

# -----------------------------------------------------------------------------

ExpitDerivative <- function(x, na_correction = T) {
  expit_deriv <- exp(x)/(1+exp(x))^2
  if(na_correction){
    expit_deriv[which(is.na(expit_deriv))] = 0
  }
  return(expit_deriv)
}

# -----------------------------------------------------------------------------

TruncatedCubic <- function(x, knot_location){
  return(((x > knot_location) * (x-knot_location))^3)
}

# -----------------------------------------------------------------------------

OneHotEncoding <- function(s){

  one_hot_encoding <- c()
  num_categories <- length(unique(s))

  for(k in 1:(num_categories - 1)){
    one_hot_encoding <- cbind(one_hot_encoding, ifelse(s == (k - 1), 1, 0))
  }

  return(one_hot_encoding)
}

# -----------------------------------------------------------------------------

InteractionBasis <- function(X, S, S_interaction = FALSE){

  basis <- X

  p <- length(X[1,])

  for (j in 2:p){
    basis <- cbind(basis, X[,1] * X[,j])
  }

  for (j in 3:p){
    basis <- cbind(basis, X[,2] * X[,j])
  }



  if(S_interaction){

      basis.S <- S
    for (k in 1:length(basis[1, ])) {
      basis.S <- cbind(basis.S, S * basis[,k])

    }}else{

      # Check with Molei about this for IC supp.
      basis.S <- OneHotEncoding(S)

    }


  basis <- cbind(basis, basis.S)
  return(basis)
}
# -----------------------------------------------------------------------------


TwoWayInteractionBasis <- function(X, S){

  basis <- c()

  for (k in 1:length(X[1, ])) {
    basis <- cbind(basis, S * X[,k])
  }

  return(basis)
}

# -----------------------------------------------------------------------------

MeanSquaredError <- function(X, beta, y, weight = NULL){

  if(is.null(weight)){weight <- rep(1, length(y))}

  pred_prob <- Expit(cbind(1, X) %*% beta)

  return(mean(((y - pred_prob)^2) * weight))
}

# -----------------------------------------------------------------------------


AbsoluteError <- function(X, beta, y, weight = NULL, threshold = NULL){

  if(is.null(weight)){weight <- rep(1, length(y))}

  if(is.null(threshold)){
    pred_prob <- Expit(cbind(1, X) %*% beta)
  }else{
    pred_prob <- I(Expit(cbind(1, X) %*% beta) > threshold)
  }

  return(mean(abs(y - pred_prob) * weight))
}

