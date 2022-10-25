# Updated: 2021-04-19

#' Computes the natural spline basis.
#'
#' @param X Covariate matrix.
#' @param S Vector with stratum id.
#' @param num_knots Number of knots.
#' @param basis_type Type of basis to use.
#' @return Matrix containing basis.
AlternativeBasis <- function(X, S, num_knots, basis_type = 'interact'){

  if (basis_type == 'interact'){

    basis <- InteractionBasis(X, S)

  }

  # Will change name later.
  if (basis_type == 'IC1'){

    basis <- InteractionBasis(X, S, S_interaction = TRUE)

  }

  # Supplemental basis generation.
  if (basis_type == 'II1'){

    natural_spline_basis <- NaturalSplineBasis(X, S, num_knots)
    interaction_basis <- TwoWayInteractionBasis(X, S)

    basis <- cbind(natural_spline_basis, interaction_basis)
  }


  return(basis)
}
