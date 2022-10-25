
test_that("The full stratifiedSSL analysis works", {
  data(simulation_data)
  n_lab <- 400
  n_unlab <- 20000
  p <- 10
  rho <- 0.4
  num_strata <- 2
  X_labeled <- simulation_data$covariates_lab
  X_unlabeled <- simulation_data$covariates_unlab
  S_unlabeled <- simulation_data$S_unlab
  S_labeled <- simulation_data$S_lab
  y <- simulation_data$Y_lab
  samp_prob <- simulation_data$samp_prob
  basis_type <- 'NS_basis'
  reps = 2
  num_perts = 2
  num_knots <- 3
  basis_type <- 'NS_basis'
  my_threshold <- 0.5
  
  stratifiedSSL_results <- RunStratifiedSSL(X_labeled, X_unlabeled, S_labeled, S_unlabeled, y, samp_prob, n_lab,
                                            num_knots, basis_type, num_folds, my_threshold, reps, num_perts)
  expect_equal(round(stratifiedSSL_results[["dr_pert_omr"]][1],1), 0.1)
})
