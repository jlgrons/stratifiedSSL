# stratifiedSSL

The stratifiedSSL R package performs semi-supervised learning under stratified sampling.  The package can fit a semi-supervised logistic regression model that has improved statistical efficiency over standard supervised logistic regression using an imputation-based approach leveraging weighted basis function regression.  The package also performs semi-supervised model evaluation with respect to the overall misclassification rate and Brier score to yield consistent estimates for the population parameters with lower variation than their supervised counterparts.  In addition to estimation procedures, cross-validation based inference derived from the influence function expansion for the regression parameter are available as well as a cross-validated perturbation resampling procedure for the model evaluation parameters.

In order to use the package, the following data is required:

* A labeled data set with information on both the covariates and the outcome of interest
* An unlabeled data set with information on only the covariates, which is of much larger size than the labeled dataset

With the required data, the user can output the following analysis:

* Regression parameter estimates from the proposed semi-supervised method, standard supervised learning, and density-ratio based estimation
* Overall misclassification rate and Brier score estimates from the proposed semi-supervised method, standard supervised learning, and density-ratio based estimation
* Standard error estimates for the regression parameter and model evaluation parameters

## Installation

```{R, eval = FALSE}
devtools::install_github(repo = "https://github.com/jlgrons/stratifiedSSL")
```

## Citation

Gronsbell J, Liu M, Tian L, and Cai T.  Efficient evaluation of prediction rules in semi-supervised settings under stratified sampling.  Accepted at the Journal of the Royal Statisical Society: Series B. [link](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12502)
