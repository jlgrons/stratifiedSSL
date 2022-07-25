# stratifiedSSL

The stratifiedSSL R package performs semi-supervised learning under stratified sampling.  The package can fit a semi-supervised logistic regression model that has improved statisical efficiency over standard supervised logistic regression.  The package also performs semi-supervised model evaluation with respect to the overall misclassification rate and Brier score to yield estimates with lower variation than their supervised counterparts.  In addition to estimation procedures, cross-validation based inference derived from influence function expansions for the regression parameter are implemented as well as a cross-validated pertrubation resampling procedure for the model evaluation matrix.

In order to use the package, the following data is required:

* A labeled data set with information on both the covariates and the outcome of interest
* An unlabeled data set with information on only the covariates which is of much larger size than the labeled dataset

With the required data, the user can output the following analysis:

* 


Citation

Gronsbell J, Liu M, Tian L, and Cai T.  Efficient evaluation of prediction rules in semi-supervised settings under stratified sampling.  Accepted at the Journal of the Royal Statisical Society: Series B. (https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12502)[link]
