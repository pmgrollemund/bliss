#' bliss: Bayesian functional Linear regression with Sparse Step functions
#'
#' A method for the Bayesian Functional Linear Regression model (functions-on-scalar),
#' including two estimators of the coefficient function and an estimator of its support.
#' A representation of the posterior distribution is also available.
#'
#' @name bliss
#' @import Rcpp RcppArmadillo
#' @importFrom Rcpp evalCpp cppFunction Rcpp.plugin.maker
#' @useDynLib bliss, .registration = TRUE
# @exportPattern "^[[:alpha:]]+"
"_PACKAGE"
NULL
