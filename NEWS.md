# bliss 1.0.5

* change @docType keyword to _PACKAGE_ in bliss.R

# bliss 1.0.4

* Use `mvrnorm()` from MASS instead of deprecated rockchalk package

# bliss 1.0.3

* Added a `NEWS.md` file to track changes to the package.
* Fix comparison problem between *unsigned* and *int* in Rcpp script
* Fix problem with progress displaying in `fit_Bliss()`
* Fix bugs from GitHub installation
* Fix potential bug with `support_estimate()` function
* Add an argument to choose if the support estimate should be computed or not
* Watch if the design matrixes contain constant columns
* Add messages indicating potential numerical problems
* Update introduction vignette
* Allow the possibility to tune *phi_l*, the hyperparameters of the exponential prior on the length of the intervals
