# bliss <img src='man/figures/logo.png' align="right" height="139" />

Bayesian functional Linear regression with Sparse Step functions (BLiSS)

A method for the Bayesian Functional Linear Regression model (functions-on-scalar),
including two estimators of the coefficient function and an estimator of its support.
A representation of the posterior distribution is also available.

https://pmgrollemund.github.io/bliss/

# Installation

To install the **bliss** package, the easiest is to install it directly from GitHub. Open an R session and run the following commands:

```R
library(remotes) 
install_github("pmgrollemund/bliss", build_vignettes=TRUE)
```

# Usage

Once the package is installed on your computer, it can be loaded into a R session:

```R
library(bliss)
help(package="bliss")
```

# Citation

As a lot of time and effort were spent in creating the **bliss** method, please cite it when using it for data analysis:

Grollemund, Paul-Marie; Abraham, Christophe; Baragatti, Meïli; Pudlo, Pierre. Bayesian Functional Linear Regression with Sparse Step Functions. Bayesian Anal. 14 (2019), no. 1, 111--135. doi:10.1214/18-BA1095. https://projecteuclid.org/euclid.ba/1524103229

You should also cite the **bliss** package:

```R
citation("bliss")
```

See also citation() for citing R itself.
