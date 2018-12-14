#########################################################
#                                                       #
#          Bliss method : simulate dataset              #
#                                                       #
#########################################################
################################# ----
#' choose_beta
################################# ----
#' @description Compute a coefficient function for the Function Linear Regression
#'              model.
#' @details Several shapes are available.
#' @return A numerical vector which corresponds to the coefficient function
#'         at given times points (\code{grid}).
#' @param param a list containing:
#' \describe{
#'  \item{grid}{a numerical vector, the time points.}
#'  \item{p}{a numerical value, the length of the vector \code{grid}.}
#'  \item{shape}{a character vector: "smooth", "random_smooth",
#'               "simple", "simple_bis", "random_simple", "sinusoid",
#'               "flat_sinusoid" and "sharp"}
#' }
#' @export
#' @examples
#' ### smooth
#' param <- list(p=100,grid=seq(0,1,length=100),shape="smooth")
#' beta_function <- choose_beta(param)
#' plot(param$grid,beta_function,type="l")
#' ### random_smooth
#' param <- list(p=100,grid=seq(0,1,length=100),shape="random_smooth")
#' beta_function <- choose_beta(param)
#' plot(param$grid,beta_function,type="l")
#' ### simple
#' param <- list(p=100,grid=seq(0,1,length=100),shape="simple")
#' beta_function <- choose_beta(param)
#' plot(param$grid,beta_function,type="s")
#' ### simple_bis
#' param <- list(p=100,grid=seq(0,1,length=100),shape="simple_bis")
#' beta_function <- choose_beta(param)
#' plot(param$grid,beta_function,type="s")
#' ### random_simple
#' param <- list(p=100,grid=seq(0,1,length=100),shape="random_simple")
#' beta_function <- choose_beta(param)
#' plot(param$grid,beta_function,type="s")
#' ### sinusoid
#' param <- list(p=100,grid=seq(0,1,length=100),shape="sinusoid")
#' beta_function <- choose_beta(param)
#' plot(param$grid,beta_function,type="l")
#' ### flat_sinusoid
#' param <- list(p=100,grid=seq(0,1,length=100),shape="flat_sinusoid")
#' beta_function <- choose_beta(param)
#' plot(param$grid,beta_function,type="l")
#' ### sharp
#' param <- list(p=100,grid=seq(0,1,length=100),shape="sharp")
#' beta_function <- choose_beta(param)
#' plot(param$grid,beta_function,type="l")
choose_beta <- function(param){
  # load objects
  p    <- param[["p"]]
  grid <- param[["grid"]]
  shape <- param[["shape"]]

  # Compute a grid on (0,1).
  grid_01 <- (grid - min(grid))/ (max(grid) - min(grid))
  # Choose a function beta
  if(shape == "smooth"){
    beta <- 5*exp(-((grid_01-0.25)*20)^2) +
      2*exp(-((grid_01-0.75)*20)^2) -
      2*exp(-((grid_01-0.5)*20)^2)
  }
  if(shape == "random_smooth"){
    beta <- runif(1,-5,5)*exp(-((grid_01-runif(1,0,1))*20)^2) +
      runif(1,-5,5)*exp(-((grid_01-runif(1,0,1))*20)^2) +
      runif(1,-5,5)*exp(-((grid_01-runif(1,0,1))*20)^2)
  }
  if(shape == "simple2"){
    beta <- rep(0,p)
    beta[round(p/10):round(3*p/10)] <- 3
    beta[round(5*p/10):round(6*p/10)] <- 4
    beta[round(8*p/10):round(9.5*p/10)] <- -1
  }
  if(shape == "simple"){
    beta <- rep(0,p)
    beta[round(p/10):round(3*p/10)] <- 3
    beta[round(5*p/10):round(6*p/10)] <- 4
    beta[round(8*p/10):round(9.5*p/10)] <- -3
  }
  if(shape == "simple_bis"){
    beta <- rep(0,p)
    beta[1:round(2*p/10)] <- 3
    beta[round(5*p/10):round(6*p/10)] <- 4
    beta[round(6*p/10):round(7.5*p/10)] <- -1
  }
  if(shape == "simple_K10"){
    beta <- rep(0,p)
    beta[round(0.5*p/10):round(2*p/10)]   <- 1 + beta[round(0.5*p/10):round(2*p/10)]
    beta[round(p/10):round(2*p/10)]       <- 2 + beta[round(p/10):round(2*p/10)]
    beta[round(0.8*p/10):round(1.7*p/10)] <- 1 + beta[round(0.8*p/10):round(1.7*p/10)]
    beta[round(4.5*p/10):round(7*p/10)]   <- 2 + beta[round(4.5*p/10):round(7*p/10)]
    beta[round(5*p/10):round(7*p/10)]     <- 1 + beta[round(5*p/10):round(7*p/10)]
    beta[round(5*p/10):round(6*p/10)]     <- 2 + beta[round(5*p/10):round(6*p/10)]
    beta[round(8*p/10):round(9.5*p/10)]   <- -0.5 + beta[round(8*p/10):round(9.5*p/10)]
    beta[round(8*p/10):round(10*p/10)]    <- -1 + beta[round(8*p/10):round(10*p/10)]
    beta[round(8*p/10):round(9.5*p/10)]   <- -1 + beta[round(8*p/10):round(9.5*p/10)]
    beta[round(8.8*p/10):round(9.5*p/10)] <- -0.5 + beta[round(8.8*p/10):round(9.5*p/10)]
  }
  if(shape == "random_simple"){
    beta <- rep(0,p)
    boundaries <- sort(sample(1:p,6))
    beta[boundaries[1]:boundaries[2]] <- runif(1,-5,5)
    beta[boundaries[3]:boundaries[4]] <- runif(1,-5,5)
    beta[boundaries[5]:boundaries[6]] <- runif(1,-5,5)
  }
  if(shape == "sinusoid"){
    beta <- sin(grid_01 * 2* pi)
  }
  if(shape == "flat_sinusoid"){
    beta <- rep(0,p)
    flat          <- round(p/3)
    beta[1:flat]  <- sin(10/(p-flat+10)               * 2* pi)
    beta[flat:p]  <- sin((10:(p-flat+10))/(p-flat+10) * 2* pi)
    beta          <- beta * sigmoid(1:p-flat)
  }
  if(shape == "sharp"){
    beta  <- rep(0,p)
    shift <- max(grid_01) - min(grid_01)
    beta  <- beta +
      2 * sigmoid_sharp(grid_01,min(grid_01) + 0.2 * shift,v=100,asym=1) -
      3 * sigmoid_sharp(grid_01,min(grid_01) + 0.6 * shift,v=100,asym=1)
  }

  # Return the chosen function
  return(beta)
}

################################# ----
#' sim
################################# ----
#' @description Simulate a dataset for the Function Linear Regression model.
#' @return a list containing:
#' \describe{
#'  \item{Q}{an integer, the number of functional covariates.}
#'  \item{y}{a numerical vector, the outcome observations.}
#' \item{x}{a list of matrices, the qth matrix contains the observations of the
#'       qth functional covariate at time points given by \code{grids}.}
#'  \item{grids}{a list of numerical vectors, the qth vector is the grid of
#'        time points for the qth functional covariate.}
#'  \item{betas}{a list of numerical vectors, the qth vector is the 'true' coefficient
#'        function associated to the qth covariate on a grid of time points
#'        given with \code{grids}.}
#' }
#' @param param a list containing:
#' \describe{
#'  \item{beta_shapes}{a character vector. The qth item indicates the shape of
#'        the coefficient function associated to the qth functional covariate.}
#'  \item{n}{an integer, the sample size.}
#'  \item{p}{a vector of integers, the qth component is the number of
#'        times for the qth covariate.}
#'  \item{Q}{an integer, the number of functional covariates.}
#'  \item{autocorr_diag}{a list of numerical vectors (optional), the qth vector is the
#'        diagonal of the autocorrelation matrix of the qth functional
#'        covariate.}
#'  \item{autocorr_spread}{a vector of numerical values (optional) which are related to the
#'        autocorrelation of the functional covariates.}
#'  \item{grids}{a list of numerical vectors (optional), the qth vector is the grid
#'        of time points for the qth functional covariate.}
#'  \item{grids_lim}{a list of numerical vectors  (optional), the qth item is the lower
#'        and upper boundaries of the domain for the qth functional covariate.}
#'  \item{link}{a function (optional) to simulate data from the Generalized Functional
#'        Linear Regression model.}
#'  \item{mu}{a numerical value (optional), the 'true' intercept of the model.}
#'  \item{r}{a nonnegative value (optional), the signal to noise ratio.}
#'  \item{x_shapes}{a character vector (optional). The qth item indicates the shape of the
#'        functional covariate observations.}
#' }
#' @param verbose write stuff if TRUE.
#' @export
#' @examples
#' library(RColorBrewer)
#' param <- list(Q=2,n=25,p=c(50,50),grids_lim=list(c(0,1),c(-1,2)))
#' data <- sim(param)
#' data$y
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(10)
#' q=2
#' par(mfrow=c(2,1))
#' matplot(data$grids[[q]],t(data$x[[q]]),type="l",lty=1,col=cols)
#' plot(data$grids[[q]],data$betas[[q]],type="l")
#' abline(h=0,lty=2,col="gray")
#' par(mfrow=c(1,1))
sim <- function(param,verbose=FALSE){
  if(verbose) cat("Simulation of the data.\n")
  # load objects
  Q <- param[['Q']]
  n <- param[['n']]
  p <- param[['p']]

  # load optional objects
  grids_lim <- param[['grids_lim']]
  mu     <- param[['mu']]
  r      <- param[['r']]
  link   <- param[['link']]
  grids  <- param[['grids']]
  beta_shapes <- param[['beta_shapes']]
  x_shapes    <- param[['x_shapes']]
  autocorr_spread <- param[['autocorr_spread']]
  autocorr_diag   <- param[['autocorr_diag']]

  # Initialize the required unspecified objects
  if(is.null(grids_lim)){
    grids_lim <- list()
    for(q in 1:Q) grids_lim[[q]] <- c(0,1)
  }
  if(is.null(mu))   mu   <- 1
  if(is.null(r))    r    <- 5
  if(is.null(link)) link <- function(expectation) expectation
  if(is.null(x_shapes))    x_shapes    <- rep(NULL,Q)
  if(is.null(beta_shapes)) beta_shapes <- rep("simple",Q)

  # Derive some objects
  if(is.null(grids)) {
    grids  <- list()
    for (q in 1:Q) grids[[q]] <- seq(grids_lim[[q]][1],grids_lim[[q]][2],length=p[q])
    param[['grids']] <- grids
  }
  if(!is.null(grids)) {
    check <- TRUE
    for (q in 1:Q) check <- check & length(grids[[q]])==p[q]
    if(check == FALSE) stop("The length of each grid (parameter grids) should correspond to the number of observation times (parameter p).")
  }

  # Simulate the functional covariate observed on the grids.
  if(verbose) cat("\t Simulate functional covariate observations.\n")
  x <- list()
  for (q in 1:Q){
    param_sim_x <- list(n=n,p=p[q],grid=grids[[q]],shape=x_shapes[q],
                        autocorr_spread=autocorr_spread[q],
                        autocorr_diag=autocorr_diag[[q]])
    x[[q]] <- sim_x(param_sim_x)
  }

  # Choose a coefficient function beta
  if(verbose) cat("\t Choose a coefficient function.\n")
  betas <- list()
  for (q in 1:Q){
    param_choose_beta <- list(p=p[q],grid=grids[[q]],shape=beta_shapes[q])
    betas[[q]] <- choose_beta(param_choose_beta)
  }

  if(verbose) cat("\t Compute the outcome values.\n")
  # Compute the expectation of the outcome
  y_expe <- rep(mu,n)
  for(i in 1:n){
    for(q in 1:Q){
      x_beta    <- x[[q]][i,] * betas[[q]]
      y_expe[i] <- y_expe[i] + integrate_trapeze(grids[[q]],x_beta)
    }
  }

  # Compute the error
  err <- rnorm(n,0,1)
  err <- sd(y_expe) * err / (sd(err) * sqrt(r))

  # Compute the outcome values
  y <- link(y_expe) + err

  # Return the data.
  return(list(Q     = Q,
              y     = y,
              x     = x,
              betas = betas,
              grids = grids))
}

################################# ----
#' sim_x
################################# ----
#' @description Simulate functional covariate observations.
#' @details Several shape are available for the observations: "Fourier",
#'          "Fourier2", "random_walk", "random_sharp", "uniform", "gaussian",
#'          "mvgauss", "mvgauss_different_scale", "mvgauss_different_scale2",
#'          "mvgauss_different_scale3" and "mvgauss_different_scale4".
#' @return a matrix which contains the functional covariate observations at time
#'         points given by \code{grid}.
#' @param param a list containing :
#' \describe{
#'  \item{grid}{a numerical vector, the observation times.}
#'  \item{n}{an integer, the sample size.}
#'  \item{p}{an integer, the number of observation times.}
#'  \item{diagVar}{a numerical vector (optional), the diagonal of the autocorrelation matrix.}
#'  \item{dim}{a numerical value (optional), the dimension of the Fourier basis,
#'             if "shape" is "Fourier" or "Fourier2". }
#'  \item{ksi}{a numerical value (optional) related to the observations correlation.}
#'  \item{x_shape}{a character vector (optional), the shape of the observations. }
#' }
#' @importFrom rockchalk mvrnorm
#' @importFrom stats pexp runif
#' @export
#' @examples
#' library(RColorBrewer)
#' ### Fourier
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_shape="Fourier")
#' x <- sim_x(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### Fourier2
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="Fourier2")
#' x <- sim_x(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### random_walk
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="random_walk")
#' x <- sim_x(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### random_sharp
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="random_sharp")
#' x <- sim_x(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### uniform
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="uniform")
#' x <- sim_x(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### gaussian
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="gaussian")
#' x <- sim_x(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### mvgauss
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="mvgauss")
#' x <- sim_x(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
sim_x <- function(param){
  # load objects
  n <- param$n
  p <- param$p
  grid <- param$grid

  # load optional objects
  shape    <- param$shape
  dim     <- param$dim
  ksi     <- param$ksi
  diagVar <- param$diagVar

  # Initialize the necessary unspecified objects
  if(is.null(shape))   shape   <- "mvgauss"
  if(is.null(dim))     dim     <- 4
  if(is.null(ksi))     ksi     <- 1
  if(is.null(diagVar)) diagVar <- abs(rnorm(p,1,1/10))

  # Deduce objects
  by <- diff(grid)[1]

  # Simulate the functions x_i(t)
  if(shape == "Fourier"){
    # Set a Fourier basis
    Fourier_basis <- build_Fourier_basis(grid = grid,
                                         dim  = dim,
                                         per  = 1.5*(max(grid)-min(grid)))
    # Choose the coefficients
    a_n <- rockchalk::mvrnorm(n,(dim:1)/dim, diag((dim:1)/(50*dim)))
    b_n <- rockchalk::mvrnorm(n,(dim:1)/dim, diag((dim:1)/(50*dim)))

    # Compute the functions x_i(t)
    x <- a_n %*% Fourier_basis[1:dim,] + b_n %*% Fourier_basis[(dim+1):(2*dim),]
  }
  if(shape == "Fourier2"){
    # Set a Fourier basis
    Fourier_basis <- build_Fourier_basis(grid = grid,
                                         dim  = dim,
                                         per  = 1.5*(max(grid)-min(grid)))
    # Choose the coefficients
    a_n <- runif(n*dim,-3,3)
    dim(a_n) <- c(n,dim)
    b_n <- runif(n*dim,-3,3)
    dim(b_n) <- c(n,dim)

    # Determiner les courbes
    x <- a_n %*% Fourier_basis[1:dim,] + b_n %*% Fourier_basis[(dim+1):(2*dim),]
  }
  if(shape == "random_walk"){
    start <- rnorm(n,0,2)
    x <- compute_random_walk(n,p,0,1,start)
  }
  if(shape == "random_sharp"){
    locs <- runif(n*2,grid[1],tail(grid,1))
    dim(locs) <- c(n,2)

    asyms <- runif(n*2,1,5)
    dim(asyms) <- c(n,2)

    vs <- runif(n*2, 1/(4*by), 1/(3*by) )
    dim(vs) <- c(n,2)

    s <- sample(c(-1,1),2*n,replace=T)
    dim(s) <- c(n,2)

    x <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- s[i,1] * sigmoid_sharp(grid,locs[i,1],asyms[i,1],vs[i,1]) +
        s[i,2] * sigmoid_sharp(grid,locs[i,2],asyms[i,2],vs[i,2])
    }
  }
  if(shape == "uniform"){
    x <- matrix(0,n,p)
    for(j in 1:p){
      x[,j] <- runif(n,-5,5)
    }
  }
  if(shape == "gaussian"){
    x <- matrix(0,n,p)
    for(j in 1:p){
      x[,j] <- rnorm(n,0,4)
    }
  }
  if(shape == "mvgauss"){
    mu      <- (1:p-p/2)^2/(p^2/4)
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- rockchalk::mvrnorm(1,mu,Sigma)
    }
  }
  if(shape == "mvgauss_different_scale"){
    mu      <- (1:p-p/2)^2/(p^2/4)
    diagVar[1:floor(p/3)] <- 10 * diagVar[1:floor(p/3)]
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- rockchalk::mvrnorm(1,mu,Sigma)
    }
  }
  if(shape == "mvgauss_different_scale2"){
    mu      <- (1:p-p/2)^2/(p^2/4)
    diagVar[1:floor(p/3)] <- 100 * diagVar[1:floor(p/3)]
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- rockchalk::mvrnorm(1,mu,Sigma)
    }
  }
  if(shape == "mvgauss_different_scale3"){
    mu      <- (1:p-p/2)^2/(p^2/4)
    diagVar[1:floor(p/3)] <- 1000 * diagVar[1:floor(p/3)]
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- rockchalk::mvrnorm(1,mu,Sigma)
    }
  }
  if(shape == "mvgauss_different_scale4"){
    mu      <- (1:p-p/2)^2/(p^2/4)
    diagVar[floor(2*p/3):p] <- 100 * diagVar[floor(2*p/3):p]
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- rockchalk::mvrnorm(1,mu,Sigma)
    }

  }

  # Return the functions
  return(x)
}

################################# ----
#' build_Fourier_basis
################################# ----
#' @description Define a Fourier basis to simulate functional covariate observations.
#' @return a matrix. Each row is an functional observation evaluated on the
#'         \code{grid} time points.
#' @param grid a numerical vector.
#' @param dim a numerical value. It corresponds to \code{dim(basis)/2}.
#' @param per a numerical value which corresponds to the period of the sine and
#'        cosine functions.
#' @details See the \code{\link[=sim_x]{sim_x}} function.
#' @export
#' @examples
#' # See the function \code{sim_x}.
build_Fourier_basis <- function(grid,dim,per=2*pi){
  sapply(grid,function(x) c(cos(2*pi*x*(1:dim)/per),sin(2*pi*x*(1:dim)/per) )  )
}

################################# ----
#' compute_random_walk
################################# ----
#' @description Compute a (Gaussian) random walk.
#' @return a matrix where each row is a random walk.
#' @param n an integer, the number of random walks.
#' @param p an integer, the length of the random walks.
#' @param mu a numerical vector, the mean of the random walks.
#' @param sigma a numerical value which is the standard deviation of the
#'                 gaussian distribution used to compute the random walks.
#' @param start a numerical vector (optional) which is the initial value of
#'                 the random walks.
#' @details See the \code{\link[=sim_x]{sim_x}} function.
#' @importFrom stats rnorm
#' @export
#' @examples
#' # see the sim_x() function.
compute_random_walk <- function(n,p,mu,sigma,start=rep(0,n)){
  res <- matrix(0,n,p)
  for(i in 1:n){
    add     <- rnorm(p,mu,sigma)
    res[i,] <- cumsum(add)
    res[i,] <- start[i] + res[i,]
  }
  res <- start + res
  return(res)
}

################################# ----
#' sigmoid
################################# ----
#' @description Compute a sigmoid function.
#' @return a numerical vector.
#' @param x a numerical vector, time points.
#' @param asym a numerical value (optional), the asymptote of the sigmoid function.
#' @param v a numerical value (optional), related to the slope at the origin.
#' @details see the function \code{\link[=sim_x]{sim_x}}.
#' @export
#' @examples
#' ## Test 1 :
#' x <- seq(-7,7,0.1)
#' y <- sigmoid(x)
#' plot(x,y,type="l",main="Sigmoid function")
#' ## Test 2 :
#' x  <- seq(-7,7,0.1)
#' y  <- sigmoid(x)
#' y2 <- sigmoid(x,asym=0.5)
#' y3 <- sigmoid(x,v   =  5)
#' plot(x,y,type="l",main="Other sigmoid functions")
#' lines(x,y2,col=2)
#' lines(x,y3,col=3)
sigmoid <- function(x,asym=1,v=1){
  (asym^-1 + exp(-v*x))^-1
}

################################# ----
#' sigmoid_sharp
################################# ----
#' @description Compute a sharp sigmoid function.
#' @return a numerical vector.
#' @param x a numerical vector, time points.
#' @param loc a numerical value (optional), the time of the sharp.
#' @param ... Arguments (optional) for the function sigmoid.
#' @details see the function \code{\link[=sim_x]{sim_x}}.
#' @export
#' @examples
#' ## Test 1 :
#' x <- seq(-7,7,0.1)
#' y <- sigmoid_sharp(x)
#' plot(x,y,type="l",main="Sharp sigmoid")
#' ## Test 2 :
#' x  <- seq(-7,7,0.1)
#' y  <- sigmoid_sharp(x,loc=3)
#' y2 <- sigmoid_sharp(x,loc=3,asym=0.5)
#' y3 <- sigmoid_sharp(x,loc=3,v   =  5)
#' plot(x,y,type="l",main="Other sharp sigmoids")
#' lines(x,y2,col=2)
#' lines(x,y3,col=3)
sigmoid_sharp <- function(x,loc=0,...){
  # 4 should be replace by (a+1)^2 such that the maximum of the curve
  # provided by sigmoid_sharp is 1.
  4*(sigmoid(x-loc,...) * sigmoid(-x+loc,...))
}


################################# ----
#' corr_matrix
################################# ----
#' @description Compute an autocorrelation matrix.
#' @return a symmetric matrix.
#' @param diagonal a numerical vector corresponding to the diagonal.
#' @param ksi a numerical value, related to the correlation.
#' @export
#' @examples
#' ### Test 1 : weak autocorrelation
#' ksi     <- 1
#' diagVar <- abs(rnorm(100,50,5))
#' Sigma   <- corr_matrix(diagVar,ksi^2)
#' persp(Sigma)
#' ### Test 2 : strong autocorrelation
#' ksi     <- 0.2
#' diagVar <- abs(rnorm(100,50,5))
#' Sigma   <- corr_matrix(diagVar,ksi^2)
#' persp(Sigma)
corr_matrix <- function(diagonal,ksi){
  # Initialize
  p <- length(diagonal)
  res <- diag(diagonal)

  # Compute the correlation matrix
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      res[i,j] <- exp(-ksi*(i-j)^2/p)*sqrt(res[i,i]*res[j,j])
      res[j,i] <- exp(-ksi*(i-j)^2/p)*sqrt(res[i,i]*res[j,j])
    }
  }

  # return the matrix
  return(res)
}
