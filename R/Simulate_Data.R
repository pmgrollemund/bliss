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
#' @return A numerical vector which is the observation of the coefficient function
#'         at the given grid of times points (\code{grid}).
#' @param param a list containing: 
#' \describe{
#'  \item{grid}{a numerical vector, the observation time points.}
#'  \item{p}{a numerical value, the length of the vector \code{grid}.}
#'  \item{shape}{a character vector to choose between "smooth", "random_smooth",
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
 if(shape == "simple"){
  beta <- rep(0,p)
  beta[round(p/10):round(3*p/10)] <- 3
  beta[round(5*p/10):round(6*p/10)] <- 4
  beta[round(8*p/10):round(9.5*p/10)] <- -1
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
#'  \item{x}{a list of matrices, the qth matrix contains the observation of the 
#'        qth functional covariate on a grid of time points given with \code{grids}.}
#'  \item{grids}{a list of numerical vectors, the qth vector is the grid of 
#'        time points for the qth functional covariate.}
#'  \item{betas}{a list of numerical vectors, the qth vector is the coefficient 
#'        function associated to the qth covariate on a grid of time points 
#'        given with \code{grids}.}
#' }
#' @param param a list containing:
#' \describe{
#'  \item{n}{an integer, the sample size.}
#'  \item{p}{a vector of integers, the qth component is the number of 
#'        times for the qth covariate.}
#'  \item{grids}{a list of numerical vectors, the qth vector is the grid 
#'        of time points for the qth functional covariate (optional, can be 
#'        deduced from grid_lim and p if it is \code{NULL}).}
#'  \item{grids_lim}{a list of numerical vectors, the qth item is the lower
#'        and upper boundaries of the domain for the qth functional covariate. (optional)}        
#'  \item{mu}{a numerical value, the intercept of the model. (optional)}
#'  \item{r}{a nonnegative value, the signal to noise ratio. (optional)}
#'  \item{link}{a function. XXXX (optional)}
#'  \item{beta_shapes}{a character vector. The qth item indicates the shape of 
#'        the coefficient function associated to the qth functional covariate.}
#'  \item{x_shapes}{a character vector. The qth item indicates the shape of the 
#'        functional covariate observations. (optional)}
#'  \item{autocorr_diag}{a list of numerical vectors, the qth vector is the 
#'        diagonal of the autocorrelation matrix of the qth functional 
#'        covariate (optional).}
#'  \item{autocorr_spread}{a vector of numerical values which are related to the 
#'        autocorrelation of the functional covariates. (optional)}
#'  \item{correlation}{a matrix which gives the correlation structure between the 
#'        functional covariates (optional).}
#' }
#' @param progress a logical value. If TRUE, the algorithm progress is displayed.
#'         (optional)
#' @export
#' @examples
#' library(RColorBrewer)
#' param <- list(n=10,p=c(200,100),shapes=c("simple","smooth"),grid_min=c(0,-2.12),grid_max=c(1,3.14))
#' data <- sim_multiple(param)
#' data$y
#' data$expe
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(10)
#' q=2
#' par(mfrow=c(2,1))
#' matplot(data$grids[[q]],t(data$x_mult[[q]]),type="l",lty=1,col=cols)
#' plot(data$grids[[q]],data$beta_function_mult[[q]],type="l")
#' abline(h=0,lty=2,col="gray")
#' par(mfrow=c(1,1))

sim <- function(param,progress=FALSE){
 if(progress) cat("Simulation of the data.\n")
 # load objects
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
 correlation     <- param[['correlation']]
 
 # Initialize the required unspecified objects
 Q <- length(p)
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
 if(progress) cat("\t Simulate functional covariate observations.\n")
 if( (Q == 1) && !(is.null(correlation)) )
  stop("'Correlation' is a correlation structure between different functional covariates.")
 if( (Q > 1) && !(is.null(correlation)) ){
  x <- sim_correlated_x(param)
 }
 if( is.null(correlation)){
  x <- list()
  for (q in 1:Q){
   param.sim_x <- list(n=n,p=p[q],grid=grids[[q]],shape=x_shapes[q],
                       autocorr_spread=autocorr_spread[q],
                       autocorr_diag=autocorr_diag[[q]])
   x[[q]] <- sim_x(param.sim_x)
  }
 }
 
 # Choose a coefficient function beta
 if(progress) cat("\t Choose a coefficient function.\n")
 betas <- list()
 for (q in 1:Q){
  param.choose_beta <- list(p=p[q],grid=grids[[q]],shape=beta_shapes[q])
  betas[[q]] <- choose_beta(param.choose_beta)
 }
 
 if(progress) cat("\t Compute the outcome values.\n")
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
#' @description simulate the functions x_qi(t) for a given q.
#' @details Several shape are available for the functions x_qi(t): "Fourier", "Fourier2", "random_walk", "random_sharp",
#'          "uniform", "gaussian", "mvgauss", "mvgauss_different_scale", "mvgauss_different_scale2", "mvgauss_different_scale3",
#'          "mvgauss_different_scale4"
#' @return a matrix. Each row is a function x_qi(t) (i=1,...,n and q is fixed).
#' @param param a list containing :
#' \describe{
#'  \item{n}{an integer, the number of functions.}
#'  \item{p}{an integer, the number of observation times.}
#'  \item{grid}{a numerical vector, the grid of observation times.}
#'  \item{x_shape}{a character vector, the shape of the functions x_i(t). (optional)}
#'  \item{dim}{a numerical value, the dimension of the Fourier basis, if "shape" is "Fourier" or "Fourier2". (optional)}
#'  \item{ksi}{a numerical value, a "coefficient of correlation", see the Bliss article Section 3.1 for more details.}
#'  \item{diagVar}{a numerical vector, the diagonal of the autocorrelation matrix of the functions x_i(t).}
#' }
#' @importFrom rockchalk mvrnorm
#' @export
#' @examples
#' library(RColorBrewer)
#' ### Fourier
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_shape="Fourier")
#' x <- sim_functions(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' \donttest{
#' ### Fourier2
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="Fourier2")
#' x <- sim_functions(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### random_walk
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="random_walk")
#' x <- sim_functions(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### random_sharp
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="random_sharp")
#' x <- sim_functions(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### uniform
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="uniform")
#' x <- sim_functions(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### gaussian
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="gaussian")
#' x <- sim_functions(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' ### mvgauss
#' param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="mvgauss")
#' x <- sim_functions(param)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
#' matplot(param$grid,t(x),type="l",lty=1,col=cols)
#' }
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
  x <- random_walk(n,p,0,1,start)
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
#' sim_correlated_x
################################# ----
#' @description todo
#' @return todo
#' @param param a list containing :
#' \describe{
#'  \item{n}{an integer, the number of functions.}
#'  \item{p}{an integer, the number of observation times.}
#'  \item{grid}{a numerical vector, the grid of observation times.}
#'  \item{x_shape}{a character vector, the shape of the functions x_i(t). (optional)}
#'  \item{dim}{a numerical value, the dimension of the Fourier basis, if "shape" is "Fourier" or "Fourier2". (optional)}
#'  \item{ksi}{a numerical value, a "coefficient of correlation", see the Bliss article Section 3.1 for more details.}
#'  \item{diagVar}{a numerical vector, the diagonal of the autocorrelation matrix of the functions x_i(t).}
#' }
#' @export
#' @examples
#' #todo
sim_correlated_functions <- function(param){
 # load objects
 n <- param$n
 p <- param$p
 grids <- param$grids
 
 # load optional objects
 shape    <- param$shape
 dim     <- param$dim
 ksi     <- param$ksi
 diagVar <- param$diagVar
 correlation <- param$correlation # IS 16/02/2018
 
 # Initialize the necessary unspecified objects
 Q <- length(p) # IS 06/02/2018
 if(is.null(correlation))  correlation <- diag(Q)  # IS 16/02/2018
 if(is.null(shape))    shape    <- "mvgauss"
 if(is.null(ksi))     ksi     <- c(1,1)
 diagVar <- list()
 for(q in 1:Q){
  diagVar[[q]] <- abs(rnorm(p[q],1,1/10))
 }
 
 # Simulate the functions x_i(t)
 if(shape == "mvgauss"){
  mu <- NULL
  for(q in 1:Q){
   mu <- c(mu,q*(1:p[q]-p[q]/2)^2/(p[q]^2/4))
   # mu <- c(mu,rnorm(1,0,0.001)+(1:p[q]/p[q])^rnorm(1,3,1)/(rnorm(1,1,1))^2/4)
  }
  Sigma <- matrix(0,sum(p),sum(p))
  
  count <- 0
  for(q in 1:Q){
   count2 <- 0
   for(qq in 1:Q){
    if(q == qq){
     tmp <- diag(diagVar[[q]])
    }else{
     tmp <- matrix(0,p[q],p[qq])
    }
    
    for(i in 1:p[q]){
     for(j in 1:p[qq]){
      if( q==qq ){
       tmp[i,j] <- exp(-ksi[q]*abs(grids[[q]][i]-grids[[q]][j]))*
        sqrt(diagVar[[q]][i]*diagVar[[qq]][j])
      }else{
       # IS 06/02/2018: d'ou vient correlation???
       tmp[i,j] <- correlation[q,qq] *
        exp(-sqrt(ksi[q]*ksi[qq])*abs(grids[[q]][i]-grids[[qq]][j]))*
        sqrt(diagVar[[q]][i]*diagVar[[qq]][j])
      }
     }
    }
    Sigma[ (1+count):(p[q]+count) , (1+count2):(p[qq]+count2)] <- tmp
    count2 <- count2 + p[qq]
   }
   count <- count + p[q]
  }
  
  x       <- matrix(0,n,sum(p))
  for(i in 1:n){
   x[i,] <- rockchalk::mvrnorm(1,mu,Sigma)
  }
 }
 
 x_mult <- list()
 count <- 0
 for( q in 1:Q){
  x_mult[[q]] <- x[,(1+count):(p[q]+count)]
  count <- count + p[q]
 }
 
 # Return the functions
 return(x_mult)
}


################################# ----
#' build_Fourier_basis
################################# ----
#' @description Define a Fourier basis to simulate functional covariate observation.
#' @return a matrix. Each row is an functional observation evaluated on the 
#'         \code{grid} time points.
#' @param grid a numerical vector.
#' @param dim a numerical value. It corresponds to dim(basis)/2.
#' @param per a numerical value which corresponds to the period of the sine and 
#'        cosine functions.
#' @details See the \code{\link[=sim_x]{sim_x}} function.
#' @export
#' @examples
#' # See the sim_x() function.
build_Fourier_basis <- function(grid,dim,per=2*pi){
 sapply(grid,function(x) c(cos(2*pi*x*(1:dim)/per),sin(2*pi*x*(1:dim)/per) )  )
}
################################# ----
#' random_walk
################################# ----
#' @description compute a random walk. (gaussian)
#' @return a matrix where each row is a random walk.
#' @param n an integer, the number of random walk.
#' @param p an integer, the length of the random walks.
#' @param mu a numerical vector, the average random walk.
#' @param sigma a numerical value which is the standard deviation of the
#'                 gaussian distribution used to compute the random walk.
#' @param start a numerical vector which is the initial value of
#'                 the random walks. (optional)
#' @details see the \code{\link[=sim_functions]{sim_functions}} function.
#' @importFrom stats rnorm
#' @export
#' @examples
#' # see the sim_functions() function.
random_walk <- function(n,p,mu,sigma,start=rep(0,n)){
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
#' @description compute a sigmoid function.
#' @details used to simulate a coefficient function or functions x_i(t).
#' @return a numerical vector.
#' @param x a numerical vector, a grid of points.
#' @param asym the value of the asymptote of the sigmoid function. (optional)
#' @param v a numerical value which is related to the slope at the origin. (optional)
#' @details see the function \code{\link[=sim_functions]{sim_functions}}.
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
#' @description compute a sharp function from the sigmoid function
#' @details used to simulate a coefficient function or functions x_i(t).
#' @return a numerical vector.
#' @param x a numerical vector, a grid of points.
#' @param loc a numerical value, the instant of the sharp. (optional)
#' @param ... Arguments to be passed to the function sigmoid. (optional)
#' @details see the function \code{\link[=sim_functions]{sim_functions}}.
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
#' @description compute an autocorrelation matrix to simulate functions x_i(t).
#' @return a symmetric matrix.
#' @param diagonal a numerical vector corresponding to the diagonal of the final matrix.
#' @param ksi a "coefficient of correlation". See the article Bliss, Section 3.1 for more details.
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