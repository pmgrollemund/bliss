################################# ----
#' Bliss
################################# ----
#' @description perform the Bliss method to fit the Bayesian Functional
#' Linear Regression model (with Q functional covariates).
#' @return return a list containing :
#' \describe{
#'  \item{alpha}
#'  \item{beta_sample}{a list of matrices. For the qth covariate, beta_functions[[q]] is a matrix where each row is a
#'                 function beta_qi associated to the iteration i of the Gibbs sampler.}
#'  \item{bliss_estimate}{a list of numerical vectors corresponding to the
#'  Bliss estimates of each coefficient functions.}
#'  \item{chain}
#'  \item{param_density}
#'  \item{posterior_density_estimate}{a list of Q levels. Each level contains:
#'  Firts level: a matrix containing a posterior sample of coefficient function.
#'  Second level: for each covariate, there is a list of three levels containing
#'  1) a grid x, 2) a grid y and 3) matrix z containing approximations of the
#'  posterior density of the coefficient function on the 2-dimensional grid
#'  defined from x and y.}
#'
#'  \item{res.Simulated_Annealing}{a list of lists: res.Simulated_Annealing[[q]] is the result of the
#'                 function Bliss_Simulated_Annealing applied for the qth covariate.}
#'  \item{res.Gibbs_Sampler}{a list of lists: res.Gibbs_Sampler[[q]] is the the result of the function
#'                 Bliss_Gibbs_Sampler for the qth covariate.}
#' }
#'
#' @param data a list containing
#' \describe{
#' \item{Q}{an integer, the number of covariates,}
#' \item{x}{a list of matrices, the functions x_qi(t) observed at time points given by grids,}
#' \item{y}{a numerical vector, the outcome values y_i.}
#' }
#'
#' @param param a list
#' \describe{
#' \item{K}{a vector of integers, hyperparameters of the Bliss model corresponding
#' to the number of intervals for the Q covariates.}
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of
#' observation points for the qth covariate.}
#' \item{iter}{an integer, the number of iterations of the Gibbs sampler algorithm.}
#' \item{basis}{a vector of characters among : "uniform" (default),
#'                 "epanechnikov", "gauss" and "triangular" which correspond to
#'                 different basis expansion of coefficient function and functional data (optional)}
#' \item{g}{hyperparameter of the Bliss model,  a nonnegative value,
#'                the coefficient of the Zellner prior.}
#' \item{burnin}{an integer, the number of iteration to drop of the Gibbs sampler. (optional)}
#' \item{thin}{an integer, used to thin the Gibbs sample to compute an
#'                 estimate of the posterior density of beta(t). (optional)}
#' \item{iter_sann}{an integer, the number of iteration of the Simulated Annealing. (optional)}
#' \item{n_chains}{number of chains to do in the Gibbs sampler.}
#'
#'
#' \item{theta_posterior_density}{}
#' \item{Temp}{a vector of nonnegative values, the Q initial temperatures for the
#'                 cooling function of the Q Simulated Annealings. (optional)}
#' }
#' @param plot a logical value. If TRUE, the Bliss estimates and
#'                 the estimates of the posterior density are plotted. (optional)
#' @param density a logical value. If TRUE, the posterior density is estimated. (optional)
#' @param sann a logical value. If TRUE, the bliss estimate is computed. (optional)
#' @param display a logical value. If TRUE, the algorithm progress is displayed.  (optional)
#' @importFrom utils tail
#' @export
#' @examples
#' # see the vignette BlissIntro.
Bliss   <- function(data,param,plot=FALSE,density=FALSE,sann=FALSE,display=FALSE){
 # define Q
 Q <- length(data$x)
 # centering the data and define p
 p  <- numeric()
 for (q in 1:Q){
  data$x[[q]]   <- apply(data$x[[q]],2,function(v) v - mean(v))
  p[q] <- length(param$grids[[q]])
 }
 param$p <- p
 
 # How many chains i have to do ?
 n_chains <- param[["n_chains"]]
 if(is.null(n_chains)){
  n_chains <- 1
  param[["n_chains"]] <- n_chains
 }
 
 # Initialize the list "chains"
 chains <- list()
 # For each chain :
 for(j in 1:n_chains){
  if(display) cat("Chain ",j,": \n",sep="")
  chains[[j]] <- list()
  
  # Execute the Gibbs Sampler algorithm to sample the posterior distribution
  chains[[j]]$posterior_sample <- Bliss_Gibbs_Sampler(data,param)
  
  # Compute a posterior sample of coefficient function
  chains[[j]]$beta_sample <- compute_beta(chains[[j]]$posterior_sample,param)
  
  # Compute the posterior density of the coefficient function
  chains[[j]]$beta_posterior_density <- list()
  param_density <- list()
  if (density){
   for(q in 1:Q){
    diff_grid <- diff(param$grids[[q]])[1]
    param$grids2[[q]] <- c(param$grids[[q]]-diff_grid/2,
                           tail(param$grids[[q]],1)+diff_grid/2)
    param$xlim[[q]] <- range(param$grids2[[q]])
    
    param_density <- list(grid= param$grids[[q]],
                          iter= param$iter,
                          p   = param[["p"]][q],
                          n        = param[["n"]],
                          thin     = param$thin,
                          burnin   = param[["burnin"]],
                          lims.kde = param$lims.kde[[q]],
                          h1       = param$h1,
                          new_grid = param[["new_grid"]],
                          xlim     = range(param$grids[[q]]) + c(-diff_grid,
                                                                 diff_grid),
                          display = display
    )
    chains[[j]]$beta_posterior_density[[q]] <-
     compute_beta_posterior_density(chains[[j]]$beta_sample[[q]],
                                    param_density)
   }
  }
 }
 
 # Choose a chain for inference
 j <- sample(n_chains,1)
 
 posterior_sample       <- chains[[j]]$posterior_sample
 beta_sample            <- chains[[j]]$beta_sample
 beta_posterior_density <- chains[[j]]$beta_posterior_density
 
 # Execute the Simulated Annealing algorithm to estimate the coefficient function
 Bliss_estimation <- list()
 for(q in 1:Q){
  param_sann <- list( grid = param$grids[[q]],
                      iter = param[["iter"]],
                      p    = param$p[q],
                      Temp = param$Temp[q],
                      k_max = param$k_max[q],
                      iter_sann = param[["iter_sann"]],
                      burnin    = param[["burnin"]],
                      l_max     = param[["l_max_sann"]][q],
                      basis     = param[["basis"]][q])
  
  Bliss_estimation[[q]] <- Bliss_Simulated_Annealing(beta_sample[[q]],param_sann,
                                                     posterior_sample$param$scale_ml[[q]])
 }
 
 # Compute the support estimate
 if(display) cat("Support estimation.\n")
 support_estimate <- list()
 alpha <- list()
 for(q in 1:Q){
  res_support <- support_estimation(beta_sample[[q]])
  support_estimate[[q]] <- res_support$estimate
  alpha[[q]]            <- res_support$alpha
 }
 
 # Plot the Bliss estimate and the posterior density estimate
 if(density & plot){
  for(q in 1:Q){
   image_Bliss(beta_posterior_density[[q]],param)
   lines(param$grid,Bliss_estimation[[q]]$Bliss_estimate,type="s",lwd=2)
   lines(param$grid,Bliss_estimation[[q]]$posterior_expe,type="l",lty=2)
  }
 }
 
 # Do not return the list "chains" if n_chains is 1.
 if(n_chains == 1) chains <- NULL
 
 # compute the posterior density of theta=(mu,beta,sigma)
 if(length(param$K)==1)
  if(param$theta_posterior_density){
   if(display) cat("Compute the (log) densities. \n")
   posterior_sample$posterior_density <- dposterior(posterior_sample,data)
  }
 
 # The object to return
 Bliss_estimate <- list()
 for(q in 1:Q){
  Bliss_estimate[[q]] <- Bliss_estimation[[q]]$Bliss_estimate
 }
 
 res <- list(Bliss_estimate             = Bliss_estimate,
             support_estimate           = support_estimate,
             alpha                    = alpha,
             beta_posterior_density = beta_posterior_density,
             beta_sample             = beta_sample,
             Bliss_estimation    = Bliss_estimation,
             posterior_sample          = posterior_sample,
             chains                     = chains,
             param                      = param,
             param_density              = param_density
 )
 class(res) = c("bliss")
 return(invisible(res))
}

#' Print a bliss Object
#'
#' @param x input bliss Object
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom utils str
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1,density=TRUE)
#' printbliss(res_Bliss_mult)
printbliss<-function(x,...){
 if(!any(class(x) == "bliss"))
  stop("Input must have class \"bliss\".")
 
 cat("This is a bliss x\n")
 print(str(x))
}
