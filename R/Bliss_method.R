################################# ----
#' Bliss method
################################# ----
#' @description Perform the Bliss method to fit the Bayesian Functional
#' Linear Regression model (with Q functional covariates).
#' @return return a list containing :
#' \describe{
#'  \item{alpha}{a numerical vector. XXXXXXX}
#'  \item{beta_sample}{a list of matrices. For the qth covariate, beta_sample[[q]] 
#'        is a matrix where each row is a coefficient function associated  to 
#'        the iteration i of the Gibbs sampler.}
#'  \item{bliss_estimate}{a list of numerical vectors corresponding to the
#'  Bliss estimates of each coefficient function.}
#'  \item{chain}{a list. XXXXX}
#'  \item{param.beta_density}{a list. XXXX} 
#'  \item{beta_posterior_density}{a list of Q items. Each item contains: XXXX des trucs a changer ? XXXX
#'  Firts level: a matrix containing a posterior sample of coefficient function.
#'  Second level: for each covariate, there is a list of three levels containing
#'  1) a grid x, 2) a grid y and 3) matrix z containing approximations of the
#'  posterior density of the coefficient function on the 2-dimensional grid
#'  defined from x and y.}#'
#'  \item{res.Simulated_Annealing}{a list: the qth item is the result of the 
#'        function Bliss_Simulated_Annealing applied for the qth functional 
#'        covariate.}
#'  \item{res.Gibbs_Sampler}{a list: the qth item is the result of the function 
#'        Bliss_Gibbs_Sampler for the qth covariate.}
#'  \item{support_estimate}{XXXXX}
#' }
#' @param data a list containing:
#' \describe{
#' \item{Q}{an integer, the number of covariates,}
#' \item{x}{a list of matrices, the qth matrix contains the observation of the 
#'       qth functional covariate at time points given by grids,}
#' \item{y}{a numerical vector, the outcome values y_i.}
#' }
#' @param param a list containing: XXXXXX
#' \describe{ 
#' \item{iter}{an integer, the number of iterations of the Gibbs sampler algorithm.}
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for 
#'       each covariate.}
#' \item{basis}{a vector of characters among : "uniform" (default),
#'       "epanechnikov", "gauss" and "triangular" which correspond to
#'       different basis functions to expand the coefficient function and the 
#'       functional covariates (optional)}
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of
#'       observation points of the qth covariate.}
#' \item{g}{a nonnegative value, hyperparameter of the Bliss model which is the
#'          coefficient of the Ridge Zellner prior. (optional)}
#' \item{burnin}{an integer, the number of iteration to drop from the Gibbs 
#'       sample. (optional)}
#' \item{thin}{an integer, used to thin the Gibbs sample to compute an
#'       approximation of the posterior density of beta(t). (optional)}
#' \item{iter_sann}{an integer, the number of iteration of the Simulated 
#'       Annealing algorithm. (optional)}
#' \item{n_chains}{number of chains to do in the Gibbs sampler.}
#' \item{theta_posterior_density}{a logical value. XXXXX}
#' \item{Temp_init}{a nonnegative value, the initial temperatures for the
#'                 cooling function of the Q Simulated Annealings. (optional)}
#' }
#' @param do_beta_posterior_density a logical value. If TRUE, the posterior density 
#'         of the coefficient function is computed. (optional)
#' @param sann a logical value. If TRUE, the Bliss estimate is computed. (optional)
#' @param progress a logical value. If TRUE, the algorithm progress is displayed.
#'         (optional)
#' @importFrom utils tail
#' @export
#' @examples
#' # see the vignette BlissIntro.
Bliss   <- function(data,param,beta_posterior_density=FALSE,sann=FALSE,
                    progress=FALSE){
 # define Q
 Q <- param[["Q"]]
 if(is.null(Q))
  stop("Please specify Q: the number of functional covariates.")
 # Define p
 param$p  <- sapply(param$grids,length)
 # Centering the data 
 data$x <- sapply(data$x,function(m) scale(m,scale=F))
 
 # How many chains i have to do ?
 if(is.null(param[["n_chains"]])){
  param[["n_chains"]] <- 1
 }
 n_chains <- param[["n_chains"]]
 # Initialize the list "chains"
 chains <- list()
 
 # For each chain :
 for(j in 1:n_chains){
  if(progress & n_chains > 1) cat("Chain ",j,": \n",sep="")
  chains[[j]] <- list()
  
  # Execute the Gibbs Sampler algorithm to sample the posterior distribution
  chains[[j]]$posterior_sample <- Bliss_Gibbs_Sampler(data,param)
  
  # Compute a posterior sample of coefficient function
  chains[[j]]$beta_sample <- compute_beta(chains[[j]]$posterior_sample,param)
 }
 
 # Choose a chain for inference
 j <- sample(n_chains,1)
 posterior_sample <- chains[[j]]$posterior_sample
 beta_sample      <- chains[[j]]$beta_sample
 
 # Compute an approximation of the posterior density of the coefficient function
 beta_posterior_density <- list()
 if (do_beta_posterior_density){
  for(q in 1:Q){
   diff_grid <- diff(param$grids[[q]])[1]
   param$grids2[[q]] <- c(param$grids[[q]]-diff_grid/2,
                          tail(param$grids[[q]],1)+diff_grid/2) #XXXXXXXX
   param$xlim[[q]] <- range(param$grids2[[q]])
   
   param_density <- list(grid= param$grids[[q]], #XXXXXXXX
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
   beta_posterior_density[[q]] <-
    compute_beta_posterior_density(beta_sample[[q]],param_density) #XXXXXXXX
  }
 }
 
 # Execute the Simulated Annealing algorithm to estimate the coefficient function
 Bliss_estimation <- list()
 Bliss_estimate <- list()
 for(q in 1:Q){
  param$Simulated_Annealing <- list( grid = param$grids[[q]], #XXXXXXXX
                      iter = param[["iter"]],
                      p    = param$p[q],
                      Temp = param$Temp[q],
                      k_max = param$k_max[q],
                      iter_sann = param[["iter_sann"]],
                      burnin    = param[["burnin"]],
                      l_max     = param[["l_max_sann"]][q],
                      basis     = param[["basis"]][q])
  
  Bliss_estimation[[q]] <- Bliss_Simulated_Annealing(beta_sample[[q]],
                                                     param$Simulated_Annealing,
                                                     posterior_sample$param$scale_ml[[q]])  #XXXXXXXX
  Bliss_estimate[[q]] <- Bliss_estimation[[q]]$Bliss_estimate
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
 
 # Do not return the list "chains" if n_chains is 1.
 if(n_chains == 1) chains <- NULL
 
 # compute the posterior density of theta=(mu,(b)_k,(m)_k,(l)_k,sigma)
 if(length(param$K)==1) #XXXXXXXX
  if(param$theta_posterior_density){
   if(display) cat("Compute the (log) densities. \n")  #XXXXXXXX mettrte lkh et prior
   posterior_sample$posterior_density <- dposterior(posterior_sample,data)
  }
 
 # The object to return
 res <- list(alpha                  = alpha,
             beta_posterior_density = beta_posterior_density,
             beta_sample            = beta_sample,
             Bliss_estimate         = Bliss_estimate,
             Bliss_estimation       = Bliss_estimation, # A enlever ?
             chains                 = chains,
             param                  = param, # voir si c'est utile de le renvoyer et s'il faut en renvoyer d'autres aussi
             posterior_sample       = posterior_sample,
             support_estimate       = support_estimate
 )
 class(res) = c("bliss")
 return(invisible(res))
}
### XXXXXXXXXXX : changer ce qui suit en fonction du resultat final

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
