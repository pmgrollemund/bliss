################################# ----
#' fit_Bliss
################################# ----
#' @description Fit the Bayesian Functional
#' Linear Regression model (with Q functional covariates).
#' @return return a list containing:
#' \describe{
#'  \item{alpha}{a list of Q numerical vector. Each vector is the function
#'        alpha(t) associated to a functional covariate. For each t, alpha(t)
#'        is the posterior probabilities of the event "the support covers t".}
#'  \item{beta_posterior_density}{a list of Q items. Each item contains a list
#'        containing information to plot the posterior density of the
#'        coefficient function with the \code{image} function.
#'        \describe{
#'        \item{\code{grid_t}}{a numerical vector: the x-axis.}
#'        \item{\code{grid_beta_t}}{a numerical vector: the y-axis.}
#'        \item{\code{density}}{a matrix: the z values.}
#'        \item{\code{new_beta_sample}}{a matrix: beta sample used to compute
#'              the posterior densities.}
#'        }
#'        }
#'  \item{beta_sample}{a list of Q matrices. The qth matrix is a posterior
#'        sample of the qth functional covariates.}
#'  \item{Bliss_estimate}{a list of numerical vectors corresponding to the
#'        Bliss estimates of each functional covariates.}
#'  \item{chains}{a list of \code{posterior_sample}. \code{chains} is \code{NULL} if
#'        \code{n_chains}=1.}
#'  \item{chains_info}{a list for each chain providing: a mu estimate, a sigma_sq estimate,
#'  the Smooth estimate of the coefficient function and the autocorrelation of the
#'  Markov Chain.}
#'  \item{data}{a list containing the data.}
#'  \item{posterior_sample}{a list of information about the posterior sample:
#'        the trace matrix of the Gibbs sampler, a list of Gibbs sampler parameters
#'        and the posterior densities.}
#'  \item{support_estimate}{a list of support estimates of each functional covariate.}
#'  \item{support_estimate_fct}{another version of the support estimates.}
#'  \item{trace_sann}{a list of Q matrices which are the trace of the
#'        Simulated Annealing algorithm.}
#' }
#' @param data a list containing:
#' \describe{
#' \item{Q}{an integer, the number of functional covariates.}
#' \item{y}{a numerical vector, the outcomes.}
#' \item{x}{a list of matrices, the qth matrix contains the observations of the
#'       qth functional covariate at time points given by \code{grids}.}
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of
#'        time points for the qth functional covariate.}
#' }
#' @param param a list containing:
#' \describe{
#' \item{iter}{an integer, the number of iterations of the Gibbs sampler algorithm.}
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for
#'       each covariate.}
#' \item{basis}{a character vector (optional). The possible values are "uniform" (default),
#'       "epanechnikov", "gauss" and "triangular" which correspond to
#'       different basis functions to expand the coefficient function and the
#'       functional covariates}
#' \item{burnin}{an integer (optional), the number of iteration to drop from the
#'       posterior sample.}
#' \item{iter_sann}{an integer (optional), the number of iteration of the Simulated
#'       Annealing algorithm.}
#' \item{k_max}{an integer (optional), the maximal number of intervals for the
#'       Simulated Annealing algorithm.}
#' \item{l_max}{an integer (optional), the maximal interval length for the
#'       Simulated Annealing algorithm.}
#' \item{lims_kde}{an integer (optional), correspond to the \code{lims} option
#'       of the \code{kde2d} funtion.}
#' \item{n_chains}{an integer (optional) which corresponds to the number of
#'       Gibbs sampler runs.}
#' \item{new_grids}{a list of Q vectors (optional) to compute beta samples on
#'       different grids.}
#' \item{Temp_init}{a nonnegative value (optional), the initial temperature for
#'      the cooling function of the Simulated Annealing algorithm.}
#' \item{thin}{an integer (optional) to thin the posterior sample.}
#' \item{times_sann}{an integer (optional), the number of times the algorithm
#'       will be executed}
#' }
#' @param compute_density a logical value. If TRUE, the posterior density
#'         of the coefficient function is computed. (optional)
#' @param sann a logical value. If TRUE, the Bliss estimate is
#'         computed with a Simulated Annealing Algorithm. (optional)
#' @param support_estimate a logical value. If TRUE, the estimate of the
#' coefficient function support is computed. (optional)
#' @param verbose write stuff if TRUE (optional).
#' @importFrom utils tail
#' @export
#' @examples
#' # see the vignette BlissIntro.
fit_Bliss <- function(data,param,compute_density=TRUE,sann=TRUE,
                      support_estimate=TRUE,
                      verbose=FALSE){
  # Define Q
  Q <- data[["Q"]]
  if(is.null(Q))
    stop("Please specify Q: the number of functional covariates (in the 'data' object).")
  # Define p
  param$p  <- sapply(data$grids,length)
  # Centering the data
  data$x_save <- data$x
  for(q in 1:Q){
    data$x[[q]] <- scale(data$x[[q]],scale=F)
  }

  # How many chains i have to do ?
  if(is.null(param[["n_chains"]])){
    param[["n_chains"]] <- 1
  }
  n_chains <- param[["n_chains"]]
  # Initialize the list "chains"
  chains <- list()
  chains_info <- list()

  if(verbose) cat("Sample from the posterior distribution.\n")
  # For each chain :
  for(j in 1:n_chains){
    if(verbose & n_chains > 1) cat("Chain ",j,": \n",sep="")
    chains[[j]] <- list()

    # Execute the Gibbs Sampler algorithm to sample the posterior distribution
    param_Gibbs_Sampler <- list(iter  = param[["iter"]],
                                K     = param[["K"]],
                                basis = param[["basis"]],
                                p     = param[["p"]],
                                phi_l = param[["phi_l"]],
                                grids = data[["grids"]])
    chains[[j]] <- Bliss_Gibbs_Sampler(data,param_Gibbs_Sampler,verbose)

    chains_info[[j]] <- compute_chains_info(chains[[j]],param_Gibbs_Sampler)
  }

  # Choose a chain for inference
  j <- sample(n_chains,1)
  posterior_sample <- chains[[j]]

  # Compute a posterior sample of coefficient function
  if(verbose) cat("Coefficient function: smooth estimate.\n")
  beta_sample <- compute_beta_sample(posterior_sample,param_Gibbs_Sampler,Q)

  # Execute the Simulated Annealing algorithm to estimate the coefficient function
  Bliss_estimation <- list()
  Bliss_estimate   <- list()
  trace_sann       <- list()
  Smooth_estimate  <- list()
  if(sann){
    if(verbose) cat("Coefficient function: Bliss estimate.\n")
    for(q in 1:Q){
      param_Simulated_Annealing <- list( grid = data[["grids"]][[q]],
                                         iter = param[["iter"]],
                                         p    = param[["p"]][q],
                                         Temp_init = param[["Temp_init"]],
                                         K    = param[["K"]][q],
                                         k_max = param[["k_max"]][q],
                                         iter_sann = param[["iter_sann"]],
                                         times_sann= param[["times_sann"]],
                                         burnin    = param[["burnin"]],
                                         l_max     = param[["l_max"]][q],
                                         basis     = param[["basis"]][q])

      Bliss_estimation[[q]] <- Bliss_Simulated_Annealing(beta_sample[[q]],
                                                         posterior_sample$param$normalization_values[[q]],
                                                         param_Simulated_Annealing)
      Bliss_estimate[[q]]  <- Bliss_estimation[[q]]$Bliss_estimate
      trace_sann[[q]]      <- Bliss_estimation[[q]]$trace
      Smooth_estimate[[q]] <- Bliss_estimation[[q]]$Smooth_estimate
    }
    rm(Bliss_estimation)
  }

  # Compute an approximation of the posterior density of the coefficient function
  beta_posterior_density <- list()
  if (compute_density){
    if(verbose) cat("Compute the approximation of the posterior distribution.\n")
    for(q in 1:Q){
      param_beta_density <- list(grid= data[["grids"]][[q]],
                                 iter= param[["iter"]],
                                 p   = param[["p"]][q],
                                 n        = length(data[["y"]]),
                                 thin     = param[["thin"]],
                                 burnin   = param[["burnin"]],
                                 lims_kde = param[["lims_kde"]][[q]],
                                 new_grid = param[["new_grids"]][[q]],
                                 lims_estimate = range(Smooth_estimate[[q]]),
                                 verbose = verbose)

      beta_posterior_density[[q]] <-
        compute_beta_posterior_density(beta_sample[[q]],param_beta_density)
    }
  }

  # Compute the support estimate
  if(support_estimate){
    if(verbose) cat("Support estimation.\n")
    support_estimate <- list()
    support_estimate_fct <- list()
    alpha <- list()
    for(q in 1:Q){
      res_support <- support_estimation(beta_sample[[q]])
      support_estimate[[q]]     <- res_support$estimate
      support_estimate_fct[[q]] <- res_support$estimate_fct
      alpha[[q]]                <- res_support$alpha
    }
    rm(res_support)
  }else{
    support_estimate <- list()
    support_estimate_fct <- list()
    alpha <- list()
  }

  # Do not return the list "chains" if n_chains is 1.
  if(n_chains == 1) chains <- NULL

  if(verbose) cat("Compute the (log) densities of the posterior sample. \n")
  posterior_sample$posterior_density <- dposterior(posterior_sample,data)

  # The object to return
  res <- list(alpha                  = alpha,
              beta_posterior_density = beta_posterior_density,
              beta_sample            = beta_sample,
              Bliss_estimate         = Bliss_estimate,
              chains                 = chains,
              chains_info            = chains_info,
              data                   = data,
              posterior_sample       = posterior_sample,
              Smooth_estimate        = Smooth_estimate,
              support_estimate       = support_estimate,
              support_estimate_fct   = support_estimate_fct,
              trace_sann             = trace_sann
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
#' # See fit_Bliss() function
printbliss<-function(x,...){
  if(!any(class(x) == "bliss"))
    stop("Input must have class \"bliss\".")

  cat("This is a bliss object:\n")
  cat("----- \n")
  cat("\n")
  print(str(x))
}

