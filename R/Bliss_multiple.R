################################# ----
#' Bliss_multiple
################################# ----
#' @description perform the Bliss method to obtain an estimate
#'               of the Q coefficient functions of the functional linear
#'               regression model with Q functional covariates.
#' @return return a list containing :
#' \describe{
#'  \item{Bliss_estimate}{a list of numerical vectors, the Bliss estimates : we obtain Q estimates, one for each covariable.}
#'  \item{posterior_density_estimate}{a list of lists containing the estimates of the posterior density.
#'            (obtained with the kde2d function). Firts level of the list: Q lists, for the Q covariates. Second
#'            level: for each covariate we have a list of three components, that is the  x and y coordinates of the grid points, and z
#'            a matrix of the estimated density: rows correspond to the value of x, columns to the value of y.}
#'  \item{beta_functions}{a list of matrices. For the qth covariate, beta_functions[[q]] is a matrix where each row is a
#'                 function beta_qi associated to the iteration i of the Gibbs sampler.}
#'  \item{res.Simulated_Annealing}{a list of lists: res.Simulated_Annealing[[q]] is the result of the
#'                 function Bliss_Simulated_Annealing applied for the qth covariate.}
#'  \item{res.Gibbs_Sampler}{a list of lists: res.Gibbs_Sampler[[q]] is the the result of the function
#'                 Bliss_Gibbs_Sampler for the qth covariate.}
#' }
#' @param data a list containing
#' \describe{
#' \item{1)}{the number of covariates}
#' \item{2)}{the functions x_qi(t) observed at grids of time points}
#' \item{3)}{the outcome values y_i}
#' }
#' @param param a list
#' \describe{
#' \item{iter}{an integer, the number of iterations of the Gibbs sampler algorithm.}
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of observation points for the qth covariate.}
#' \item{K}{a vector of integers, hyperparameters of the Bliss model, the number of intervals for the Q covariates.}
#' \item{l_max}{a vector of integers, hyperparameters of the Bliss model. Beware, give the index
#'              corresponding to the value of the maximum width of the intervals (optional)}
#' \item{basis}{a vector of characters among : "uniform" (default),
#'                 "epanechnikov", "gauss" and "triangular". This indicates the
#'                 shapes of the Q coefficient functions on the intervals. (optional)}
#' \item{eta_tilde}{a numerical vector of length (1+sum(K)), hyperparameter of the Bliss model.
#'                 By default, eta is the vector (0,...,0). (optional)}
#' \item{V_tilde}{a matrix of dimension (1+sum(K))*(1+sum(K)), hyperparameter of the Bliss model. (optional) (nasty code)}
#' \item{a}{a nonnegative value, hyperparameter of the Bliss model. Default, a = 0.1. (optional)}
#' \item{b}{a nonnegative value, hyperparameter of the Bliss model. Default, b = 0.1. (optional)}
#' \item{g}{hyperparameter of the Bliss model,  a nonnegative value,
#'                the coefficient of the Zellner prior.}
#' \item{phi_m}{a list of numerical vectors. The priors of the mq, q=1,...,Q. If not specified, a uniform
#'                 distribution is used for phi_m[[q]].}
#' \item{phi_l}{a list of numerical vectors and/or characters. The priors of the lq, q=1,...,Q. If not specified, a uniform
#'                 distribution is used for phi_l[[q]].}
#' \item{phi_l_mean}{a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#'                corresponds to the mean of the Gamma prior of lq.}
#' \item{phi_l_sd}{a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#'                corresponds to the standard deviation of the Gamma prior of lq.}
#' \item{prior_beta}{a character string, which indicates the prior on the
#'                beta_star[q]. The possible values are :
#'                1) "diag" (default) for a diagonal matrix prior,
#'                2) "Ridge_Zellner" for the Ridge Zellner prior}
#' \item{burnin}{an integer, the number of iteration to drop of the Gibbs sampler. (optional)}
#' \item{thin}{an integer, used to thin the Gibbs sample to compute an
#'                 estimate of the posterior density of beta(t). (optional)}
#' \item{lims.kde}{a Qx2 matrix. lims.kde[q,] corresponds to the parameter (yl,yu) for the representation
#'                 of the posterior density of the qth covariable. There are limits of the y-axis for the function kde2d. (optional)}
#' \item{n}{an integer, the number of grid points in each direction, for the kde2d function.
#'                 See function kde2d. (optional)}
#' \item{iter_sann}{an integer, the number of iteration of the Simulated Annealing. (optional)}
#' \item{Temp}{a vector of nonnegative values, the Q initial temperatures for the
#'                 cooling function of the Q Simulated Annealings. (optional)}
#' \item{k_max}{a vector of integers, k_max[q] is the maximum number of intervals for the
#'                 function beta_q(t) at each iteration. (optional)}
#' \item{cols}{a vector of colors for the function image. Only if plot=TRUE. (optional)}
#' \item{new_grids}{a list of Q numerical vectors. If new_grids is not NULL, the
#'                 coefficient functions beta_q(t) at each iteration are computed on these grids
#'                 (only) to plot a graphical representation of the posterior distribution.}
#' \item{h1}{a vector of numerical values which are the Q bandwidths of the kernel
#'                 density estimation for the t-axis (optional). h1[q] the bandwidth of the kernel density estimation for the qth covariate.}
#' \item{ylim}{a Qx2 matrix, the qth line gives the limits for the y-axis for the plotting function
#'                 "image_Bliss" for the qth covariate.}
#' \item{main}{a vector of characters, main[q] is for the plotting function
#'                 "image_Bliss" for the qth covariate.}
#' \item{n_chains}{number of chains to do in the Gibbs sampler.}
#' }
#' @param plot a logical value. If TRUE, the Bliss estimates and
#'                 the estimates of the posterior density are plotted. (optional)
#' @param density a logical value. If TRUE, the posterior density is estimated. (optional)
#' @importFrom utils tail
#' @export
#' @examples
#' # see the vignette BlissIntro.
Bliss_multiple   <- function(data,param,plot=FALSE,density=FALSE){
  display <- param$display
  if(is.null(display)) display <- TRUE
  # preprocessing
  Q <- length(data$x_mult)
  p  <- numeric()
  for (q in 1:Q){
    data$x_mult[[q]]   <- apply(data$x_mult[[q]],2,function(vect) vect - mean(vect))
    p[q] <- length(param$grids[[q]])
  }
  param$p <- p

  if(is.null(param$density))           param$density <- TRUE
  if(is.null(param$sann))              param$sann <- TRUE
  if(is.null(param$compute_posterior)) param$compute_posterior <- TRUE

  # How many chains i have to do ?
  n_chains <- param[["n_chains"]]
  if(is.null(n_chains)){
    n_chains <- 1
    param[["n_chains"]] <- n_chains
  }

  # Initialize the list "chains"
  chains <- list()
  # Each chain :
  for(j in 1:n_chains){
    if(display) cat("Chain ",j,": \n",sep="")
    cat("Chain ",j,": \n",sep="")
    # Initialize the list "chains[[j]]"
    chains[[j]] <- list()
    # Execute the Gibbs Sampler algorithm to obtain a posterior sample
    chains[[j]]$res.Gibbs_Sampler <- Bliss_Gibbs_Sampler_multiple(data,param)

    # Compute the functions beta_i for each iteration i of the Gibbs sample.
    chains[[j]]$beta_functions <- compute_beta_functions_mult(chains[[j]]$res.Gibbs_Sampler,param)

    # Estimate the density of the posterior sample of functions beta_i
    chains[[j]]$posterior_density_estimate <- list()
    param_density <- list()
    # IS 16/02/2018: condition sur calcul density
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

        chains[[j]]$posterior_density_estimate[[q]] <- density_estimation(chains[[j]]$beta_functions[[q]],param_density)
      }
    }

  }


  # Choose a chain
  j <- sample(n_chains,1)

  res.Gibbs_Sampler          <- chains[[j]]$res.Gibbs_Sampler
  beta_functions             <- chains[[j]]$beta_functions
  posterior_density_estimate <- chains[[j]]$posterior_density_estimate

  # Execute the Simulated Annealing algorithm to obtain an estimate
  res.Simulated_Annealing <- list()
  for(q in 1:Q){
    beta_functions_tmp <- beta_functions[[q]]
    scale_tmp <- res.Gibbs_Sampler$param$scale_ml[[q]]
    param_sann <- list( grid = param$grids[[q]],
                        iter = param[["iter"]],
                        p    = param$p[q],
                        Temp = param$Temp[q],
                        k_max = param$k_max[q],
                        iter_sann = param[["iter_sann"]],
                        burnin    = param[["burnin"]],
                        l_max     = param[["l_max_sann"]][q],
                        basis     = param[["basis"]][q])

    res.Simulated_Annealing[[q]] <- Bliss_Simulated_Annealing(beta_functions_tmp,param_sann,scale_tmp)
  }

  # Compute the estimate of the support
  if(display) cat("Estimation of the support. \n")
  res.support <- list()
  res.support_estimate <- list()
  res.alpha_t <- list()
  for(q in 1:Q){
    res.support[[q]] <- support_estimation(beta_functions[[q]])
    res.support_estimate[[q]] <- res.support[[q]]$estimate
    res.alpha_t[[q]] <- res.support[[q]]$alpha_t
  }

  # Plot the Bliss estimate and the posterior density estimate
  if(density & plot){
    for(q in 1:Q){
      image_Bliss(posterior_density_estimate[[q]],param)
      lines(param$grid,res.Simulated_Annealing[[q]]$Bliss_estimate,type="s",lwd=2)
      #lines(param$grid,res.Simulated_Annealing$posterior_expe,type="l",lty=2)
    }
  }

  # Do not return the list "chains" if n_chains is 1.
  if(n_chains == 1) chains <- NULL

  # compute the posterior density on the MCMC sample
  if(length(param$K)==1)
    if(param$compute_posterior){
      if(display) cat("Compute the (log) densities. \n")
      res.Gibbs_Sampler$posterior_density <-
        dposterior(res.Gibbs_Sampler,data)
    }

  # The object to return
  Bliss_estimate <- list()
  for(q in 1:Q){
    Bliss_estimate[[q]] <- res.Simulated_Annealing[[q]]$Bliss_estimate
  }

  res <- list(Bliss_estimate             = Bliss_estimate,
              support_estimate           = res.support_estimate,
              alpha_t                    = res.alpha_t,
              posterior_density_estimate = posterior_density_estimate,
              beta_functions             = beta_functions,
              res.Simulated_Annealing    = res.Simulated_Annealing,
              res.Gibbs_Sampler          = res.Gibbs_Sampler,
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
