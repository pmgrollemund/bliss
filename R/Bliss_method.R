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
#' \item{new_grids}{a list of Q vectors (optional) to compute beta samples on
#'       different grids.}
#' \item{Temp_init}{a nonnegative value (optional), the initial temperature for
#'      the cooling function of the Simulated Annealing algorithm.}
#' \item{thin}{an integer (optional) to thin the posterior sample.}
#' \item{times_sann}{an integer (optional), the number of times the algorithm
#'       will be executed}
#' \item{times_sann}{an integer (optional), the number of times the algorithm
#'       will be executed}
#' \item{allow_reducing}{a boolean value (optional), indicate if the function is
#'       allowed to reduce the number of sample times of each functional covariate.}
#' \item{verbose_cpp}{a boolean value (optional). Write stuff from the Rcpp scripts
#'       if TRUE.}
#' }
#' @param compute_density a logical value. If TRUE, the posterior density
#'         of the coefficient function is computed. (optional)
#' @param sann a logical value. If TRUE, the Bliss estimate is
#'         computed with a Simulated Annealing Algorithm. (optional)
#' @param sann_trace a logical value. If TRUE, the trace of the Simulated
#'         Annealing algorithm is included into the result object. (optional)
#' @param support_estimate a logical value. If TRUE, the estimate of the
#' coefficient function support is computed. (optional)
#' @param verbose write stuff if TRUE (optional).
#' @importFrom utils tail
#' @export
#' @examples
#' # see the vignette BlissIntro.
fit_Bliss <- function(data,param,
                      sann=TRUE,compute_density=TRUE,support_estimate=TRUE,
                      sann_trace=FALSE,
                      verbose=TRUE){
  if(verbose) cat("===========================================================\n")
  ###### Initialisation
  param$n <- length(data$y)
  param$p <- sapply(data$grids,length)
  param$Q <- length(data$x)
  if(is.null(param$grids)) param$grids <- data$grids
  param$sann_trace <- sann_trace
  param$verbose <- verbose
  if(is.null(param[["allow_reducing"]])) param[["allow_reducing"]] <- TRUE

  ###### Alerte
  if(!is.list(data$x))
    stop("data$x must be a list.")
  if(!is.list(data$grids))
    stop("data$grids must be a list.")
  if(is.null(data$grids))
    stop("data must contain an element 'grids'.")
  if(length(data$x) != length(data$grids))
    stop("data$x and data$grids must be related to the same amount of functional covariate. \n\t length(data$x) == length(data$grids)")

  ##### Pretreatmentt
  # If p is too large for some functional covariate, it is better
  # to reduce this number in order to get faster computation.c
  if(verbose) cat("- Pretreatment.\n")
  reduction_required <- do_need_to_reduce(param)
  if(reduction_required & param[["allow_reducing"]]){
    warning(paste("We need to reduce the dimension of at least one functional ",
                  "covariate, otherwise computational time would be too large. ",
                  "If you do want to overpass this, please indicate: \n\t",
                  "param$allow_reducing <- FALSE\n",
                  sep="")
    )
    data <- reduce_x(data,param)
    param$p <- sapply(data$grids,length)
    param$grids <-data$grids
  }

  ###### Fit the Bliss model
  if(verbose) cat("- Model fiting.\n")
  # Sample from the posterior distribution
  if(verbose) cat("   Gibbs Sampler: Sample from the posterior distribution.\n")
  posterior_sample <- Bliss_Gibbs_Sampler(data,param,verbose)

  # Compute a posterior sample of coefficient function
  if(verbose) cat("   Compute beta functions related to the posterior sample.\n")
  beta_sample <- compute_beta_sample(posterior_sample,param)

  ###### Compute estimates
  if(verbose) cat("- Derive estimates.\n")
  # Run the Simulated Annealing algorithm to estimate the coefficient function
  if(sann){
    if(verbose) cat("   Simulated Annealing (to get the Bliss and smooth estimates):\n")
    res_sann <- Bliss_Simulated_Annealing(beta_sample,posterior_sample,param)

    Bliss_estimate  <- res_sann$Bliss_estimate
    trace_sann      <- res_sann$trace
    Smooth_estimate <- res_sann$Smooth_estimate
  }else{
    Bliss_estimate   <- list()
    trace_sann       <- list()
    Smooth_estimate  <- list()
  }

  # Compute an approximation of the posterior density of the coefficient function
  if (compute_density){
    if(verbose) cat("   Compute the approximation of the posterior distribution.\n")
    beta_posterior_density <- compute_beta_posterior_density(beta_sample,param)
  }else{
    beta_posterior_density <- list()
  }

  # Compute the support estimate
  if(support_estimate){
    if(verbose) cat("   Support estimation.\n")
    res_support <- support_estimation(beta_sample,param)

    support_estimate     <- res_support$estimate
    support_estimate_fct <- res_support$estimate_fct
    alpha                <- res_support$alpha

    rm(res_support)
  }else{
    support_estimate <- list()
    support_estimate_fct <- list()
    alpha <- list()
  }

  ###### Compute posterior quantities
  if(verbose) cat("- Posterior quantities.\n")

  # Densities
  if(verbose) cat("   Compute the (log) densities of the posterior sample. \n")
  posterior_sample$posterior_density <- dposterior(posterior_sample,data)

  if(sann){
    if(verbose) cat("   Fitting posterior y values. \n")
    # Fitted values
    y_fitted <- predict_bliss(data$x,data$grids,param$burnin,posterior_sample,Smooth_estimate)
  }else{
    y_fitted <- NULL
  }

  # Fitted value distribution
  y_fitted_distribution <-
    predict_bliss_distribution(data$x,data$grids,param$burnin,posterior_sample,
                               beta_sample)

  # Usual quantities (freq and bayes)
  if(verbose) cat("   Usual quantities. \n")
  res_post_treatment <- post_treatment_bliss(posterior_sample,param,data)
  BIC_value <- res_post_treatment$BIC
  loglik <- res_post_treatment$loglik
  nb_param <- res_post_treatment$nb_param
  MSE <- mean( (data$y - y_fitted)^2)

  ###### Output
  # The object to return
  res <- list(alpha                  = alpha,
              beta_posterior_density = beta_posterior_density,
              beta_sample            = beta_sample,
              Bliss_estimate         = Bliss_estimate,
              data                   = data,
              y_fitted               = y_fitted,
              y_fitted_distribution  = y_fitted_distribution,
              MSE                    = MSE,
              posterior_sample       = posterior_sample,
              smooth_estimate        = Smooth_estimate,
              support_estimate       = support_estimate,
              support_estimate_fct   = support_estimate_fct,
              trace_sann             = trace_sann,
              BIC                    = BIC_value,
              loglik                 = loglik,
              nb_param               = nb_param
  )
  class(res) = c("bliss")

  if(verbose)
    printbliss(res)

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
  to_print <- c(nb_param=x$nb_param,
                loglik=x$loglik,
                BIC=x$BIC,
                MSE=x$MSE
  )

  cat("===========================================================\n")
  print(to_print)
  cat("===========================================================\n")

  cat("* Useful fields \n")
  cat("   $Bliss_estime, $smooth_estimate, $support_estimate, $alpha,\n")
  cat("   $beta_posterior_density, $posterior_sample, $beta_sample \n")
  cat("* Useful post-treatment functions \n")
  cat("   image_Bliss(), interpretation_plot() \n")

}

