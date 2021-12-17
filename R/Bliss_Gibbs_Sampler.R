################################# ----
#' Bliss_Gibbs_Sampler
################################# ----
#' @description A Gibbs Sampler algorithm to sample the posterior distribution of
#'              the Bliss model.
#' @return a list containing :
#' \describe{
#'  \item{trace}{a matrix, the trace of the Gibbs Sampler.}
#'  \item{param}{a list containing parameters used to run the function. }
#' }
#' @param data a list containing:
#' \describe{
#' \item{Q}{an integer, the number of functional covariates.}
#' \item{y}{a numerical vector, the outcome values \code{y_i}.}
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
#' \item{p}{an integer, the number of time points.}
#' \item{basis}{a character vector (optional). The possible values are "uniform" (default),
#'       "epanechnikov", "gauss" and "triangular" which correspond to
#'       different basis functions to expand the coefficient function and the
#'       functional covariates}
#' \item{phi_l}{a numerical (optional). An hyperparameters related to the exponential prior
#' on the length of the intervals. Lower values promotes wider intervals.}
#' }
#' @param verbose write stuff if TRUE (optional).
#' @importFrom stats var
#' @export
#' @examples
#' \donttest{
#' # May take a while
#' param_sim <- list(Q=1,n=25,p=50,grids_lim=list(c(0,1)),iter=1e4,K=2)
#' data_sim <- sim(param_sim,verbose=FALSE)
#' res_Bliss_Gibbs_Sampler <- Bliss_Gibbs_Sampler(data_sim,param_sim)
#' theta_1 <- res_Bliss_Gibbs_Sampler$trace[1,]
#' theta_1
#' # Resultat for few iterations
#' param_sim <- list(Q=1,n=25,p=50,grids_lim=list(c(0,1)),iter=5e2,K=2)
#' data_sim <- sim(param_sim,verbose=FALSE)
#' res_Bliss_Gibbs_Sampler <- Bliss_Gibbs_Sampler(data_sim,param_sim)
#' theta_1 <- res_Bliss_Gibbs_Sampler$trace[1,]
#' theta_1
#' }
Bliss_Gibbs_Sampler <- function(data,param,verbose=FALSE){
  # load objects
  x     <- data[["x"]]
  y     <- data[["y"]]
  Q     <- data[["Q"]]
  grids <- data[["grids"]]
  basis <- param[["basis"]]
  iter  <- param[["iter"]]
  K     <- param[["K"]]
  p     <- param[["p"]]
  phi_l <- param[["phi_l"]]

  # Initialize the required unspecified objects
  if(is.null(K)) stop("Please specify a value for the vector K.")
  if(!is.null(K)){
    for (q in 1:Q){
      if(is.na(K[q]))
        stop("Please specify a value for all components of the vector K.")
    }
  }
  if(is.null(basis)){
    basis <- character()
    for (q in 1:Q){
      basis[q] <- "Uniform"
    }
  }
  if(is.null(phi_l)){
    phi_l <- 5
  }
  if(!is.null(basis)){
    for (q in 1:Q){
      if(is.na(basis[q])){basis[q] <- "Uniform"}
    }
  }
  lambda <- 5 # a mettre dans le cpp ?
  V_tilde <- diag(1+sum(K)) # a mettre dans le cpp ?
  g <- length(y) # a mettre dans le cpp ?

  # Determine the possible values of l
  l_values <- list()
  for (q in 1:Q){
    l_values[[q]] <- (grids[[q]] - grids[[q]][1])[-1]
  }
  l_values_length <- sapply(l_values,length)

  # Determine the prior distribution of l
  Probs_l <- l_values
  for(q in 1:Q){
    # Probs_l[[q]] <- pdexp( 5*K[q] , l_values[[q]] )
    Probs_l[[q]] <- pdexp( phi_l*K[q] , l_values[[q]] )
  }

  if(verbose){
    verbose_cpp <- TRUE
  }else{
    verbose_cpp <- FALSE
  }

  # Determine if there are columns for which there is no variation
  pb_count <- 0
  for(q in 1:length(x)){
    index <- which(apply(x[[q]],2, function(v) length(unique(v)) == 1))
    if(length(index) > 0){
      pb_count <- pb_count + 1
      if(q == 1) number_exposant <- "st"
      if(q == 2) number_exposant <- "nd"
      if(q >  3) number_exposant <- "th"
      if(verbose) cat(paste("\tFor the ",q,number_exposant," functional covariate, ",
                      "all the individus have the same values for the columns:\n\t  ",
                      paste(index,collapse = " "),".\n\tThis would",
                      " eventually leads to non-invertible partial design ",
                      "matrixes.\n",sep=""))
    }
  }
  if(pb_count > 0){
    stop(paste("Please provide design columns with non-unique values.",
    "You could jitter the columns with unique values, or remove them."))
  }
  rm(pb_count)

  # Perfome the Gibbs Sampler and return the result.
  res <- Bliss_Gibbs_Sampler_cpp(Q,y,x,grids,
                                 iter,K,basis,
                                 g,lambda,V_tilde, l_values_length,Probs_l,
                                 verbose_cpp,tol=sqrt(.Machine$double.eps))

  trace_names <- NULL
  if(Q == 1){
    for(k in 1:K[q]){
      trace_names <- c(trace_names,paste("b",k,sep="_"))
    }
    for(k in 1:K[q]){
      trace_names <- c(trace_names,paste("m",k,sep="_"))
    }
    for(k in 1:K[q]){
      trace_names <- c(trace_names,paste("l",k,sep="_"))
    }
  }else{
    for(q in 1:Q){
      for(k in 1:K[q]){
        trace_names <- c(trace_names,paste("b",k,"q",q,sep="_"))
      }
      for(k in 1:K[q]){
        trace_names <- c(trace_names,paste("m",k,"q",q,sep="_"))
      }
      for(k in 1:K[q]){
        trace_names <- c(trace_names,paste("l",k,"q",q,sep="_"))
      }
    }
  }
  colnames(res$trace) <- c(trace_names ,"mu","sigma_sq")
  return(res)
}
