#########################################################
#                                                       #
#    Bliss method : MCMC and optimization algorithms    #
#                                                       #
#########################################################
################################# ----
#' Bliss_Gibbs_Sampler_multiple
################################# ----
#' @description Gibbs Sampler to sample from the posterior distribution of the Bliss model.
#' @return a list containing :
#' \describe{
#'  \item{trace}{a matrix. Each row is an iteration of the Gibbs Sampler.}
#'  \item{param}{a list containing a, b, V_tilde, K, eta_tilde, l_max, grids and scale_ml}
#' }
#' @param data a list containing
#' \describe{
#'  \item{x_mult}{a list containing the functions x_qi(t) observed at grids of time points}
#'  \item{y}{the outcome values y_i}
#'  \item{grids}{a list containing the observation time points of the covariates functions.}
#' }
#' @param param a list which have to contain :
#' \describe{
#' \item{iter}{an integer, the number of iterations of the Gibbs sampler algorithm.}
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of observation time points for the qth covariate.}
#' \item{K}{a vector of integers, hyperparameters of the Bliss model, the number of intervals for the Q covariates.}
#' \item{l_max}{a vector of integers, hyperparameters of the Bliss model. }
#' \item{basis}{a vector of characters among : "uniform" (default),
#'                 "epanechnikov", "gauss" and "triangular". This indicates the
#'                 shapes of the Q coefficient functions on the intervals. (optional)}
#' \item{eta_tilde}{a numerical vector of length (1+sum(K)), hyperparameter of the Bliss model.
#'                 By default, eta is the vector (0,...,0). (optional)}
#' \item{V_tilde}{a matrix of dimension (1+sum(K))*(1+sum(K)), hyperparameter of the Bliss model. (optional)}
#' \item{a}{a nonnegative value, hyperparameter of the Bliss model. By default, a = 0.1. (optional)}
#' \item{b}{a nonnegative value, hyperparameter of the Bliss model. By default, b = 0.1. (optional)}
#' \item{g}{hyperparameter of the Bliss model,  a nonnegative value, the coefficient of the Zellner prior.}
#' \item{phi_m}{a list of numerical vectors. The priors of the mq, q=1,...,Q. If not specified, a uniform
#'                 distribution is used for phi_m[[q]].}
#' \item{phi_l}{a list of numerical vectors. The priors of the lq, q=1,...,Q. If not specified, a uniform
#'                 distribution is used for phi_l[[q]].}
#' \item{phi_l_mean}{a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#'               corresponds to the mean of the Gamma prior of lq.}
#' \item{phi_l_sd}{a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#'                corresponds to the standard deviation of the Gamma prior of lq.}
#' \item{prior_beta}{a character string, which indicates the prior on the
#'                beta_star[q]. The possible values are :
#'                1) "diag" (default) for a diagonal matrix prior,
#'                2) "Ridge_Zellner" for the Ridge Zellner prior}
#' \item{lambda}{a numerical value relative to the Ridge Zellner prior.}
#' \item{display}{a logical value. If FALSE nothing is printed.}
#' }
#' @importFrom stats var
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' res_Bliss_Gibbs_Sampler <- Bliss_Gibbs_Sampler_multiple(data1,param1)
#' K       <- param1$K
#' theta_1 <- res_Bliss_Gibbs_Sampler$trace[1,]
#' theta_1
Bliss_Gibbs_Sampler_multiple <- function(data,param){
  display <- param$display
  if(is.null(display)) display <- TRUE
  # Initialize
  x_mult <- data$x_mult
  y <- data$y
  Q <- length(x_mult)
  # load objects
  iter   <- param[["iter"]]
  grids   <- param$grids
  grids_l <- list()
  for (q in 1:Q){
    grids_l[[q]] <- (grids[[q]] - grids[[q]][1])[-1]
  }
  K      <- param$K
  p      <- param[["p"]]

  # load optional objects
  l_max   <- param[["l_max"]]
  eta_tilde     <- param$eta_tilde
  V_tilde <- param[["V_tilde"]]
  lambda  <- param[["lambda"]]
  g       <- param[["g"]]
  a       <- param$a
  b       <- param[["b"]]
  phi_m   <- param[["phi_m"]]
  phi_l   <- param[["phi_l"]]
  basis   <- param[["basis"]]
  phi_l_mean <- param[["phi_l_mean"]]
  phi_l_sd   <- param[["phi_l_sd"]]
  prior_beta <- param$prior_beta

  # Initialize the necessary unspecified objects
  p <- numeric()
  for(q in 1:Q){
    p[q] <- length(grids[[q]])
  }
  if(is.null(prior_beta))
    stop("Please specify a value for the vector prior_beta.")
  if(is.null(K)) stop("Please specify a value for the vector K.")
  if(!is.null(K)){
    for (q in 1:Q){
      if(is.na(K[q]))
        stop("Please specify a value for all components of the vector K.")
    }
  }
  if(is.null(eta_tilde))     eta_tilde <- rep(0,1+sum(K))
  if(!is.null(eta_tilde)){
    for (q in 1:Q){
      if(is.na(eta_tilde[q]))
        stop("Please specify a value for all components of the vector K.")
    }
  }
  if(is.null(a))       a       <- 1e-1
  if(is.null(b))       b       <- 1e-1
  if(is.null(basis)){
    basis <- character()
    for (q in 1:Q){
      basis[q] <- "uniform"
    }
  }
  if(!is.null(basis)){
    for (q in 1:Q){
      if(is.na(basis[q])){basis[q] <- "uniform"}
    }
  }
  if (is.null(phi_l_mean)) phi_l_mean <- rep(NA,Q)
  if (is.null(phi_l_sd))   phi_l_sd <- rep(NA,Q)
  for(q in 1:Q){
    if(!is.null(phi_l[[q]]) &&
       is.character(phi_l[[q]]) && phi_l[[q]] != "Gamma")
      stop("The qth component of phi_l should be a numeric vector or 'Gamma'.")
    if(!is.null(phi_l[[q]]) &&
       is.character(phi_l[[q]]) && phi_l[[q]] == "Gamma"){

      if(is.na(phi_l_mean[q])) phi_l_mean[q] <-
          diff(range(grids[[q]]))/5 + grids[[q]][1]
      if(is.na(phi_l_sd[q]))   phi_l_sd[q]   <-
          diff(range(grids[[q]]))/5
      phi_l[[q]] <- prior_l(phi_l_mean[q]/K[q],phi_l_sd[q]/K[q],grids_l[[q]])
    }
  }
  if(is.null(phi_l)){
    phi_l <- list()
    for (q in 1:Q){
      if(is.null(l_max)) l_max <- floor(p/5)
      if(!is.null(l_max) & is.na(l_max[q])){l_max[q] <- floor(p[q]/5)}
      phi_l[[q]] <- rep(1/l_max[q],l_max[q])
    }
  }

  for (q in 1:Q){
    l_max[q] <- length(phi_l[[q]])
  }
  if(!is.null(phi_m)){
    for (q in 1:Q){
      if(is.na(phi_m[[q]])){phi_m[[q]] <- rep(1/p[q],p[q])}
    }
  }
  if(is.null(phi_m)){
    phi_m <- list()
    for (q in 1:Q){
      phi_m[[q]] <- rep(1/p[q],p[q])
    }
  }
  if(is.null(prior_beta)) prior_beta <- "diag"
  if(is.character(prior_beta) & prior_beta != "diag" &
     prior_beta != "Ridge_Zellner")
    stop("The prior for beta should be 'diag' or 'Ridge_Zellner'.")
  if (prior_beta == "diag"){
    V_tilde <- diag(1+sum(K))
    V_tilde[1,1] <- 1e2 * var(y)
    indice_q <- 1
    for(q in 1:Q){
      var_x <- min(apply(x_mult[[q]],2,var))
      V_tilde[(indice_q+1):(indice_q+K[q]),(indice_q+1):(indice_q+K[q])] <-
        diag(K[q]) * 1e2 * var(y)/var_x
      indice_q <- indice_q+K[q]
    }
    if(is.null(g))      g <- length(y)
    if(is.null(lambda)) lambda <- 5
  }
  if (prior_beta == "Ridge_Zellner"){
    if(is.null(g))       g <- length(y)
    if(is.null(lambda))  lambda <- 5
    if(is.null(V_tilde)) V_tilde <- diag(1+sum(K))
    #     V_tilde[1,1] <- 1e2 * var(y)
    #     V_tilde[-1,-1] <- diag(sum(K))  * 5

  }

  # Perfome the Gibbs Sampler and return the result.
  res <- Bliss_Gibbs_Sampler_multiple_cpp(Q,y,x_mult,iter,grids,K,l_max,
                                          eta_tilde,a,b,phi_m,phi_l,prior_beta,
                                          g,lambda,V_tilde,
                                          tol=sqrt(.Machine$double.eps),basis)

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
        trace_names <- c(trace_names,paste("b",k," q",q,sep="_"))
      }
      for(k in 1:K[q]){
        trace_names <- c(trace_names,paste("m",k," q",q,sep="_"))
      }
      for(k in 1:K[q]){
        trace_names <- c(trace_names,paste("l",k," q",q,sep="_"))
      }
    }
  }
  colnames(res$trace) <- c(trace_names ,"mu","sigma_sq")
  return(res)
}
