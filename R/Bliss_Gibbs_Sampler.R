################################# ----
#' Bliss_Gibbs_Sampler
################################# ----
#' @description A Gibbs Sampler algorithm to sample posterior distribution of the Bliss model.
#' @return a list containing :
#' \describe{
#'  \item{trace}{a matrix. Each row is a draw from the posterior.}
#'  \item{param}{a list containing K and scale_ml}
#' }
#' @param data a list containing
#' \describe{
#'  \item{x}{a list containing the functions x_qi(t) observed at time points given by grids,}
#'  \item{y}{the outcome values y_i,}
#'  \item{grids}{a list containing the observation time points of the covariates functions.}
#' }
#' @param param a list containing :
#' \describe{
#' \item{iter}{an integer, the number of iterations of the Gibbs sampler algorithm.}
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for each covariate.}
#' \item{basis}{a vector of characters among : "uniform" (default),
#'                 "epanechnikov", "gauss" and "triangular" which correspond to
#'                 different basis expansion of coefficient function and functional data (optional)}
#' \item{g}{hyperparameter of the Bliss model,  a nonnegative value, the coefficient of the Zellner prior.}
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
Bliss_Gibbs_Sampler <- function(data,param){
 display <- param$display
 if(is.null(display)) display <- TRUE
 # Initialize
 x <- data$x
 y <- data$y
 Q <- length(x)
 # load objects
 basis <- param[["basis"]]
 grids <- param[["grids"]]
 iter  <- param[["iter"]]
 K     <- param[["K"]]
 g     <- param[["g"]]
 p     <- param[["p"]]
 
 # Initialize the necessary unspecified objects
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
   basis[q] <- "uniform"
  }
 }
 if(!is.null(basis)){
  for (q in 1:Q){
   if(is.na(basis[q])){basis[q] <- "uniform"}
  }
 }
 
 grids_l <- list()
 for (q in 1:Q){
  grids_l[[q]] <- (grids[[q]] - grids[[q]][1])[-1] # Utilite ? a diviser par max aussi ?
 }
 
 
 #### Passage a bien changer !
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
 ######
 
 if(is.null(g))       g <- length(y)
 if(is.null(lambda))  lambda <- 5
 if(is.null(V_tilde)) V_tilde <- diag(1+sum(K))
 
 
 # Perfome the Gibbs Sampler and return the result.
 res <- Bliss_Gibbs_Sampler_cpp(Q,y,x_mult,iter,grids,K,l_max,
                                eta_tilde,a,b,phi_m,phi_l,prior_beta,
                                g,lambda,V_tilde,
                                tol=sqrt(.Machine$double.eps),basis)
 # option a changer ?
 
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
