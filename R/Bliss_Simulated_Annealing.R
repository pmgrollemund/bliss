################################# ----
#' Bliss_Simulated_Annealing
################################# ----
#' @description A Simulated Annealing algorithm to determine the Bliss estimate
#'              (for only functional covariate), i.e. the minimum of the 
#'              posterior expectation of the Bliss loss.
#' @return a list containing:
#' \describe{
#'  \item{Bliss_estimate}{a numerical vector, corresponding to the Bliss estimate
#'        of the coefficient function.}
#'  \item{Smooth_estimate}{a numerical vector, which is the posterior expectation
#'        of beta(t), for each t in the grid ot time points.}
#'  \item{trace}{Each row is an iteration of the Simulated Annealing algorithm.}
#'  \item{argmin}{an integer, which is the index of the iteration minimizing
#'        the Bliss loss.}
#' }
#' @param beta_sample a matrix. Each row is a coefficient function computed from the
#'        posterior sample.
#' @param normalization_values a matrix given by the function "Bliss_Gibbs_Sampler". XXXXXXX
#' @param param a list containing:
#' \describe{
#' \item{grid}{a numerical vector, the observation time points.}
#' \item{burnin}{an integer, the number of iterations to drop from the Gibbs
#'       sample. (optional)}
#' \item{iter}{an integer, the number of iteration of the Gibbs Sampler algorithm.}
#' \item{iter_sann}{an integer, the number of iteration of the Simulated
#'       Annealing algorithm. (optional)} XXXXXXX
#' \item{Temp_init}{a non negative value, the initial temperature for the
#'       cooling function of the Simulated Annealing. (optional)}
#' \item{k_max}{an integer, the maximum number of intervals. (optional)} XXXXXXX
#' \item{l_max_sann}{an integer, the maximum value for the parameter l. (optional)} XXXXXXX
#' \item{basis}{a character vectors, used to compute the coefficient function,
#'       see the function beta_build. (optional)} XXXXXXX
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for
#'       each covariate.}
#' }
#' @param progress a logical value. If TRUE, the algorithm progress is displayed.
#'         (optional)
#' @importFrom stats median
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' #res_Simulated_Annealing <- Bliss_Simulated_Annealing(beta_functions,param1)
#' #ylim <- range(c(res_Simulated_Annealing$Bliss_estimate,
#'  #                res_Simulated_Annealing$posterior_expe))
#' #plot(param$grid,res_Simulated_Annealing$Bliss_estimate,type="l",ylim=ylim)
#' #lines(param$grid,res_Simulated_Annealing$posterior_expe,lty=2)
Bliss_Simulated_Annealing <- function(beta_sample,normalization_values,param,progress=FALSE){
 # load optional objects
 grid <- param[["grid"]]
 iter <- param[["iter"]]
 Temp_init <- param[["Temp_init"]]
 K         <- param[["K"]]
 k_max     <- param[["k_max"]]
 iter_sann <- param[["iter_sann"]]
 times_sann<- param[["times_sann"]]
 burnin    <- param[["burnin"]]
 l_max     <- param[["l_max"]] # pas bon : changement de nom
 basis     <- param[["basis"]]
 p    <- length(grid)
 
 # Initialize the necessary unspecified objects
 if(is.null(Temp_init)) Temp_init <- 1000
 if(is.null(k_max))     k_max     <- K  # PMG 22/06/18 
 if(is.null(iter_sann)) iter_sann <- 5e4
 if(is.null(times_sann))times_sann<- 100
 if(is.null(burnin))    burnin    <- floor(iter/5) 
 if(is.null(l_max))     l_max     <- floor(p/5) # XXX a changer ?
 if(is.null(basis))     basis     <- "Uniform"
 
 # Check if the burnin isn't too large.
 if(2*burnin > iter){
  burnin <- floor(iter/5)
  if(progress) 
   cat("\t Burnin is too large. New burnin : ",burnin,".\n")
 }
 
 # dm and dl are related to the random walk.
 dm <- floor(p/5)+1
 dl <- floor(l_max/2)+1
 
 # Compute the starting point
 starting_point = compute_starting_point(apply(beta_sample,2,mean))
 if(k_max < nrow(starting_point)) k_max <- nrow(starting_point)  # PMG 04/08/18
  
 # Compute the Simulated Annealing algorithm (3 times to find a suitable value
 # for the initial temperature)
 res_sann <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_sample,
                                           grid,burnin,Temp_init,k_max,
                                           l_max,dm,dl,p,basis,
                                           normalization_values,
                                           progress,starting_point)
 min_loss <- min(res_sann$trace[,3*k_max+4])
 for(times in 1:times_sann ){
  if(progress) 
   cat("\t ",times,".\n")
  res_sann_tmp <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_sample,
                                                grid,burnin,Temp_init,k_max,
                                                l_max,dm,dl,p,basis,
                                                normalization_values,
                                                progress,starting_point)
  min_loss_tmp <- min(res_sann_tmp$trace[,3*k_max+4])
  if(min_loss_tmp < min_loss){
   res_sann <- res_sann_tmp
   min_loss <- min_loss_tmp
  }
 }
 
 # Determine the estimate
 index <- which(res_sann$trace[,ncol(res_sann$trace)] %in%
                 min(res_sann$trace[,ncol(res_sann$trace)]))[1]
 
 argmin <- res_sann$trace[index,]
 res_k  <- argmin[length(argmin)-1] # XXXXX selection avec nom "K"
 
 # Compute the estimate
 b        <- argmin[1:res_k]
 m        <- argmin[1:res_k+  k_max]
 l        <- argmin[1:res_k+2*k_max]
 k        <- argmin[3*k_max+3]
 estimate <- compute_beta_cpp(b,m,l,grid,p,k,basis,normalization_values)
 
 difference  = abs(estimate-res_sann$posterior_expe)
 sdifference = moving_average_cpp(difference,floor(p/10))
 
 trace_names <- NULL
 for(k in 1:k_max){
  trace_names <- c(trace_names,paste("b",k,sep="_"))
 }
 for(k in 1:k_max){
  trace_names <- c(trace_names,paste("m",k,sep="_"))
 }
 for(k in 1:k_max){
  trace_names <- c(trace_names,paste("l",k,sep="_"))
 }
 colnames(res_sann$trace) <- c(trace_names ,"accepted","choice","k","loss")
 
 return(list(Bliss_estimate  = estimate,
             Smooth_estimate = res_sann$posterior_expe,
             trace           = res_sann$trace,
             argmin          = argmin,
             difference      = difference,
             sdifference     = sdifference))
}
