################################# ----
#' Bliss_Simulated_Annealing
################################# ----
#' @description A Simulated Annealing algorithm to determine the Bliss estimate,
#'              i.e. the minimum of the posterior expectation of the Bliss loss.
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
#' @param scale_ml a matrix given by the function "Bliss_Gibbs_Sampler". XXXXXXX
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
#' #res.Simulated_Annealing <- Bliss_Simulated_Annealing(beta_functions,param1)
#' #ylim <- range(c(res.Simulated_Annealing$Bliss_estimate,
#'  #                res.Simulated_Annealing$posterior_expe))
#' #plot(param$grid,res.Simulated_Annealing$Bliss_estimate,type="l",ylim=ylim)
#' #lines(param$grid,res.Simulated_Annealing$posterior_expe,lty=2)
Bliss_Simulated_Annealing <- function(beta_sample,scale_ml,param,progress=FALSE){
 # load optional objects
 grid <- param[["grid"]]
 iter <- param[["iter"]]
 p    <- length(grid)
 Temp_init <- param[["Temp_init"]]
 k_max     <- param[["k_max"]]
 iter_sann <- param[["iter_sann"]]
 burnin    <- param[["burnin"]]
 l_max     <- param[["l_max_sann"]] # pas bon : changement de nom
 basis     <- param[["basis"]]
 
 # Initialize the necessary unspecified objects
 if(is.null(Temp_init)) Temp_init <- 1000
 # if(is.null(k_max))     k_max     <- 5
 if(is.null(k_max))     k_max     <- K  # PMG 22/06/18
 if(is.null(iter_sann)) iter_sann <- 1e5
 if(is.null(burnin))    burnin    <- floor(iter/5) # utile ? XXXXXXX
 if(is.null(l_max))     l_max     <- floor(p/5)
 if(is.null(basis))     basis     <- "uniform"
 
 # Check if the burnin value is correct.
 if(iter <= burnin+1){
  burnin <- floor(iter/5)
  if(progress) cat("Burnin value is too large. New burnin value: ",burnin,"\n")
 }
 
 # dm and dl are related to the random walk.
 dm <- floor(p/5)+1
 dl <- floor(l_max/2)+1
 
 # Compute the Simulated Annealing algorithm (3 times to find a suitable value
 # for the initial temperature)
 res_sann_list <- list()
 
 # first time
 if(progress) cat("First Simulated Annealing. \n")
 res_sann_list[[1]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_sample,
                                                     grid,burnin,Temp_init,k_max,
                                                     l_max,dm,dl,p,basis,
                                                     scale_ml,
                                                     progress)
 # Derive a new initial temperature
 Temp_init <- min(abs(range(res_sann_list[[1]]$trace[,ncol(res_sann_list[[1]]$trace)])
                      - median(res_sann_list[[1]]$trace[,ncol(res_sann_list[[1]]$trace)])))
 # second time
 if(progress) cat("Second Simulated Annealing. \n")
 res_sann_list[[2]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_sample,
                                                     grid,burnin,Temp_init,k_max,
                                                     l_max,dm,dl,p,basis,
                                                     scale_ml,
                                                     progress)
 # Derive a new initial temperature
 Temp_init <- min(abs(range(res_sann_list[[2]]$trace[,ncol(res_sann_list[[2]]$trace)])
                      - median(res_sann_list[[2]]$trace[,ncol(res_sann_list[[2]]$trace)])))
 # third time
 if(progress) cat("Third Simulated Annealing. \n")
 res_sann_list[[3]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_sample,
                                                     grid,burnin,Temp_init,k_max,
                                                     l_max,dm,dl,p,basis,
                                                     scale_ml,
                                                     progress)
 # Comparison and selection
 mins      <- c(min(res_sann_list[[1]]$trace[,ncol(res_sann_list[[1]]$trace)]),
                min(res_sann_list[[2]]$trace[,ncol(res_sann_list[[2]]$trace)]),
                min(res_sann_list[[3]]$trace[,ncol(res_sann_list[[3]]$trace)]))
 index <- which(mins==min(mins))
 res_sann <- res_sann_list[[index]]
 
 # Determine the estimate
 index <- which(res_sann$trace[,ncol(res_sann$trace)] %in%
                 min(res_sann$trace[,ncol(res_sann$trace)]))[1]
 
 argmin <- res_sann$trace[index,]
 res_k  <- argmin[length(argmin)-1] # XXXXX selection avec nom "K"
 
 # Compute the estimate
 b        <- argmin[1:res_k]
 m        <- argmin[1:res_k+  k_max]
 l        <- argmin[1:res_k+2*k_max]
 k        <- argmin[3*k_max+2]
 estimate <- beta_build_cpp(b,m,l,grid,p,k,basis,scale_ml)
 
 return(list(Bliss_estimate  = estimate,
             Smooth_estimate = res_sann$posterior_expe,
             trace           = res_sann$trace,
             argmin          = argmin))
}
