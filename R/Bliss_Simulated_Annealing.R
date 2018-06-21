################################# ----
#' Bliss_Simulated_Annealing
################################# ----
#' @description Simulated Annealing algorithm to determine the Bliss estimate,
#' i.e. the minimum of the posterior expectation of the Bliss loss.
#' @return a list containing :
#' \describe{
#'  \item{Bliss_estimate}{a numerical vector, corresponding to the Bliss estimate
#'  of the coefficient function.}
#'  \item{Smooth_estimate}{a numerical vector, which is the posterior expectation
#'   of beta(t), for each t in the grid ot time points.}
#'  \item{trace}{Each row is an iteration of the Simulated Annealing algorithm.}
#'  \item{argmin}{an integer, which is the index of the iteration minimizing
#'  the loss.}
#' }
#' @param beta a matrix. Each row is a coefficient function computed from the
#' posterior sample.
#' @param param a list containing :
#' \describe{
#' \item{grid}{a numerical vector, the observation time points.}
#' \item{burnin}{an integer, the number of iterations to drop from the Gibbs sample. (optional)}
#' \item{iter}{an integer, the number of iteration of the Gibbs Sampler algorithm.}
#' \item{iter_sann}{an integer, the number of iteration of the Simulated Annealing algorithm. (optional)}
#' \item{Temp_init}{a non negative value, the initial temperature for the cooling function of the Simulated Annealing. (optional)}
#' \item{k_max}{an integer, the maximum number of intervals. (optional)}
#' \item{l_max_sann}{an integer, the maximum value for the parameter l. (optional)}
#' \item{basis}{a character vectors, used to compute the coefficient function, see the function beta_build. (optional)}
#' }
#' @param scale_ml a matrix given by the function "Bliss_Gibbs_Sampler".
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
Bliss_Simulated_Annealing <- function(beta,param,scale_ml){
 display <- param$display
 if(is.null(display)) display <- TRUE
 
 # Initialize
 grid <- param$grid
 iter <- param[["iter"]]
 #p <- param$p         ## PMG - 01/03/18
 p    <- length(grid)  ## PMG - 01/03/18
 
 # load optional objects
 Temp      <- param[["Temp"]] # Commence a mettre du Temp_init partout
 k_max     <- param[["k_max"]]
 iter_sann <- param[["iter_sann"]]
 burnin    <- param[["burnin"]]
 l_max     <- param[["l_max_sann"]] # pas bon : changement de nom
 basis     <- param[["basis"]]
 
 # Initialize the necessary unspecified objects
 if(is.null(Temp))      Temp      <- 1000
 if(is.null(k_max))     k_max     <- 5
 if(is.null(iter_sann)) iter_sann <- 1e5
 if(is.null(burnin))    burnin    <- floor(iter/5)
 if(is.null(l_max))     l_max     <- floor(p/5)
 if(is.null(basis))     basis     <- "uniform"
 
 # Check if the burnin value is correct.
 if(iter <= burnin+1){
  burnin <- floor(iter/5)
  if(display) cat("Burnin value is too large. New burnin value: ",burnin,"\n")
 }
 
 # dm is related to the random walk. A new m_k' is chosen in [m_k - dm , m_k + dm].
 dm <- floor(p/5)+1
 # dl is related to the random walk. A new l_k' is chosen in [l_k - dl , l_k + dl].
 dl <- floor(l_max/2)+1
 
 # Compute the Simulated Annealing algorithm (3 times to find a suitable value
 # for the initial temperature)
 res_sann_list <- list()
 
 # first time
 if(display) cat("First Simulated Annealing. \n")
 res_sann_list[[1]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_functions,
                                                     grid,burnin,Temp,k_max,
                                                     l_max,dm,dl,p,basis,
                                                     scale_ml)
 # Derive a new initial temperature
 Temp <- min(abs(range(res_sann_list[[1]]$trace[,ncol(res_sann_list[[1]]$trace)])
                 - median(res_sann_list[[1]]$trace[,ncol(res_sann_list[[1]]$trace)])))
 # second time
 if(display) cat("Second Simulated Annealing. \n")
 res_sann_list[[2]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_functions,
                                                     grid,burnin,Temp,k_max,
                                                     l_max,dm,dl,p,basis,
                                                     scale_ml)
 # Derive a new initial temperature
 Temp <- min(abs(range(res_sann_list[[2]]$trace[,ncol(res_sann_list[[2]]$trace)])
                 - median(res_sann_list[[2]]$trace[,ncol(res_sann_list[[2]]$trace)])))
 # third time
 if(display) cat("Third Simulated Annealing. \n")
 res_sann_list[[3]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_functions,
                                                     grid,burnin,Temp,k_max,
                                                     l_max,dm,dl,p,basis,
                                                     scale_ml)
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
 res_k  <- argmin[length(argmin)-1]
 
 # Compute the estimate
 b <- argmin[1:res_k]
 m         <- argmin[1:res_k+  k_max]
 l         <- argmin[1:res_k+2*k_max]
 k         <- argmin[3*k_max+2]
 estimate  <- beta_build_cpp(b,m,l,grid,p,k,basis,scale_ml)
 
 return(list(Bliss_estimate = estimate,
             posterior_expe = res_sann$posterior_expe,
             posterior_var  = res_sann$posterior_var,
             trace          = res_sann$trace,
             argmin         = argmin))
}
