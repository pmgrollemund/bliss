################################# ----
#' Bliss_Simulated_Annealing
################################# ----
#' @description A Simulated Annealing algorithm to compute the Bliss estimate.
#' @return a list containing:
#' \describe{
#'  \item{Bliss_estimate}{a numerical vector, corresponding to the Bliss estimate
#'        of the coefficient function.}
#'  \item{Smooth_estimate}{a numerical vector, which is the posterior expectation
#'        of the coefficient function for each time points.}
#'  \item{trace}{a matrix, the trace of the algorithm.}
#' }
#' @param beta_sample a matrix. Each row is a coefficient function computed from the
#'        posterior sample.
#' @param posterior_sample a list resulting from the \code{Bliss_Gibbs_Sampler}
#'        function.
#' @param param a list containing:
#' \describe{
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of
#'        time points for the qth functional covariate.}
#' \item{basis}{a character (optional). The possible values are "uniform"
#'       (default), "epanechnikov", "gauss" and "triangular" which correspond to
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
#' \item{Temp_init}{a nonnegative value (optional), the initial temperature for
#'      the cooling function of the Simulated Annealing algorithm.}
#' \item{Q}{an integer, the number of functional covariates.}
#' \item{p}{a vector of integers, the numbers of time point of each functional
#'      covariate.}
#' \item{verbose}{write stuff if TRUE (optional).}
#' }
#' @param verbose_cpp Rcpp writes stuff if TRUE (optional).
#' @importFrom stats median
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#' param1$Q <- length(data1$x)
#' param1$grids <- data1$grids
#' param1$p <- sapply(data1$grids,length)
#'
#' posterior_sample <- res_bliss1$posterior_sample
#' beta_sample <- compute_beta_sample(posterior_sample,param1)
#'
#' res_sann <- Bliss_Simulated_Annealing(beta_sample,posterior_sample,param1)
#' }
Bliss_Simulated_Annealing <- function(beta_sample,posterior_sample,param,
                                      verbose_cpp=FALSE){

  ###### Initialisation - Load objects
  Temp_init <- param[["Temp_init"]]
  k_maxs    <- param[["k_max"]]
  iter_sann <- param[["iter_sann"]]
  times_sann<- param[["times_sann"]]
  burnin    <- param[["burnin"]]
  l_maxs    <- param[["l_max"]]
  basis     <- param[["basis"]]
  ps     <- param[["p"]]
  grids  <- param[["grids"]]
  Q      <- param[["Q"]]
  verbose      <- param[["verbose"]]

  ###### Initialisation - Define values
  if(is.null(Temp_init)) Temp_init <- 1000
  if(length(k_maxs) == 1)k_maxs    <- rep(k_maxs,Q)  # PMG 24/04/26
  if(is.null(k_maxs))    k_maxs    <- param[["K"]]  # PMG 22/06/18
  if(is.null(iter_sann)) iter_sann <- 5e4
  if(is.null(times_sann))times_sann<- 5 #50 # PMG 24/04/26
  if(is.null(burnin))    burnin    <- ceiling(nrow(posterior_sample$trace) / 10) #floor(iter_sann/5) # PMG 24/05/02
  if(is.null(l_maxs))    l_maxs    <- floor(ps/5)
  if(is.null(basis))     basis     <- "Uniform"
  if(is.null(verbose))   verbose   <- FALSE
  if(is.null(param[["sann_trace"]])){
    sann_trace <- FALSE
  }else{
    sann_trace <- param[["sann_trace"]]
  }

  ###### Alerte
  # Check if the burnin isn't too large.
  iter <- nrow(beta_sample[[1]]) -1 # PMG 11/11/18
  if(2*burnin > iter){
    burnin <- floor(iter/5)
    cat("\t SANN -- Burnin is too large. New burnin : ",burnin,".\n")
  }

  ###### Initialisation - Output object
  res_bliss_sann <- list()
  res_bliss_sann$Bliss_estimate <- list() ; length(res_bliss_sann$Bliss_estimate) <- Q
  res_bliss_sann$trace <- list() ; length(res_bliss_sann$trace) <- Q
  res_bliss_sann$Smooth_estimate <- list() ; length(res_bliss_sann$Smooth_estimate) <- Q

  ###### Run the Simulated Annealing algorithm
  # for each functional covariates
  for(q in 1:Q){
    if(verbose){
      message_to_print <- paste("\r \t Functional Dimension ",q,"/",Q,sep="")
      message(message_to_print, appendLF = FALSE)
    }
    # get values related to the qth functional covariates
    grid <- grids[[q]]
    p <- ps[[q]]
    l_max <- l_maxs[[q]]
    k_max <- k_maxs[[q]]
    normalization_values <- posterior_sample$param$normalization_values[[q]]
    dm <- floor(p/5)+1
    dl <- floor(l_max/2)+1

    # beta sample after burnin
    beta_sample_tmp <- beta_sample[[q]][-(1:burnin),] # PMG 2024-04-29

    # Compute the starting point
    starting_point = compute_starting_point_sann(apply(beta_sample_tmp,2,mean))
    if(k_max < nrow(starting_point)) k_max <- nrow(starting_point)  # PMG 04/08/18

    # Compute the Simulated Annealing algorithm
    #   (3 times to find a suitable value for the initial temperature)
    if(verbose){
      message_to_print_2 <- c(message_to_print,paste("  Repeat 1/",times_sann,sep=""))
      message(message_to_print_2, appendLF = FALSE)
    }
    res_sann <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_sample_tmp,
                                              grid,Temp_init,k_max,
                                              l_max,dm,dl,p,basis,
                                              normalization_values,
                                              verbose_cpp,starting_point)

    # Determine the argmin
    min_loss <- min(res_sann$trace[,3*k_max+4])

    # Repeat
    if(times_sann >= 2){
      for(times in 2:times_sann ){
        if(verbose){
          message_to_print_2 <- c(message_to_print,paste("  Repeat ",times,"/",times_sann,sep=""))
          message(message_to_print_2, appendLF = FALSE)
        }
        res_sann_tmp <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_sample_tmp,
                                                      grid,Temp_init,k_max,
                                                      l_max,dm,dl,p,basis,
                                                      normalization_values,
                                                      verbose_cpp,starting_point)

        # Determine the argmin (and check if this is better of the first ones)
        min_loss_tmp <- min(res_sann_tmp$trace[,3*k_max+4])
        if(min_loss_tmp < min_loss){
          res_sann <- res_sann_tmp
          min_loss <- min_loss_tmp
        }
      }
    }

    # Determine the estimate
    index <- which(res_sann$trace[,ncol(res_sann$trace)] %in%
                     min(res_sann$trace[,ncol(res_sann$trace)]))[1]

    argmin <- res_sann$trace[index,]
    res_k  <- argmin[length(argmin)-1]

    # Compute the estimate
    b        <- argmin[1:res_k]
    m        <- argmin[1:res_k+  k_max]
    l        <- argmin[1:res_k+2*k_max]
    k        <- argmin[3*k_max+3]
    estimate <- compute_beta_cpp(b,m,l,grid,p,k,basis,normalization_values)

    # Output of the qth functional covariate
    if(sann_trace){
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

      res_bliss_sann$trace[[q]] <- res_sann$trace
      colnames(res_bliss_sann$trace[[q]]) <- c(trace_names ,"accepted","choice","k","loss")
    }

    res_bliss_sann$Bliss_estimate[[q]] <- estimate
    res_bliss_sann$Smooth_estimate[[q]] <- res_sann$posterior_expe
  }
  if(verbose) cat("\n")

  # Output
  return(res_bliss_sann)
}
