################################# ----
#' BIC_model_choice
################################# ----
#' @description Model selection with BIC criterion.
#' @return A numerical vector, the BIC values for the Bliss model for different
#'         K value.
#' @param Ks a numerical vector containing the K values.
#' @param iter an integer, the number of iteration for each run of \code{fit_Bliss}.
#' @param data a list containing required options to run the function
#'        \code{fit_Bliss}.
#' @param verbose write stuff if TRUE (optional).
#' @export
#' @examples
#' \donttest{
#' param_sim <- list(Q=1,n=100,p=c(50),grids_lim=list(c(0,1)))
#' data      <- sim(param_sim,verbose=TRUE)
#' iter = 1e2
#' Ks <- 1:5
#'
#' res_BIC <- BIC_model_choice(Ks,iter,data)
#' plot(res_BIC,xlab="K",ylab="BIC")
#' }
BIC_model_choice <- function(Ks,iter,data,verbose=T){
  BIC <- rep(0,length(Ks))
  for(i in 1:length(Ks)){
    if(verbose) cat("K = ", Ks[i], "\n",sep="")
    param_BIC <- list(iter=iter,K=Ks[i])

    res_bliss <- fit_Bliss(data,param_BIC,verbose=F,sann=F,compute_density=F,
                           support_estimate=F)

    llkh <- dposterior(res_bliss$posterior_sample,data)
    llkh <- llkh[,'log likelihood']
    lML <- llkh[which.max(llkh)]
    BIC[i] <- (3*Ks[i]+2) * log(length(data$y)) - 2 * lML
  }
  return(BIC)
}

################################# ----
#' post_treatment_bliss
################################# ----
#' @description Compute the post treatment values.
#' @return A list of important post treatment value: BIC, the maximum of the log
#' likelihood and the numbre of parameters.
#' @param posterior_sample a list provided by the function \code{Bliss_Gibbs_Sampler}.
#' @param param a list containing:
#' \describe{
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for
#'       each covariate.}
#' }
#' @param data a list containing required options to run the function
#'        \code{dposterior}.
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#'
#' post_treatment_bliss(res_bliss1$posterior_sample,param1,data1)
#' }
post_treatment_bliss <- function(posterior_sample,param,data){
  ###### Initialisation
  K <- param[["K"]]

  ###### Compute posterior values
  llkh <- dposterior(posterior_sample,data)
  llkh <- llkh[,'log likelihood']
  lML <- llkh[which.max(llkh)][1]
  nb_param <- (3*K+2)
  BIC <- nb_param * log(length(data$y)) - 2 * lML

  ###### Output
  res <- list(BIC=BIC,loglik=lML,nb_param=nb_param)
  return(res)
}

################################# ----
#' predict_bliss_distribution
################################# ----
#' @description Compute the distribution of the predictions.
#' @return A matrix containing predictions for each individual data \code{x}.
#' @param x a list containing the design matrices related to the functional
#'        covariates. Must be similar to the result of the function \code{sim_x}.
#' @param grids a list of numerical vectors, the qth vector is the grid of
#'        time points for the qth functional covariate.
#' @param burnin an integer (optional), the number of iteration to drop from the
#'       posterior sample.
#' @param posterior_sample a list provided by the function \code{Bliss_Gibbs_Sampler}.
#' @param beta_sample a list provided by the function \code{compute_beta_sample}.
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#'
#' predict_bliss_distribution(data1$x,data1$grids,50,res_bliss1$posterior_sample,
#'    res_bliss1$beta_sample)
#' }
predict_bliss_distribution <- function(x,grids,burnin,posterior_sample,beta_sample){
  ###### Initialisation
  Q <- length(x)
  n <- nrow(x[[1]])
  n_iter <- nrow(posterior_sample$trace)
  if(is.null(burnin)) burnin <- ceiling(n_iter / 10)
  n_iter <- n_iter - burnin
  y_fitted <- matrix(NA,nrow = n,ncol = n_iter)

  ###### Pretreatment
  for(q in 1:Q){
    x[[q]] <- scale(x[[q]],scale=F)
  }

  ###### Compute fitted values

  for(i in 1:n){
    for(j in 1:n_iter){
      y_tmp <- posterior_sample$trace$mu[j+burnin]

      for(q in 1:Q){
        y_tmp <- y_tmp +
          integrate_trapeze_cpp(grids[[q]],beta_sample[[q]][j+burnin,] * x[[q]][i,])
      }
      y_fitted[i,j] <- y_tmp
    }
  }

  ###### Output
  return(y_fitted)
}
################################# ----
#' predict_bliss
################################# ----
#' @description Compute predictions.
#' @return A vector of predictions for each individual data \code{x}.
#' @param x a list containing the design matrices related to the functional
#'        covariates. Must be similar to the result of the function \code{sim_x}.
#' @param grids a list of numerical vectors, the qth vector is the grid of
#'        time points for the qth functional covariate.
#' @param burnin an integer (optional), the number of iteration to drop from the
#'       posterior sample.
#' @param posterior_sample a list provided by the function \code{Bliss_Gibbs_Sampler}.
#' @param Smooth_estimate one of the objects resulting from \code{Bliss_Simulated_Annealing}.
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#'
#' predict_bliss(data1$x,data1$grids,50,res_bliss1$posterior_sample,res_bliss1$smooth_estimate)
#' }
predict_bliss <- function(x,grids,burnin,posterior_sample,Smooth_estimate){
  ###### Initialisation
  Q <- length(x)
  n <- nrow(x[[1]])
  y_fitted <- rep(NA,n)
  n_iter <- nrow(posterior_sample$trace)
  if(is.null(burnin)) burnin <- ceiling(n_iter / 10)

  ###### Pretreatment
  for(q in 1:Q){
    x[[q]] <- scale(x[[q]],scale=F)
  }

  ###### Compute fitted values
  mu_hat <- mean(posterior_sample$trace$mu[-(1:burnin)])

  for(i in 1:n){
    y_tmp <- mu_hat
    for(q in 1:Q){
      y_tmp <- y_tmp +
        integrate_trapeze_cpp(grids[[q]],Smooth_estimate[[q]] * x[[q]][i,])
    }
    y_fitted[i] <- y_tmp
  }

  ###### Output
  return(y_fitted)
}


################################# ----
#' compute_beta_sample
################################# ----
#' @description Compute the posterior coefficient function from the posterior
#'              sample.
#' @return a matrix containing the coefficient function posterior sample.
#' @param posterior_sample a list provided by the function \code{Bliss_Gibbs_Sampler}.
#' @param param a list containing:
#' \describe{
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for
#'       each covariate.}
#' \item{grids}{a numerical vector, the observation time points.}
#' \item{basis}{a character (optional) among : "uniform" (default),
#'       "epanechnikov", "gauss" and "triangular" which correspond to
#'       different basis functions to expand the coefficient function and the
#'       functional covariates.}
#' \item{Q}{an integer, the number of functional covariates.}
#' \item{p}{a vector of integers, the numbers of time points of each functional
#'      covariate.}
#' }
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#' param1$grids <- data1$grids
#' param1$p <- sapply(data1$grids,length)
#' param1$Q <- length(data1$x)
#' beta_sample <- compute_beta_sample(posterior_sample=res_bliss1$posterior_sample,
#'                                    param=param1)
compute_beta_sample <- function(posterior_sample,param){
  ###### Initialisation - Load objects
  Ks     <- param[["K"]]
  grids  <- param[["grids"]]
  basis  <- param[["basis"]]
  Q      <- param[["Q"]]
  ps     <- param[["p"]]

  ####### Initialization - Define values
  if(is.null(basis)) basis <- rep("Uniform",length(Ks))

  ####### Compute the coefficient function
  # for each iteration of the Gibbs Sampler and for each covariable
  beta_sample <- list()
  count <- 0
  for(q in 1:Q){
    K     <- Ks[q]
    grid  <- grids[[q]]
    p     <- ps[q]
    trace_tmp <- posterior_sample$trace[,(1+count):(count+3*K)]
    trace_tmp <- as.matrix(trace_tmp)
    norm_val <- posterior_sample$param$normalization_values[[q]]

    beta_sample[[q]] <- compute_beta_sample_cpp(trace_tmp,
                                                K,grid,p,basis[q],norm_val)
    count <- count + 3*K[q]
  }

  ####### Output
  return(beta_sample)
}

################################# ----
#' compute_beta_posterior_density
################################# ----
#' @description Compute the posterior density of the coefficient function.
#' @details The posterior densities correponds to approximations of the marginal
#'          posterior distribitions (of beta(t) for each t).
#'          The sample is thinned in order to reduce the correlation and the
#'          computational time of the function \code{\link[=kde2d]{kde2d}}.
#' @return An approximation of the posterior density on a two-dimensional grid
#'         (corresponds to the result of the \code{\link[=kde2d]{kde2d}} function).
#' @param beta_sample a matrix. Each row is a coefficient function computed from the
#'        posterior sample.
#' @param param a list containing:
#' \describe{
#' \item{grid}{a numerical vector, the time points.}
#' \item{lims_estimate}{a numerical vector, the time points.}
#' \item{burnin}{an integer (optional), the number of iteration to drop from the Gibbs
#'       sample.}
#' \item{lims_kde}{an integer (optional), correspond to the \code{lims} option
#'       of the \code{kde2d} funtion.}
#' \item{new_grid}{a numerical vector (optional) to compute beta sample on a
#'       different grid.}
#' \item{thin}{an integer (optional) to thin the posterior sample.}
#' }
#' @importFrom MASS bandwidth.nrd kde2d
#' @export
#' @examples
#' \donttest{
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#' param1$grids <- data1$grids
#' param1$p <- sapply(data1$grids,length)
#' param1$Q <- length(data1$x)
#'
#' density_estimate <- compute_beta_posterior_density(res_bliss1$beta_sample,param1)
#' }
compute_beta_posterior_density <- function(beta_sample,param){

  ####### Initialization - Load objects
  grids <- param[["grids"]]
  thin     <- param[["thin"]]
  burnin   <- param[["burnin"]]
  lims_kdes <- param[["lims_kdes"]]
  new_grids <- param[["new_grids"]]
  ps <- param[["p"]]
  n <- param[["n"]]
  Q <- param[["Q"]]

  ####### Initialization - Define values
  N <- 512          # PMG 11/11/18

  iter <- nrow(beta_sample[[1]]) -1 # PMG 11/11/18
  max_points <- 1e5
  if(is.null(burnin)) burnin <- floor(iter/5)
  if(is.null(thin))   thin   <- floor((iter-burnin)*min(ps)/max_points)


  # Compute the coefficient function on the new grid (if required).
  if(!is.null(new_grids)){
    old_beta_sample <- beta_sample
    for(q in 1:Q){
      beta_sample[[q]] <- matrix(0,nrow(beta_sample[[q]]),ps[q])
      for(i in 1:nrow(beta_sample[[q]])){
        beta_sample[[q]][i,] <- change_grid(old_beta_sample[[q]][i,],grids[[q]],
                                            new_grids[[q]])
      }
    }
    param$old_grids <- grids
    param$grids     <- new_grids
    param$new_grids <- NULL # PMG 22/06/18
    grids           <- new_grids
  }

  # Output object
  res <- list() ; length(res) <- Q

  ####### Pretreatment
  if(2*burnin > iter){
    burnin <- floor(iter/5)
  }


  ####### Compute the posterior density for each functional covariates
  for(q in 1:Q){
    # Get objet related to the qth functional covariate
    lims_estimate <- range(apply(beta_sample[[q]],2,mean)) # PMG 11/11/18

    # Thin the posterior sample
    thin_min   <- max(1,floor((iter-burnin)*ps[q]/max_points))
    if(thin <  thin_min){
      thin <- thin_min
    }
    beta_sample[[q]] <- beta_sample[[q]][seq(1+burnin,iter,by=thin),]

    # Perform the kde2d function
    beta_x <- rep(grids[[q]],nrow(beta_sample[[q]]))
    beta_y <- as.vector(t(beta_sample[[q]]))

    h1 <- MASS::bandwidth.nrd(beta_x)
    h2 <- MASS::bandwidth.nrd(beta_y)
    if(h2 == 0){
      h2 <- 4 * 1.06 * sd(beta_y) * length(beta_y)^(-1/5)
    }
    points <- cbind(beta_x,beta_y)

    lims_kde <- c(range(beta_x), stats::quantile(beta_y,c(0.025,0.975))) # PMG 04/07/18
    if(lims_kde[3] >= 0 ) lims_kde[3] <- -h2/2 # PMG 04/07/18
    if(lims_kde[4] <= 0 ) lims_kde[4] <- -h2/2 # PMG 04/07/18
    if(lims_kde[3] >= lims_estimate[1] ) lims_kde[3] <- lims_estimate[1]-h2/2 # PMG 04/07/18
    if(lims_kde[4] <= lims_estimate[2] ) lims_kde[4] <- lims_estimate[2]+h2/2 # PMG 04/07/18
    res_kde2d <- MASS::kde2d(x=beta_x,y=beta_y,lims=lims_kde,
                       n=N,h=c(h1,h2))

    # Update output result
    res[[q]]$grid_t <- res_kde2d$x
    res[[q]]$grid_beta_t <- res_kde2d$y
    res[[q]]$density <- res_kde2d$z
  }

  ####### Output
  # update beta_sample
  if(!is.null(param$old_grids)){
    res$new_beta_sample <- beta_sample
  }else{
    res$new_beta_sample <- NULL
  }

  return(res)
}

################################# ----
#' between
################################# ----
#' @description Check if a number belong to a given interval.
#' @return a logical value.
#' @param value a numerical value.
#' @param interval a numerical vector: (lower,upper).
#' @export
#' @examples
#' 1 %between% c(0,2)
#' 2 %between% c(0,2)
#' 3 %between% c(0,2)
"%between%" <- function(value,interval){
 (value >= interval[1]) & (value <= interval[2])
}

################################# ----
#' support_estimation
################################# ----
#' @description Compute the support estimate.
#' @return a list containing:
#' \describe{
#'  \item{alpha}{a numerical vector. The approximated posterior probabilities
#'        that the coefficient function support covers \code{t} for each time
#'        points \code{t}.}
#'  \item{estimate}{a numerical vector, the support estimate.}
#'  \item{estimate_fct}{a numerical vector, another version of the support
#'        estimate.}
#' }
#' @param beta_sample the result of the function \code{compute_beta_sample}.
#' @param param a list containing the value \code{Q} and an optional
#'        parameter \code{gamma}.
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#' param1$Q <- length(data1$x)
#'
#' res_support <- support_estimation(res_bliss1$beta_sample,param1)
support_estimation <- function(beta_sample,param){

  ####### Initialization - Load objects
  Q     <- param[["Q"]]
  gamma <- param[["gamma"]]

  ####### Initialization - Define values
  if(is.null(gamma)) gamma <- 0.5

  res_support <- list()
  res_support$estimate <- list() ; length(res_support$estimate) <- Q
  res_support$support_estimate_fct <- list() ; length(res_support$support_estimate_fct) <- Q
  res_support$alpha <- list() ; length(res_support$alpha) <- Q

  ####### Run the support estimation
  # for each functional covariates
  for(q in 1:Q){

    # alpha: posterior probabilities
    alpha <- apply(beta_sample[[q]],2, function(vec) sum(vec != 0)/length(vec))

    # Support estimate
    tmp   <- rep(0,ncol(beta_sample[[q]]))
    tmp2  <- which(alpha >= gamma)
    tmp[tmp2] <- 1
    estimate <- determine_intervals(tmp)
    estimate <- estimate[estimate$value != 0,]
    estimate$value <- NULL

    # Support estimate (vectorial version)
    estimate_fct <- rep(0,ncol(beta_sample[[q]]))
    if(nrow(estimate) > 0){
      for(i in 1:nrow(estimate)){
        estimate_fct[estimate$begin[i]:estimate$end[i]] <- 1
      }
    }

    # Update output object
    res_support$estimate[[q]] <- estimate
    res_support$estimate_fct[[q]] <- estimate_fct
    res_support$alpha[[q]] <- alpha
  }

  ####### Output
  return(res_support)
}

################################# ----
#' determine_intervals
################################# ----
#' @description Determine for which intervals a function is nonnull.
#' @return a matrix with 3 columns : "begin", "end" and "value". The two first
#'         columns define the begin and the end of the intervals and the third
#'         gives the mean values of each interval.
#' @param beta_fct a numerical vector.
#' @importFrom stats qnorm sd
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' # result of res_bliss1<-fit_Bliss(data=data1,param=param1)
#' data(res_bliss1)
#' intervals <- determine_intervals(res_bliss1$Bliss_estimate[[1]])
#' plot(data1$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s")
#' for(k in 1:nrow(intervals)){
#'    segments(data1$grids[[1]][intervals[k,1]],intervals[k,3],
#'            data1$grids[[1]][intervals[k,2]],intervals[k,3],col=2,lwd=4)
#' }
determine_intervals <- function(beta_fct){
  intervals <- data.frame()
  begin <- 1
  for (i in 2:length(beta_fct)){
    if (beta_fct[i] != beta_fct[i-1]) {
      end <- i - 1
      intervals <- rbind(intervals,
                         c(begin,
                           i-1,
                           beta_fct[i-1]))
      begin <- i
    }
  }
  intervals <- rbind(intervals,c(begin,i,beta_fct[i]))

  # IS 06/09/2018
  if(nrow(intervals) != 0) names(intervals) <- c("begin", "end", "value")
  return(intervals)
}

################################# ----
#' compute_starting_point_sann
################################# ----
#' @description Compute a starting point for the Simulated Annealing algorithm.
#' @return a matrix with 3 columns : "m", "l" and "b". The two first
#'         columns define the begin and the end of the intervals and the third
#'         gives the mean values of each interval.
#' @param beta_expe a numerical vector, the expectation of the coefficient
#'        function posterior sample.
#' @importFrom stats qnorm sd
#' @export
#' @examples
#' data(res_bliss1)
#' mystart<-compute_starting_point_sann(apply(res_bliss1$beta_sample[[1]],2,mean))
compute_starting_point_sann <- function(beta_expe){
  positive_vec <- unlist(sapply(beta_expe,function(value) if(value<0) 0  else value))
  negative_vec <- unlist(sapply(beta_expe,function(value) if(value>0) 0  else value))

  for(i in 1:length(beta_expe)){
    positive_vec[positive_vec < max(positive_vec)/20] <- 0
    negative_vec[negative_vec > min(negative_vec)/20] <- 0
  }

  tmp <- NULL
  count_p = 0;
  count_n = 0;
  for(i in 1:length(beta_expe)){
    # positive
    if(positive_vec[i] != 0){
      count_p = count_p +1
    }else{
      if(count_p > 0){
        upper = i-1
        lower = upper - count_p +1
        value = mean(beta_expe[lower:upper])
        tmp <- rbind(tmp,c(lower,upper,value))
        count_p = 0
      }
    }

    # negative
    if(negative_vec[i] != 0){
      count_n = count_n +1
    }else{
      if(count_n > 0){
        upper = i-1
        lower = upper - count_n +1
        value = mean(beta_expe[lower:upper])
        tmp <- rbind(tmp,c(lower,upper,value))
        count_n = 0
      }
    }
  }
  if(count_p > 0){
    upper = i
    lower = upper - count_p +1
    value = mean(beta_expe[lower:upper])
    tmp <- rbind(tmp,c(lower,upper,value))
  }
  if(count_n > 0){
    upper = i
    lower = upper - count_n +1
    value = mean(beta_expe[lower:upper])
    tmp <- rbind(tmp,c(lower,upper,value))
  }

  m <- NULL
  l <- NULL
  b <- NULL
  for(j in 1:nrow(tmp)){
    m <- c(m,floor(mean(tmp[j,1:2])))
    l_tmp <- m[j]-tmp[j,1]
    if(l_tmp == 0) l_tmp <- 1
    l <- c(l,l_tmp)
    b <- c(b,tmp[j,3])
  }

  res <- cbind(m,l,b)
  res
}

################################# ----
#' change_grid
################################# ----
#' @description Compute a function (evaluated on a grid) on a given (finer) grid.
#' @return a numerical vector, the approximation of the function on the new grid.
#' @param fct a numerical vector, the function to evaluate on the new grid.
#' @param grid a numerical vector, the initial grid.
#' @param new_grid a numerical vector, the new grid.
#' @export
#' @examples
#' grid <- seq(0,1,l=1e1)
#' new_grid <- seq(0,1,l=1e2)
#' fct <- 3*grid^2 + sin(grid*2*pi)
#' plot(grid,fct,type="o",lwd=2,cex=1.5)
#' lines(new_grid,change_grid(fct,grid,new_grid),type="o",col="red",cex=0.8)
change_grid <- function(fct,grid,new_grid){
  res <- rep(0,length(new_grid))
  for(i in 1:(length(grid)-1)){
    index <- new_grid %between% grid[i:(i+1)]

    res[index] = fct[i] + (fct[i+1]-fct[i]) *
      abs(new_grid[index] - grid[i]) / abs(grid[i] - grid[i+1])
  }

  index <- new_grid < min(grid)
  if(sum(index) == 1 ) res[index] = fct[1]
  if(sum(index) > 1  ) stop("The range of 'new_grid' is too large." )

  index <- new_grid > max(grid)
  if(sum(index) == 1 ) res[index] = fct[length(fct)]
  if(sum(index) > 1  ) stop("The range of 'new_grid' is too large." )

  return(res)
}


################################# ----
#' pdexp
################################# ----
#' @description Probability function of a discretized Exponentiel distribution.
#' @return a numerical vector, which is the prability function on \code{l_values}.
#' @param a a positive value, the mean of the Exponential prior.
#' @param l_values a numerical value, the discrete support of the parameter l.
#' @importFrom stats pgamma
#' @export
#' @examples
#' pdexp(10,seq(0,1,1))
#'
#' x <- seq(0,10,le=1e3)
#' plot(x,dexp(x,0.5),lty=2,type="l")
#' lines(pdexp(0.5,1:10),type="p")
pdexp <- function(a,l_values){
  step <- diff(l_values)[1] / 2
  probs <- pexp(l_values + step ,a) -
    pexp(l_values - step ,a)

  return(probs)
}

################################# ----
#' integrate_trapeze
################################# ----
#' @description Trapezoidal rule to approximate an integral.
#' @return a numerical value, the approximation.
#' @param x a numerical vector, the discretization of the domain.
#' @param y a numerical value, the discretization of the function to integrate.
#' @importFrom stats pgamma
#' @export
#' @examples
#' x <- seq(0,1,le=1e2)
#' integrate_trapeze(x,x^2)
#'
#' integrate_trapeze(data1$grids[[1]],t(data1$x[[1]]))
integrate_trapeze <- function(x,y){
  apply(as.matrix(y),2,function(vect)
    sum(diff(x)*(vect[-1]+vect[-length(vect)]))/2)
}

################################# ----
#' do_need_to_reduce
################################# ----
#' @description Determine if it is required to reduce the size of the grid time
#'        points for each functional covariate.
#' @return a boolean value.
#' @param param a list containing  p_threshold the maximum number of time points
#'        and p the actual number of time points for each functional covariate.
#' @export
#' @examples
#' data(param1)
#' param1$p <- sapply(data1$grids,length)
#'
#' do_need_to_reduce(param1)
do_need_to_reduce <- function(param){
  ###### Initialisation - Load objects
  p_threshold <- param[["p_threshold"]]

  ###### Initialisation - Define objects
  if(is.null(p_threshold)) p_threshold <- 100

  ###### Determine the condition
  cond <- any( param$p > p_threshold)

  ###### Output
  return(cond)
}

################################# ----
#' reduce_x
################################# ----
#' @description Reduce the number of time points.
#' @return a numerical value, the approximation.
#' @param data similar to \code{fit_Bliss}.
#' @param param a list containing values \code{Q}, \code{p} and \code{p}
#' @export
#' @examples
#' param <- list(Q=1,n=10,p=c(150),grids_lim=list(c(0,1)))
#' data <- sim(param)
#'
#' data(param1)
#' param1$n <- nrow(data$x[[1]])
#' param1$p <- sapply(data$grids,length)
#' param1$Q <- length(data$x)
#'
#' data <- reduce_x(data,param1)
reduce_x <- function(data,param){
  ###### Initialisation - Load objects
  p_threshold <- param[["p_threshold"]]
  Q <- param[["Q"]]
  p <- param[["p"]]
  n <- param[["n"]]

  ###### Initialisation - Define objects
  if(is.null(p_threshold)) p_threshold <- 100

  ###### Change the grid for the required functional covariates
  for(q in 1:Q){
    if(p[q] > p_threshold){
      new_x <- matrix(NA,n,p_threshold)

      # Define the new grid
      new_grid <- seq(min(data$grids[[q]]),
                      max(data$grids[[q]]),le=p_threshold)

      # Change all individual functional covariate
      for(i in 1:n){
        new_x[i,] <- change_grid(fct = data$x[[q]][i,],grid = data$grids[[q]],new_grid=new_grid)
      }

      # Change the daa object
      data$x[[q]] <- new_x
      data$grids[[q]] <- new_grid
    }
  }

  ###### Ouput
  return(data)
}

