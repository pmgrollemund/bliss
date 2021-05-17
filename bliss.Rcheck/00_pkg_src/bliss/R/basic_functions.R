################################# ----
#' BIC_model_choice
################################# ----
#' @description Model selection with BIC criterion.
#' @return A numerical vector, the BIC values for the Bliss model for different
#'         K value.
#' @param Ks a numerical vector containing the K values.
#' @param iter an integer, the number of iteration for each run of \code{fit_Bliss}.
#' @param data a list containing:
#' \describe{
#' \item{Q}{an integer, the number of functional covariates.}
#' \item{y}{a numerical vector, the outcomes.}
#' \item{x}{a list of matrices, the qth matrix contains the observations of the
#'       qth functional covariate at time points given by \code{grids}.}
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of
#'        time points for the qth functional covariate.}
#' }
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

    res_bliss <- fit_Bliss(data,param_BIC,verbose=F,sann=F,compute_density=F)

    llkh <- dposterior(res_bliss$posterior_sample,data)
    llkh <- llkh[,'log likelihood']
    lML <- llkh[which.max(llkh)]
    BIC[i] <- (3*Ks[i]+2) * log(length(data$y)) - 2 * lML
  }
  return(BIC)
}


################################# ----
#' compute_beta_sample
################################# ----
#' @description Compute the posterior coefficient function from the posterior
#'              sample.
#' @return return a matrix containing the coefficient function posterior sample.
#' @param posterior_sample a list provided by the function \code{Bliss_Gibbs_Sampler}.
#' @param param a list containing:
#' \describe{
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for
#'       each covariate.}
#' \item{grids}{a numerical vector, the observation time points.}
#' \item{basis}{a vector of characters (optional) among : "uniform" (default),
#'       "epanechnikov", "gauss" and "triangular" which correspond to
#'       different basis functions to expand the coefficient function and the
#'       functional covariates.}
#' }
#' @param Q numeric
#' @param verbose write stuff if TRUE (optional).
#' @export
#' @examples
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' param1$grids<-data1$grids
#' # result of res_bliss1<-fit_Bliss(data=data1,param=param1)
#' data(res_bliss1)
#' beta_sample <- compute_beta_sample(posterior_sample=res_bliss1$posterior_sample,
#'                                    param=param1,Q=1)
#' indexes <- sample(nrow(beta_sample[[1]]),1e2,replace=FALSE)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
#' matplot(param1$grids[[1]],t(beta_sample[[1]][indexes,]),type="l",lty=1,col=cols,
#' xlab="grid",ylab="")
compute_beta_sample <- function(posterior_sample,param,Q,verbose=FALSE){
  if(verbose) cat("Compute the coefficient function posterior sample. \n")
  # Initialize parameters
  Ks     <- param[["K"]]
  grids  <- param[["grids"]]
  basiss <- param[["basis"]]
  ps <- sapply(grids,length)
  if(is.null(basiss)) basiss <- rep("Uniform",length(Ks))

  # Compute the coefficient function for each iteration of the Gibbs Sampler
  # and for each covariable
  beta_sample <- list()
  count <- 0
  for(q in 1:Q){
    K     <- Ks[q]
    grid  <- grids[[q]]
    p     <- ps[q]
    basis <- basiss[q]
    trace_tmp <- posterior_sample$trace[,(1+count):(count+3*K)]
    norm_val <- posterior_sample$param$normalization_values[[q]]

    beta_sample[[q]] <- compute_beta_sample_cpp(trace_tmp,
                                                K,grid,p,basis,norm_val)
    count <- count + 3*K[q]
  }

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
#' @param verbose write stuff if TRUE (optional).
#' @importFrom MASS bandwidth.nrd kde2d
#' @export
#' @examples
#' \donttest{
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' # result of res_bliss1<-fit_Bliss(data=data1,param=param1)
#' data(res_bliss1)
#' q <- 1
#' param_beta_density <- list(grid= data1[["grids"]][[q]],
#'                            iter= param1[["iter"]],
#'                            p   = param1[["p"]][q],
#'                            n        = length(data1[["y"]]),
#'                            thin     = param1[["thin"]],
#'                            burnin   = param1[["burnin"]],
#'                            lims_kde = param1[["lims_kde"]][[q]],
#'                            new_grid = param1[["new_grids"]][[q]],
#'                            lims_estimate = range(res_bliss1$Smooth_estimate[[q]]))
#' density_estimate <- compute_beta_posterior_density(res_bliss1$beta_sample[[q]],param_beta_density)
#' image(density_estimate$grid_t,
#'       density_estimate$grid_beta_t,
#'       density_estimate$density,col=rev(heat.colors(100)))
#' }
compute_beta_posterior_density <- function(beta_sample,param,verbose=FALSE){
  if(verbose)
    cat("Compute an approximation of the posterior density of the coefficient function.\n")
  # load optional objects
  grid <- param[["grid"]]
  thin     <- param[["thin"]]
  burnin   <- param[["burnin"]]
  lims_kde <- param[["lims_kde"]]
  new_grid <- param[["new_grid"]]

  p <- length(grid) # PMG 11/11/18
  N <- 512          # PMG 11/11/18
  iter <- nrow(beta_sample) -1 # PMG 11/11/18
  lims_estimate <- range(apply(beta_sample,2,mean)) # PMG 11/11/18

  # Initialize the necessary unspecified objects
  max_points <- 1e5
  if(is.null(burnin)) burnin <- floor(iter/5)
  if(is.null(thin))   thin   <- floor((iter-burnin)*p/max_points)

  # Check if the burnin isn't too large.
  if(2*burnin > iter){
    burnin <- floor(iter/5)
    if(verbose)
      cat("\t Burnin is too large. New burnin : ",burnin,".\n")
  }

  # Compute the coefficient function on the new grid (if required).
  if(!is.null(new_grid)){
    old_beta_sample <- beta_sample
    beta_sample <- matrix(0,nrow(beta_sample),p)
    if(verbose)
      cat("Compute the coefficient functions on the new grid.\n")
    for(i in 1:nrow(beta_sample)){
      beta_sample[i,] <- change_grid(old_beta_sample[i,],grid,new_grid)
    }
    param$old_grid <- grid
    param$grid     <- new_grid
    param$new_grid <- NULL # PMG 22/06/18
    grid           <- new_grid
  }

  # Thin the posterior sample
  thin_min   <- max(1,floor((iter-burnin)*p/max_points))
  if(thin <  thin_min){
    if(verbose)
      cat("\t 'thin = ",thin,"' is too small. Now, thin = ",
          thin_min,".\n",sep="")
    thin <- thin_min
  }
  if(verbose)
    cat("\t Thin the sample.\n")
  beta_sample <- beta_sample[seq(1+burnin,iter,by=thin),]

  # Perform the kde2d function
  if(verbose)
    cat("\t Perform the 'kde2d' function.\n")
  beta_x <- rep(grid,nrow(beta_sample))
  beta_y <- as.vector(t(beta_sample))

  h1 <- bandwidth.nrd(beta_x)
  h2 <- bandwidth.nrd(beta_y)
  if(h2 == 0){
    h2 <- 4 * 1.06 * sd(beta_y) * length(beta_y)^(-1/5)
  }
  points <- cbind(beta_x,beta_y)

  lims_kde <- c(range(beta_x), quantile(beta_y,c(0.025,0.975))) # PMG 04/07/18
  if(lims_kde[3] >= 0 ) lims_kde[3] <- -h2/2 # PMG 04/07/18
  if(lims_kde[4] <= 0 ) lims_kde[4] <- -h2/2 # PMG 04/07/18
  if(lims_kde[3] >= lims_estimate[1] ) lims_kde[3] <- lims_estimate[1]-h2/2 # PMG 04/07/18
  if(lims_kde[4] <= lims_estimate[2] ) lims_kde[4] <- lims_estimate[2]+h2/2 # PMG 04/07/18
  res_kde2d <- kde2d(x=beta_x,y=beta_y,lims=lims_kde,
                     n=N,h=c(h1,h2))

  # What to return ? # PMG 22/06/18
  if(!is.null(param$old_grid)){
    new_beta_sample <- beta_sample
  }else{
    new_beta_sample <- NULL
  }

  return(list(grid_t          = res_kde2d$x,
              grid_beta_t     = res_kde2d$y,
              density         = res_kde2d$z,
              new_beta_sample = beta_sample
  ))
  if(verbose)
    cat("\t Done.\n")
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
#' @param beta_sample_q a matrix. Each row is a coefficient function computed from the
#'        posterior sample.
#' @param gamma a numeric value, the default value is \code{0.5}.
#' @return a list containing:
#' \describe{
#'  \item{alpha}{a numerical vector. The approximated posterior probabilities
#'        that the coefficient function support covers \code{t} for each time
#'        points \code{t}.}
#'  \item{estimate}{a numerical vector, the support estimate.}
#'  \item{estimate_fct}{a numerical vector, another version of the support
#'        estimate.}
#' }
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' # result of res_bliss1<-fit_Bliss(data=data1,param=param1)
#' data(res_bliss1)
#' res_support <- support_estimation(res_bliss1$beta_sample[[1]])
#'
#' ### The estimate
#' res_support$estimate
#' ### Plot the result
#' grid <- res_bliss1$data$grids[[1]]
#' plot(grid,res_support$alpha,ylim=c(0,1),type="l",xlab="",ylab="")
#' for(k in 1:nrow(res_support$estimate)){
#'     segments(grid[res_support$estimate[k,1]],0.5,
#'              grid[res_support$estimate[k,2]],0.5,lwd=2,col=2)
#'     points(grid[res_support$estimate[k,1]],0.5,pch="|",lwd=2,col=2)
#'     points(grid[res_support$estimate[k,2]],0.5,pch="|",lwd=2,col=2)
#' }
#' abline(h=0.5,col=2,lty=2)
support_estimation <- function(beta_sample_q,gamma=0.5){
  # alpha: posterior probabilities
  alpha <- apply(beta_sample_q,2, function(vec) sum(vec != 0)/length(vec))

  # Support estimate
  tmp   <- rep(0,ncol(beta_sample_q))
  tmp2  <- which(alpha >= gamma)
  tmp[tmp2] <- 1
  estimate <- determine_intervals(tmp)
  estimate <- estimate[estimate$value != 0,]
  estimate$value <- NULL

  # Support estimate (vectorial version)
  estimate_fct <- rep(0,ncol(beta_sample_q))
  for(i in 1:nrow(estimate)){
    estimate_fct[estimate$begin[i]:estimate$end[i]] <- 1
  }

  return(list(alpha        = alpha,
              estimate     = estimate,
              estimate_fct = estimate_fct))
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

