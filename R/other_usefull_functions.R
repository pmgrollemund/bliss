################################# ----
#' dposterior
################################# ----
#' @description Compute the posterior density for a given parameter set.
#' @param res.Gibbs_Sampler a list given by the "Bliss_Gibbs_Sampler" function.
#' @param data a list containing
#' \describe{
#' \item{Q}{an integer, the number of covariates,}
#' \item{x}{a list of matrices, the functions x_qi(t) observed at time points given by grids,}
#' \item{y}{a numerical vector, the outcome values y_i.}
#' }
#' @param theta a matrix or a vector which contains the parameter set.
#' @details If the option theta is NULL, the posterior density is computed for
#' the MCMC sample given in the "res.Gibbs_Sampler" object.
#'
#' @return return the posterior density and the posterior log density for
#' the given parameter set.
#' @useDynLib bliss
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' # Compute the posterior density of the MCMC sample :
#' res_poste <- dposterior(res_Bliss_mult$res.Gibbs_Sampler,data1)
#' ### Compute the posterior density of the MCMC sample :
#' theta <- list( beta= rbind(c(3,4,-1), c(3,2,1)),
#' m = rbind(c(10,50,80),c(50,10,10)),
#' l = rbind(c(10,10,5),c(10,5,5)),
#' mu=c(1,1.5),sigma_sq=c(1,0.5))
#' dposterior(res_Bliss_mult$res.Gibbs_Sampler,data1,theta=theta)
dposterior <- function(res.Gibbs_Sampler,data,theta=NULL){
  if(!is.null(theta)){
    if( length(theta$mu) != 1){
      rposterior <- cbind(theta$beta,theta$m,theta$l,theta$mu,theta$sigma_sq)
    }else{
      rposterior <- t(matrix(
        c(theta$beta,theta$m,theta$l,theta$mu,theta$sigma_sq)))
    }
    K <- ncol(theta$beta)
  }else{
    rposterior <- res.Gibbs_Sampler$trace
    K <- res.Gibbs_Sampler$param$K[1]
  }
  N <- nrow(rposterior)

  y <- data$y
  all_intervals <- res.Gibbs_Sampler$param$all_intervals
  scale_ml      <- res.Gibbs_Sampler$param$scale_ml
  all_intervals_dims <- c(ncol(data$x_mult[[1]]),
                          res.Gibbs_Sampler$param$l_max,
                          length(data$y))
  lambda <- 5
  lmax <- res.Gibbs_Sampler$param$l_max

  res <- dposterior_cpp(rposterior,y,N,K,all_intervals[[1]],all_intervals_dims,
                        lambda,lmax)
  colnames(res) <- c("posterior density","log posterior density",
                     "likelihood","log likelihood",
                     "prior density","log prior density")
  return(res)
}

################################# ----
#' interpretation_plot
################################# ----
#' @description Provide a graphical representation of the functional data
#' with a focus on the detected periods with the Bliss method.
#' @param estimate a numerical vector, the Bliss estimate.
#' @param data a list containing y, x and grid.
#' @param q an integer, the index of the functional covariate to examine.
#' @param centered a Boolean value. If TRUE, the functional data are centered.
#' @param cols a numerical vector of colours.
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' interpretation_plot(res_Bliss_mult$Bliss_estimate[[1]],data1)
#' interpretation_plot(res_Bliss_mult$Bliss_estimate[[1]],data1,centered=FALSE)
#' }
interpretation_plot <- function(estimate,data,q=1,centered=FALSE,cols=NULL){
  # load some objects
  x <- data$x[[q]]
  y <- data$y
  grid  <- data$grids[[q]]
  grid0 <- data$grids[[q]] - 0.5*diff(data$grids[[q]])[1]
  grid0 <- c(grid0 , max(data$grids[[q]]) + 0.5*diff(data$grids[[q]])[1]  )

  x_centered <- apply(x,2,function(vec) vec - mean(vec))

  new_grid <- rep(0,2*length(grid)+1)
  new_grid[seq(2,length(new_grid),by=2)] <- grid
  new_grid[seq(3,length(new_grid),by=2)] <- grid + 0.5*diff(grid)[1]
  new_grid[1] <- grid[1] - 0.5*diff(grid)[1]

  new_x <- matrix(NA,nrow=nrow(x),ncol=2*ncol(x)+1)
  new_x[,seq(2,ncol(new_x),by=2)] <- as.matrix(x)
  for(i in seq(3,ncol(new_x)-1,by=2))
    new_x[,i] <- 0.5*(new_x[,i-1]+new_x[,i+1])

  new_x_centered <- matrix(NA,nrow=nrow(x),ncol=2*ncol(x)+1)
  new_x_centered[,seq(2,ncol(new_x_centered),by=2)] <- as.matrix(x_centered)
  for(i in seq(3,ncol(new_x_centered)-1,by=2))
    new_x_centered[,i] <- 0.5*(new_x_centered[,i-1]+new_x_centered[,i+1])


  intervals <- interval_detection(estimate)
  intervals$value[3] <- 0 ## Mais qu'est-ce que c'est que ca ?
  intervals$end[4] <- 13 ## Mais qu'est-ce que c'est que ca ?
  intervals$begin[5] <- 14 ## Mais qu'est-ce que c'est que ca ?
  # Drop the intervals of length 1
  # intervals <- intervals[ intervals$begin != intervals$end,]

  if(is.null(cols)) cols  <- rev(heat.colors(length(y)))
  cols2  <- rev(gray.colors(length(y)))
  grid_y <- seq(min(y),max(y),length=length(y))
  match  <- sapply(y,function(v) order(abs(v - grid_y))[1])
  cols   <- cols[match]
  cols2  <- cols2[match]

  lwds <- seq(0.1,2,length=length(y))
  lwds <- lwds[match]

  intervals_nonnull <- intervals[intervals[,3] != 0,]
  intervals_null    <- intervals[intervals[,3] == 0,]

  index <- NULL
  if(centered){
    matplot(grid,t(x_centered) ,type="l",lty=1,col=cols,lwd=lwds+1,
            main="",xlab="",ylab="")
    for(i in 1:nrow(intervals)){
      if(intervals[i,3]!=0) text( (grid0[intervals[i,1]] + grid[intervals[i,2]])/2 ,
                                  max(x_centered), round(intervals[i,3],1),cex=1)
      # if(intervals[i,3]<0) text( (grid0[intervals[i,1]] + grid[intervals[i,2]])/2 ,
      #                            min(x_centered), round(intervals[i,3],1),cex=1)
      if(intervals[i,3] == 0) index <- c(index,
                                         (2*intervals[i,1]-1) :
                                           (2*intervals[i,2]+1))
    }
    abline(v=c(grid0[unlist(intervals[,"begin"])],
               tail(grid0,1)),lty=2,col="gray60",lwd=1)

    # if(max(index) > length(grid))
    x_centered_2 <- x_centered_new
    if(!is.null(index)){
      x_centered_2[,-index] <- NA
      matplot(grid_new,t(x_centered_2),type="l",lty=1,col=cols2,lwd=lwds+1,add=TRUE)
    }
    abline(h=0)
  }else{
    matplot(grid,t(x) ,type="l",lty=1,col=cols,lwd=lwds+1,
            main="",xlab="",ylab="")
    for(i in 1:nrow(intervals)){
      if(intervals[i,3]!=0) text( (grid0[intervals[i,1]] + grid[intervals[i,2]])/2 ,
                                  max(x), round(intervals[i,3],1),cex=1)
      # if(intervals[i,3]<0) text( (grid0[intervals[i,1]] + grid[intervals[i,2]])/2 ,
      #                            min(x_centered), round(intervals[i,3],1),cex=1)
      if(intervals[i,3] == 0) index <- c(index,
                                         (2*intervals[i,1]-1) :
                                           (2*intervals[i,2]+1))
    }
    abline(v=c(grid0[unlist(intervals[,"begin"])],
               tail(grid0,1)),lty=2,col="gray60",lwd=1)

    x_2 <- x_new
    if(!is.null(index)){
      x_2[,-index] <- NA
      matplot(grid_new,t(x_2),type="l",lty=1,col=cols2,lwd=lwds+1,add=TRUE)
    }
    x_center <- apply(x,2,mean)
    lines(grid,x_center)
  }
}
################################# ----
#' diagnostics
################################# ----
#' @description perform some diagnostics for the chains resulting of the
#'              Gibbs Sampler algorithm.
#' @return a list containing the diagnostics which can be plotted with the
#'         function "plot_diagnostics".
#' @details The implementation does not tackle the issue of several functional covariates.
#' @param chains a list
#' \describe{
#' \item{res.Gibbs_Sampler}{a list resulting of the function Bliss_Gibbs_Sampler.}
#' \item{beta_functions}{a matrix which is the result of the function compute_beta_functions.}
#' \item{posterior_density_estimate}{a list which is the result of the function density_estimation.}
#' }
#' @param param a list
#' \describe{
#' \item{iter}{the number of iterations.}
#' \item{burnin}{the number of iterations to drop.}
#' \item{K}{the number of intervals of the beta_i's functions.}
#' \item{p}{the number of time points.}
#' \item{ts}{a vector. The sample of time points t_j used for compute some diagnostics of beta(t_j). (optional)}
#' \item{l_ts}{an integer, the number of time points to sample, if ts is not specified. (optional)}
#' \item{lag_max}{an integer, the maximal lag when compute the autocorrelation of the trace.}
#' }
#' @importFrom stats cor density
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' param1$n_chains <- 3
#' param1$iter <- 1e3
#' param1$burnin <- 1e2
#' param1$display <- FALSE
#' param1$compute_posterior <- FALSE
#' res_bliss_chains <- Bliss_multiple(data1,param1)
#' res_diagnostic <- diagnostics(res_bliss_chains$chains,param1)
#' cols <- makeTransparent(2:5,0.8)
#' plot(res_diagnostic$hist_mu[[1]],border=0,col=cols[1],main="",xlab="",ylab="")
#' for(c in 2:3){
#'    plot(res_diagnostic$hist_mu[[c]],add=TRUE,border=0,col=cols[c])
#' }
#' plot(res_diagnostic$density_beta[[4]][[1]],main="",xlab="",ylab="")
#' for(c in 2:3){
#'    lines(res_diagnostic$density_beta[[4]][[c]],col=c)
#' }
diagnostics <- function(chains,param){
  display <- param$display
  if(is.null(display)) display <- TRUE

  if(display) cat("Compute some diagnostics on the posterior sample. \n")
  # load objects
  n_chains  <- param[["n_chains"]]
  burnin    <- param[["burnin"]]
  K         <- param[["K"]]
  p         <- length(param[["grids"]][[1]])
  iter      <- param[["iter"]]
  # load optional objects
  ts      <- param[["ts"]]
  l_ts    <- param[["l_ts"]]
  lag_max <- param[["lag_max"]]
  # Initialize the necessary unspecified objects
  if(is.null(l_ts))    l_ts    <- 10
  if(is.null(ts))      ts      <- floor(seq(1,p,l=l_ts))
  if(is.null(lag_max)) lag_max <- min(100,floor(iter/50))
  l_ts <- length(ts)
  # Initialize
  n_class     <- min(1e3 , floor( (iter-burnin+1)/100 ) )
  lags        <- 1:lag_max
  DF_mu       <-  data.frame()
  DF_sigma    <-  data.frame()
  DF_beta     <-  list()
  trace_mu    <- NULL
  trace_sigma <- NULL
  trace_beta  <- NULL
  autocorr_lag_mu    <- NULL
  autocorr_lag_sigma <- NULL
  autocorr_lag_beta  <- NULL
  ylim_mu        <- NULL
  ylim_sigma     <- NULL
  ylim_beta      <- rep(NA,2*l_ts)
  dim(ylim_beta) <- c(2,l_ts)
  hist_mu       <- list()
  hist_sigma    <- list()
  hist_beta     <- list()
  density_mu    <- list()
  density_sigma <- list()
  density_beta  <- list()
  length(DF_beta)       <- l_ts
  length(hist_mu)       <- length(n_chains)
  length(hist_sigma)    <- length(n_chains)
  length(hist_beta)     <- l_ts
  length(density_mu)    <- length(n_chains)
  length(density_sigma) <- length(n_chains)
  length(density_beta)  <- l_ts

  #### diagnostic for mu and sigma_sq
  # Gather the different chains
  for(j in 1:n_chains){
    trace_chain <- chains[[j]]$res.Gibbs_Sampler$trace[-(1:burnin),]
    trace_mu    <- cbind(trace_mu   ,trace_chain[,1])
    trace_sigma <- cbind(trace_sigma,trace_chain[,1+K+1])

    DF_mu <- rbind(DF_mu,data.frame(chain = paste("chain_",j,sep=""),
                                    obs   =  trace_chain[,1]))
    DF_sigma <- rbind(DF_sigma,data.frame(chain = paste("chain_",j,sep=""),
                                          obs   =  trace_chain[,1+K+1]))
  }

  # Autocorr mu
  for(j in 1:n_chains){
    n_iter <- nrow(trace_mu)
    autocorr_lag_chain <- NULL
    for(l in lags){
      indice     <- 1:(n_iter-l)
      indice_lag <- 1:(n_iter-l) + l

      autocorr_lag_chain <- c(autocorr_lag_chain,
                              cor(trace_mu[indice,j],
                                  trace_mu[indice_lag,j]))

    }
    autocorr_lag_mu <- cbind(autocorr_lag_mu,autocorr_lag_chain)
  }

  # Autocorr sigma
  for(j in 1:n_chains){
    n_iter <- nrow(trace_sigma)
    autocorr_lag_chain <- NULL
    for(l in lags){
      indice     <- 1:(n_iter-l)
      indice_lag <- 1:(n_iter-l) + l

      autocorr_lag_chain <- c(autocorr_lag_chain,
                              cor(trace_sigma[indice,j],
                                  trace_sigma[indice_lag,j]))

    }

    autocorr_lag_sigma <- cbind(autocorr_lag_sigma,autocorr_lag_chain)
  }

  # Compute the histograms of mu and sigma_sq
  breaks_mu    <- seq(min(DF_mu$obs),max(DF_mu$obs),l=n_class+1)
  breaks_sigma <- seq(min(DF_sigma$obs),max(DF_sigma$obs),l=n_class+1)
  step_mu      <- diff(breaks_mu)[1]
  step_sigma   <- diff(breaks_sigma)[1]
  for(j in 1:n_chains){
    hist_mu[[j]] <- hist(DF_mu$obs[DF_mu$chain==paste("chain_",j,sep="")],
                         plot = FALSE,breaks=breaks_mu)
    density_mu[[j]] <- density(DF_mu$obs[DF_mu$chain==paste("chain_",j,sep="")])
    ylim_mu <- range(ylim_mu, hist_mu[[j]]$density,density_mu[[j]]$y)

    hist_sigma[[j]] <- hist(DF_sigma$obs[DF_sigma$chain==paste("chain_",j,sep="")],
                            plot = FALSE,breaks=breaks_sigma)
    density_sigma[[j]] <- density(DF_sigma$obs[DF_sigma$chain==paste("chain_",j,sep="")])
    ylim_sigma <- range(ylim_sigma, hist_sigma[[j]]$density,density_sigma[[j]]$y)
  }

  #### diagnostic for beta(t) (for some t)
  breaks_beta <- rep(0,(n_class+1)*l_ts)
  dim(breaks_beta) <- c(n_class+1, l_ts)
  trace_beta <- rep(0,l_ts*(iter-burnin+1)*n_chains)
  dim(trace_beta) <- c(l_ts,iter-burnin+1,n_chains)
  autocorr_lag_beta      <- rep(0,l_ts*lag_max*n_chains)
  dim(autocorr_lag_beta) <- c(l_ts,n_chains,lag_max)
  for(o in 1:l_ts){
    for(j in 1:n_chains){
      trace_chain <- chains[[j]]$beta_functions[[1]][-(1:burnin),ts[o]]
      trace_beta[o,,j] <- trace_chain

      DF_beta[[o]] <- rbind(DF_beta[[o]],data.frame(chain = paste("chain_",j,sep=""),
                                                    obs   =  trace_chain))
    }

    #  autocorrelation of beta
    for(j in 1:n_chains){
      n_iter <- nrow(as.matrix(trace_beta[o,,]))
      autocorr_lag_chain <- NULL
      for(l in lags){
        indice     <- 1:(n_iter-l)
        indice_lag <- 1:(n_iter-l) + l

        autocorr_lag_chain <- c(autocorr_lag_chain,
                                cor(trace_beta[o,indice,j],
                                    trace_beta[o,indice_lag,j]))

      }
      autocorr_lag_beta[o,j,] <- autocorr_lag_chain
    }
    # Compute the histograms and the kernel density estimations
    breaks_beta[,o] <- seq(min(DF_beta[[o]]$obs),max(DF_beta[[o]]$obs),l=n_class+1)
    step_beta   <- diff(breaks_beta[,o])[1]
    for(j in 1:n_chains){
      hist_beta[[o]][[j]] <- hist(DF_beta[[o]]$obs[DF_beta[[o]]$chain==paste("chain_",j,sep="")],
                                  plot = FALSE,breaks=breaks_beta[,o])
      density_beta[[o]][[j]] <- density(DF_beta[[o]]$obs[DF_beta[[o]]$chain==paste("chain_",j,sep="")])
      ylim_beta[,o] <- range(ylim_beta[,o], hist_beta[[o]][[j]]$density,
                             density_beta[[o]][[j]]$y,na.rm = T)
    }
  }

  result<-list(DF_mu=DF_mu, trace_mu=trace_mu, autocorr_lag_mu=autocorr_lag_mu,
               hist_mu=hist_mu, density_mu=density_mu, ylim_mu=ylim_mu,

               DF_sigma=DF_sigma, trace_sigma=trace_sigma,
               autocorr_lag_sigma=autocorr_lag_sigma,hist_sigma=hist_sigma,
               density_sigma=density_sigma, ylim_sigma=ylim_sigma,

               DF_beta=DF_beta, trace_beta=trace_beta,
               autocorr_lag_beta=autocorr_lag_beta,hist_beta=hist_beta,
               density_beta=density_beta, ylim_beta=ylim_beta,
               breaks_beta=breaks_beta,

               lags=lags, ts=ts)
  class(result) = c("blissDiag")
  return(invisible(result))
}

################################# ----
#' autocorr
################################# ----
#' @description compute the autocorrelation of the sample ( x_i(t) )_{i=1,...,n}.
#' @return  a symmetric matrix.
#' @param data a list containing 1) x, the functions x_i(t) and 2) grid, the grid of time points.
#' @param plot a logical value. If it is TRUE, an image (heat map) of the
#'                 autocorrelation matrix is plotted. (optional)
#' @importFrom graphics image
#' @export
#' @examples
#' library(RColorBrewer)
#' ### Autocorrelation of the function x_i(t)
#' param <- list(n=50,p=100,beta_type="smooth")
#' data <- sim(param)
#' res_autocorr <- autocorr(data)
#' cols <- rev(colorRampPalette(brewer.pal(9,"YlOrRd"))(50))
#' image(res_autocorr,col=cols)
#' \donttest{
#' ### Autocorrelation of the function beta_j(t).
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' beta_functions_autocorr <- autocorr(list(grid = data1$grids[[1]],
#'                                          x = res_Bliss_mult$beta_functions[[1]]))
#' image(beta_functions_autocorr)
#' }
autocorr <- function(data,plot=F){
  # Initialize
  x     <- data$x
  x.cor <- matrix(0,ncol(x),ncol(x))

  # Compute the correlations
  for(i in 1:ncol(x)){
    for(j in i:ncol(x)){
      x.cor[i,j] <- cor(x[,i],x[,j])
      x.cor[j,i] <- cor(x[,i],x[,j])
    }
  }

  # Plot the autocorrelation ?
  if(plot) image(x.cor)

  # Return the result
  return(x.cor)
}


################################# ----
#' corr_matrix
################################# ----
#' @description compute an autocorrelation matrix to simulate functions x_i(t).
#' @return a symmetric matrix.
#' @param diagonal a numerical vector corresponding to the diagonal of the final matrix.
#' @param ksi a "coefficient of correlation". See the article Bliss, Section 3.1 for more details.
#' @export
#' @examples
#' ### Test 1 : weak autocorrelation
#' ksi     <- 1
#' diagVar <- abs(rnorm(100,50,5))
#' Sigma   <- corr_matrix(diagVar,ksi^2)
#' persp(Sigma)
#' ### Test 2 : strong autocorrelation
#' ksi     <- 0.2
#' diagVar <- abs(rnorm(100,50,5))
#' Sigma   <- corr_matrix(diagVar,ksi^2)
#' persp(Sigma)
corr_matrix <- function(diagonal,ksi){
  # Initialize
  p <- length(diagonal)
  res <- diag(diagonal)

  # Compute the correlation matrix
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      res[i,j] <- exp(-ksi*(i-j)^2/p)*sqrt(res[i,i]*res[j,j])
      res[j,i] <- exp(-ksi*(i-j)^2/p)*sqrt(res[i,i]*res[j,j])
    }
  }

  # return the matrix
  return(res)
}
