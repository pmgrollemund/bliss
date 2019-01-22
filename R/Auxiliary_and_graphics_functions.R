################################# ----
#' image_Bliss
################################# ----
#' @description Plot an approximation of the posterior density.
#' @param beta_posterior_density a list. The result of the function
#'                 \code{compute_beta_posterior_density}.
#' @param param a list containing: (optional)
#' \describe{
#' \item{cols}{a vector of colors for the function image.}
#' \item{main}{an overall title for the plot.}
#' \item{xlab}{a title for the x axis.}
#' \item{ylab}{a title for the y axis.}
#' \item{ylim}{a numeric vectors of length 2, giving the y coordinate range.}
#' }
#' @param q an integer (optional), the index of the functional covariate to plot.
#' @importFrom stats quantile
#' @importFrom grDevices heat.colors
#' @export
#' @examples
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#' param1$cols <- colorRampPalette(brewer.pal(9,"Reds"))(1e2)
#' image_Bliss(res_bliss1$beta_posterior_density,param1,q=1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=3,lwd=2,type="s")
#'
#' # ---- not run
#' param1$cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
#' image_Bliss(res_bliss1$beta_posterior_density,param1,q=1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=3,lwd=2,type="s")
#'
#' param1$cols <- rev(heat.colors(12))
#' param1$col_scale <- "quantile"
#' image_Bliss(res_bliss1$beta_posterior_density,param1,q=1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=3,lwd=2,type="s")
#'
#' param1$cols <- rev(terrain.colors(12))
#' image_Bliss(res_bliss1$beta_posterior_density,param1,q=1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=2,lwd=2,type="s")
#'
#' param1$cols <- rev(topo.colors(12))
#' image_Bliss(res_bliss1$beta_posterior_density,param1,q=1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=2,lwd=2,type="s")
image_Bliss <- function(beta_posterior_density,param=list(),q=1){
  if(length(param) != 0){ # PMG 2018-11-13
    cols      <- param[["cols"]] #Ceci n'est pas une modification pmg 08-03-18
    main      <- param[["main"]]
    ylab      <- param[["ylab"]]
    xlab      <- param[["xlab"]]
    ylim      <- param[["ylim"]]
  }else{
    cols      <- NULL
    main      <- NULL
    ylab      <- NULL
    xlab      <- NULL
    ylim      <- NULL
  }
 if(is.null(cols)){
   cols <- rev(heat.colors(100))
 }
 if(is.null(main)){
   main <- ""
 }
  if(is.null(ylab)){
    ylab <- ""
  }
  if(is.null(xlab)){
    xlab <- ""
  }
  if(is.null(ylim)){
    ylim  <- range(ylim , beta_posterior_density[[q]]$grid_beta_t)
  }

 image(beta_posterior_density[[q]]$grid_t,
       beta_posterior_density[[q]]$grid_beta_t,
       beta_posterior_density[[q]]$density,
       col=cols,main=main,xlab=xlab,ylab=ylab,ylim=ylim)

 # Fix problems about the axis (PMG 2018-11-11)
 x_tmp <- max(abs(beta_posterior_density[[q]]$grid_t))
 y_tmp <- max(abs(beta_posterior_density[[q]]$grid_beta_t))
 axis(1,at=ylim+c(-x_tmp,x_tmp) )
 axis(2,at=ylim+c(-y_tmp,y_tmp) )
 axis(3,at=ylim+c(-x_tmp,x_tmp) )
 axis(4,at=ylim+c(-y_tmp,y_tmp) )
}

################################# ----
#' plot_bliss
################################# ----
#' @description A suitable representation of the Bliss estimate.
#' @param x the coordinates of points in the plot.
#' @param y the y coordinates of points in the plot.
#' @param connect a logical value (optional), to handle discontinuous function.
#'        If \code{connect} is TRUE, the plot is one line. Otherwise, several
#'        lines are used.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param ylim a numeric vectors of length 2, giving the y coordinate range.
#' @param ... Arguments to be passed to methods, such as graphical parameters
#'        (see \code{par}).
#' @importFrom grDevices gray.colors
#' @importFrom graphics abline segments axis image lines matplot plot text
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' # res_bliss1 <- fit_Bliss(data=data1,param=param1,verbose=TRUE)
#' }
#' data(res_bliss1)
#' ### Plot the BLiss estimate on a suitable grid
#' plot_bliss(res_bliss1$data$grids[[1]],
#'            res_bliss1$Bliss_estimate[[1]],lwd=2,bound=FALSE)
plot_bliss <- function(x,y,connect=FALSE,xlab="",ylab="",ylim=NULL,...){
 if(is.null(ylim)) ylim <- range(y)
 plot(x,x,type="n",ylim=ylim,ylab=ylab,xlab=xlab,...)
 lines_bliss(x,y,connect=connect,...)
}

################################# ----
#' lines_bliss
################################# ----
#' @description A suitable representation of the Bliss estimate.
#' @param x the coordinates of points in the plot.
#' @param y the y coordinates of points in the plot.
#' @param connect a logical value (optional), to handle discontinuous function.
#'        If \code{connect} is TRUE, the plot is one line. Otherwise, several
#'        lines are used.
#' @param ... Arguments to be passed to methods, such as graphical parameters
#'        (see \code{par}).
#' @export
#' @examples
#' ### Plot the BLiss estimate on a suitable grid
#' \donttest{
#' data(data1)
#' data(param1)
#' # res_bliss1 <- fit_Bliss(data=data1,param=param1,verbose=TRUE)
#' }
#' data(res_bliss1)
#' ### Plot the BLiss estimate on a suitable grid
#' plot_bliss(res_bliss1$data$grids[[1]],
#'            res_bliss1$Bliss_estimate[[1]],lwd=2,bound=FALSE)
#' lines_bliss(res_bliss1$data$grids[[1]],
#'             res_bliss1$Smooth_estimate[[1]],lty=2)
lines_bliss <- function(x,y,connect=FALSE,...){
  # Compute a more 'interpretable grid'
  extended_grid <- x - 0.5*diff(x)[1]
  extended_grid <- c(extended_grid,max(extended_grid) + diff(extended_grid)[1]) # PMG 11-11-18

 for(i in 1:length(extended_grid)){
  segments(extended_grid[i],y[i],
           extended_grid[i+1],y[i],
           ...)
  if(connect & i > 1)
   segments(extended_grid[i],y[i],
            extended_grid[i],y[i-1],
            ...)
 }
}



################################# ----
#' interpretation_plot
################################# ----
#' @description Provide a graphical representation of the functional data
#'              with a focus on the detected periods with the Bliss method.
#' @param data a list containing:
#' \describe{
#' \item{y}{a numerical vector, the outcomes.}
#' \item{x}{a list of matrices, the qth matrix contains the observations of the
#'       qth functional covariate at time points given by \code{grids}.}
#' \item{grids}{a list of numerical vectors, the qth vector is the grid of
#'        time points for the qth functional covariate.}
#' }
#' @param Bliss_estimate a numerical vector, the Bliss estimate.
#' @param q an integer (optional), the index of the functional covariate to plot.
#' @param centered a logical value (optional), If TRUE, the functional data are centered.
#' @param cols a numerical vector of colours (optional).
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' # result of res_bliss1 <- fit_Bliss(data=data1,param=param1,verbose=TRUE)
#' data(res_bliss1)
#' interpretation_plot(data=data1,Bliss_estimate=res_bliss1$Bliss_estimate,q=1)
#' interpretation_plot(data=data1,Bliss_estimate=res_bliss1$Bliss_estimate,q=1,centered=TRUE)
interpretation_plot <- function(data,Bliss_estimate,q=1,centered=FALSE,cols=NULL){
  # load some objects
  x <- data$x[[q]]
  y <- data$y
  grid  <- data$grids[[q]]

  # Graphical options
  if(is.null(cols)) cols  <- rev(heat.colors(length(y)))
  bg_col  <- rev(gray.colors(length(y)))
  grid_y <- seq(min(y),max(y),length=length(y))
  match  <- sapply(y,function(v) order(abs(v - grid_y))[1])
  cols   <- cols[match]
  bg_col  <- bg_col[match]

  lwds <- seq(0.1,2,length=length(y))
  lwds <- lwds[match]

  # Compute an extended grid which is simpler to understand graphical results
  extended_grid <- data$grids[[q]] - 0.5*diff(data$grids[[q]])[1]
  extended_grid <- c(extended_grid,max(extended_grid) + diff(extended_grid)[1]) # PMG 11-11-18

  # Need a new grid to plots
  new_grid <- rep(0,2*length(grid)+1)
  new_grid[seq(2,length(new_grid),by=2)] <- grid
  new_grid[seq(3,length(new_grid),by=2)] <- grid + 0.5*diff(grid)[1]
  new_grid[1] <- grid[1] - 0.5*diff(grid)[1]

  # Scale the data
  scaled_x <- scale(x,scale=F)

  # Compute intervals
  intervals <- determine_intervals(Bliss_estimate[[q]])
  intervals_nonnull <- intervals[intervals[,3] != 0,]
  intervals_null    <- intervals[intervals[,3] == 0,]

  # Drop the intervals of length 1 ##### XXXX to keep or not ?
  # intervals <- intervals[ intervals$begin != intervals$end,]

  index <- NULL
  if(centered){
    new_scaled_x <- matrix(NA,nrow=nrow(x),ncol=2*ncol(x)+1)
    new_scaled_x[,seq(2,ncol(new_scaled_x),by=2)] <- as.matrix(scaled_x)
    for(i in seq(3,ncol(new_scaled_x)-1,by=2))
      new_scaled_x[,i] <- 0.5*(new_scaled_x[,i-1]+new_scaled_x[,i+1])

    matplot(grid,t(scaled_x) ,type="l",lty=1,col=cols,lwd=lwds+1,
            main="",xlab="",ylab="")
    for(i in 1:nrow(intervals)){
      if(intervals[i,3]!=0) text( (extended_grid[intervals[i,1]] + grid[intervals[i,2]])/2 ,
                                  max(scaled_x), round(intervals[i,3],1),cex=1)
      if(intervals[i,3] == 0) index <- c(index,
                                         (2*intervals[i,1]-1) :
                                           (2*intervals[i,2]+1))
    }
    abline(v=c(extended_grid[unlist(intervals[,"begin"])],
               tail(extended_grid,1)),lty=2,col="gray60",lwd=1)

    # if(max(index) > length(grid))
    if(!is.null(index)){
      new_scaled_x[,-index] <- NA
      matplot(new_grid,t(new_scaled_x),type="l",lty=1,col=bg_col,lwd=lwds+1,add=TRUE)
    }
    abline(h=0)
  }else{
    new_x <- matrix(NA,nrow=nrow(x),ncol=2*ncol(x)+1)
    new_x[,seq(2,ncol(new_x),by=2)] <- as.matrix(x)
    for(i in seq(3,ncol(new_x)-1,by=2))
      new_x[,i] <- 0.5*(new_x[,i-1]+new_x[,i+1])

    matplot(grid,t(x) ,type="l",lty=1,col=cols,lwd=lwds+1,
            main="",xlab="",ylab="")
    for(i in 1:nrow(intervals)){
      if(intervals[i,3]!=0) text( (extended_grid[intervals[i,1]] + grid[intervals[i,2]])/2 ,
                                  max(x), round(intervals[i,3],1),cex=1)
      if(intervals[i,3] == 0) index <- c(index,
                                         (2*intervals[i,1]-1) :
                                           (2*intervals[i,2]+1))
    }
    abline(v=c(extended_grid[unlist(intervals[,"begin"])],
               tail(extended_grid,1)),lty=2,col="gray60",lwd=1)

    if(!is.null(index)){
      new_x[,-index] <- NA
      matplot(new_grid,t(new_x),type="l",lty=1,col=bg_col,lwd=lwds+1,add=TRUE)
    }
    x_center <- apply(x,2,mean)
    lines(grid,x_center)
  }
}

################################# ----
#' dposterior
################################# ----
#' @description Compute (non-normalized) posterior densities for a given parameter set.
#' @param posterior_sample a list given by the \code{Bliss_Gibbs_Sampler} function.
#' @param data a list containing
#' \describe{
#' \item{y}{a numerical vector, the outcomes.}
#' \item{x}{a list of matrices, the qth matrix contains the observations of the
#'       qth functional covariate at time points given by \code{grids}.}
#' }
#' @param theta a matrix or a vector which contains the parameter set.
#' @details If the \code{theta} is NULL, the posterior density is computed from
#'          the MCMC sample given in the \code{posterior_sample}.
#' @return Return the (log) posterior density, the (log) likelihood and the
#'         (log) prior density for the given parameter set.
#' @useDynLib bliss
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' # result of res_bliss1<-fit_Bliss(data=data1,param=param1)
#' data(res_bliss1)
#' # Compute the posterior density of the MCMC sample :
#' res_poste <- dposterior(res_bliss1$posterior_sample,data1)
dposterior <- function(posterior_sample,data,theta=NULL){
 if(!is.null(theta)){
  if(is.null(dim(theta))){
   rposterior <- as.matrix(t(theta))
  }
  K <- (ncol(theta)-2)/3
 }else{
  rposterior <- posterior_sample$trace
  K <- posterior_sample$param$K
 }
 N <- nrow(rposterior)
 Q <- length(as.vector(K))

 y <- data$y
 potential_intervals  <- posterior_sample$param$potential_intervals
 potential_intervals_dims <- list()
 for(q in 1:Q){
  potential_intervals_dims[[q]] <- c(ncol(data$x[[q]]),
                                     posterior_sample$param$l_values_length[[q]],
                                     length(data$y))
 }

 res <- dposterior_cpp(rposterior,y,N,as.vector(K),potential_intervals,potential_intervals_dims,
                       as.vector(posterior_sample$param$l_values_length),Q)
 colnames(res) <- c("posterior density","log posterior density",
                    "likelihood","log likelihood",
                    "prior density","log prior density")
 return(res)
}

################################# ----
#' compute_chains_info
################################# ----
#' @description Compute summaries of Gibbs Sampler chains.
#' @param chain a list given by the \code{Bliss_Gibbs_Sampler} function.
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
#' @return Return a list containing the estimates of \code{mu} and \code{sigma_sq}, the
#'         Smooth estimate and the chain autocorrelation for \code{mu}, \code{sigma_sq} and \code{beta}.
#' @useDynLib bliss
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cor
#' @export
#' @examples
#' \donttest{
#' param_sim <- list(Q=1,
#'                   n=100,
#'                   p=c(50),
#'                   grids_lim=list(c(0,1)))
#' data <- sim(param_sim,verbose=TRUE)
#'
#' param <- list(iter=5e2,
#'               K=c(3),
#'               n_chains = 3)
#' res_bliss <- fit_Bliss(data,param,verbose=TRUE,compute_density=FALSE,sann=FALSE)
#'
#' param$grids <- data$grids
#' chains_info1 <- compute_chains_info(res_bliss$chains[[1]],param)
#' chains_info2 <- compute_chains_info(res_bliss$chains[[2]],param)
#' chains_info3 <- compute_chains_info(res_bliss$chains[[3]],param)
#'
#' # Smooth estimates
#' ylim <- range(range(chains_info1$estimates$Smooth_estimate),
#' range(chains_info2$estimates$Smooth_estimate),
#' range(chains_info3$estimates$Smooth_estimate))
#' plot(data$grids[[1]],chains_info1$estimates$Smooth_estimate,type="l",ylim=ylim,
#' xlab="grid",ylab="")
#' lines(data$grids[[1]],chains_info2$estimates$Smooth_estimate,col=2)
#' lines(data$grids[[1]],chains_info3$estimates$Smooth_estimate,col=3)
#'
#' # Autocorrelation
#' plot(chains_info1$autocorr_lag[,1],type="h")
#' }
compute_chains_info <- function(chain,param){
  # Trace
  trace <- chain$trace
  # Beta sample
  Q <- length(param$K)
  beta_sample <- compute_beta_sample(chain,param,Q)

  # Estimate mu beta sigma
  mu_hat <- mean(trace[,'mu'])
  sigma_sq_hat <- mean(trace[,'sigma_sq'])
  Smooth_estimate <- apply(beta_sample[[1]],2,mean)

  # Autocorrelation
  lags <- 1:50
  n_iter <- nrow(trace)
  autocorr_lag <- NULL
  for(l in lags){
    indice     <- 1:(n_iter-l)
    indice_lag <- 1:(n_iter-l) + l

    options(warn = -1)
    cor_beta <- max(apply(beta_sample[[1]],2,function(v) {
      cor(v[indice],
          v[indice_lag])
    }),na.rm=T)
    options(warn = 0)

    autocorr_lag <- rbind(autocorr_lag,
                                c(cor(trace[indice,'mu'],
                                      trace[indice_lag,'mu']),
                                cor(trace[indice,'sigma_sq'],
                                    trace[indice_lag,'sigma_sq']),
                                cor_beta))
  }
  colnames(autocorr_lag) <- c("mu","sigma_sq","beta")

  # Result
  return(list(estimates = list(mu_hat          = mu_hat,
                               sigma_sq_hat    = sigma_sq_hat,
                               Smooth_estimate = Smooth_estimate),
              autocorr_lag = autocorr_lag))
}
