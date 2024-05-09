################################# ----
#' image_Bliss
################################# ----
#' @description Plot an approximation of the posterior density.
#' @param beta_posterior_density a list. The result of the function
#'                 \code{compute_beta_posterior_density}.
#' @param param an optional  list containing arguments: col_low, col_mid, col_high,
#'          ylim, xlab, ylab, title.
#' @param q an integer (optional), the index of the functional covariate to plot.
#' @param to_print display the plot if TRUE.
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#'
#' image_Bliss(res_bliss1$beta_posterior_density,param1,q=1)
image_Bliss <- function(beta_posterior_density,param=list(),q=1,to_print=TRUE){
  ########## Initialization - Load objects
  col_low <- param[["col_low"]]
  col_mid <- param[["col_mid"]]
  col_high <- param[["col_high"]]

  if(is.null(col_low)) col_low <- "white"
  if(is.null(col_mid)) col_mid <- "yellow"
  if(is.null(col_high)) col_high <- "red"

  grid_t <- beta_posterior_density[[q]][["grid_t"]]
  grid_beta_t <- beta_posterior_density[[q]][["grid_beta_t"]]
  density <- beta_posterior_density[[q]][["density"]]

  ########## Pretreatment
  df_density <- expand.grid(grid_t,grid_beta_t)
  names(df_density) <- c("grid_t","grid_beta_t")

  df_density$density = as.vector(density)
  density_max <- max(df_density$density)

  ########## Do the plot
  p <- ggplot2::ggplot(df_density) + ggplot2::aes(x=grid_t,y=grid_beta_t) +
    ggplot2::geom_tile(ggplot2::aes(fill=density)) +
    ggplot2::scale_fill_gradient2(low = col_low,mid = col_mid,high = col_high,
                         midpoint = density_max/2,guide = "none") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme(panel.background = ggplot2::element_blank()) # theme_classic()

  ########## Add options
  if(!is.null(param[["ylim"]]))
    p <- p + ggplot2::ylim(param[["ylim"]][1],param[["ylim"]][2])
  if(!is.null(param[["xlab"]]))
    p <- p + ggplot2::xlab(param[["xlab"]])
  if(!is.null(param[["ylab"]]))
    p <- p + ggplot2::ylab(param[["ylab"]])
  if(!is.null(param[["title"]]))
    p <- p + graphics::title(param[["title"]])

  ########## Output
  if(to_print) print(p)

  ########## Output
  return(p)
}

################################# ----
#' lines_bliss
################################# ----
#' @description Add a line to a plot obtained with \code{image_Bliss}.
#' @param x the coordinates of points in the plot.
#' @param y the y coordinates of points in the plot.
#' @param col a color.
#' @param lty option corresponding to "linetype" of \code{geom_line}.
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#'
#' image_Bliss(res_bliss1$beta_posterior_density,param1,q=1) +
#' lines_bliss(res_bliss1$data$grids[[1]],res_bliss1$smooth_estimate[[1]])+
#' lines_bliss(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],col="purple")
#'
lines_bliss <- function(x,y,col="black",lty="solid"){
  ########## Initialization - Pretreatment
  df <- data.frame(x=as.vector(x),y=as.vector(y))

  ########## Do the plot
  word <- ggplot2::geom_line(data=df,ggplot2::aes(x=x,y=y),col=col,linetype=lty)

  ########## Output
  return(word)
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
  if(is.null(cols)) cols  <- rev(grDevices::heat.colors(length(y)))
  bg_col  <- rev(grDevices::gray.colors(length(y)))
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

    graphics::matplot(grid,t(scaled_x) ,type="l",lty=1,col=cols,lwd=lwds+1,
            main="",xlab="",ylab="")
    for(i in 1:nrow(intervals)){
      if(intervals[i,3]!=0) graphics::text( (extended_grid[intervals[i,1]] + grid[intervals[i,2]])/2 ,
                                  max(scaled_x), round(intervals[i,3],1),cex=1)
      if(intervals[i,3] == 0) index <- c(index,
                                         (2*intervals[i,1]-1) :
                                           (2*intervals[i,2]+1))
    }
    graphics::abline(v=c(extended_grid[unlist(intervals[,"begin"])],
               tail(extended_grid,1)),lty=2,col="gray60",lwd=1)

    # if(max(index) > length(grid))
    if(!is.null(index)){
      new_scaled_x[,-index] <- NA
      graphics::matplot(new_grid,t(new_scaled_x),type="l",lty=1,col=bg_col,lwd=lwds+1,add=TRUE)
    }
    graphics::abline(h=0)
  }else{
    new_x <- matrix(NA,nrow=nrow(x),ncol=2*ncol(x)+1)
    new_x[,seq(2,ncol(new_x),by=2)] <- as.matrix(x)
    for(i in seq(3,ncol(new_x)-1,by=2))
      new_x[,i] <- 0.5*(new_x[,i-1]+new_x[,i+1])

    graphics::matplot(grid,t(x) ,type="l",lty=1,col=cols,lwd=lwds+1,
            main="",xlab="",ylab="")
    for(i in 1:nrow(intervals)){
      if(intervals[i,3]!=0) graphics::text( (extended_grid[intervals[i,1]] + grid[intervals[i,2]])/2 ,
                                  max(x), round(intervals[i,3],1),cex=1)
      if(intervals[i,3] == 0) index <- c(index,
                                         (2*intervals[i,1]-1) :
                                           (2*intervals[i,2]+1))
    }
    graphics::abline(v=c(extended_grid[unlist(intervals[,"begin"])],
               tail(extended_grid,1)),lty=2,col="gray60",lwd=1)

    if(!is.null(index)){
      new_x[,-index] <- NA
      graphics::matplot(new_grid,t(new_x),type="l",lty=1,col=bg_col,lwd=lwds+1,add=TRUE)
    }
    x_center <- apply(x,2,mean)
    graphics::lines(grid,x_center)
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
    rposterior <- as.matrix(posterior_sample$trace)
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
#' @examples
#' \donttest{
#' a=1
#' }
compute_chains_info <- function(chain,param){
  # Trace
  trace <- chain$trace
  # Beta sample
  Q <- length(param$K)
  beta_sample <- compute_beta_sample(chain,param)

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

    # 07/07/2020 IS: replace warnings options (warn=-1) by suppressWarnings()
    suppressWarnings(
    cor_beta <- max(apply(beta_sample[[1]],2,function(v) {
                  cor(v[indice],
                  v[indice_lag])
                  }),na.rm=T)
    )

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
