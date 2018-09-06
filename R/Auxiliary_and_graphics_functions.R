################################# ----
#' image_Bliss
################################# ----
#' @description Plot an approximation of the posterior density.
#' @param beta_posterior_density a list. The result of the function
#'                 \code{compute_beta_posterior_density}.
#' @param param a list containing:
#' \describe{
#' \item{cols}{a vector of colors for the function image (optional).}
#' \item{col_scale}{a character vector.} XXXXX
#' \item{ylim}{a numerical two-vector (optional).}
#' \item{main}{a character string.}
#' }
#' @importFrom stats quantile
#' @importFrom grDevices heat.colors
#' @export
#' @examples
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' data(res_bliss1)
#' param1$cols <- colorRampPalette(brewer.pal(9,"Reds"))(1e2)
#' image_Bliss(res_bliss1$beta_posterior_density[[1]],param1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=3,lwd=2,type="s")
#'
#' \donttest{
#' # ---- not run
#' param1$cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
#' image_Bliss(res_bliss1$beta_posterior_density[[1]],param1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=3,lwd=2,type="s")
#'
#' param1$cols <- rev(heat.colors(12))
#' param1$col_scale <- "quantile"
#' image_Bliss(res_bliss1$beta_posterior_density[[1]],param1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=3,lwd=2,type="s")
#'
#' param1$cols <- rev(terrain.colors(12))
#' image_Bliss(res_bliss1$beta_posterior_density[[1]],param1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=2,lwd=2,type="s")
#'
#' param1$cols <- rev(topo.colors(12))
#' image_Bliss(res_bliss1$beta_posterior_density[[1]],param1)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(res_bliss1$data$grids[[1]],res_bliss1$data$betas[[1]],col=2,lwd=2,type="s")
#' }
image_Bliss <- function(beta_posterior_density,param){
 cols      <- param[["cols"]] #Ceci n'est pas une modification pmg 08-03-18
 col_scale <- param[["col_scale"]] # XXX
 ylim      <- param[["ylim"]]
 main      <- param[["main"]]
 if(is.null(cols)){
  cols <- rev(heat.colors(100))
 }

 # image() needs x, y and z for inputs
 # x is objectBliss$beta_posterior_density$grid_t
 # x is objectBliss$beta_posterior_density$grid_beta_t
 # Z is objectBliss$beta_posterior_density$density (not sure!!)

 breaks <- seq(min(as.vector(beta_posterior_density$density)),
                max(as.vector(beta_posterior_density$density)),
                length=length(cols)+1)

 xlim <- range(beta_posterior_density$grid_t)
 if(is.null(ylim)) ylim <- range(beta_posterior_density$grid_beta_t)
 if(is.null(main)) main <- ""
 image(beta_posterior_density$density,
       col=cols,breaks = breaks,main=main,ylim=ylim,xlim=xlim)
}


################################# ----
#' plot_bliss
################################# ----
#' @description A suitable representation of the Bliss estimate.
#' @param extended_grid a numerical vector, a new grid to provide a suitable
#'        representation.
#' @param fct a numerical vector, the function to plot.
#' @param bound a logical value. If bound is TRUE, the plot of fct is one line.
#' Otherwise, several lines are used to plot fct.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @importFrom grDevices gray.colors
#' @importFrom graphics abline segments
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' res_bliss1 <- fit_Bliss(data=data1,param=param1,progress=TRUE)
#' }
#' data(res_bliss1)
#' ### Plot the BLiss estimate on a suitable grid
#' plot_bliss(res_bliss1$data$grids[[1]],
#'            res_bliss1$Bliss_estimate[[1]],lwd=2,bound=FALSE)
plot_bliss <- function(extended_grid,fct,bound=FALSE,...){
 ylim <- range(fct)
 plot(extended_grid,extended_grid,type="n",ylim=ylim,...)
 lines_bliss(extended_grid,fct,bound=bound,...)
}




################################# ----
#' lines_bliss
################################# ----
#' @description A suitable representation of the Bliss estimate.
#' @param extended_grid a numerical vector, a new grid to provide a suitable
#'        representation.
#' @param fct a numerical vector, the function to plot.
#' @param bound a logical value. If bound is TRUE, the plot of fct is one line.
#' Otherwise, several lines are used to plot fct.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @export
#' @examples
#' ### Plot the BLiss estimate on a suitable grid
#' \donttest{
#' data(data1)
#' data(param1)
#' res_bliss1 <- fit_Bliss(data=data1,param=param1,progress=TRUE)
#' }
#' data(res_bliss1)
#' ### Plot the BLiss estimate on a suitable grid
#' plot_bliss(res_bliss1$data$grids[[1]],
#'            res_bliss1$Bliss_estimate[[1]],lwd=2,bound=FALSE)
#' lines_bliss(res_bliss1$data$grids[[1]],
#'             res_bliss1$Smooth_estimate[[1]],lty=2)
lines_bliss <- function(extended_grid,fct,bound=FALSE,...){
 for(i in 1:length(extended_grid)){
  segments(extended_grid[i],fct[i],
           extended_grid[i+1],fct[i],
           ...)
  if(bound & i > 1)
   segments(extended_grid[i],fct[i],
            extended_grid[i],fct[i-1],
            ...)
 }
}



################################# ----
#' interpretation_plot
################################# ----
#' @description Provide a graphical representation of the functional data
#'              with a focus on the detected periods with the Bliss method.
#' @param estimate a numerical vector, the Bliss estimate.
#' @param data a list containing:
#' \describe{
#' \item{Q}{an integer, the number of covariates,}
#' \item{x}{a list of matrices, the qth matrix contains the observation of the
#'       qth functional covariate at time points given by grids,}
#' \item{y}{a numerical vector, the outcome values y_i.}
#' }
#' @param q an integer, the index of the functional covariate to plot.
#' @param centered a logical value. If TRUE, the functional data are centered.
#' @param cols a numerical vector of colours.
#' @export
#' @examples
#' \donttest{
#' # Not run!!
#' data(data1)
#' data(param1)
#' res_bliss1 <- fit_Bliss(data=data1,param=param1,progress=TRUE)
#' data(res_bliss1)
#' interpretation_plot(res_bliss1$Bliss_estimate[[1]],data1)
#' interpretation_plot(res_bliss1$Bliss_estimate[[1]],data1,centered=FALSE)
#' }
interpretation_plot <- function(estimate,data,q=1,centered=FALSE,cols=NULL){
  # load some objects
  x <- data$x[[q]]
  y <- data$y
  grid  <- data$grids[[q]]
  grid0 <- data$grids[[q]] - 0.5*diff(data$grids[[q]])[1]
  grid0 <- c(grid0 , max(data$grids[[q]]) + 0.5*diff(data$grids[[q]])[1]  )

  x_centered <- scale(x,scale=F)

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
#' autocorr
################################# ----
#' @description Compute the autocorrelation of the functional covariate sample.
#' @return A symmetric matrix.
#' @param data a list containing:
#' \describe{
#' \item{Q}{an integer, the number of covariates,}
#' \item{x}{a list of matrices, the qth matrix contains the observation of the
#'       qth functional covariate at time points given by grids,}
#' \item{y}{a numerical vector, the outcome values y_i.}
#' }
#' @param plot a logical value. If it is TRUE, an image (heat map) of the
#'                 autocorrelation matrix is plotted. (optional)
#' @param q XXXX
#' @importFrom graphics image
#' @export
#' @examples
#' library(RColorBrewer)
#' ### Autocorrelation of the function x_i(t)
#' param <- list(n=50,p=100,beta_type="smooth",Q=1)
#' data <- sim(param)
#' res_autocorr <- autocorr(data)
#' cols <- rev(colorRampPalette(brewer.pal(9,"YlOrRd"))(50))
#' image(res_autocorr,col=cols)
#' ### Autocorrelation of the function beta_j(t).
#' data(data1)
#' data(param1)
#' # result of res_bliss1<-fit_Bliss(data=data1,param=param1)
#' data(res_bliss1)
#' \donttest{
#' # Example to modify!!
#' beta_functions_autocorr <- autocorr(list(grid = data1$grids[[1]],
#'                                          x = res_bliss1[[1]]))
#' image(beta_functions_autocorr)
#' }
autocorr <- function(data,plot=F,q=1){
 # Initialize
 x     <- data$x[[q]]
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
#' diagnostics
################################# ----
#' @description Compute some diagnostics for the chains resulting of the
#'              function \code{Bliss_Gibbs_Sampler}.
#' @return A list containing the diagnostics which can be plotted with the
#'         function "plot_diagnostics".
#' @details The implementation does not tackle the issue of several
#'          functional covariates.
#' @param chains a list. Each element contains:
#' \describe{
#' \item{res.Gibbs_Sampler}{a list resulting of the function Bliss_Gibbs_Sampler.}
#' }
#' @param param a list
#' \describe{
#' \item{iter}{an integer, the number of iterations of the Gibbs sampler algorithm.}
#' \item{burnin}{an integer, the number of iteration to drop from the Gibbs
#'       sample. (optional)}
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for
#'       each covariate.}
#' \item{p}{XXXXXX}
#' \item{ts}{a vector. The sample of time points t_j used for compute some diagnostics of beta(t_j). (optional)}
#' \item{l_ts}{an integer, the number of time points to sample, if ts is not specified. (optional)}
#' \item{lag_max}{an integer, the maximal lag when compute the autocorrelation of the trace.}
#' }
#' @param progress a logical value. If TRUE, the algorithm progress is progressed.
#'         (optional)
#' @importFrom stats cor density
#' @export
#' @examples
#' \donttest{
#' # Not run!!!
#' data(data1)
#' data(param1)
#' param1$n_chains <- 3
#' param1$iter <- 1e3
#' param1$burnin <- 1e2
#' param1$progress <- FALSE
#' param1$compute_posterior <- FALSE
#' res_bliss_chains <- fit_Bliss(data1,param1)
#' res_diagnostic <- diagnostics(res_bliss_chains$chains,param1)
#' cols <- c(2:5)
#' plot(res_diagnostic$hist_mu[[1]],border=0,col=cols[1],main="",xlab="",ylab="")
#' for(c in 2:3){
#'    plot(res_diagnostic$hist_mu[[c]],add=TRUE,border=0,col=cols[c])
#' }
#' plot(res_diagnostic$density_beta[[4]][[1]],main="",xlab="",ylab="")
#' for(c in 2:3){
#'    lines(res_diagnostic$density_beta[[4]][[c]],col=c)
#' }
#' }
diagnostics <- function(chains,param,progress=FALSE){
  if(progress) cat("Compute some diagnostics on the posterior sample. \n")
  # load objects
  n_chains  <- param[["n_chains"]]
  burnin    <- param[["burnin"]]
  K         <- param[["K"]]
  p         <- param[["p"]]
  iter      <- param[["iter"]]
  # load optional objects
  ts      <- param[["ts"]]
  l_ts    <- param[["l_ts"]]
  lag_max <- param[["lag_max"]]
  # Initialize the necessary unspecified objects XXXXX
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
#' plot_diagnostics
################################# ----
#' @description : Plot the diagnostics computed by \code{\link[=diagnostics]{diagnostics}}
#' @return a graph
#' @param res_diagnostics a list containing the result of \code{\link[=diagnostics]{diagnostics}} function.
#' @param param a list containing n_chains, the number of chains.
#' @param chain an integer, corresponding to the index of the chain to those diagnositcs have to be plotted.
#'                  If chain is NULL (defaut), the diagnostics of all the chains are plotted on
#'                 the same plots. (optional)
#' @param which_plot a string character indicating which parameter is of interest. The possible value : 'mu', 'sigma_sq' or 'beta'.
#' @param time a numerical value belonging to the "ts" vector (option
#'                 of \code{\link[=diagnostics]{diagnostics}} function). (optional)
#' @importFrom graphics hist axis layout lines matplot plot points text
#' @importFrom grDevices colorRampPalette
#' @importFrom utils data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 aes aes_string geom_histogram ggplot scale_fill_manual
#' @export
#' @examples
#' \donttest{
#' # Not run!
#' data(data1)
#' data(param1)
#' param1$n_chains <- 3
#' param1$iter <- 1e3
#' param1$burnin <- 1e2
#' param1$progress <- FALSE
#' param1$compute_posterior <- FALSE
#' res_bliss_chains <- fit_Bliss(data1,param1)
#' res_diagnostic <- diagnostics(res_bliss_chains$chains,param1)
#' plot_diagnostics(res_diagnostic,param1,which_plot="mu")
#' }
plot_diagnostics <- function(res_diagnostics,param,chain=NULL,which_plot=NULL,time=NULL){
 if(is.null(which_plot)) stop("Please specify which diagnostics you want : 'mu', 'sigma_sq' or 'beta' ?")
 if(which_plot == "beta"){
  if(is.null(time)){
   time <- sample(res_diagnostics$ts,1)
   cat("You didn't specify which time point have to be taken into account ",
       "for the diagnostics. As an example, here is the diagnostics for time = ",
       time,". \n",sep="")
  }
  if(!( time %in% res_diagnostics$ts)) stop("'time' must belong to 'ts'.")
 }

 # load objects
 if(which_plot == "mu"){
  DF_mu           <- res_diagnostics$DF_mu
  trace_mu        <- res_diagnostics$trace_mu
  autocorr_lag_mu <- res_diagnostics$autocorr_lag_mu
  hist_mu         <- res_diagnostics$hist_mu
  density_mu      <- res_diagnostics$density_mu
  ylim_mu         <- res_diagnostics$ylim_mu
  n_class <- length(hist_mu[[1]]$breaks)-1
 }
 if(which_plot == "sigma_sq"){
  DF_sigma           <- res_diagnostics$DF_sigma
  trace_sigma        <- res_diagnostics$trace_sigma
  autocorr_lag_sigma <- res_diagnostics$autocorr_lag_sigma
  hist_sigma         <- res_diagnostics$hist_sigma
  density_sigma      <- res_diagnostics$density_sigma
  ylim_sigma         <- res_diagnostics$ylim_sigma
  n_class <- length(hist_sigma[[1]]$breaks)-1
 }
 if(which_plot == "beta"){
  DF_beta           <- res_diagnostics$DF_beta
  trace_beta        <- res_diagnostics$trace_beta
  autocorr_lag_beta <- res_diagnostics$autocorr_lag_beta
  hist_beta         <- res_diagnostics$hist_beta
  density_beta      <- res_diagnostics$density_beta
  ylim_beta         <- res_diagnostics$ylim_beta
  breaks_beta       <- res_diagnostics$breaks_beta
  n_class <- length(hist_beta[[1]][[1]]$breaks)-1
 }

 # Initialize
 n_chains     <- param[["n_chains"]]
 if(n_chains < 5){
  cols_chains <- 1:n_chains+1
 } else cols_chains  <- colorRampPalette(brewer.pal(n_chains,"Spectral"))(n_chains)
 lags    <- res_diagnostics$lags
 ts      <- res_diagnostics$ts
 l_ts    <- length(ts)

 if(is.numeric(chain)){ # For a specific chain or ...
  if(which_plot == "mu"){
   # Trace of mu
   readline(paste("Press [enter] to plot the trace of mu for the chain ",chain,": ",sep=""))
   plot(trace_mu[,chain],type="l",lty=1,col=cols_chains[chain],xlab="Iterations",
        ylab="",main=paste("Trace of mu for the Chain ",chain,sep=""))
   # Autocorrelation of mu
   readline(paste("Press [enter] to plot the autocorrelation of mu for the chain ",chain,":",sep=""))
   plot(lags,autocorr_lag_mu[,chain],type="h",xlab="lag",ylab="correlation",
        col=cols_chains[chain],
        main=paste("Autocorrelation of mu for the Chain ",chain,sep=""))
   # Histogram of the empirical posterior sample of mu
   readline(paste("Press [enter] to plot the histogram of the posterior sample of mu for the chain",chain,":",sep=""))
   DF_mu_tmp <- DF_mu$obs[DF_mu$chain == paste("chain",1,sep="_")]
   hist(DF_mu_tmp,nclass=n_class,col=cols_chains[chain],border=0,xlab="",ylab="density",
        main=paste("Histogram of the posterior sample of mu \n for the chain ",chain,sep="" ))
  }
  if(which_plot == "sigma_sq"){
   # Trace of sigma_sq
   readline(paste("Press [enter] to plot the trace of sigma_sq for the chain" ,chain,":",sep=""))
   plot(trace_sigma[,chain],type="l",lty=1,col=cols_chains[chain],xlab="Iterations",
        ylab="",main=paste("Trace of sigma for the Chain ",chain,sep=""))
   # Autocorrelation of sigma_sq
   readline(paste("Press [enter] to plot the autocorrelation of sigma_sq for the chain",chain,":",sep=""))
   plot(lags,autocorr_lag_sigma[,chain],type="h",xlab="lag",ylab="correlation",
        col=cols_chains[chain],
        main=paste("Autocorrelation of sigma for the Chain ",chain,sep=""))
   # Histogram of the empirical posterior distribution of sigma_sq
   readline(paste("Press [enter] to plot the histogram of the posterior sample of sigma_sq for the chain",chain,":",sep=""))
   DF_sigma_tmp <- DF_sigma$obs[DF_sigma$chain == paste("chain",1,sep="_")]
   hist(DF_sigma_tmp,nclass=n_class,col=cols_chains[chain],border=0,xlab="",ylab="density",
        main=paste("Histogram of the posterior sample of sigma \n for the chain ",chain,sep="" ))
  }
  if(which_plot == "beta"){
   ts_index <- which(ts == time)
   # Trace of beta( time )
   readline(paste("Press [enter] to plot the trace of beta(",round(data$grid[time],2),") for the chain",chain,":",sep=""))
   plot(trace_beta[ts_index,,chain],type="l",lty=1,col=cols_chains[chain],xlab="Iterations",ylab="",
        main=paste("Trace of beta(",round(data$grid[time],2),") for the Chain ",chain,sep=""))
   # Autocorrelation of beta( time )
   readline(paste("Press [enter] to plot the autocorrelation of beta(",round(data$grid[time],2),") for the chain",chain,":",sep=""))
   plot(lags,autocorr_lag_beta[ts_index,chain,],type="h",xlab="lag",ylab="correlation",col=cols_chains[chain],
        main=paste("Autocorrelation of beta(",round(data$grid[time],2),") for the Chain ",chain,sep=""))
   # Histogram of the empirical posterior distribution of beta( time )
   readline(paste("Press [enter] to plot the histogram of the posterior sample of beta(",round(data$grid[time],2),") for the chain",chain,":",sep=""))
   DF_beta_tmp <- DF_beta[[ts_index]]$obs[DF_beta[[ts_index]]$chain == paste("chain",1,sep="_")]
   hist(DF_beta_tmp,nclass=n_class,col=cols_chains[chain],border=0,xlab="",ylab="density",
        main=paste("Histogram of the posterior sample of beta(",round(data$grid[time],2),") \n for the chain",chain,sep=""))
  }
 }else{ # ... or all the chain on the same plot.
  matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
  layout(matrix_layout)
  if(which_plot == "mu"){
   # Trace of mu
   readline(paste("Press [enter] to plot the trace of mu :",sep=""))
   matplot(trace_mu,type="l",lty=1,col=cols_chains,xlab="Iterations",ylab="",
           main="Trace of mu")
   plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
   points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
   text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],
        labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
   # Autocorrelation of mu
   readline(paste("Press [enter] to plot the autocorrelation of mu :",sep=""))
   plot(lags,type="n",ylim=range(autocorr_lag_mu),xlab="lag",ylab="correlation",
        main="Autocorrelation of mu ")
   for(j in 1:n_chains){
    lines(lags+j/(n_chains+1),autocorr_lag_mu[,j],type="h",col=cols_chains[j])
   }

   breaks_mu <- hist_mu[[1]]$breaks
   step_mu   <- diff( breaks_mu )[1]
   # Histograms of the empirical posterior distribution of mu
   readline(paste("Press [enter] to plot the histogram of the posterior sample of mu :",sep=""))
   # IS 19/02/2018
   print(ggplot(data=DF_mu, aes_string(x="obs", fill="chain")) +
          geom_histogram(binwidth=5*step_mu, colour="black", position="dodge") +
          scale_fill_manual(breaks=paste("chain_",1:n_chains,sep=""), values=cols_chains))

   readline(paste("Press [enter] to plot another histogram of the posterior sample of mu :",sep=""))
   matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
   layout(matrix_layout)
   plot(1,type="n",xlab="",main="mu for different chains",axes = FALSE,ylab="",
        xlim=range(breaks_mu),ylim=ylim_mu)
   axis(1) ; axis(2)
   for(j in 1:n_chains){
    hist_mu[[j]]$counts <- hist_mu[[j]]$density
    lines(hist_mu[[j]],border = 0,col = cols_chains[j],lty=2)
    lines(density_mu[[j]],col=cols_chains[j])
   }
   plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
   points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
   text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],
        labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
   layout(1)
  }
  if(which_plot == "sigma_sq"){
   # Trace of sigma_sq
   readline(paste("Press [enter] to plot the trace of sigma_sq :",sep=""))
   matplot(trace_sigma,type="l",lty=1,col=cols_chains,xlab="Iterations",ylab="",
           main="Trace of sigma")
   plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
   points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
   text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],
        labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
   layout(1)
   # Autocorrelation of sigma_sq
   readline(paste("Press [enter] to plot the autocorrelation of sigma_sq :",sep=""))
   plot(lags,type="n",ylim=range(autocorr_lag_sigma),xlab="lag",ylab="correlation",
        main="Autocorrelation of mu ")
   for(j in 1:n_chains){
    lines(lags+j/(n_chains+1),autocorr_lag_sigma[,j],type="h",col=cols_chains[j])
   }

   breaks_sigma <- hist_sigma[[1]]$breaks
   step_sigma   <- diff( breaks_sigma )[1]
   # Histograms of the empirical posterior distribution of sigma_sq
   readline(paste("Press [enter] to plot the histogram of the posterior sample of sigma_sq :",sep=""))
   # IS 19/02/2018
   print(ggplot(DF_sigma, aes_string(x="obs", fill="chain")) +
          geom_histogram(binwidth=5*step_sigma, colour="black", position="dodge") +
          scale_fill_manual(breaks=paste("chain_",1:n_chains,sep=""), values=cols_chains))

   readline(paste("Press [enter] to plot another histogram of the posterior sample of sigma_sq :",sep=""))
   matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
   layout(matrix_layout)
   plot(1,type="n",xlab="",main="sigma for different chains",axes = FALSE,ylab="",
        xlim=range(breaks_sigma),ylim=ylim_sigma)
   axis(1) ; axis(2)
   for(j in 1:n_chains){
    hist_sigma[[j]]$counts <- hist_sigma[[j]]$density
    lines(hist_sigma[[j]],border = 0,col = cols_chains[j],lty=2)
    lines(density_sigma[[j]],col=cols_chains[j])
   }
   plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
   points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
   text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],
        labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
   layout(1)

  }
  if(which_plot == "beta"){
   ts_index <- which(ts == time)
   # Trace of beta( time )
   readline(paste("Press [enter] to plot the trace of beta(",round(data$grid[time],2),"):",sep=""))
   matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
   layout(matrix_layout)
   matplot(trace_beta[ts_index,,],type="l",lty=1,col=cols_chains,xlab="Iterations",ylab="",
           main=paste("Trace of beta(",round(data$grid[time],2),")",sep=""))
   plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
   points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
   text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],
        labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
   layout(1)

   # Autocorrelation of beta( time )
   readline(paste("Press [enter] to plot the autocorrelation of beta(",round(data$grid[time],2),"):",sep=""))
   matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
   layout(matrix_layout)
   plot(lags,type="n",ylim=range(autocorr_lag_beta[ts_index,,]),xlab="lag",ylab="correlation",
        main=paste("Autocorrelation of beta(",round(data$grid[time],2),") ",sep=""))
   for(j in 1:n_chains){
    lines(lags+j/(n_chains+1),autocorr_lag_beta[ts_index,j,],type="h",col=cols_chains[j])
   }
   plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
   points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
   text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],
        labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
   layout(1)

   step_beta   <- diff( breaks_beta[,ts_index] )[1]
   # Histograms of the empirical posterior distribution of beta( time )
   readline(paste("Press [enter] to plot the histogram of the posterior sample of beta(",round(data$grid[time],2),"):",sep=""))
   # IS 19/02/2018
   print(ggplot(DF_beta[[ts_index]], aes_string(x="obs", fill="chain")) +
          geom_histogram(binwidth=5*step_beta, colour="black", position="dodge") +
          scale_fill_manual(breaks=paste("chain_",1:n_chains,sep=""), values=cols_chains))

   readline(paste("Press [enter] to plot another trace of beta(",round(data$grid[time],2),"):",sep=""))
   matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
   layout(matrix_layout)
   plot(1,type="n",xlab="",main=paste("beta(",round(data$beta_function[time],2),
                                      ") for different chains",sep=""),axes = FALSE,ylab="",
        xlim=range(breaks_beta[,ts_index]),ylim=ylim_beta[,ts_index])
   axis(1) ; axis(2)
   for(j in 1:n_chains){
    hist_beta[[ts_index]][[j]]$counts <- hist_beta[[ts_index]][[j]]$density
    lines(hist_beta[[ts_index]][[j]],border = 0,col = cols_chains[j],lty=2)
    lines(density_beta[[ts_index]][[j]],col=cols_chains[j])
   }
   plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
   points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
   text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],
        labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
   layout(1)

  }
 }
}
################################# ----
#' dposterior
################################# ----
#' @description Compute the posterior density for a given parameter set. XXX (non normalisé)
#' @param posterior_sample a list given by the "Bliss_Gibbs_Sampler" function.
#' @param data a list containing
#' \describe{
#' \item{x}{a list of matrices, the functions x_qi(t) observed at time points given by grids,}
#' \item{y}{a numerical vector, the outcome values y_i.}
#' }
#' @param Q an integer, the number of covariates
#' @param theta a matrix or a vector which contains the parameter set. XXXX
#' @details If the option theta is NULL, the posterior density is computed for
#'          the MCMC sample given in the "res.Gibbs_Sampler" object.
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
#' res_poste <- dposterior(res_bliss1$posterior_sample,data1,Q=2)
dposterior <- function(posterior_sample,data,Q,theta=NULL){
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

 y <- data$y
 potential_intervals  <- posterior_sample$param$potential_intervals
 potential_intervals_dims <- list()
 for(q in 1:Q){
  potential_intervals_dims[[q]] <- c(ncol(data$x[[q]]),
                                     posterior_sample$param$l_values_length[[q]],
                                     length(data$y))
 }

 res <- dposterior_cpp(rposterior,y,N,K,potential_intervals,potential_intervals_dims,
                       posterior_sample$param$l_values_length,Q)
 colnames(res) <- c("posterior density","log posterior density",
                    "likelihood","log likelihood",
                    "prior density","log prior density")
 return(res)
}

################################# ----
#' compute_chains_info
################################# ----
#' compute_chains_info
#'
#' @param x numeric
#'
#' @return numeric
#' @export
#'
#' @examples
#' # compute x
compute_chains_info <- function(x){
  x
}
