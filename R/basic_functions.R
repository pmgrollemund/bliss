################################# ----
#' dposterior
################################# ----
#' @description Compute the posterior distribution of the model parameters
#' @param res.Gibbs_Sampler a list corresponding to the object returns by the "Bliss_Gibbs_Sampler_multiple" function.
#' @param data a list containing
#' \describe{
#' \item{1)}{the number of covariates,}
#' \item{2)}{the functions x_qi(t) observed at grids of time points,}
#' \item{3)}{the outcome values y_i.}
#' }
#' @param theta a matrix or a vector which contains the parameter value(s) for which the posterior density is claimed.
#' @details If the option theta is null, the posterior density is computed for the MCMC sample given in the "res.Gibbs_Sampler" object.
#'
#' @return return the posterior density and the posterior log density for a given parameters set.
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
#' plot_step_function
################################# ----
#' @description a suitable representation of the Bliss estimate of the coefficient function.
#' @param extended_grid a vector. Similar to the grid of the functional
#'                 covariate but it is adapted to provide suitable representation.
#' @param fct the numerical vector to plot.
#' @param bound a boolean value. If bound is TRUE, the plot of fct is
#'                 one line. Otherwise, several lines are used to plot fct.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @importFrom grDevices gray.colors
#' @importFrom graphics abline segments
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' ### Plot the BLiss estimate on a suitable grid
#' plot_step_function(res_Bliss_mult$param$grids[[1]],
#'                    res_Bliss_mult$Bliss_estimate[[1]],lwd=2,bound=FALSE)
plot_step_function <- function(extended_grid,fct,bound=FALSE,...){
  ylim <- range(fct)
  plot(extended_grid,extended_grid,type="n",ylim=ylim,...)
  lines_step_function(extended_grid,fct,bound=bound,...)
}

################################# ----
#' lines_step_function
################################# ----
#' @description a suitable representation of the Bliss estimate of the coefficient function.
#' @param extended_grid a vector. Similar to the grid of the functional
#'                 covariate but it is adapted to provide suitable representation.
#' @param fct the numerical vector to plot.
#' @param bound a boolean value. If bound is TRUE, the plot of fct is
#'                 one line. Otherwise, several lines are used to plot fct.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @export
#' @examples
#' ### Plot the BLiss estimate on a suitable grid
lines_step_function <- function(extended_grid,fct,bound=FALSE,...){
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
#' support_estimation
################################# ----
#' @description Compute the estimate of the support of the coefficient function.
#' @param beta_functions a matrix for which each row is an evaluation of a function on a grid.
#' @return a list
#' \describe{
#' \item{alpha_t}{a numerical vector. For each point t, alpha(t) is the probability that beta(t) is non-null}
#' \item{estimate}{a matrix given by the function "interval_detection", which give the estimation of
#' the support of the coefficient function}
#' }
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' res.support <- support_estimation(res_Bliss_mult$beta_functions[[1]])
#'  ### The estimate
#'  res.support$estimate
#'  ### Plot the result
#'  grid <- res_Bliss_mult$param$grids[[1]]
#'  plot(grid,res.support$alpha_t,ylim=c(0,1),type="l",xlab="",ylab="")
#'  for(k in 1:nrow(res.support$estimate)){
#'     segments(grid[res.support$estimate[k,1]],0.5,
#'            grid[res.support$estimate[k,2]],0.5,lwd=2,col=2)
#'    points(grid[res.support$estimate[k,1]],0.5,pch="|",lwd=2,col=2)
#'    points(grid[res.support$estimate[k,2]],0.5,pch="|",lwd=2,col=2)
#'  }
#'  abline(h=0.5,col=2,lty=2)
support_estimation <- function(beta_functions){
  alpha_t <- apply(beta_functions,2, function(vec) sum(vec != 0)/length(vec))
  tmp  <- rep(0,ncol(beta_functions))
  tmp2 <- which(alpha_t >= 0.5)
  tmp[tmp2] <- 1
  estimate <- interval_detection(tmp)
  estimate <- estimate[ estimate[,3]==1 , -3]

  return(list(alpha_t=alpha_t,estimate=estimate))
}


################################# ----
#' plot_impact_mult
################################# ----
#' @description provide a graphical representation of the data with a focus on the detected periods with the Bliss method.
#' @param estimate a numerical vector, the estimate of the coefficient function.
#' @param data a list containing y, x and grid.
#' @param q an integer, the index of the functional covariate to consider.
#' @param centered a logical value. If TRUE, the curve x are centered.
#' @param cols a numerical vector of colors.
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' plot_impact_mult(res_Bliss_mult$Bliss_estimate[[1]],data1)
#' plot_impact_mult(res_Bliss_mult$Bliss_estimate[[1]],data1,centered=FALSE)
#' }
plot_impact_mult <- function(estimate,data,q=1,centered=TRUE,cols=NULL){
  # load some objects
  x <- data$x_mult[[q]]
  y <- data$y
  grid <- data$grids[[q]]
  grid0 <- data$grids[[q]] - 0.5*diff(data$grids[[q]])[1]
  grid0 <- c(grid0 , max(data$grids[[q]]) + 0.5*diff(data$grids[[q]])[1]  )

  x_centered <- apply(x,2,function(vec) vec - mean(vec))

  grid_new <- rep(0,2*length(grid)+1)
  grid_new[seq(2,length(grid_new),by=2)] <- grid
  grid_new[seq(3,length(grid_new),by=2)] <- grid + 0.5*diff(grid)[1]
  grid_new[1] <- grid[1] - 0.5*diff(grid)[1]

  x_new <- matrix(NA,nrow=nrow(x),ncol=2*ncol(x)+1)
  x_new[,seq(2,ncol(x_new),by=2)] <- as.matrix(x)
  for(i in seq(3,ncol(x_new)-1,by=2))
    x_new[,i] <- 0.5*(x_new[,i-1]+x_new[,i+1])

  x_centered_new <- matrix(NA,nrow=nrow(x),ncol=2*ncol(x)+1)
  x_centered_new[,seq(2,ncol(x_centered_new),by=2)] <- as.matrix(x_centered)
  for(i in seq(3,ncol(x_centered_new)-1,by=2))
    x_centered_new[,i] <- 0.5*(x_centered_new[,i-1]+x_centered_new[,i+1])


  intervals <- interval_detection(estimate)
  intervals$value[3] <- 0
  intervals$end[4] <- 13
  intervals$begin[5] <- 14
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
#' interval_detection
################################# ----
#' @description Cut a numerical vector in several homogenous pieces.
#' @return a matrix with 3 columns : "begin", "end" and "value". The two first
#'           columns define the begin and the end of the pieces and the third
#'           gives the mean values of each piece.
#' @param vec a numerical vector.
#' @param smooth a logical value, indicates if vec is smooth or not, default FALSE.
#' @param q a two-vector, if smooth is TRUE, it is used to detect a rupture.
#' @importFrom stats qnorm sd
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' intervals <- interval_detection(res_Bliss_mult$Bliss_estimate[[1]])
#' par(mfrow=c(2,1))
#' plot(data1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="s")
#' plot(data1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="n")
#' for(k in 1:nrow(intervals)){
#'   segments(intervals[k,1],intervals[k,3],
#'           intervals[k,2],intervals[k,3],col=2,lwd=2)
#' }
#' par(mfrow=c(1,1))
#' intervals <- interval_detection(res_Bliss_mult$res.Simulated_Annealing[[1]]$posterior_expe,
#'                                smooth=TRUE)
#' plot(data1$grids[[1]],res_Bliss_mult$res.Simulated_Annealing[[1]]$posterior_expe,type="l")
#' for(k in 1:nrow(intervals)){
#'    segments(intervals[k,1],intervals[k,3],
#'             intervals[k,2],intervals[k,3],col=2,lwd=2)
#' }
interval_detection <- function(vec, smooth=FALSE, q = c(0.1, 0.9)){
  intervals <- data.frame()
  begin <- 1
  value <- vec[1]
  add   <- diff(vec)
  tresh <- qnorm(q, mean(add), sd(add))
  for (i in 2:length(vec)) {
    if(smooth){
      if( !(add[i - 1] %between% tresh) ){
        end <- i - 1
        intervals <- rbind(intervals,
                           c(begin, end, mean(vec[begin:end])))
        begin <- i
      }
    }
    if(!smooth){
      if (vec[i] != value) {
        end <- i - 1
        intervals <- rbind(intervals,
                           c(begin, end, value))
        begin <- i
        value <- vec[i]
      }
    }

  }
  end <- i
  if(smooth){
    intervals <- rbind(intervals, c(begin, end, mean(vec[begin:end])))
    colnames(intervals) <- c("begin", "end", "mean")
  }
  if(!smooth){
    intervals <- rbind(intervals, c(begin, end, value))
    colnames(intervals) <- c("begin", "end", "value")
  }
  return(intervals)
}

################################# ----
#' finer_grid
################################# ----
#' @description Compute a curve on a new finer gird.
#' @return a numerical vector, the curve evaluated on the new grid.
#' @param curve a numerical vector.
#' @param grid a numerical vector, the former grid.
#' @param new_grid a numerical vector.
#' @details This is nasty code.
#' @export
#' @examples
#' grid <- seq(0,1,l=1e1)
#' new_grid <- seq(0,1,l=1e2)
#' curve <- 3*grid^2 + sin(grid*2*pi)
#' plot(grid,curve,type="o")
#' lines(new_grid,finer_grid(curve,grid,new_grid),type="o",col=makeTransparent(2))
finer_grid <- function(curve,grid,new_grid){
  res <- rep(0,length(new_grid))
  for(i in 1:(length(grid)-1)){
    index <- new_grid %between% grid[i:(i+1)]
    res[index] <- curve[i] + 0:(sum(index)-1) /
      (sum(index)-1) * (curve[i+1]-curve[i])
  }
  return(res)
}

################################# ----
#' prior_l
################################# ----
#' @description Compute the probability function of the Gamma prior on l.
#' @return a numerical vector, which is the prability function on "grid_l".
#' @param m a positive value, the mean of the Gamma prior.
#' @param s a nonnegative value, the standard deviation of the Gamma prior.
#' @param grid_l a numerical value, the discrete support of the parameter l.
#' @importFrom stats pgamma
#' @export
#' @examples
#' grid_l <- seq(0,5,l=100)
#' f <- prior_l(3,1,grid_l)
#' plot(grid_l,f,type="h",xlab="",ylab="")
#' f <- prior_l(1,0.5,grid_l)
#' plot(grid_l,f,type="h",xlab="",ylab="")
#' f <- prior_l(1,5,grid_l)
#' plot(grid_l,f,type="h",xlab="",ylab="")
prior_l <- function(m,s,grid_l){
  # Compute the scale and the rate
  alpha <- m^2/s^2
  beta  <- m/s^2

  # Compute the probability function on the grid
  step <- diff(grid_l)[1] / 2
  probs <- pgamma(grid_l + step ,shape=alpha,rate=beta) -
           pgamma(grid_l - step ,shape=alpha,rate=beta)

  return(probs)
}

################################# ----
#' makeTransparent
################################# ----
#' @description Make transparent some color.
#' @return a character string coding for a color.
#' @param ... must be a color. Numerical or character string.
#' @param alpha a numerical value between 0 and 1, corresponding to the
#'               transparency of the returned color.
#' @seealso Thanks to Ricardo Oliveros-Ramos for the function on the web :
#'               \url{http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color}
#' @importFrom grDevices col2rgb rgb
#' @export
#' @examples
#' cols <- makeTransparent(2:6,0.3)
#' res_hist <- hist(rnorm(1e4,0,1),nclass=2e2,border = 0,col = cols[1],
#'                  xlim=c(-2,6),xlab="",main = "")
#' hist(rnorm(1e4,1,1),border = 0,col = cols[2],add=TRUE,nclass=2e2)
#' hist(rnorm(1e4,2,1),border = 0,col = cols[3],add=TRUE,nclass=2e2)
#' hist(rnorm(1e4,3,1),border = 0,col = cols[4],add=TRUE,nclass=2e2)
#' hist(rnorm(1e4,4,1),border = 0,col = cols[5],add=TRUE,nclass=2e2)
makeTransparent <- function(..., alpha=0.5) {

  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")

  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)

  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }

  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)

  return(newColor)
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
#' between
################################# ----
#' @description check if a value is in an interval.
#' @return a logical value.
#' @param value a numerical value.
#' @param interval a numerical vector of lenght 2 : (low,up).
#' @export
#' @examples
#' 1 %between% c(0,2)
#' 2 %between% c(0,2)
#' 3 %between% c(0,2)
"%between%" <- function(value,interval){
  (value >= interval[1]) & (value <= interval[2])
}

################################# ----
#' compute_beta_functions
################################# ----
#' @description Compute the function beta_i(t) for each iteration i of the Gibbs sample.
#' @return return a matrix. Each row is a function observed on the grid of time points.
#' @param res.Gibbs_Sampler a list (provided by the function Bliss_Gibbs_Sampler).
#' @param param (optional) a list containing
#' \describe{
#' \item{p}{the number of time points}
#' \item{K}{the number of intervals}
#' \item{grid}{the grid of time points}
#' \item{basis}{a character vector}
#' }
#' @export
#' @examples
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' beta_functions <- compute_beta_functions_mult(res_Bliss_mult$res.Gibbs_Sampler,param1)
#' indexes <- sample(nrow(beta_functions[[1]]),1e2,replace=FALSE)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
#' matplot(param1$grids[[1]],t(beta_functions[[1]][indexes,]),type="l",lty=1,col=cols)
compute_beta_functions_mult <- function(res.Gibbs_Sampler,param){
  cat("Compute the functions beta_i. \n")
  # Initialize parameters
  K     <- param[["K"]]
  grids <- param$grids
  p     <- param[["p"]]

  if(is.null(p)) p <- sapply(grids,length)

  # Initialize parameters
  basis <- param[["basis"]]
  if(is.null(basis)) basis <- rep("uniform",length(K))

  # Compute the functions beta_i for each iteration i of the Gibbs Sampler and for each covariable
  beta_functions <- list()
  count <- 0
  for(q in 1:length(K)){
    trace_tmp <- res.Gibbs_Sampler$trace[,(1+count):(count+3*K[q])]

    beta_functions[[q]] <- compute_beta_functions_cpp(trace_tmp,p[q],K[q],grids[[q]],
                                                      basis[q],res.Gibbs_Sampler$param$scale_ml[[q]])
    count <- count + 3*K[q]
  }
  # Return the functions
  return(beta_functions)
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

################################# ----
#' density_estimation
################################# ----
#' @description compute a graphical representation of the posterior distribution of beta.
#' @details The sample is thinned in order to reduce the number of points and
#'           so the time of the computation of the function \code{\link[=kde2d]{kde2d}}.
#' @return the result of the \code{\link[=kde2d]{kde2d}} function, i.e. an estimate of the
#'         posterior density on a two-dimensional grid.
#' @param beta_functions a list (provided by the function compute_beta_functions).
#' @param param (optional), a list containing
#' \describe{
#'  \item{burnin}{an interger which is the number of iteration to drop (optional)}
#'  \item{thin}{a numerical value (optional)}
#'  \item{lims.kde}{a numerical vector (yl,yu) (optional)}
#'  \item{n}{an integer related to the precision of the heat map (optional)}
#'  \item{h1}{a numerical value which is the bandwidth of the kernel density estimation for the t-axis (optional)}
#'  \item{new_grid}{a numerical vector.}
#' }
#' @importFrom MASS bandwidth.nrd kde2d
#' @export
#' @examples
#' \donttest{
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1,density=TRUE)
#' q <- 1
#' diff_grid <- diff(param1$grids[[q]])[1]
#' param1$grids2[[q]] <- c(param1$grids[[q]]-diff_grid/2,
#'                        tail(param1$grids[[q]],1)+diff_grid/2)
#' param1$xlim[[q]] <- range(param1$grids2[[q]])
#' param_density<-list(grid= param1$grids[[q]],
#'                     iter= param1$iter,
#'                     p   = param1[["p"]][q],
#'                     n        = param1[["n"]],
#'                     thin     = 10,
#'                     burnin   = param1[["burnin"]],
#'                     lims.kde = param1$lims.kde[[q]],
#'                     h1       = param1$h1,
#'                     new_grid = param1[["new_grid"]],
#'                     xlim = range(param1$grids[[q]]) + c(-diff_grid,diff_grid),
#'                     display = FALSE
#' )
#' density_estimate <- density_estimation(res_Bliss_mult$beta_functions[[1]],param_density)
#' image(density_estimate$res.kde2d,col=rev(heat.colors(100)))
#' }
density_estimation <- function(beta_functions,param){
  cat("Compute the estimation of the posterior density.\n")
  # Initialize
  grid <- param$grid
  iter <- param[["iter"]]
  p    <- param$p

  # load optional objects
  n        <- param[["n"]]
  thin     <- param$thin
  burnin   <- param[["burnin"]]
  lims.kde <- param$lims.kde
  h1       <- param$h1
  new_grid <- param[["new_grid"]]

  # Initialize the necessary unspecified objects
  max_points <- 1e5
  if(!is.null(new_grid)) p      <- length(new_grid)
  if(is.null(n))         n      <- 512
  if(is.null(burnin))    burnin <- floor(iter/5)
  if(is.null(thin))      thin   <- floor((iter-burnin)*p/max_points)

  # Check if the burnin isn't too large.
  if(2*burnin > iter){
    burnin <- floor(iter/5)
    cat("\t Burnin is too large. New burnin : ",burnin,".\n")
  }

  # Thin the sample of beta functions
  thin_min   <- max(1,floor((iter-burnin)*p/max_points))
  if(thin <  thin_min){
    cat("\t 'thin = ",thin,"' is too small. Now, thin = ",
        thin_min,".\n",sep="")
    thin <- thin_min
  }
  cat("\t Thin the sample.\n")
  beta_functions <- beta_functions[seq(1+burnin,iter,by=thin),]


  # Compute the functions beta_i on the new grid (if claimed).
  if(!is.null(new_grid)){
    beta_functions_save <- beta_functions
    beta_functions <- matrix(0,nrow(beta_functions),p)
    cat("Compute the coefficient functions on the new grid.\n")
    for(i in 1:nrow(beta_functions)){
      beta_functions[i,] <- finer_grid(beta_functions_save[i,],
                                       grid,new_grid)
    }
    param$old_grid <- grid
    param$grid     <- new_grid
    grid           <- new_grid
  }

  # Filter the functions with a range too large which make a
  # noninterpretable graphical representation
  cat("\t Drop the extreme values.\n")
  t_beta <- rep(grid,nrow(beta_functions))
  y_beta <- as.vector(t(beta_functions))

  t_beta <- t_beta[y_beta %between% quantile(y_beta,c(0.01,0.99))]
  y_beta <- y_beta[y_beta %between% quantile(y_beta,c(0.01,0.99))]

  # If a window is given to plot the posterior density.
  if(!is.null(lims.kde)){
    index   <- which(y_beta %between% lims.kde)
    t_beta  <- t_beta[index]
    y_beta  <- y_beta[index]
  }

  # If there are too much of y_beta=0 (more than 50%), the default bandwidth of the
  # density estimation of the kde function is 0, using the defaut function
  # bw.nrd0. Indeed, with this function the bandwidth is :
  #     bw = 0.9 * min( \hat{\sigma} , IQR /1.34 ) * n^{-1/5}
  # So, if there are to much of y_beta=0, we choose to reduce them to 25%.
  if(sum(y_beta==0)/length(y_beta) > 0.25){
    cat("\t Bandwidth untractrable. Some points are dropped.\n")
    res.table <- table(t_beta[y_beta==0])
    # remove all the points which y=0
    t_beta <- t_beta[y_beta!=0]
    y_beta <- y_beta[y_beta!=0]
    Toadd <- length(y_beta)*0.25/(1-0.25)

    # Add 25% of them for some correct t_beta
    t_beta_0 <- sample(as.numeric(names(res.table)),Toadd,prob=res.table,replace=T)
    y_beta_0 <- rep(0,Toadd)
    t_beta <- c(t_beta,t_beta_0)
    y_beta <- c(y_beta,y_beta_0)
  }

  cat("\t Perform the 'kde2d' function.\n")
  h2 <- bandwidth.nrd(y_beta)
  if(is.null(h1)) h1 <- bandwidth.nrd(t_beta)
  points <- cbind(t_beta,y_beta)

  lims.kde <- c(range(t_beta),range(y_beta))
  if(!is.null(param$xlim))
    lims.kde[1:2] <- param$xlim
  if(!is.null(param$ylim))
    lims.kde[3:4] <- param$ylim

  res.kde2d <- kde2d(x=t_beta,y=y_beta,lims=lims.kde,
                     n=n,h=c(h1,h2))
  return(list(beta_functions = beta_functions,
              points         = points,
              thin           = thin,
              res.kde2d      = res.kde2d))
  cat("\t Done.\n")
}

################################# ----
#' Fourier_basis_build
################################# ----
#' @description define a Fourier basis to simulate functions x_i(t).
#' @return a matrix.
#' @param grid a numerical vector.
#' @param dim a numerical value. It corresponds to dim(basis)/2.
#' @param per a numerical value which corresponds to the period of the sine and cosine functions.
#' @details see the \code{\link[=sim_functions]{sim_functions}} function.
#' @export
#' @examples
#' # see the sim_functions() function.
Fourier_basis_build <- function(grid,dim,per=2*pi){
  sapply(grid,function(x) c(cos(2*pi*x*(1:dim)/per),sin(2*pi*x*(1:dim)/per) )  )
}

################################# ----
#' random_walk
################################# ----
#' @description compute a random walk. (gaussian)
#' @return a matrix where each row is a random walk.
#' @param n an integer, the number of random walk.
#' @param p an integer, the length of the random walks.
#' @param mu a numerical vector, the average random walk.
#' @param sigma a numerical value which is the standard deviation of the
#'                 gaussian distribution used to compute the random walk.
#' @param start a numerical vector which is the initial value of
#'                 the random walks. (optional)
#' @details see the \code{\link[=sim_functions]{sim_functions}} function.
#' @importFrom stats rnorm
#' @export
#' @examples
#' # see the sim_functions() function.
random_walk <- function(n,p,mu,sigma,start=rep(0,n)){
  res <- matrix(0,n,p)
  for(i in 1:n){
    add     <- rnorm(p,mu,sigma)
    res[i,] <- cumsum(add)
    res[i,] <- start[i] + res[i,]
  }
  res <- start + res
  return(res)
}

################################# ----
#' sigmoid
################################# ----
#' @description compute a sigmoid function.
#' @details used to simulate a coefficient function or functions x_i(t).
#' @return a numerical vector.
#' @param x a numerical vector, a grid of points.
#' @param asym the value of the asymptote of the sigmoid function. (optional)
#' @param v a numerical value which is related to the slope at the origin. (optional)
#' @details see the function \code{\link[=sim_functions]{sim_functions}}.
#' @export
#' @examples
#' ## Test 1 :
#' x <- seq(-7,7,0.1)
#' y <- sigmoid(x)
#' plot(x,y,type="l",main="Sigmoid function")
#' ## Test 2 :
#' x  <- seq(-7,7,0.1)
#' y  <- sigmoid(x)
#' y2 <- sigmoid(x,asym=0.5)
#' y3 <- sigmoid(x,v   =  5)
#' plot(x,y,type="l",main="Other sigmoid functions")
#' lines(x,y2,col=2)
#' lines(x,y3,col=3)
sigmoid <- function(x,asym=1,v=1){
  (asym^-1 + exp(-v*x))^-1
}

################################# ----
#' sigmoid_sharp
################################# ----
#' @description compute a sharp function from the sigmoid function
#' @details used to simulate a coefficient function or functions x_i(t).
#' @return a numerical vector.
#' @param x a numerical vector, a grid of points.
#' @param loc a numerical value, the instant of the sharp. (optional)
#' @param ... Arguments to be passed to the function sigmoid. (optional)
#' @details see the function \code{\link[=sim_functions]{sim_functions}}.
#' @export
#' @examples
#' ## Test 1 :
#' x <- seq(-7,7,0.1)
#' y <- sigmoid_sharp(x)
#' plot(x,y,type="l",main="Sharp sigmoid")
#' ## Test 2 :
#' x  <- seq(-7,7,0.1)
#' y  <- sigmoid_sharp(x,loc=3)
#' y2 <- sigmoid_sharp(x,loc=3,asym=0.5)
#' y3 <- sigmoid_sharp(x,loc=3,v   =  5)
#' plot(x,y,type="l",main="Other sharp sigmoids")
#' lines(x,y2,col=2)
#' lines(x,y3,col=3)
sigmoid_sharp <- function(x,loc=0,...){
  # 4 should be replace by (a+1)^2 such that the maximum of the curve
  # provided by sigmoid_sharp is 1.
  4*(sigmoid(x-loc,...) * sigmoid(-x+loc,...))
}
