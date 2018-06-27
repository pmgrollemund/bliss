################################# ----
#' compute_beta_sample
################################# ----
#' @description Compute the coefficient function for each iteration of the 
#'              Gibbs sample.
#' @return return a matrix. Each row is a function evaluated on the grid of 
#'         time points.
#' @param res.Gibbs_Sampler a list provided by the function Bliss_Gibbs_Sampler.
#' @param param a list containing: (optional)
#' \describe{
#' \item{p}{XXXXXX}
#' \item{K}{a vector of integers, corresponding to the numbers of intervals for
#'       each covariate.}
#' \item{grids}{a numerical vector, the observation time points.}
#' \item{basis}{a vector of characters among : "uniform" (default),
#'       "epanechnikov", "gauss" and "triangular" which correspond to
#'       different basis functions to expand the coefficient function and the
#'       functional covariates (optional)}
#' }
#' @param progress a logical value. If TRUE, the algorithm progress is displayed.
#'         (optional)
#' @export
#' @examples
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' beta_sample <- compute_beta_sample_mult(res_Bliss_mult$res.Gibbs_Sampler,param1)
#' indexes <- sample(nrow(beta_sample[[1]]),1e2,replace=FALSE)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
#' matplot(param1$grids[[1]],t(beta_sample[[1]][indexes,]),type="l",lty=1,col=cols)
compute_beta_sample <- function(res.Gibbs_Sampler,param,progress=FALSE){
 if(progress) cat("Compute the coefficient function posterior sample. \n")
 # Initialize parameters
 K     <- param[["K"]]
 grids <- param[["grids"]]
 p     <- param[["p"]]
 if(is.null(p)) p <- sapply(grids,length)
 basis <- param[["basis"]]
 if(is.null(basis)) basis <- rep("uniform",length(K))
 
 # Compute the coefficient function for each iteration of the Gibbs Sampler 
 # and for each covariable
 beta_sample <- list()
 count <- 0
 for(q in 1:length(K)){
  trace_tmp <- res.Gibbs_Sampler$trace[,(1+count):(count+3*K[q])]
  
  beta_sample[[q]] <- compute_beta_sample_cpp(trace_tmp,p[q],K[q],grids[[q]],
                                              basis[q],
                                              res.Gibbs_Sampler$param$normalization_values[[q]]) # XXXXX
  count <- count + 3*K[q]
 }
 
 return(beta_sample)
}
################################# ----
#' compute_beta_posterior_density
################################# ----
#' @description Compute a graphical representation of the marginal posterior
#'              distributions of beta(t) for each t.
#' @details The sample is thinned in order to reduce the correlation and
#'           so the time of the computation of the function \code{\link[=kde2d]{kde2d}}.
#' @return An approximation of the posterior density on a two-dimensional grid.
#'         (corresponds to the result of the \code{\link[=kde2d]{kde2d}} function)
#' @param beta_sample a list (given by the compute_beta function).
#' @param param (optional), a list containing:
#' \describe{
#' \item{grid}{a numerical vector, the observation time points.}
#' \item{burnin}{an integer, the number of iteration to drop from the Gibbs 
#'       sample. (optional)}
#' \item{thin}{an integer, used to thin the Gibbs sample to compute an
#'       approximation of the posterior density of beta(t). (optional)}
#'  \item{lims.kde}{a numerical vector (yl,yu) XXXXXXXXXXX (optional)}
#'  \item{N}{an integer related to the precision of the approximation. (optional)}
#'  \item{new_grid}{a numerical vector.} 
#' }
#' @param progress a logical value. If TRUE, the algorithm progress is displayed.
#'         (optional)
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
#' density_estimate <- density_estimation(res_Bliss_mult$beta_sample[[1]],param_density)
#' image(density_estimate$res.kde2d,col=rev(heat.colors(100)))
#' }
compute_beta_posterior_density <- function(beta_sample,param,progress=FALSE){
 if(progress) 
  cat("Compute an approximation of the posterior density of the coefficient function.\n")
 # load optional objects
 grid <- param[["grid"]]
 iter <- param[["iter"]]
 p    <- param[["p"]]
 N        <- param[["N"]]
 thin     <- param[["thin"]]
 burnin   <- param[["burnin"]]
 lims.kde <- param[["lims.kde"]]
 new_grid <- param[["new_grid"]] # XXXXX
 
 # Initialize the necessary unspecified objects
 max_points <- 1e5 
 if(!is.null(new_grid)) p      <- length(new_grid)  # XXXXX
 if(is.null(N))         N      <- 512
 if(is.null(burnin))    burnin <- floor(iter/5)  # XXXXX
 if(is.null(thin))      thin   <- floor((iter-burnin)*p/max_points)
 
 # Check if the burnin isn't too large.
 if(2*burnin > iter){
  burnin <- floor(iter/5)
  if(progress) 
   cat("\t Burnin is too large. New burnin : ",burnin,".\n")
 }
 
 # Thin the posterior sample 
 thin_min   <- max(1,floor((iter-burnin)*p/max_points))
 if(thin <  thin_min){
  if(progress) 
   cat("\t 'thin = ",thin,"' is too small. Now, thin = ",
      thin_min,".\n",sep="")
  thin <- thin_min
 }
 if(progress) 
  cat("\t Thin the sample.\n")
 beta_sample <- beta_sample[seq(1+burnin,iter,by=thin),]
 
 
 # Compute the coefficient function on the new grid (if required).
 if(!is.null(new_grid)){
  beta_sample_save <- beta_sample
  beta_sample <- matrix(0,nrow(beta_sample),p)
  if(progress) 
   cat("Compute the coefficient functions on the new grid.\n")
  for(i in 1:nrow(beta_sample)){
   beta_sample[i,] <- finer_grid(beta_sample_save[i,],grid,new_grid)
  }
  param$old_grid <- grid
  param$grid     <- new_grid
  param$new_grid <- NULL # PMG 22/06/18
  grid           <- new_grid
 }
 
 # Remove the functions with a range too large which leads to a noninterpretable
 # graphical representation
 if(progress) 
  cat("\t Drop the extreme values.\n") # XXXXXXX
 beta_x <- rep(grid,nrow(beta_sample))
 beta_y <- as.vector(t(beta_sample))
 
 beta_x <- beta_x[beta_y %between% quantile(beta_y,c(0.01,0.99))]
 beta_y <- beta_y[beta_y %between% quantile(beta_y,c(0.01,0.99))]
 
 # If a window is given to plot the posterior density.
 if(!is.null(lims.kde)){ # XXXXXXX # XXXXXXX
  index   <- which(beta_y %between% lims.kde)
  beta_x  <- beta_x[index]
  beta_y  <- beta_y[index]
 }
 
 ############ est-ce vraiment utile ? (surtout avec bandwidth.nrd ?)
 # If there are too much of beta_y=0 (more than 50%), the default bandwidth of the
 # density estimation of the kde function is 0, using the defaut function
 # bw.nrd0. Indeed, with this function the bandwidth is :
 #     bw = 0.9 * min( \hat{\sigma} , IQR /1.34 ) * n^{-1/5}
 # So, if there are to much of beta_y=0, we choose to reduce them to 25%.
 if(sum(beta_y==0)/length(beta_y) > 0.25){
  if(progress) 
   cat("\t Bandwidth untractrable. Some points are dropped.\n")
  res.table <- table(beta_x[beta_y==0])
  # remove all the points which y=0
  beta_x <- beta_x[beta_y!=0]
  beta_y <- beta_y[beta_y!=0]
  Toadd <- length(beta_y)*0.25/(1-0.25)
  
  # Add 25% of them for some correct beta_x
  beta_x_0 <- sample(as.numeric(names(res.table)),Toadd,prob=res.table,replace=T)
  beta_y_0 <- rep(0,Toadd)
  beta_x <- c(beta_x,beta_x_0)
  beta_y <- c(beta_y,beta_y_0)
 }
 ############
 
 # Perform the kde2d function
 if(progress) 
  cat("\t Perform the 'kde2d' function.\n")
 h1 <- bandwidth.nrd(beta_x)
 h2 <- bandwidth.nrd(beta_y)
 points <- cbind(beta_x,beta_y)
 
 lims.kde <- c(range(beta_x),range(beta_y)) # XXXX
 if(!is.null(param$xlim)) )) # XXXX
  lims.kde[1:2] <- param$xlim
 if(!is.null(param$ylim)) )) # XXXX
  lims.kde[3:4] <- param$ylim
 
 res.kde2d <- kde2d(x=beta_x,y=beta_y,lims=lims.kde,
                    n=N,h=c(h1,h2))
 
 # What to return ? # PMG 22/06/18
 if(!is.null(param$old_grid)){
  new_beta_sample <- beta_sample
 }else{
  new_beta_sample <- NULL
 }
 if(param[["thin"]] != thin){
  new_thin <- thin
 }else{
  new_thin <- NULL
 }
 
 return(list(grid_t          = res.kde2d$x,
             grid_beta_t     = res.kde2d$y,
             density         = res.kde2d$z,
             new_beta_sample = beta_sample,
             new_thin        = new_thin
             ))
 if(progress) 
  cat("\t Done.\n")
}
################################# ----
#' between
################################# ----
#' @description Check if a number belong to an interval.
#' @return a logical value.
#' @param value a numerical value.
#' @param interval a numerical vector of lenght 2 : (lower,upper).
#' @export
#' @examples
#' 1 %between% c(0,2)
#' 2 %between% c(0,2)
#' 3 %between% c(0,2)
"%between%" <- function(value,interval){
 (value >= interval[1]) & (value <= interval[2])
}

################################# ----
#' support_estimation XXXXXXXXX
################################# ----
#' @description Compute the support estimate.
#' @param beta_sample a matrix. Each row is a coefficient function computed from the
#'        posterior sample.
#' @return a list containing: 
#' \describe{
#'  \item{alpha}{a numerical vector. XXXXXXX}
#'  \item{estimate}{a numerical vector. XXXXXXX} 
#' }
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' res.support <- support_estimation(res_Bliss_mult$beta_sample[[1]])
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
support_estimation <- function(beta_sample){
 alpha <- apply(beta_sample,2, function(vec) sum(vec != 0)/length(vec))
 tmp   <- rep(0,ncol(beta_sample))
 tmp2  <- which(alpha >= 0.5)
 tmp[tmp2] <- 1
 estimate  <- interval_detection(tmp)
 estimate  <- estimate[ estimate[,3]==1 , -3]
 # Rajouter d'autres calculs
 return(list(alpha=alpha,estimate=estimate))
}
################################# ----
#' interval_detection XXXXXXXXX
################################# ----
#' @description Determine for which intervals a functions ???
#' @return a matrix with 3 columns : "begin", "end" and "value". The two first
#'           columns define the begin and the end of the pieces and the third
#'           gives the mean values of each piece.
#'           Or something else.
#' @param vec a numerical vector.
#' @param smooth a Boolean value, indicates if vec is a smooth function or a
#' stepwise function.
#' @param q a two-vector, if smooth is TRUE, it is used to detect a rupture.
#' The numerical values have to belong in [0,1].
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
#' change_grid XXXXXXXXXXXXx
################################# ----
#' @description Compute an approximation of a (discretized) function on a new 
#'              finer grid.
#' @return a numerical vector, the approximation of the function on the new grid.
#' @param fct a numerical vector, the function to evaluate on the new grid.
#' @param grid a numerical vector, the initial grid.
#' @param new_grid a numerical vector.
#' @details This is nasty code. XXXXXXXXXXXXx
#' @export
#' @examples
#' grid <- seq(0,1,l=1e1)
#' new_grid <- seq(0,1,l=1e2)
#' fct <- 3*grid^2 + sin(grid*2*pi)
#' plot(grid,fct,type="o")
#' lines(new_grid,finer_grid(fct,grid,new_grid),type="o",col=makeTransparent(2))
change_grid <- function(fct,grid,new_grid){
 res <- rep(0,length(new_grid))
 for(i in 1:(length(grid)-1)){
  index <- new_grid %between% grid[i:(i+1)]
  res[index] <- fct[i] + 0:(sum(index)-1) /
   (sum(index)-1) * (fct[i+1]-fct[i])
 }
 return(res)
}

################################# ----
#' prior_l XXXXXXXXXXXXXx
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
#' integrate_trapeze XXXXXXXXXXX
################################# ----


integrate_trapeze <- function(x,y){
 apply(as.matrix(y),2,function(vect) 
  sum(diff(x)*(vect[-1]+vect[-length(vect)]))/2)
} 


