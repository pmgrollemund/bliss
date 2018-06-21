################################# ----
#' support_estimation
################################# ----
#' @description Compute the support estimate.
#' @param beta a matrix. Each row is a coefficient function computed from the
#' posterior sample.
#' @return a list containing
#' \describe{
#' \item{alpha}
#' \item{estimate}
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
support_estimation <- function(beta){
 alpha <- apply(beta,2, function(vec) sum(vec != 0)/length(vec))
 tmp  <- rep(0,ncol(beta))
 tmp2 <- which(alpha >= 0.5)
 tmp[tmp2] <- 1
 estimate <- interval_detection(tmp)
 estimate <- estimate[ estimate[,3]==1 , -3]
 # Rajouter d'autres calculs
 return(list(alpha=alpha,estimate=estimate))
}
################################# ----
#' interval_detection
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
#' finer_grid
################################# ----
#' @description Compute an approximaton of a discrete function on a new finer grid.
#' @return a numerical vector, the approximation of the function on the new grid.
#' @param fct a numerical vector, the function to evaluate on the new grid.
#' @param grid a numerical vector, the initial grid.
#' @param new_grid a numerical vector.
#' @details This is nasty code.
#' @export
#' @examples
#' grid <- seq(0,1,l=1e1)
#' new_grid <- seq(0,1,l=1e2)
#' fct <- 3*grid^2 + sin(grid*2*pi)
#' plot(grid,fct,type="o")
#' lines(new_grid,finer_grid(fct,grid,new_grid),type="o",col=makeTransparent(2))
finer_grid <- function(fct,grid,new_grid){
 res <- rep(0,length(new_grid))
 for(i in 1:(length(grid)-1)){
  index <- new_grid %between% grid[i:(i+1)]
  res[index] <- fct[i] + 0:(sum(index)-1) /
   (sum(index)-1) * (fct[i+1]-fct[i])
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
#' between
################################# ----
#' @description Check if a number belong to an interval.
#' @return a Boolan value.
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
#' compute_beta
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
compute_beta <- function(res.Gibbs_Sampler,param){
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
#' compute_beta_posterior_density
################################# ----
#' @description Compute a graphical representation of the marginal posterior
#' distributions of beta(t) for each t.
#' @details The sample is thinned in order to reduce the correlation and
#'           so the time of the computation of the function \code{\link[=kde2d]{kde2d}}.
#' @return an approximation of the posterior density on a two-dimensional grid.
#' (corresponds to the result of the \code{\link[=kde2d]{kde2d}} function)
#' @param beta a list (given by the compute_beta function).
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
compute_beta_posterior_density <- function(beta,param){
 cat("Compute an approximation of the posterior density of the coefficient function.\n")
 # load optional objects
 grid <- param[["grid"]]
 iter <- param[["iter"]]
 p    <- param[["p"]]
 n        <- param[["n"]]
 thin     <- param[["thin"]]
 burnin   <- param[["burnin"]]
 lims.kde <- param[["lims.kde"]]
 h1       <- param[["h1"]]
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
#' plot.bliss
################################# ----
#' @description A suitable representation of the Bliss estimate.
#' @param extended_grid a numerical vector, a new grid to provide a suitable representation.
#' @param fct a numerical vector, the function to plot.
#' @param bound a boolean value. If bound is TRUE, the plot of fct is one line.
#' Otherwise, several lines are used to plot fct.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @importFrom grDevices gray.colors
#' @importFrom graphics abline segments
#' @export
#' @examples
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1)
#' ### Plot the BLiss estimate on a suitable grid
#' plot.bliss(res_Bliss_mult$param$grids[[1]],
#'                    res_Bliss_mult$Bliss_estimate[[1]],lwd=2,bound=FALSE)
plot.bliss <- function(extended_grid,fct,bound=FALSE,...){
 ylim <- range(fct)
 plot(extended_grid,extended_grid,type="n",ylim=ylim,...)
 lines.bliss(extended_grid,fct,bound=bound,...)
}


################################# ----
#' lines.bliss
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
lines.bliss <- function(extended_grid,fct,bound=FALSE,...){
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
