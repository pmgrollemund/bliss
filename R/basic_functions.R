################################# ----
#' compute_beta_sample
################################# ----
#' @description Compute the coefficient function for each iteration of the
#'              Gibbs sample.
#' @return return a matrix. Each row is a function evaluated on the grid of
#'         time points.
#' @param posterior_sample a list provided by the function Bliss_Gibbs_Sampler.
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
#' @param Q numeric
#' @param progress a logical value. If TRUE, the algorithm progress is displayed.
#'         (optional)
#' @export
#' @examples
#' \donttest{
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- fit_Bliss(data=data1,param=param1)
#' beta_sample <- compute_beta_sample_mult(res_Bliss_mult$res_Gibbs_Sampler,param1)
#' indexes <- sample(nrow(beta_sample[[1]]),1e2,replace=FALSE)
#' cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
#' matplot(param1$grids[[1]],t(beta_sample[[1]][indexes,]),type="l",lty=1,col=cols)
#' }
compute_beta_sample <- function(posterior_sample,param,Q,progress=FALSE){
 if(progress) cat("Compute the coefficient function posterior sample. \n")
 # Initialize parameters
 Ks     <- param[["K"]]
 grids  <- param[["grids"]]
 ps     <- param[["p"]]
 basiss <- param[["basis"]]
 if(is.null(ps)) ps <- sapply(grids,length)
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
#'  \item{lims_kde}{a numerical vector (yl,yu) XXXXXXXXXXX (optional)}
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
#' res_Bliss_mult <- fit_Bliss(data=data1,param=param1)
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
#'                     lims_kde = param1$lims_kde[[q]],
#'                     h1       = param1$h1,
#'                     new_grid = param1[["new_grid"]],
#'                     xlim = range(param1$grids[[q]]) + c(-diff_grid,diff_grid),
#'                     progress = FALSE
#' )
#' density_estimate <- density_estimation(res_Bliss_mult$beta_sample[[1]],param_density)
#' image(density_estimate$res_kde2d,col=rev(heat.colors(100)))
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
 lims_kde <- param[["lims_kde"]]
 new_grid <- param[["new_grid"]]
 lims_estimate <- param[["lims_estimate"]]


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


 # Compute the coefficient function on the new grid (if required).
 if(!is.null(new_grid)){
  old_beta_sample <- beta_sample
  beta_sample <- matrix(0,nrow(beta_sample),p)
  if(progress)
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
  if(progress)
   cat("\t 'thin = ",thin,"' is too small. Now, thin = ",
      thin_min,".\n",sep="")
  thin <- thin_min
 }
 if(progress)
  cat("\t Thin the sample.\n")
 beta_sample <- beta_sample[seq(1+burnin,iter,by=thin),]

 # Perform the kde2d function
 if(progress)
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
#' @param gamma a numeric, default 0.5
#' @return a list containing:
#' \describe{
#'  \item{alpha}{a numerical vector. XXXXXXX}
#'  \item{estimate}{a numerical vector. XXXXXXX}
#' }
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- fit_Bliss(data=data1,param=param1)
#' res_support <- support_estimation(res_Bliss_mult$beta_sample[[1]])
#'  ### The estimate
#'  res_support$estimate
#'  ### Plot the result
#'  grid <- res_Bliss_mult$param$grids[[1]]
#'  plot(grid,res_support$alpha_t,ylim=c(0,1),type="l",xlab="",ylab="")
#'  for(k in 1:nrow(res_support$estimate)){
#'     segments(grid[res_support$estimate[k,1]],0.5,
#'            grid[res_support$estimate[k,2]],0.5,lwd=2,col=2)
#'    points(grid[res_support$estimate[k,1]],0.5,pch="|",lwd=2,col=2)
#'    points(grid[res_support$estimate[k,2]],0.5,pch="|",lwd=2,col=2)
#'  }
#'  abline(h=0.5,col=2,lty=2)
#' }
support_estimation <- function(beta_sample,gamma=0.5){
 alpha <- apply(beta_sample,2, function(vec) sum(vec != 0)/length(vec))
 tmp   <- rep(0,ncol(beta_sample))
 tmp2  <- which(alpha >= gamma)
 tmp[tmp2] <- 1
 estimate  <- interval_detection(tmp)
 estimate$value <- NULL
 # Rajouter d'autres calculs
 estimate2 <- NULL
 return(list(alpha=alpha,estimate=estimate,estimate2=estimate2))
}

################################# ----
#' interval_detection XXXXXXXXX
################################# ----
#' @description Determine for which intervals a functions ???
#' @return a matrix with 3 columns : "begin", "end" and "value". The two first
#'           columns define the begin and the end of the pieces and the third
#'           gives the mean values of each piece.
#'           Or something else.
#' @param beta_sample a numerical vector.
#' @importFrom stats qnorm sd
#' @export
#' @examples
#' \donttest{
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- fit_Bliss(data=data1,param=param1)
#' intervals <- interval_detection(res_Bliss_mult$Bliss_estimate[[1]])
#' par(mfrow=c(2,1))
#' plot(data1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="s")
#' plot(data1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="n")
#' for(k in 1:nrow(intervals)){
#'   segments(intervals[k,1],intervals[k,3],
#'           intervals[k,2],intervals[k,3],col=2,lwd=2)
#' }
#' par(mfrow=c(1,1))
#' intervals <- interval_detection(res_Bliss_mult$res_Simulated_Annealing[[1]]$posterior_expe,
#'                                smooth=TRUE)
#' plot(data1$grids[[1]],res_Bliss_mult$res_Simulated_Annealing[[1]]$posterior_expe,type="l")
#' for(k in 1:nrow(intervals)){
#'    segments(intervals[k,1],intervals[k,3],
#'             intervals[k,2],intervals[k,3],col=2,lwd=2)
#' }
#' }
interval_detection <- function(beta_sample){
 intervals <- data.frame()
 begin <- 1
 for (i in 2:length(beta_sample)){
  if (beta_sample[i] != beta_sample[i-1]) {
   end <- i - 1
   intervals <- rbind(intervals,
                      c(begin,
                        i-1,
                        beta_sample[i-1]))
   begin <- i
  }
 }
 if(begin == length(beta_sample))
  intervals <- rbind(intervals,
                     c(begin,
                       i,
                       beta_sample[i]))
 # temp IS 06/09/2018
 #names(intervals) <- c("begin", "end", "value")
 return(intervals)
}

#' compute_starting_point
#'
#' @param beta_expe numeric
#'
#' @return a value
#' @export
#' @examples
#' data(data1)
#' data(param1)
compute_starting_point <- function(beta_expe){
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
#' lines(new_grid,change_grid(fct,grid,new_grid),type="o",col="red")
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
#' dexp_grid
################################# ----
#' @description Compute the probability function of the Exponential prior on l. #XXXXXXXX
#' @return a numerical vector, which is the prability function on "l_values".
#' @param a a positive value, the mean of the Exponential prior.
#' @param l_values a numerical value, the discrete support of the parameter l.
#' @importFrom stats pgamma
#' @export
#' @examples
#' data(data1)
#' data(param1)
pexp_grid <- function(a,l_values){
 step <- diff(l_values)[1] / 2
 probs <- pexp(l_values + step ,a) -
  pexp(l_values - step ,a)

 return(probs)
}



################################# ----
#' integrate_trapeze XXXXXXXXXXX
################################# ----
#' integrate_trapeze
#'
#' @param x numeric
#' @param y numeric
#'
#' @return a graphic
#' @export
#' @examples
#' data(data1)
#' data(param1)
integrate_trapeze <- function(x,y){
 apply(as.matrix(y),2,function(vect)
  sum(diff(x)*(vect[-1]+vect[-length(vect)]))/2)
}


