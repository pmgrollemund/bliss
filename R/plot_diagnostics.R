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
#' data(data1)
#' data(param1)
#' param1$n_chains <- 3
#' param1$iter <- 1e3
#' param1$burnin <- 1e2
#' param1$display <- FALSE
#' param1$compute_posterior <- FALSE
#' res_bliss_chains <- Bliss_multiple(data1,param1)
#' res_diagnostic <- diagnostics(res_bliss_chains$chains,param1)
#' plot_diagnostics(res_diagnostic,param1,which_plot="mu")
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
  cols_chains2 <- makeTransparent(cols_chains,alpha=1/(max(5,n_chains)))
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
        lines(hist_mu[[j]],border = 0,col = cols_chains2[j])
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
        lines(hist_sigma[[j]],border = 0,col = cols_chains2[j])
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
        lines(hist_beta[[ts_index]][[j]],border = 0,col = cols_chains2[j])
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
