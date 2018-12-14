## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE,message=FALSE, warning=FALSE------------------------------
  library(bliss)

## ----eval=TRUE,include = TRUE--------------------------------------------
  param <- list(                        # define the "param" to simulate data
                Q=1,                    # the number of functional covariate
                n=100,                  # n is the sample size and p is the
                p=c(50),                # number of time observations of the curves
                beta_types=c("smooth"), # define the shape of the "true" coefficient function
                grids_lim=list(c(0,1))) # Give the beginning and the end of the observation's domain of the functions.

  data <- sim(param) # Simulate the data

## ----eval=TRUE, include = TRUE-------------------------------------------
  param <- list(                        # define the required values of the Bliss method.
                iter=1e3,               # The number of iteration of the main numerical algorithm of Bliss.
                K=c(3))                 # The number of intervals of the beta

  res_bliss<-fit_Bliss(data=data,param=param,verbose=T)

  # Structure of a Bliss object
  str(res_bliss)

## ----eval=TRUE, include = TRUE,fig.height=5,fig.width=7------------------
  param$cols <- rev(heat.colors(100))
  image_Bliss(res_bliss$beta_posterior_density,param,q=1)
  lines(res_bliss$data$grids[[1]],res_bliss$Bliss_estimate[[1]],type="s",lwd=2)
  lines(res_bliss$data$grids[[1]],res_bliss$data$betas[[1]],col=2,lwd=2,type="s")

## ----eval=TRUE, include = TRUE,fig.height=5,fig.width=7------------------
  plot_bliss(res_bliss$data$grids[[1]],
             res_bliss$Bliss_estimate[[1]],lwd=2)
  lines_bliss(res_bliss$data$grids[[1]],
             res_bliss$Smooth_estimate[[1]],lty=2)

## ----eval=FALSE, include = TRUE------------------------------------------
#    param <- list(Q=2,n=50,p=c(15,12),beta_types=c("simple","smooth"),
#                 b_inf=c(0,0),b_sup=c(1,1),
#                 iter=1e3,burnin=2e2,K=c(3,3),grids=data$grids,
#                 prior_beta="Ridge_Zellner",phi_l=list("Gamma","Gamma"))
#    data <- sim(param)

## ----eval=FALSE, include = TRUE------------------------------------------
#    res_Bliss_mult <- fit_Bliss(data=data,param=param)

## ----eval=FALSE, include = TRUE,fig.height=5,fig.width=7-----------------
#    # for(q in 1:Q){
#    #   ylim <- range(c(res_Bliss_mult$posterior_density_estimate[[q]]$res.kde2d$y,
#    #                   data$beta[[q]]))
#    #   ylim <- ylim +  (ylim[2] - ylim[1])/20 * c(-1,1)
#    #   param$ylim <- ylim
#    #
#    #   image_Bliss(res_Bliss_mult$posterior_density_estimate[[q]],param)
#    #   # the Bliss estimate
#    #   lines_step_function(res_Bliss_mult$param$grids2[[q]],
#    #                       res_Bliss_mult$Bliss_estimate[[q]],lwd=2,bound=F)
#    #   # the posterior expection of beta(t)
#    #   lines_step_function(res_Bliss_mult$param$grids2[[q]],
#    #                       res_Bliss_mult$res.Simulated_Annealing[[q]]$posterior_expe,
#    #                       lty=2,bound=T)
#    #   # plot the true coefficient function
#    #   lines(data$grids[[q]],data$beta[[q]],col=3)
#    # }

## ----session,echo=FALSE,message=FALSE, warning=FALSE---------------------
  sessionInfo()

