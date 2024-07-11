## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE,message=FALSE, warning=FALSE-----------------------------------
  library(bliss)

## ----eval=TRUE,include = TRUE-------------------------------------------------
  set.seed(1)
  param <- list(                        # define the "param" to simulate data
                Q=1,                    # the number of functional covariate
                n=50,                  # n is the sample size and p is the
                p=c(20),                # number of time observations of the curves
                beta_types=c("smooth"), # define the shape of the "true" coefficient function
                grids_lim=list(c(0,1))) # Give the beginning and the end of the observation's domain of the functions.

  data <- sim(param) # Simulate the data

## ----eval=TRUE, include = TRUE------------------------------------------------
  param <- list(            # define the required values of the Bliss method.
                iter=5e2,   # The number of iteration of the main numerical algorithm of Bliss.
                burnin=2e2, # The number of burnin iteration for the Gibbs Sampler
                K=c(3))     # The number of intervals of the beta

   
  res_bliss<-fit_Bliss(data=data,param=param)
  
  # Structure of a Bliss object
  # str(res_bliss)

## ----eval=TRUE, include = TRUE,fig.height=5,fig.width=7-----------------------
  library(ggplot2)
  image_Bliss(res_bliss$beta_posterior_density,param,q=1) 

## ----eval=TRUE, include = TRUE,fig.height=5,fig.width=7-----------------------
  image_Bliss(res_bliss$beta_posterior_density,param,q=1) + 
    lines_bliss(res_bliss$data$grids[[1]],res_bliss$Bliss_estimate[[1]]) + 
    lines_bliss(res_bliss$data$grids[[1]],res_bliss$smooth_estimate[[1]],lty = "dashed")+ 
    lines_bliss(res_bliss$data$grids[[1]],data$betas[[1]],col="purple")

## ----eval=TRUE, include = TRUE,fig.height=5,fig.width=7-----------------------
  plot(res_bliss$alpha[[1]],type="o",xlab="time",ylab="posterior probabilities")

## ----eval=TRUE, include = TRUE,fig.height=5,fig.width=7-----------------------
  plot(res_bliss$alpha[[1]],type="o",xlab="time",ylab="posterior probabilities")
  abline(h=0.5,col=2,lty=2)
  
  for(i in 1:nrow(res_bliss$support_estimate[[1]])){
  segments(res_bliss$support_estimate[[1]]$begin[i],0.05,
           res_bliss$support_estimate[[1]]$end[i],0.05,col="red"
           )
  points(res_bliss$support_estimate[[1]]$begin[i],0.05,col="red",pch="|",lwd=2)
  points(res_bliss$support_estimate[[1]]$end[i],0.05,col="red",pch="|",lwd=2)
  }

## ----eval=TRUE, include = TRUE,fig.height=5,fig.width=7-----------------------
res_bliss$support_estimate[[1]]

## ----eval=FALSE, include = TRUE-----------------------------------------------
#     param <- list(Q=2,
#                   n=50,
#                   p=c(40,10),
#                   beta_shapes=c("simple","smooth"),
#                   grids_lim=list(c(0,1),c(0,2)))
#  
#    data <- sim(param)

## ----eval=FALSE, include = TRUE-----------------------------------------------
#    param <- list(       # define the required values of the Bliss method.
#       iter=1e3,         # The number of iteration of the main numerical algorithm of Bliss.
#       burnin=2e2,       # The number of burnin iteration for the Gibbs Sampler
#       K=c(3,3))         # The number of intervals of the beta
#  
#    res_Bliss_mult <- fit_Bliss(data=data,param=param)

## ----eval=FALSE, include = TRUE,fig.height=5,fig.width=7----------------------
#    image_Bliss(res_Bliss_mult$beta_posterior_density,param,q=1) +
#      lines_bliss(res_Bliss_mult$data$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]]) +
#      lines_bliss(res_Bliss_mult$data$grids[[1]],res_Bliss_mult$smooth_estimate[[1]],lty = "dashed")+
#      lines_bliss(res_Bliss_mult$data$grids[[1]],data$betas[[1]],col="purple")
#  
#    image_Bliss(res_Bliss_mult$beta_posterior_density,param,q=2) +
#      lines_bliss(res_Bliss_mult$data$grids[[2]],res_Bliss_mult$Bliss_estimate[[2]]) +
#      lines_bliss(res_Bliss_mult$data$grids[[2]],res_Bliss_mult$smooth_estimate[[2]],lty = "dashed")+
#      lines_bliss(res_Bliss_mult$data$grids[[2]],data$betas[[2]],col="purple")

## ----session,echo=FALSE,message=FALSE, warning=FALSE--------------------------
  sessionInfo()

