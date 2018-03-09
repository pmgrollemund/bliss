################################# ----
#' image_Bliss
################################# ----
#' @description plot the representation of the estimated posterior density.
#' @param posterior_density_estimate a list. The result of the function
#'                 density_estimation.
#' @param param a list containing
#' \describe{
#' \item{cols}{a vector of colors for the function image (optional)}
#' \item{col_scale}{a character vector}
#' \item{ylim}{a numerical two-vector (optional)}
#' \item{main}{a character string.}
#' }
#' @importFrom stats quantile
#' @importFrom grDevices heat.colors
#' @export
#' @examples
#' library(RColorBrewer)
#' data(data1)
#' data(param1)
#' res_Bliss_mult <- Bliss_multiple(data=data1,param=param1,density=TRUE)
#' param1$cols <- colorRampPalette(brewer.pal(9,"Reds"))(1e2)
#' image_Bliss(res_Bliss_mult$posterior_density_estimate[[1]],param1)
#' lines(param1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(param1$grids[[1]],data1$beta_function[[1]],col=3,lwd=2,type="s")
#'
#' param1$cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
#' image_Bliss(res_Bliss_mult$posterior_density_estimate[[1]],param1)
#' lines(param1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(param1$grids[[1]],data1$beta_function[[1]],col=3,lwd=2,type="s")
#'
#' param1$cols <- rev(heat.colors(12))
#' param1$col_scale <- "quantile"
#' image_Bliss(res_Bliss_mult$posterior_density_estimate[[1]],param1)
#' lines(param1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(param1$grids[[1]],data1$beta_function[[1]],col=3,lwd=2,type="s")
#'
#' param1$cols <- rev(terrain.colors(12))
#' image_Bliss(res_Bliss_mult$posterior_density_estimate[[1]],param1)
#' lines(param1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(param1$grids[[1]],data1$beta_function[[1]],col=2,lwd=2,type="s")
#'
#' param1$cols <- rev(topo.colors(12))
#' image_Bliss(res_Bliss_mult$posterior_density_estimate[[1]],param1)
#' lines(param1$grids[[1]],res_Bliss_mult$Bliss_estimate[[1]],type="s",lwd=2)
#' lines(param1$grids[[1]],data1$beta_function[[1]],col=2,lwd=2,type="s")
image_Bliss <- function(posterior_density_estimate,param){
  cols      <- param$cols #Ceci n'est pas une modification pmg 08-03-18
  col_scale <- param$col_scale
  ylim      <- param$ylim
  xlim      <- param$xlim
  main      <- param$main
  if(is.null(cols)){
    #cols <- colorRampPalette(brewer.pal(9,"Spectral"))(1e2)
    cols <- rev(heat.colors(100))
  }
  nbre_cols <-length(cols)
  if(is.null(col_scale)) col_scale <- "regular"
  if(col_scale == "regular"){
    breaks <- seq(min(as.vector(posterior_density_estimate$res.kde2d$z)),
                  max(as.vector(posterior_density_estimate$res.kde2d$z)),
                  length=nbre_cols+1)
  }
  if(col_scale == "quantile"){
    breaks <- quantile(as.vector(posterior_density_estimate$res.kde2d$z),
                       seq(0,1,length=nbre_cols+1))
  }
  if(is.null(ylim)) ylim <- range(posterior_density_estimate$res.kde2d$y)
  if(is.null(xlim)) xlim <- range(posterior_density_estimate$res.kde2d$x)
  if(is.null(main)) main <- "Bliss estimates"
  image(posterior_density_estimate$res.kde2d,col=cols,breaks = breaks,
        main=main,ylim=ylim,xlim=xlim,useRaster = TRUE)
}

# image_Bliss2 <- function(posterior_density_estimate,param,...){
#   cols      <- param$cols
#   col_scale <- param$col_scale
#   ylim      <- param$ylim
#   xlim      <- param$xlim
#   main      <- param$main
#   if(is.null(cols)){
#     #cols <- colorRampPalette(brewer.pal(9,"Spectral"))(1e2)
#     cols <- rev(heat.colors(100))
#   }
#   nbre_cols <-length(cols)
#   if(is.null(col_scale)) col_scale <- "regular"
#   if(col_scale == "regular"){
#     breaks <- seq(min(as.vector(posterior_density_estimate$res.kde2d$z)),
#                   max(as.vector(posterior_density_estimate$res.kde2d$z)),
#                   length=nbre_cols+1)
#   }
#   if(col_scale == "quantile"){
#     breaks <- quantile(as.vector(posterior_density_estimate$res.kde2d$z),
#                        seq(0,1,length=nbre_cols+1))
#   }
#   if(is.null(ylim)) ylim <- range(posterior_density_estimate$res.kde2d$y)
#   if(is.null(xlim)) xlim <- range(posterior_density_estimate$res.kde2d$x)
#   if(is.null(main)) main <- "BLiSS estimates"
#   image(posterior_density_estimate$res.kde2d,col=cols,breaks = breaks,
#         main=main,ylim=ylim,xlim=xlim,useRaster = TRUE,...)
# }
