#' a list of data
#'
#' A data object for bliss model
#' @format a list of data
#' \describe{
#'   \item{y}{y coordinate}
#'   \item{x}{x coordinate}
#'   \item{betas}{the coefficient function used to generate the data}
#'   \item{grids}{the grid of the observation times}
#' }
"data1"

#'
#'
#' A list of param for bliss model
#' @format a list of param for bliss model
#' \describe{
#'   \item{Q}{the number of functional covariates}
#'   \item{n}{the sample size}
#'   \item{p}{the number of observation times}
#'   \item{beta_shapes}{the shapes of the coefficient functions}
#'   \item{grids_lim}{the range of the observation times}
#'   \item{grids}{the grids of the observation times}
#'   \item{K}{the number of intervals for the coefficient function}
#' }
"param1"

#'
#'
#' A result of the BliSS method
#' @format a Bliss object (list)
#' \describe{
#'   \item{alpha}{a list of Q numerical vector. Each vector is the function
#'        alpha(t) associated to a functional covariate. For each t, alpha(t)
#'        is the posterior probabilities of the event "the support covers t".}
#'   \item{beta_posterior_density}{a list of Q items. Each item contains a list
#'        containing information to plot the posterior density of the
#'        coefficient function with the \code{image} function.
#'        \describe{
#'        \item{\code{grid_t}}{a numerical vector: the x-axis.}
#'        \item{\code{grid_beta_t}}{a numerical vector: the y-axis.}
#'        \item{\code{density}}{a matrix: the z values.}
#'        \item{\code{new_beta_sample}}{a matrix: beta sample used to compute
#'              the posterior densities.}
#'        }}
#'   \item{beta_sample}{a list of Q matrices. The qth matrix is a posterior
#'        sample of the qth functional covariates.}
#'   \item{Bliss_estimate}{a list of numerical vectors corresponding to the
#'        Bliss estimates of each functional covariates.}
#'   \item{data}{see the description of the object \code{data1}.}
#'   \item{posterior_sample}{a list containing (for each chain) the result of the
#'   \code{Bliss_Gibbs_Sampler} function.}
#'   \item{Smooth_estimate}{a list containing the Smooth estimates of the
#'   coefficient functions.}
#'   \item{support_estimate}{a list containing the estimations of the support.}
#'   \item{support_estimate_fct}{a list containing the estimation of the support.}
#'   \item{trace_sann}{a list containing (for each chain) the trace of the
#'   Simulated Annealing algorithm.}
#' }
"res_bliss1"

