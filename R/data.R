#' a list of data
#'
#' A data object for bliss model
#' @format a list of data
#' \describe{
#'   \item{Q}{the number of functional covariate}
#'   \item{y}{y coordinate}
#'   \item{x}{x coordinate}
#'   \item{betas}{}
#'   \item{grids}{The grid of the time observations}
#' }
"data1"

#'
#'
#' A list of param for bliss model
#' @format a list of param for bliss model
#' \describe{
#'   \item{Q}{the number of functional covariate}
#'   \item{n}{n is the sample size}
#'   \item{p}{p is the number of time observations of the curves}
#'   \item{iter}{The number of iteration of the main numerical algorithm of Bliss}
#'   \item{grids_lim}{The grid of the time observations}
#'   \item{K}{The number of intervals of the beta}
#' }
"param1"

#'
#'
#' A result of the BliSS method
#' @format a Bliss object (list)
#' \describe{
#'   \item{alpha}{a list}
#'   \item{beta_posterior_density}{a list containing grid_t, grid_beta_t, density and new_beta_sample}
#'   \item{beta_sample}{a list of beta sample}
#'   \item{Bliss_estimate}{a list of bliss estimate}
#'   \item{chains_info}{a list of information of computed chains}
#'   \item{param}{a list of the parameters used in the method}
#'   \item{data}{a list of the used data}
#'   \item{posterior_sample}{a list posterior_sample}
#'   \item{support_estimate}{a list support_estimate}
#'   \item{trace_sann}{a list of trace_sann if simulated annealing algo was used}
#' }
"res_bliss1"

