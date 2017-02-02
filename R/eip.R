#' @name eip
#' @title eip (estimate inclusion probabilties) creates a probability surface based on
#'  the location of known survey sites
#' @export
#' @param known.sites a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in sampling (of size Nxdimension).
#'   Note: SpatialPoints* will only be suitiable for \code{dimension=2}.
#' @param study.area an optional extent, SpatialPolygons* or Raster* object giving the
#'   area over which to generate background points. If ignored, a rectangle
#'   defining the extent of \code{known.sites} will be used instead.
#' @param params needs to be changed to actual parameters used to estimate the density layer.... more to follow.
#' @description Estimating the probabilty of presence from a series of spatial points.
#' The probability of *absence in an area of size A* according to the poisson distribution is
#'  \deqn{pr(y=0) = exp(-\lambda(u)*A)}
#'  The prob of *presence* is then
#'  \deqn{pr( y=1) = 1-pr(y=0)
#'                 = 1-exp(-\lambda(u)*A)}
#' \eqn{\lambda(u)} = the intensity value at point \eqn{u}. This is estimated using \code{\link[spatstat]{density}}

eip <- function(known.sites,study.area=NULL,params){

}

