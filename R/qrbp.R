#' @title qrbp package
#' @name qrbp-package
#' @description qrbp is a package to create quasi-random background points for use in
#'   Poisson point process modelling in R. The main function is \code{qrbp}, which takes some
#'   coordinates, an area of interest and covariates and produces a \code{qrbp}
#'   dataset. It is essentially a wrapper around the existing \code{\link[MBHdesign]} package which alreay implements
#'   quasi-random sampling within a basic domain.
#' @docType package
#' @import raster
#' @import sp
#' @import MBHdesign

NULL

#' @name qrbp
#' @title create quasi-random background points from study area or point data
#' @export
#' @param n integer The number of background point to produce.
#' @param dimension The number of dimensions to sample the study extent.
#' \code{dimension=2} is the default and samples an area from an areal perspective.
#' @param study.area an optional extent, SpatialPolygons* or Raster* object giving the
#'   area over which to generate background points. If ignored, a rectangle
#'   defining the extent of \code{coords} will be used instead.
#' @param inclusion.probs a vector specifying the inclusion probability for each of the
#' N potential sampling sites. This is the probability that each site will be included
#' in the final sample. Locations are ordered the same as the potential.sites argument.
#' If NULL (default) equal inclusion probabilities are specified.
#'


qrbp <- function(n,
                 dimension = 2,
                 coords = NULL,
                 study.area = NULL,
                 inclusion.probs = NULL,
                 covariates = NULL){


}

