#' @title qrbp package
#' @name qrbp-package
#' @description qrbp is a package to create quasi-random background points for use in
#'   Poisson point process modelling in R. The main function is \code{qrbp}, which takes some
#'   coordinates, an area of interest and covariates (and other parameters) and produces a \code{qrbp}
#'   dataset. It is essentially a wrapper around the existing \code{\link[MBHdesign]{quasiSamp}}
#'   function in the MBHdesign package which alreay implements quasi-random sampling within a
#'   basic domain.
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

     coords <- coords_match_dim(coords,dimension)

}

# convert coordinates (in whatever format they arrive in) into a dataframe that matches the dimensions
# used for sampling.
coords_match_dim <- function (coords,dimension){

  # check object classes
  expectClasses(coords,
                c('matrix',
                  'data.frame',
                  'SpatialPoints',
                  'SpatialPointsDataFrame'),
                name = 'coords')

  if (is.matrix(coords) | is.data.frame(coords)) {

    # for matrices/dataframes, make sure there are only two columns
    if (ncol(coords) != dimension ) {
      stop (sprintf('coords should match the number of dimensions used in quasi-random sampling,
                    giving the horizontal (x/longitude) then vertical (y/latitude) coordinates.
                    The object passed had %i columns, while the sampling dimensions are %i',NCOL(coords),dimension))
    }

    # otherwise, coerce into a data.frame and rename the columns
    df <- data.frame(coords)

  } else {
    # otherwise, for SpatialPoints* objects, just grab the coordinates
    df <- data.frame(coords@coords)
  }

  # set column names
  if (ncol(coords)>2) colnames(df) <- c("x","y",paste0("var",seq_len(ncol(coords)-2)))
  else colnames(df) <- c('x', 'y')
  return (df)
}

