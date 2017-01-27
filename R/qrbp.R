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
#' @param coords a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in sampling (of size Nxdimension).
#'   Note: SpatialPoints* will only be suitiable for \code{dimension=2}.
#'   If NULL (default) N=10000 samples are placed on a regular grid.
#' @param study.area an optional extent, SpatialPolygons* or Raster* object giving the
#'   area over which to generate background points. If ignored, a rectangle
#'   defining the extent of \code{coords} will be used instead.
#' @param inclusion.probs a vector specifying the inclusion probability for each of the
#' N potential sampling sites. This is the probability that each site will be included
#' in the final sample. Locations are ordered the same as the potential.sites argument.
#' If NULL (default) equal inclusion probabilities are specified.
#' @param covariates an optional Raster* object containing covariates for
#'   modelling the point process (best use a Raster* stack or Raster* brick)


qrbp <- function(n,
                 dimension = 2,
                 coords = NULL,
                 study.area = NULL,
                 inclusion.probs = NULL,
                 covariates = NULL,
                 res=.1){

     coords <- coords_match_dim(coords,dimension)

     #if no study area is provided create a polygon around coordinates.
     if(is.null(study.area)) study.area <- default_study_area(coords)

     X <- studyarea_to_gridded_points(study.area, res = res)

     if(!is.null(inclusion.probs) & !is.null(coords)){
       sprintf('coords included in qrbp and so are inclusion.probs calling on
alterInclProbs to adjust sampling probabilities with legacy site information')
       # p2 <- alterInclProbs( X[legacySites,], X, p)
       }

     bg_points <- quasiSamp(n=n,dimension = dimension, study.area = NULL,
                            potential.sites = X, inclusion.probs = inclusion.probs)

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
      stop (sprintf('coords should match the number of dimensions used in quasi-random sampling.
                    The object passed had %i columns, while the sampling dimensions are %i',NCOL(coords),dimension))
    }

    # otherwise, coerce into a data.frame and rename the columns
    df <- data.frame(coords)

  } else {
    # otherwise, for SpatialPoints* objects, just grab the coordinates
    df <- data.frame(coords@coords)
  }

  # set column names
  if (ncol(coords)>2) {
    if(!is.null(colnames(coords))) colnames(df) <- c("x","y",colnames(coords)[-1:-2])
    colnames(df) <- c("x","y",paste0("var",seq_len(ncol(coords)-2)))
  } else {
    colnames(df) <- c("x","y")
  }
  return (df)
}

## new function: create_grid
## This function will create a grid from scratch, polygon or raster.
default_study_area <- function (coords) {
  # get limits
  xlim <- range(coords$x)
  ylim <- range(coords$y)

    # make a SpatialPolygons object
  p <- Polygon(cbind(x = xlim[c(1, 1, 2, 2)],
                     y = ylim[c(1, 2, 2, 1)]))
  ps <- Polygons(list(p), 1)
  sp <- SpatialPolygons(list(ps))
  return (sp)
}

#study.area needs to be a shapefile atm. I need to add an option for raster.
studyarea_to_gridded_points <- function(study.area,res=1,dimension=2){

  if(inherits(study.area, c('SpatialPolygons', 'SpatialPolygonsDataFrame'))){
    #create a bounding box around polygons/shapefile
    bb <- bbox(study.area)
    #clean up the bounding box
    bb <- res*round(bb/res)
    #create a grid to generate points
    gt <- GridTopology(cellcentre.offset = bb[,1],
                       cellsize = c(res, res),
                       cells.dim = (c(diff(bb[1,]), diff(bb[2,]))/res) + 1)

    bg_pts <- SpatialPoints(gt, proj4string = CRS(proj4string(study.area)))

    #mask out points outside shapefile
    vals <- over(bg_pts, study.area)
    bg_pts_vals <- cbind(coordinates(bg_pts), vals)
    x <- bg_pts_vals[!is.na(bg_pts_vals[,3]),]
    # do I want a spatial points data frame?
    # x2 <- SpatialPoints(x[,1:2], proj4string = CRS(proj4string(study.area)))
    # x2
  }
  if(inherits(study.area, c('RasterLayer','RasterStack','RasterBrick'))){
    x <- rasterToPoints(study.area)#if you want a SpatialPointsDataFrame object spatial=TRUE
  }
  x<-x[,1:dimension]
  rownames(x)<-colnames(x)<-NULL
  return(x)
}

