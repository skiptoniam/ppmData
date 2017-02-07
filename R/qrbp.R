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
#' @param known.sites a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in sampling (of size Nxdimension).
#'   Note: SpatialPoints* will only be suitiable for \code{dimension=2}.
#' @param include.known.sites logical If set to TRUE (default = FALSE), then will attempt
#'  to use known.sites as legacy.sites. See \code{\link[MBHdesign]{alterInclProbs}} for details.
#' @param study.area an optional extent, SpatialPolygons* or Raster* object giving the
#'   area over which to generate background points. If ignored, a rectangle
#'   defining the extent of \code{known.sites} will be used instead.
#' @param inclusion.probs a vector specifying the inclusion probability for each of the
#' N potential sampling sites. This is the probability that each site will be included
#' in the final sample. Locations are ordered the same as the potential.sites argument.
#' If NULL (default) equal inclusion probabilities are specified.
#' @param covariates an optional Raster* object containing covariates for
#'   modelling the point process (best use a Raster* stack or Raster* brick)
#' @param reso resolution to convert polygon to grid (default is .1).
#' This call is needed if no Raster* is provided as the study.area - make sure that the resolution matches
#' the coordinate system. i.e if lat/lon = .1, while if equal area (meters), reso = 10000 (~.1).
#' Default is NULL, and will setup resolution on 100X100 grid.
#' @param sigma Shape parameter for gaussian kernel. See \code{\link[MBHdesign]{alterInclProbs}} for details.
#' @param plot.prbs logical if TRUE plots the underlying probability layer an background points on top.
#' @examples
#' # Generate some random points and a raster to represent study area.
#' N <- 100
#' ks <- as.data.frame(cbind(x1=runif(N, min=-10, max=10),x2=runif(N, min=-10, max=10)))
#' sa <- raster(nrows=100, ncols=100, xmn=-10, xmx=10,ymn=-10,ymx=10)
#' sa[]<-rnorm(10000)
#' projection(sa) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#'
#' # How many qausiRandom background points do we need?
#' set.seed(123)
#' n <- 200
#' bkpts <- qrbp(n,dimension = 2,known.sites=ks,include.known.sites=TRUE,
#'               study.area = sa,inclusion.probs = NULL,sigma=1,plot.prbs=TRUE)
#'

qrbp <- function(n,
                 dimension = 2,
                 known.sites,
                 include.known.sites=FALSE,
                 study.area = NULL,
                 inclusion.probs = NULL,
                 sigma=NULL,
                 covariates = NULL,
                 reso=NULL,
                 plot.prbs=TRUE){

     known.sites <- coords_match_dim(known.sites,dimension)

     #if no study area is provided create a polygon around coordinates.
     if(is.null(study.area)) study.area <- default_study_area(known.sites)

     # create a gridded area based on a raster or polygon - resolution for polygon must be defined.
     # while resolution for raster in based on crs projection.
     if(is.null(reso))reso <- null_reso(study.area)
     X <- studyarea_to_gridded_points(study.area, reso = reso)

     #extract potential sites from spatial data frame
     potential.sites <- X@coords

     #make a nameless matrix otherwise Scott's function flips out.
     rownames(potential.sites) <- colnames(potential.sites) <- NULL

     if(include.known.sites){
       sprintf('known.sites included in qrbp to alter inclusion.probs. Calling on
                eip to adjust sampling probabilities with legacy site information')
       if(is.null(sigma)) sprintf('you need a sigma - go on have a guess or estiamte it using a
                                  spatstat function - see density.ppp for details.')
       legacy.sites <- find_known_sites(study.area = study.area, known.sites = known.sites)
       inclusion.probs <- eip(known.sites = known.sites, study.area = study.area,
                              sigma = sigma,plot.prbs = plot.prbs)
       inclusion.probs <- n * inclusion.probs / sum(inclusion.probs)
       }
      bg_points <- MBHdesign::quasiSamp(n=n,dimension = dimension, study.area = NULL,
                            potential.sites = potential.sites, inclusion.probs = inclusion.probs)
      if(plot.prbs)points(bg_points[,1:2],col='springgreen',pch=16)
      return(bg_points)

}

# convert coordinates (in whatever format they arrive in) into a dataframe that matches the dimensions
# used for sampling.
coords_match_dim <- function (known.sites,dimension){

  # check object classes
  expectClasses(known.sites,
                c('matrix',
                  'data.frame',
                  'SpatialPoints',
                  'SpatialPointsDataFrame'),
                name = 'known.sites')

  if (is.matrix(known.sites) | is.data.frame(known.sites)) {

    # for matrices/dataframes, make sure there are only two columns
    if (ncol(known.sites) != dimension ) {
      stop (sprintf('known.sites should match the number of dimensions used in quasi-random sampling.
                    The object passed had %i columns, while the sampling dimensions are %i',NCOL(known.sites),dimension))
    }

    # otherwise, coerce into a data.frame and rename the columns
    df <- data.frame(known.sites)

  } else {
    # otherwise, for SpatialPoints* objects, just grab the coordinates
    df <- data.frame(known.sites@coords)
  }

  # set column names
  if (ncol(known.sites)>2) {
    if(!is.null(colnames(known.sites))) colnames(df) <- c("x","y",colnames(known.sites)[-1:-2])
    colnames(df) <- c("x","y",paste0("var",seq_len(ncol(known.sites)-2)))
  } else {
    colnames(df) <- c("x","y")
  }
  return (df)
}

#create a default study area based on known.sites - produces a simple box.
default_study_area <- function (known.sites) {
  # get limits
  xlim <- range(known.sites$x)
  ylim <- range(known.sites$y)

  # add on 5%
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.05
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.05

    # make a SpatialPolygons object
  p <- Polygon(cbind(x = xlim[c(1, 1, 2, 2)],
                     y = ylim[c(1, 2, 2, 1)]))
  ps <- Polygons(list(p), 1)
  sp <- SpatialPolygons(list(ps))
  return (sp)
}

#study.area needs to be a shapefile atm. I need to add an option for raster.
studyarea_to_gridded_points <- function(study.area,reso=reso){

  if(inherits(study.area, c('SpatialPolygons', 'SpatialPolygonsDataFrame'))){
    #create a bounding box around polygons/shapefile
    bb <- bbox(study.area)
    #clean up the bounding box
    bb <- reso*round(bb/reso)
    #create a grid to generate points
    gt <- GridTopology(cellcentre.offset = bb[,1],
                       cellsize = c(reso, reso),
                       cells.dim = (c(diff(bb[1,]), diff(bb[2,]))/reso) + 1)

    bg_pts <- SpatialPoints(gt,proj4string = CRS(proj4string(study.area)))
    bg_pts <- SpatialPointsDataFrame(gt,data.frame(id=1:nrow(bg_pts@coords)),
                                      proj4string = CRS(proj4string(study.area)))
    #
    #mask out points outside shapefile and keep an index of cell ids.
    vals <- over(bg_pts, study.area)
    bg_pts_vals <- cbind(coordinates(bg_pts), vals, id=1:nrow(bg_pts@coords))
    x <- bg_pts_vals[!is.na(bg_pts_vals[,3]),-3]
    # do I want a spatial points data frame? yes,
    x <- SpatialPointsDataFrame(x[,1:2],data.frame(id=x[,3]), proj4string = CRS(proj4string(study.area)))

  }
  if(inherits(study.area, c('RasterLayer','RasterStack','RasterBrick'))){
    r <- raster(study.area)
    r[] <- 1:ncell(r)
    r <- mask(r,study.area)
    x <- rasterToPoints(r,spatial = TRUE)#if you want a SpatialPointsDataFrame object spatial=TRUE
  }
  rownames(x)<-colnames(x)<-NULL
  return(x)
}

find_known_sites <- function(study.area,known.sites){
  if(inherits(study.area, c('SpatialPolygons', 'SpatialPolygonsDataFrame'))){
    bb <- bbox(study.area)
    #clean up the bounding box
    bb <- reso*round(bb/reso)
    #create a grid to generate points
    gt <- GridTopology(cellcentre.offset = bb[,1],
                       cellsize = c(reso, reso),
                       cells.dim = (c(diff(bb[1,]), diff(bb[2,]))/reso) + 1)

    # bg_pts <- SpatialPoints(gt,proj4string = CRS(proj4string(study.area)))
    bg_pts <- SpatialPointsDataFrame(gt,data.frame(id=1:nrow(bg_pts@coords)),
                                     proj4string = CRS(proj4string(study.area)))
    #
    #mask out points outside shapefile and keep an index of cell ids.
    vals <- over(bg_pts, study.area)
    bg_pts_vals <- cbind(coordinates(bg_pts), vals, id=1:nrow(bg_pts@coords))
    x <- bg_pts_vals[!is.na(bg_pts_vals[,3]),-3]
    r <- rasterFromXYZ(as.data.frame(x))

    #re-define index and get out cell numbers.
    r[]<-1:ncell(r)
    cn <- extract(r,known.sites)
    cn_clean <- cn[!is.na(cn)]
  }

  if(inherits(study.area, c('RasterLayer','RasterStack','RasterBrick'))){

  r <- study.area
  r[!is.na(r[])] <- 1:length(r[!is.na(r[])])
  cn <- cellFromXY(r,known.sites)
  na_sites <- extract(r,known.sites)
  cn_na <- cbind(cn,na_sites)
  cn_clean <- cn_na[!is.na(cn_na[,2]),2]
  }
  return(cn_clean)
}

null_reso <- function (study.area) {
  ext <- extent(study.area)
  height <- abs(diff(ext[1:2]))
  width <-  abs(diff(ext[3:4]))
  reso <- round(abs(seq(min(ext[1]),max(ext[2]),length.out = 100)[1]-seq(min(ext[1]),max(ext[2]),length.out = 100)[2]))
  return (reso/2)
}
