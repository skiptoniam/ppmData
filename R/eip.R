#' @name eip
#' @title eip (estimate inclusion probabilties) creates a probability surface based on
#'  the location of known survey sites
#' @export
#' @param known_sites a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in sampling (of size Nxdimension).
#'   Note: SpatialPoints* will only be suitiable for \code{dimension=2}.
#' @param study_area an optional extent Raster* object giving the
#'   area over which to generate background points. If ignored, a rectangle
#'   defining the extent of \code{known_sites} will be used instead.
#' @param pfr (projection for raster) a equal area projection appropriate for the region. This will be
#' used when setting up the ppp density layer.
#' @param sigma A numeric value or a function to estimate sigma see \code{\link[spatstat]{density.ppp}}
#' for details.
#' @param plot_prb plot the return probabiltiy surface.
#' @description Estimating the probabilty of presence from a series of spatial points.
#' The probability of *absence in an area of size A* according to the poisson distribution is
#'  \deqn{pr(y=0) = exp(-\lambda(u)*A)}
#'  The prob of *presence* is then
#'  \deqn{pr( y=1) = 1-pr(y=0)
#'                 = 1-exp(-\lambda(u)*A)}
#' where \eqn{\lambda(u)} = the intensity value at point \eqn{u} and A is the area of the sampling unit (cell size).
#' This is estimated using \code{\link[spatstat]{density}}
#' @importFrom maptools as.im.SpatialGridDataFrame
#' @importFrom maptools as.SpatialGridDataFrame.im

#' @examples
#' #generate some random points and a raster to represent study area.
#' n <- 100
#' ks <- as.data.frame( cbind( x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10)))
#' sa <- raster(nrows=100, ncols=100, xmn=-10, xmx=10,ymn=-10,ymx=10)
#' sa[]<-rnorm(10000)
#' projection(sa) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#' # run eip sigma in this case is chosen as ~5 times cell size.
#' inclusion.probs <- eip(known_sites = ks, study_area = sa,sigma = 1)
#'

eip <- function(known_sites,
                study_area = NULL,
                pfr = NULL,
                sigma = NULL,
                plot_prbs=TRUE){

  if (is.null(study_area)) {
    message("No study area supplied generating bounding box")
    dimension <- dim(known_sites)[2]
    known_sites <- coords_match_dim(known_sites,dimension)
    study_area <- default_study_area(known_sites)
    reso <- null_reso(study_area)
    bb <- bbox(study_area)
    bb <- reso*round(bb/reso)
    gt <- GridTopology(cellcentre.offset = bb[,1],
                       cellsize = c(reso, reso),
                       cells.dim = (c(diff(bb[1,]), diff(bb[2,]))/reso) + 1)


    study_area.grid <- sp::SpatialGridDataFrame(gt,data.frame(id=1:(gt@cells.dim[1]*gt@cells.dim[2])),
                              proj4string = CRS(proj4string(study_area)))

    #let's create a gridded window.
    study_area.im <- maptools::as.im.SpatialGridDataFrame(study_area.grid)
    w <- spatstat::as.owin(study_area.im)

    surveys<- coords_match_dim(known_sites,dimension)
    coordinates(surveys) <- ~x+y
    proj4string(surveys) <- pfr #need to change this to actual crs.

    # create a ppp object
    surveys.ppp <- spatstat::as.ppp(coordinates(surveys),w)

    #create a density layer
    if(is.null(sigma))sigma <- reso
    dpp <- spatstat::density.ppp(surveys.ppp,sigma=sigma,w=w)

    #work out the area per unit - this could be an issue if non-equal area. I'm sure there is away to work this out.
    unit_area <- dpp$xstep*dpp$ystep

    #this will produce the probabiltiy intensity for qausi-random sampling.
    # prb_ppp <- 1-exp(-dpp*unit_area)
    prb_ppp <- exp(-dpp*unit_area)
    if(plot_prbs==TRUE)plot(prb_ppp,main='probability of sampling intensity')
    # convert to a spatial grid data frame and export results.
    prb_sgdf <- maptools::as.SpatialGridDataFrame.im(prb_ppp)

    #currently uses a vector - some just extract the vector for now.
    inclus.probs <- prb_sgdf@data[!is.na(prb_sgdf@data)]

    }  else {

    #check that study area is a raster.
    if(!is(study_area, 'Raster')) stop("you have included a study area - but it's not a raster - start again.")
    if(!is(pfr, 'CRS')) pfr <- crs(study_area)
    if(is.na(proj4string(study_area))) stop(substitute(study_area), ' lacks a CRS.')

  # Then read in known sites
  dimension <- dim(known_sites)[2]
  surveys<- coords_match_dim(known_sites,dimension)
  coordinates(surveys) <- ~x+y
  proj4string(surveys) <- pfr #need to change this to actual crs.

  # now les't create a owin from a image
  study_area.img <- maptools::as.im.RasterLayer(study_area)
  w <- spatstat::as.owin(study_area.img)

  # create a mask from ppp
  wm <- spatstat::as.mask(w)

  # create a ppp object
  surveys.ppp <- spatstat::as.ppp(coordinates(surveys),wm)

  #create a density layer
  dpp <- spatstat::density.ppp(surveys.ppp,sigma=sigma,w=w)

  #work out the area per unit - this could be an issue if non-equal area. I'm sure there is away to work this out.
  unit_area <- dpp$xstep*dpp$ystep

  #this will produce the probabiltiy intensity for qausi-random sampling.
  prb_ppp <- 1-exp(-dpp*unit_area)
  # prb_ppp <- exp(-dpp*unit_area)
  if(plot_prbs==TRUE)plot(prb_ppp,main='probability of sampling intensity')
  # convert to a spatial grid data frame and export results.
  prb_sgdf <- maptools::as.SpatialGridDataFrame.im(prb_ppp)

  #currently uses a vector - some just extract the vector for now.
  inclus.probs <- prb_sgdf@data[!is.na(prb_sgdf@data)]
  r <- study_area
  r[!is.na(r[])]<-inclus.probs
  }

  return(r)

}


