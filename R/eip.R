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
#' @param pfr (projection for raster) a equal area projection appropriate for the region. This will be
#' used when setting up the ppp density layer.
#' @param sigma A numeric value or a function to estimate sigma see \code{\link[spatstat]{density.ppp}}
#' for details.
#' @param plot.
#' @description Estimating the probabilty of presence from a series of spatial points.
#' The probability of *absence in an area of size A* according to the poisson distribution is
#'  \deqn{pr(y=0) = exp(-\lambda(u)*A)}
#'  The prob of *presence* is then
#'  \deqn{pr( y=1) = 1-pr(y=0)
#'                 = 1-exp(-\lambda(u)*A)}
#' \eqn{\lambda(u)} = the intensity value at point \eqn{u}. This is estimated using \code{\link[spatstat]{density}}

#' @examples
#' n <- 100
#' ks <- as.data.frame( cbind( x1=runif( n, min=-10, max=10), x2=runif( n, min=-10, max=10)))
#' sa <- raster(nrows=100, ncols=100, xmn=-10, xmx=10,ymn=-10,ymx=10)
#' sa[]<-rnorm(10000)
#' inclusion.probs <- eip(known.sites = ks, study.area = sa,sigma = 1)

eip <- function(known.sites,
                study.area = NULL,
                pfr = NULL,
                sigma = NULL,
                plot.prbs=TRUE){

  if (is.null(study.area)) {
    message("No study area supplied generating bounding box")
    known.sites <- coords_match_dim(known.sites,dimension)
    study.area <- default_study_area(known.sites)
    reso <- null_reso(study.area)
    bb <- bbox(study.area)
    bb <- reso*round(bb/reso)
    gt <- GridTopology(cellcentre.offset = bb[,1],
                       cellsize = c(reso, reso),
                       cells.dim = (c(diff(bb[1,]), diff(bb[2,]))/reso) + 1)


    study.area.grid <- SpatialGridDataFrame(gt,data.frame(id=1:(gt@cells.dim[1]*gt@cells.dim[2])),
                              proj4string = CRS(proj4string(study.area)))

    #let's create a gridded window.
    study.area.im <- as.im.SpatialGridDataFrame(study.area.grid)
    w <- as.owin(study.area.im)

    surveys<- coords_match_dim(known.sites,dimension)
    coordinates(surveys) <- ~x+y
    proj4string(surveys) <- pfr #need to change this to actual crs.

    # create a ppp object
    surveys.ppp <- as.ppp(coordinates(surveys),w)

    #create a density layer
    if(is.null(sigma))sigma <- reso
    dpp <- density.ppp(surveys.ppp,sigma=sigma,w=w)

    #work out the area per unit - this could be an issue if non-equal area. I'm sure there is away to work this out.
    unit_area <- dpp$xstep*dpp$ystep

    #this will produce the probabiltiy intensity for qausi-random sampling.
    prb_ppp <- 1-exp(-dpp*unit_area)
    if(plot.prbs==TRUE)plot(prb_ppp,main='probability of sampling intensity')
    # convert to a spatial grid data frame and export results.
    prb_sgdf <- as.SpatialGridDataFrame.im(prb_ppp)

    #currently uses a vector - some just extract the vector for now.
    inclus.probs <- prb_sgdf@data[!is.na(prb_sgdf@data)]

    }  else {

    #check that study area is a raster.
    if(!is(study.area, 'Raster')) stop("you have included a study area - but it's not a raster - start again.")
    if(!is(pfr, 'CRS')) pfr <- crs(study.area)
    if(is.na(proj4string(study.area))) stop(substitute(study.area), ' lacks a CRS.')

  # Then read in known sites
  surveys<- coords_match_dim(known.sites,dimension)
  coordinates(surveys) <- ~x+y
  proj4string(surveys) <- pfr #need to change this to actual crs.

  # now les't create a owin from a image
  study.area.img <- maptools::as.im.RasterLayer(study.area)
  w <- as.owin(study.area.img)

  # create a mask from ppp
  wm <- as.mask(w)

  # create a ppp object
  surveys.ppp <- as.ppp(coordinates(surveys),wm)

  #create a density layer
  dpp <- density.ppp(surveys.ppp,sigma=sigma,w=w)

  #work out the area per unit - this could be an issue if non-equal area. I'm sure there is away to work this out.
  unit_area <- dpp$xstep*dpp$ystep

  #this will produce the probabiltiy intensity for qausi-random sampling.
  prb_ppp <- 1-exp(-dpp*unit_area)
  if(plot.prbs==TRUE)plot(prb_ppp,main='probability of sampling intensity')
  # convert to a spatial grid data frame and export results.
  prb_sgdf <- as.SpatialGridDataFrame.im(prb_ppp)

  #currently uses a vector - some just extract the vector for now.
  inclus.probs <- prb_sgdf@data[!is.na(prb_sgdf@data)]
  }

  return(inclus.probs)

}


