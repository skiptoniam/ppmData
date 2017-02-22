#' @name qrbp2
#' @title create background points for study area or point data
#' @export
#' @param number_of_background_points number of back ground points to create - only applicable to quasirandom.
#' @param known_sites a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in sampling (of size Nxdimension).
#'   Note: SpatialPoints* will only be suitiable for \code{dimension=2}.
#' @param study_area a Raster* object giving the area over which to generate background points.
#'  If ignored, a rectangle defining the extent of \code{known_sites} will be used.
#' @param covariates an optional Raster* object containing covariates for
#'   modelling the point process (best use a Raster* stack or Raster* brick)
#' @param resolution resolution to convert polygon to grid (default is .1).
#' This call is needed if no Raster* is provided as the study.area - make sure that the resolution matches
#' the coordinate system. i.e if lat/lon = .1, while if equal area (meters), reso = 10000 (~.1).
#' Default is NULL, and will setup resolution on 100X100 grid.
#' @param multispecies_points a community matrix of speciesXsites collection methods - still need to think about how to set this up.
#' probably need to setup an index to identify collection method.
#' eg. PO abs = 0, PO pres = 1, PA abs = 2, PA pres = 3, it would be easy to then assign any new collection methods with
#' PA-trawl-presence = 4, PA-trawl-abs = 5
#'
qrbp2 <- function(number_of_background_points = 2000, # number of back ground points to create.
                  known_sites = NULL,  #a set of coordinates,
                  study_area = NULL,  #raster
                  model_covariates = NULL, # a set of covariates.
                  resolution = 1, #resolution to setup quadriature points in lat/lon
                  multispecies_points = NULL, #a community matrix of speciesXsites - still need to thing about how to set this up.
                  method = c('grid','quasirandom','multispecies')){


  eps <- sqrt(.Machine$double.eps)

  if(is.null(known_sites))stop('come on mate, you need some occurrence points.')

  #multispecies points will attempt to use other information from other species to inform sampling bias.
  #this could be a set of multi-species observations or a layer of observation bias/points.
  #hopefully in the future I can setup a fitian style offset based on collection method. e.g PO, PA, Abundance, ect.
  if(!is.null(multispecies_points) | method == 'multispecies'){
    stop("multispecies point generation isn't workingn yet \n
          Using 'multispecies_points' or method='multispecies' will only make things worse\n
          hopefully I'll have it up and running soon.")
  }

  #if there is no study_area generate potential sites for halton function and a study area.
  if (is.null(study_area)) {
      message("no study area provided, generating a regular grid based on extent of 'known_sites'.")
      N <- 100
      potential_sites <- as.matrix(expand.grid(as.data.frame(matrix(rep(1:N,
                                                                        times = 2), ncol = 2)))/N -
                                     1/(2 * N))
      study_area <- as.matrix(expand.grid(as.data.frame(matrix(c(rep(0,
                                                                     2), rep(1, 2)), nrow = 2, byrow = TRUE))))
      colnames(potential_sites) <- colnames(study_area) <- paste0("2",
                                                                  1:2)

      study_area <- study_area[c(1, 3, 4, 2), ]
      study_area_ext <- c(study_area[,1])
  }

  # create background points based on method.
  bkgrd_pts <- switch(method,
                grid = grid_method(resolution, study_area),
                quasirandom = quasirandom_method(number_of_background_points,study_area),
                multispecies = NULL)


  # could use this for bias/multispecies function.
  # region_size <- cellStats(area_study,sum, na.rm=TRUE)*1000

  #id have provided covaraites use this to extract out environmental data from rasterstack
   if (!is.null(covariates)){
     if(!inherits(study_area, c('RasterLayer','RasterStack','RasterBrick')))
       stop("covariates must be a raster, rasterstack or rasterbrick of covariates desired for modelling")
       bk_cvr <- extract(covariates,bkgrd_pts[,c("x","y")],method='simple',na.rm=TRUE)
       bk_pts <- cbind(bkgrd_pts,bk_cvr)

       po_cvr <- extract(covariates,known_sites,method='simple',na.rm=TRUE)
       po_pts <- cbind(known_sites,po_cvr)
       names(po_pts)[1:2]<-c('x','y')
       bkgrd_wts <- bkgrd_pts$weights
       covar_data <- rbind(po_pts,bk_pts[,-3:-4])
  }

  po_wts <- rep(eps,nrow(known_sites))
  bkgrd_wts <- bkgrd_pts$weights
  dat <- data.frame(presence=c(rep(1,nrow(known_sites)),rep(0,nrow(bkgrd_pts))),
                  covar_data,
                  weights=c(po_wts,bkgrd_wts))


#small function clean_pts (remove NA data)
cat(nrow(bk_pts)-nrow(bk_pts2),"background points had NA data,\n now only",
    nrow(bk_pts2),"background points from the original\n",
    nrow(bk_pts),"background points")

return(dat)
}

grid_method <- function(resolution=1,study_area){

  if(!inherits(study_area, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'grid' method currently only works a raster input as a 'study_area'")

  if(inherits(study_area, c('RasterLayer','RasterStack','RasterBrick'))){

    #set up the dissaggreation or factor
    fct <- (res(study_area)/resolution)[1]

    #if fct is >= 1 dissaggregate, else aggregate
    if(fct>=1) dd <- disaggregate(study_area, fct, na.rm=TRUE)
    else dd <- aggregate(study_area, 1/fct, na.rm=TRUE)

    #get area from new raster resolution - make it in meters^2
    areas <- area(dd)*1000

    #mask out NA data
    areas_w_data <- mask(areas,dd)

    #create a dataframe of coordinates w/ area
    grid <- rasterToPoints(areas_w_data)
  }

  return(grid)

}

quasirandom_method <- function(number_of_background_points, study_area){

  if(!inherits(study_area, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'quasirandom' method currently only works a raster input as a 'study_area'")

  #generate a set of potential sites for quasirandom generation
  potential_sites <- raster::rasterToPoints(study_area)[,-3]
  study_area_ext <- extent(study_area)[1:4]

  #setup the dimensions need to halton random numbers - let's keep it at 2 for now, could expand to alternative dimension in the future
  dimension <- dim(potential_sites)[2]

  #intialise random sequence of random numbers
  samp <- randtoolbox::halton(sample(1:10000, 1), dim = dimension +
                                1, init = TRUE)

  #generate a large number of random numbers - minimun is 10000.
  mult <- 20
  njump <- number_of_background_points * mult
  samp <- randtoolbox::halton(max(njump, 10000), dim = dimension +
                                1, init = FALSE)

  #now extract range from raster extent - could include multiple dimensions - eg - depth.
  myRange <- t(matrix(study_area_ext,2,2,byrow = TRUE))

  #expand the the random points to be on the same range as extent
  for (ii in 1:2) samp[, ii] <- myRange[1, ii] + (myRange[2,ii] - myRange[1, ii]) * samp[, ii]

  #generate sample ids
  sampIDs <- rep(NA, nrow(samp))
  sampIDs <- class::knn1(potential_sites,samp[, -(dimension + 1), drop = FALSE], 1:nrow(potential_sites))

  #remove duplicated sampIDs
  sampIDs <- sampIDs[!duplicated(sampIDs)]

  #select the number of desired background points
  if(length(sampIDs)<number_of_background_points)stop('try a few more background points.')
  sampIDs <- sampIDs[1:number_of_background_points]

  #get coordinates from potential sites
  samp <- as.data.frame(cbind(potential_sites[sampIDs, , drop = FALSE],
                              sampIDs))
  colnames(samp) <- c(colnames(potential_sites), "ID")

  #get area percell/point
  area_rast <- area(study_area)
  area_study <- mask(area_rast,study_area)

  # get the area of each pixel in meters^2 - will need to be careful if equal area map - it'll be in meters already.
  samp$weights <- extract(area_study,samp[,1:2],na.rm=TRUE)*1000

  return(samp)
}


