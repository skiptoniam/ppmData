#' @name generate_background_points
#' @title create background points for study area or point data
#' @export
#' @param number_of_background_points number of back ground points to create - only applicable to quasirandom.
#' @param known_sites a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in sampling (of size Nxdimension).
#'   Note: SpatialPoints* will only be suitiable for \code{dimension=2}.
#' @param study_area a Raster* object giving the area over which to generate background points.
#'  If ignored, a rectangle defining the extent of \code{known_sites} will be used.
#' @param model_covariates an optional Raster* object containing covariates for
#'   modelling the point process (best use a Raster* stack or Raster* brick)
#' @param resolution resolution setup grid for integration points (default is 1 deg).
#' @param multispecies_points a community matrix of speciesXsites collection methods - still need to think about how to set this up.
#' probably need to setup an index to identify collection method.
#' eg. PO abs = 0, PO pres = 1, PA abs = 2, PA pres = 3, it would be easy to then assign any new collection methods with
#' PA-trawl-presence = 4, PA-trawl-abs = 5
#' @param method the type of method that should be used to generate background points.
#' 'grid' generates a regular grid at a set resolution.
#' 'quasirandom' generates quasirandom background points based on the coordinates.
#' 'quasirandom_covariates' generates quasirandom background points across the multiple dimensions of covariate hyperdimensions.
#' 'multispecies' more to follow soon.

generate_background_points <- function(number_of_background_points = 2000, # number of back ground points to create.
                                       known_sites = NULL,  #a set of coordinates,
                                       study_area = NULL,  #raster
                                       model_covariates = NULL, # a set of covariates.
                                       resolution = 1, #resolution to setup quadriature points in lat/lon
                                       multispecies_points = NULL, #a community matrix of speciesXsites - still need to thing about how to set this up.
                                       method = c('grid','quasirandom','quasirandom_covariates','multispecies')){

  #this function was built to match the dimensions of coordinates and other dimensions
  known_sites <- coords_match_dim(known_sites,dim(known_sites)[2])
  eps <- sqrt(.Machine$double.eps)

  if(is.null(known_sites))stop('come on mate, you need some occurrence points.')

  #multispecies points will attempt to use other information from other species to inform sampling bias.
  #this could be a set of multi-species observations or a layer of observation bias/points.
  #hopefully in the future I can setup a fitian style offset based on collection method. e.g PO, PA, Abundance, ect.
  if(!is.null(multispecies_points) | method == 'multispecies'){
    stop("multispecies point generation isn't workingn yet \n Using 'multispecies_points' or method='multispecies' will only make things worse\n hopefully I'll have it up and running soon.")
  }

  #if there is no study_area generate potential sites for halton function and a study area.
  if (is.null(study_area)) {
      message("no study area provided, a raster based on the extent of 'known_sites'\n with a default resolution of 1 deg will be returned.")
      study_area <- default_study_area(known_sites)
  }

  # create background points based on method.
  bkgrd_pts <- switch(method,
                grid = grid_method(resolution, study_area),
                quasirandom = quasirandom_method(number_of_background_points,study_area),
                quasirandom_covariates = quasirandom_covariates_method(number_of_background_points,
                                                                       model_covariates),
                multispecies = NULL)


  # could use this for bias/multispecies function.
  # region_size <- cellStats(area_study,sum, na.rm=TRUE)*1000

  #id have provided covaraites use this to extract out environmental data from rasterstack
   if (!is.null(model_covariates)){
     if(!inherits(study_area, c('RasterLayer','RasterStack','RasterBrick')))
       stop("model_covariates must be a raster, rasterstack or rasterbrick of model_covariates desired for modelling")

       # extract covariate data for background points
       bk_cvr <- extract(model_covariates,bkgrd_pts[,c("x","y")],method='simple',na.rm=TRUE)
       bk_pts <- cbind(bkgrd_pts,bk_cvr)

       # extract covariate data for presence points
       po_cvr <- extract(model_covariates,known_sites,method='simple',na.rm=TRUE)
       po_pts <- cbind(known_sites,po_cvr)

       #join the two dataset together
       covar_data <- rbind(po_pts,bk_pts[,-3])

       #create near zero weights for presence points
       po_wts <- rep(eps,nrow(known_sites))

       # create an entire dataset
       dat <- data.frame(presence=c(rep(1,nrow(known_sites)),rep(0,nrow(bkgrd_pts))),
                         covar_data,
                         weights=c(po_wts,bkgrd_pts$weights))
   } else {
     message('just creating background points with no associated covariate data')

     #create near zero weights for presence points
     po_wts <- rep(eps,nrow(known_sites))

     # create an entire dataset
     dat <- data.frame(presence=c(rep(1,nrow(known_sites)),rep(0,nrow(bkgrd_pts))),
                       x=c(known_sites$x,bkgrd_pts$x),y=c(known_sites$y,bkgrd_pts$y),
                       weights=c(po_wts,bkgrd_pts$weights))
   }

  dat <- rm_na_pts(dat)
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

    #create a dataframe of coordinates w/ area
    grid <- as.data.frame(rasterToPoints(dd)[,-3])
    grid$weights <- estimate_area(dd,grid)
    colnames(grid) <- c('x','y','weights')
  }

  return(grid)

}

quasirandom_method <- function(number_of_background_points, study_area){

  if(!inherits(study_area, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'quasirandom' method currently only works a raster input as a 'study_area'")

  #generate a set of potential sites for quasirandom generation
  potential_sites <- raster::rasterToPoints(study_area)[,-3]

  #if number_of_background_points>nrow(potential_sites) recursively create a finer grid.
  if(number_of_background_points>nrow(potential_sites)){
    fct <- 2
    while(number_of_background_points>nrow(potential_sites)){
    study_area <- disaggregate(study_area, fct, na.rm=TRUE)
    potential_sites <- raster::rasterToPoints(study_area)[,-3]
    }
  }


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
  samp <- as.data.frame(potential_sites[sampIDs, , drop = FALSE])
  colnames(samp) <- colnames(potential_sites)

  #new way to estimate weights
  samp$weights <- estimate_area(study_area,samp)

  return(as.data.frame(samp))
}

# still working on this method.
quasirandom_covariates_method <- function(number_of_background_points, covariates){

  if(!inherits(covariates, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'quasirandom_covariates' method currently only works a raster input as a 'covariates'")

  #generate a set of potential sites for quasirandom generation
  potential_sites <- raster::rasterToPoints(covariates)
  potential_sites <- na.omit(potential_sites)

  #if number_of_background_points>nrow(potential_sites) recursively create a finer grid.
  if(number_of_background_points>nrow(potential_sites)){
    fct <- 2
    while(number_of_background_points>nrow(potential_sites)){
      covariates <- disaggregate(covariates, fct, na.rm=TRUE)
      potential_sites <- raster::rasterToPoints(covariates)[,-3]
    }
  }

  covariates_ext <- extent(covariates)[1:4]

  #setup the dimensions need to halton random numbers - let's keep it at 2 for now, could expand to alternative dimension in the future
  dimension <- dim(potential_sites)[2]

  N <- nrow(potential_sites)

  # add includsion probs in the future, we can use this as a bias layer.

  #if (is.null(inclusion_probs)) {
   # message("No inclusion.probs supplied, assuming uniform")
    inclusion_probs <- rep(1/N, N)
  #
  inclusion_probs1 <- inclusion_probs/max(inclusion_probs)
  mult <- 10
  samp <- randtoolbox::halton(sample(1:10000, 1), dim = dimension +
                                1, init = TRUE)
  njump <- number_of_background_points * mult
  samp <- randtoolbox::halton(max(njump, 5000), dim = dimension +
                                1, init = FALSE)
  myRange <- apply(potential_sites, -1, range)
  for (ii in 1:dimension) samp[, ii] <- myRange[1, ii] + (myRange[2,
                                                                  ii] - myRange[1, ii]) * samp[, ii]
  if (dimension == 2) {
    tmp <- mgcv::in.out(study.area, samp[, 1:dimension])
    samp <- samp[tmp, ]
  }
  sampIDs <- rep(NA, nrow(samp))
  kount <- 0
  flag <- TRUE
  while (flag & (kount < nrow(samp))) {
    if (kount == 0)
      message("Number of samples considered (number of samples found): ",
              njump, "(0) ", sep = "")
    else message(kount + njump, "(", length(sampIDs.2),
                 ") ", sep = "")
    sampIDs[kount + 1:min(njump, nrow(samp) - kount)] <- class::knn1(potential_sites,
                                                                     samp[kount + 1:min(njump, nrow(samp) - kount), -(dimension +
                                                                                                                        1), drop = FALSE], 1:nrow(potential_sites))
    sampIDs.2 <- which(samp[1:(kount + min(njump, nrow(samp) -
                                             kount)), dimension + 1] < inclusion_probs1[sampIDs[1:(kount +
                                                                                                     min(njump, nrow(samp) - kount))]])
    if (length(sampIDs.2) >= number_of_background_points) {
      sampIDs <- sampIDs[sampIDs.2][1:number_of_background_points]
      flag <- FALSE
    }
    kount <- kount + njump
  }
  message("Finished\n")
  if (kount > nrow(samp))
    stop("Failed to find a design. It is likely that the inclusion probabilities are very low and uneven. Please try again OR make inclusion probabilities more even")
  samp <- as.data.frame(cbind(potential_sites[sampIDs, , drop = FALSE],
                              inclusion_probs[sampIDs], sampIDs))
  colnames(samp) <- c(colnames(potential_sites), "inclusion.probabilities",
                      "ID")

  #new way to estimate weights
  samp$weights <- estimate_area(covariates[[1]],samp[,c('x','y')])

  return(samp[,c('x','y','weights')])
}


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

default_study_area <- function (known_sites) {
  # get limits
  xlim <- range(known_sites$x)
  ylim <- range(known_sites$y)

  # add on 5%
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.1
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.1

  xlim[1] <- floor(xlim[1])
  xlim[2] <- ceiling(xlim[2])
  ylim[1] <- floor(ylim[1])
  ylim[2] <- ceiling(ylim[2])

  e <- extent(c(xlim,ylim))
  sa <- raster(e,res=1, crs="+proj=longlat +datum=WGS84")
  sa[]<- 1:ncell(sa)
  return (sa)
}

# remove NA values from a pts object
rm_na_pts <- function (pts) {

  # remove any points that are NA and issue a warning
  if (any(is.na(pts))) {

    # copy the original data
    pts_old <- pts
    # remove NAs
    pts <- na.omit(pts)

    # report which ones were removed
    rm <- as.vector(attributes(pts)$na.action)
    which_rm <- pts_old$points[rm]

    if (any(which_rm == 0)) {
      warning(sprintf('Removed %i integration points for which covariate values could not be assigned.
                      This may affect the results of the model, please check alignment between model_covariates and area.',
                      sum(which_rm == 0)))
    }

    if (any(which_rm == 1)) {
      warning(sprintf('Removed %i observed points for which covariate values could not be assigned.
                      This may affect the results of the model, please check alignment between model_covariates and coords.',
                      sum(which_rm == 0)))
    }

  }

  return (pts)

}

estimate_area <- function(study_area,samp){
  if(raster::isLonLat(study_area)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- area(study_area)*1000
    area_study <- mask(area_rast,study_area)
    wts <- extract(area_study,samp[,1:2],na.rm=TRUE)
    } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    cell_area <- res(study_area)[1]*res(study_area)[2]
    wts <- rep(cell_area,nrow(samp))
    }
   return(wts)
}
