#' @name ppmData
#' @title Create a Point Process dataset for spatial presence only modelling.
#' @description Creates a point process data frame for modelling single species or multiple species (marked) presences. Generates a quadrature scheme based on Berman & Turner 1992; Warton & Shepard 2010. The function can generate a quadrature scheme for a regular grid, quasi-random or random points.
#' @export
#' @param npoints number of background points to create.
#' @param presences a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of each species' presence in (should be a matrix of nsites * 3) with the three columns being c("X","Y","SpeciesID"), where X is longitude, Y is latitude and SpeciesID is a integer, character or factor which assoicates each point to a species. If presences are NULL then ppmdat will return the quadrature (background) points.
#' @param window a Raster* object giving the area over which to generate background points. NA cells are ignored and masked out of returned data
#' . If ignored, a rectangle defining the extent of \code{presences} will be used.
#' @param covariates an optional Raster* object containing covariates for
#'   modelling the point process (best use a Raster* stack or Raster* brick)
#' @param resolution resolution setup grid for integration points (default is 1 deg) - but this need to be setup  with reference to original raster resolution.
#' @param method the type of method that should be used to generate background points. The options are:
#' 'grid' generates a regular grid of background points. See Berman & Turner 1992 or Warton & Shepard 2010 for details.
#' 'quasirandom' generates quasirandom background points. See Bratley & Fox 1998 or Foster etal 2015 for details.
#' 'random' generates a random set of background points - See Philips 2006 (ala MaxEnt) for details.
#' @param interpolation either 'simple' or 'bilinear' and this determines the interpolation method for interpolating data across different cell resolutions. 'simple' is nearest neighbour, 'bilinear' is bilinear interpolation.
#' @param coords is the name of site coordinates. The default is c('X','Y').
#' @importFrom mgcv in.out
#' @importFrom raster extract

ppmData <- function(npoints = 10000,
                    presences = NULL,  #a set of coordinates,
                    window = NULL,  #raster
                    covariates = NULL, # a set of covariates.
                    resolution = NULL, #resolution to setup quadriature points in lat/lon
                    method = c('grid','quasirandom','random'),
                    interpolation='bilinear',
                    coords = c('X','Y')){

  # check if there are duplicates in the presence data.
  presences <- checkDuplicates(presences)

  ## if no resolution is provided guess the nearest resolution to return npoints for grid method.
  if(is.null(resolution)) resolution <- guessResolution(npoints,window)

  if(is.null(presences)){
   message('Generating background points in the absence of species presences')
   background_sites <- switch(method,
                               grid = gridMethod(resolution, window, covariates),
                               quasirandom = quasirandomMethod(npoints,  window, covariates),
                               random = randomMethod(npoints, window, covariates))

   # dat <- data.frame(background_sites)#replace this with a function 'get_weights'
   dat <- getWeights()

  } else {

  # This function was built to match the dimensions of coordinates and other dimensions
  # Should be three check away...
  eps <- sqrt(.Machine$double.eps)
  presences <- qrbp:::coords_match_dim(presences,3)

  #multispecies points will attempt to use other information from other species to inform sampling bias.
  #this could be a set of multi-species observations or a layer of observation bias/points.
  #hopefully in the future I can setup a fitian style offset based on collection method. e.g PO, PA, Abundance, ect.
  if(method == 'multispecies'){
    if(is.null(multispecies_presences))
    stop("Using 'multispecies_presences' or method='multispecies' requires a matrix of sites x presences. \n
         Note NA in this matrix indicate a non-presence record for that species.")
    if(!is.null(multispecies_presences))
      message("'presences' must be all the unique coordinates for known species' presences")
    if(nrow(presences)!=nrow(multispecies_presences))
      stop("nrows of 'presences' must equal nrows of 'multispecies_presences'")
  }

  #if there is no window generate potential sites for halton function and a study area.
  if (is.null(window)) {
      message("no study area provided, a raster based on the extent of 'presences'\n with a default resolution of 1 deg will be returned.")
      window <- default_window(presences)
  }

  # create background points based on method.
  background_sites <- switch(method,
                grid = gridMethod(resolution, window, covariates),
                quasirandom = quasirandomMethod(npoints,  window, covariates),
                random = randomMethod(npoints,  window, covariates))

  nspp <- length(unique(presences[,"SpeciesID"]))
  if(nspp>1){
    sppweights <- lapply(1:nspp,function(x)get_weights(presences[presences$SppID==x],
                                                       background_sites[,1:2],window,coords))
  } else {
    sppweights <- get_weights(presences,background_sites[,1:2],window,coords)
  }

  #id have provided covaraites use this to extract out environmental data from rasterstack
   if (!is.null(covariates)){
     if(!inherits(window, c('RasterLayer','RasterStack','RasterBrick')))
       stop("covariates must be a raster, rasterstack or rasterbrick of covariates desired for modelling")

       # extract covariate data for background points
       if(method=='multispecies_grid'|method=='multispecies_quasi') {
         covars <- extract(covariates,site_weights$multispecies_presence[,coords],method=interpolation,na.rm=TRUE)
         NAsites <- which(!complete.cases(covars))
         if(length(NAsites)>0){
           print(paste0('A total of ',length(NAsites),' sites where removed from the background data, because they contained NAs, check environmental data and species sites data overlap.'))
            covars <- covars[-NAsites,,drop=FALSE]
            dat <- list()
            dat$model_matrix <- data.frame(site_weights$multispecies_presence[-NAsites,-1:-2],const=1,covars)
            dat$species_weights <- site_weights$multispecies_weights[-NAsites,-1]
         } else {
           dat <- list()
           dat$model_matrix <- data.frame(site_weights$multispecies_presence[,-1:-2],const=1,covars)
           dat$species_weights <- site_weights$multispecies_weights[,-1]
         }
       } else {
         # print(head(site_weights[,coords]))
         covars <- extract(covariates,site_weights[,coords],method=interpolation,na.rm=TRUE)
         NAsites <- which(!complete.cases(covars))
         if(length(NAsites)>0){
           print(paste0('A total of ',length(NAsites),' sites where removed from the background data, because they contained NAs, check environmental data and species sites data overlap.'))
           covars <- covars[-NAsites,,drop=FALSE]
           dat <- data.frame(presence=c(rep(1,nrow(presences)),rep(0,nrow(background_sites)))[-NAsites],
                             # site_weights[-NAsites,coords],
                             covars,
                             weights=site_weights$weights[-NAsites])#replace this with a function 'get_weights'
         } else {
           dat <- data.frame(presence=c(rep(1,nrow(presences)),rep(0,nrow(background_sites))),
                             # site_weights,
                             covars,
                             weights=site_weights$weights)
        }
      }
    } else {
    if(method=='multispecies_grid'|method=='multispecies_quasi'){
       dat <- list()
       dat$model_matrix <- data.frame(site_weights$multispecies_presence[,-1:-2],
                                       x=c(presences$x,background_sites$x),
                                       y=c(presences$y,background_sites$y))
       dat$species_weights <- site_weights$multispecies_weights[,-1]
     } else {
     # create an entire dataset
     dat <- data.frame(presence=c(rep(1,nrow(presences)),rep(0,nrow(background_sites))),
                       x=site_weights$x,y=site_weights$y,
                       weights=site_weights$weights)#replace this with a function 'get_weights'
   }
  }

  if(method!='multispecies_grid'|method!='multispecies_quasi') dat <- rm_na_pts(dat)
  }
  return(dat)
}

gridMethod <- function(resolution=1, window, covariates){

  if(!inherits(window, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'grid' method currently only works a raster input as a 'window'")

  if(inherits(window, c('RasterLayer','RasterStack','RasterBrick'))){

    #set up the dissaggreation or aggregate
    fct <- (res(window)/resolution)[1]

    #if fct is >= 1 dissaggregate, else aggregate
    if(fct>=1) dd <- disaggregate(window, fct, na.rm=TRUE)
    else dd <- aggregate(window, 1/fct, na.rm=TRUE)

    #create a dataframe of coordinates w/ area
    grid <- as.data.frame(rasterToPoints(dd)[,-3])
    colnames(grid) <- c('X','Y')
  }

  if(!is.null(covariates)){
    covars <- extract(covariates,grid)
    grid <- cbind(grid,covars)
  }

  return(grid)
}


# still working on this method.
quasirandomMethod <- function(npoints, window, covariates=NULL){


  if(!is.null(covariates)){
    if(!inherits(covariates, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'quasirandom_covariates' method currently only works a raster input as a 'covariates'")
  }
  #generate a set of potential sites for quasirandom generation
  if(!is.null(covariates)){
    rast_coords <- raster::xyFromCell(covariates,1:raster::ncell(covariates))
    rast_data <- raster::values(covariates)
    na_sites <- which(is.na(rast_data[,1]))
    covariates_ext <- raster::extent(covariates)[1:4]
  } else {
    rast_coords <- raster::xyFromCell(window,1:raster::ncell(window))
    rast_data <- raster::values(window)
    na_sites <- which(is.na(rast_data))
    covariates_ext <- raster::extent(window)[1:4]
    }

  potential_sites <- cbind(rast_coords,rast_data)
  potential_sites[na_sites,-1:-2] <- 0

  #if npoints>nrow(potential_sites) recursively create a finer grid.
  if(npoints>nrow(potential_sites[-na_sites,])){
    stop('more background points than cells avaliable')
  }

  #setup the dimensions need to halton random numbers - let's keep it at 2 for now, could expand to alternative dimension in the future
  dimensions <- 2#dim(potential_sites)[2]
  N <- nrow(potential_sites)
  inclusion_probs <- rep(1/N, N)
  inclusion_probs[na_sites] <- 0
  inclusion_probs1 <- inclusion_probs/max(inclusion_probs)
  mult <- 10
  samp <- randtoolbox::halton(sample(1:10000, 1), dim = dimensions +
                                1, init = TRUE)
  njump <- npoints * mult
  samp <- randtoolbox::halton(max(njump, 5000), dim = dimensions +
                                1, init = FALSE)
  myRange <- apply(potential_sites[-na_sites,], -1, range)
  for (ii in 1:dimensions) samp[, ii] <- myRange[1, ii] + (myRange[2,ii] - myRange[1, ii]) * samp[, ii]
  if (dimensions == 2) {
    tmp <- mgcv::in.out(coordinates(window), samp[, 1:dimensions])
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
                                                                     samp[kount + 1:min(njump, nrow(samp) - kount), -(dimensions + 1), drop = FALSE],
                                                                     1:nrow(potential_sites))
    sampIDs.2 <- which(samp[1:(kount + min(njump, nrow(samp) - kount)), dimensions + 1] < inclusion_probs1[sampIDs[1:(kount +
                                                                                                     min(njump, nrow(samp) - kount))]])
    if (length(sampIDs.2) >= npoints) {
      sampIDs <- sampIDs[sampIDs.2][1:npoints]
      flag <- FALSE
    }
    kount <- kount + njump
  }
  # message("Finished\n")
  if (kount > nrow(samp))
    stop("Failed to find a design. It is likely that the inclusion probabilities are very low and uneven. Please try again OR make inclusion probabilities more even")
  samp <- as.data.frame(cbind(potential_sites[sampIDs, , drop = FALSE],
                              inclusion_probs[sampIDs], sampIDs))
  colnames(samp) <- c(colnames(potential_sites), "inclusion.probabilities","ID")
  return(samp)
}

randomMethod <- function(npoints, window, covariates = NULL){

  if(!inherits(window, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'grid' method currently only works a raster input as a 'window'")

  if(inherits(window, c('RasterLayer','RasterStack','RasterBrick'))){

    if(npoints>sum(!is.na(window[]))){
      stop('Eeek! More background points than cells avaliable... maybe use a finer grid or less points')
    }

    ## use dismo to get random points.
    randpoints <- dismo::randomPoints(window,npoints)

    #create a dataframe of coordinates w/ area
    colnames(randpoints) <- c('X','Y')
  }

  if(!is.null(covariates)){
    covars <- extract(covariates,randpoints)
    randpoints <- cbind(randpoints,covars)
  }

  return(randpoints)
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
    if(!is.null(colnames(known.sites))) colnames(df) <- c("X","Y",colnames(known.sites)[-1:-2])
    colnames(df) <- c('X','Y','SppID')
  } else {
    colnames(df) <- c('X','Y')
  }
  return (df)
}

default_window <- function (presences) {
  # get limits
  xlim <- range(presences$x)
  ylim <- range(presences$y)

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
                      This may affect the results of the model, please check alignment between covariates and area.',sum(which_rm == 0)))
    }

    if (any(which_rm == 1)) {
      warning(sprintf('Removed %i observed points for which covariate values could not be assigned.
                      This may affect the results of the model, please check alignment between covariates and coords.',sum(which_rm == 0)))
    }

  }

  return (pts)

}

estimate_area <- function(window,site_coords){
  if(raster::isLonLat(window)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- raster::area(window)
    area_study <- raster::mask(area_rast,window)
    total_area <- cellStats(area_study,sum,na.rm=TRUE)
    wts <- total_area/raster::extract(area_study,site_coords,na.rm=TRUE)
    } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    cell_area <- raster::res(window)[1]*raster::res(window)[2]
    n_cell <- length(window[!is.na(window[])])
    wts <- rep((n_cell*cell_area)/nrow(site_coords),nrow(site_coords))/1000
    }
   return(wts)
}

#estimate the total area of all cells in extent in km^2
estimate_area_region <- function(window){
  if(raster::isLonLat(window)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- area(window)
    area_study <- mask(area_rast,window)
    area_of_region <- cellStats(area_study,sum, na.rm=TRUE)
  } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    area_of_region <- (ncell(window[!is.na(window)])  * xres(window) * yres(window))/1000
  }
  return(area_of_region)
}



checkPresQuadNames <- function(presences,background){

  return(all(colnames(presences)==colnames(background)))

}



# adjust the resolution to match desired number of background points
guessResolution <- function(npoints,window){
  message('Guessing resolution based on window resolution and approximately ',npoints,' background points')
  reso <- raster::res(window)
  ncello <-  sum(!is.na(window)[])
  newres <- floor(round((ncello*reso[1]))/npoints)
  newres
}

## auto guess reolution if npoints provided.
checkResolution <- function(resolution,window){
  reso <- raster::res(window)
  ncello <-  sum(!is.na(window)[])
  newncell <- round((ncello*reso[1])/resolution)
  if(newncell>500000)stop(message("Hold up, the current resolution of ",resolution," will produce a a grid of approximately ",newncell,", choose a larger resolution. Limit is currently set to 500000 quadrature points"))
  else message("Based on the provided resolution of ",resolution," a grid of approximately ",newncell," quadrature points will be produced.")
}

## check to see if there are duplicated points per species.
## duplicated points are allowed across multiple species ala marked points.
checkDuplicates <- function(presences){
  dups <- duplicated(presences)
  if(sum(dups)>0){ message("There were ",sum(dups)," duplicated points unique to X, Y & SppID, they have been removed.")
  dat <- presences[!dups,]
  } else {
  dat <- presences
  }
  dat
}

## check to see if the presences dataset is multispecies.
checkMultispecies <- function(presences){
  if(length(unique(presences[,"SppID"]))>1) mutlt <- TRUE
  else multi <- FALSE
  multi
}
