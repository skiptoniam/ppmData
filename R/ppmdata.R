#' @name ppmData
#' @title Create a Point Process dataset for spatial presence-only modelling.
#' @description Creates a point process data frame for multiple species (marked) presences. 
#' Generates a quadrature scheme based on Berman & Turner 1992; Warton & Shepard 2010. 
#' The function can generate a quadrature scheme for a regular grid, quasi-random or random points.
#' @export
#' @param npoints Approximate number of background points to generate.
#' @param presences a matrix, dataframe or SpatialPoints* object giving the coordinates of each species' presence in (should be a matrix of nsites * 3) with the three columns being c("X","Y","SpeciesID"), where X is longitude, Y is latitude and SpeciesID is a integer, character or factor which assoicates each point to a species. If presences are NULL then ppmdat will return the quadrature (background) points.
#' @param window a Raster* object giving the area over which to generate background points. NA cells are ignored and masked out of returned data.
#' If ignored, a rectangle defining the extent of \code{presences} will be used.
#' @param covariates an optional Raster object containing covariates for modelling the point process (best use a Raster stack or Raster brick).
# #' @param resolution resolution setup grid for integration points (default is 1 deg) - but this need to be setup  with reference to original raster resolution.
#' @param method the type of method that should be used to generate background points. The options are:
#' 'quasirandom' generates quasirandom background points. See Bratley & Fox 1998 or Foster etal 2015 for details.
#' @param interpolation either 'simple' or 'bilinear' and this determines the interpolation method for interpolating data across different cell resolutions.
#' 'simple' is nearest neighbour, 'bilinear' is bilinear interpolation.
#' @param coord is the name of site coordinates. The default is c('X','Y').
#' @param control \link[qrbp]{ppmData.control}.
#'

ppmData <- function(npoints = 10000,
                    presences = NULL,
                    window = NULL,
                    covariates = NULL,
                    # resolution = NULL,
                    method = c('quasirandom'),#,'grid','psuedorandom'),
                    interpolation='bilinear',
                    coord = c('X','Y'),
                    control=ppmData.control()){

  ## if no resolution is provided guess the nearest resolution to return npoints for grid method.
  # if(is.null(resolution)) resolution <- guessResolution(npoints,window)
  if(method=='quasirandom') control$quasiSamps <- ifelse(control$quasiSamps>npoints,control$quasiSamps,npoints*2)

  ## Do some checks.
  ####  SDF: Do we really need to check for duplicates?  I realise that this will be within species, but roundoff error could kill us?  Or are duplicates assumed to
  ####  the the same observation recorded multiple times in the database?
  ####  Note that this function also depends upon a very particular format for the presences data.  It might be wise to robustify?
#  presences <- checkDuplicates(presences,coord)
  ####

  ####  SDF: I don't know what all this does, so I will assume that it is sensible.
  # checkResolution(resolution,window,control,method)
  window <- checkWindow(presences,window)
  # tmp <- disaggregateWindow(window,npoints)
  # if(!is.null(tmp$ddwindow)) window <- tmp$ddwindow; newres <- tmp$newres

  if(is.null(presences)){
   message('Generating background points in the absence of species presences')
   backgroundpoints <- switch(method,
                             # grid = gridMethod(resolution, window,control),
                             quasirandom = quasirandomMethod(npoints,  window, covariates, control, coord))#,
                             # psuedorandom = randomMethod(npoints, window))

   # wts <- getTileWeights(presences,backgroundpoints[,1:2],coord)
   sitecovariates <- getCovariates(backgroundpoints$bkg_pts[,coord],covariates,
                                   interpolation=interpolation,
                                   coord=coord,control=control)

  } else {

  presences <- coordMatchDim(presences,3)

  # create background points based on method.
  backgroundpoints <- switch(method,
                            # grid = gridMethod(resolution, window,control),
                            quasirandom = quasirandomMethod(npoints,  window, covariates, control, coord))#,
                            # psuedorandom = randomMethod(npoints,  window, covariates))

  ismulti <- checkMultispecies(presences)
  if(ismulti){
      message("Developing a quadrature scheme for multiple species (marked) dataset.")
      wts <- getMultispeciesWeights(presences, backgroundpoints$bkg_pts, coord, method, areaMethod=control$weightMethod, window, epsilon=control$epsilon)
    } else {
      ####  Will need to look at areaMethod -- hard-wired for now.
      wts <- getSinglespeciesWeights(presences, backgroundpoints$bkg_pts, coord, method, window, epsilon=control$epsilon)
    }

  # pbxy <- wts[,coord],backgroundpoints[,coord])
  sitecovariates <- getCovariates(wts,covariates,interpolation=interpolation,
                                  coord=coord,control=control)
  }

  parameters <- list(npoints=npoints,#resolution=resolution,
                     newresolution=backgroundpoints$newres,method=method,
                     interpolation=interpolation,control=control)
  dat <- assembleQuadData(presences, backgroundpoints$bkg_pts, sitecovariates, wts,
                          coord, parameters, control=control)
  class(dat) <- "ppmdata"
  return(dat)
}

#'@title Controls for ppmData generation.
#'@name ppmData.control
#'@param quiet Should any reporting be performed? Default is FALSE, for reporting.
#'@param cores The number of cores to use in fitting of species_mix models. These will be largely used to model the species-specific parameteres.
#'@param maxpoints default is 250000 and sets a limit to number of integration points generated as background data.
#'@param na.rm default is TRUE and uses "na.rm=TRUE" for extracting raster covariate data.
#'@param extractBuffer default is NULL and is the amount of buffer to provide each point on extract (radius from point).
#'@param quasiSamps default is 5000 and is the default number of samples for halton random number generator.
#'@param quasiDims default 2 and is the dimension estimated the quasirandom samples over. Two is the default and results in a spatial quasirandom sample.
#'@param multispeciesFormat default "wide" and is the format required by ecomix. Alternatives are "long" which returns a long table format or "list" which returns a list per species for all species related covariates, weights and presences/background points.
#'@param \dots Other control calls.
#'@export
ppmData.control <- function(quiet = FALSE,
                            cores = 1,
                            maxpoints=250000,
                            na.rm=TRUE,
                            extractBuffer=NULL,
                            quasiSamps=5000,
                            weightMethod="quick",
                            quasiDims=2,
                            multispeciesFormat="wide",
                            epsilon = sqrt(.Machine$double.eps),
                            ...){
  #general controls
  rval <- list(maxpoints=maxpoints,
               na.rm=na.rm,
               extractBuffer=extractBuffer,
               quasiSamps=quasiSamps,
               weightMethod=weightMethod,
               quasiDims=quasiDims,
               multispeciesFormat=multispeciesFormat,
               epsilon = epsilon)
  rval <- c(rval, list(...))
  rval
}



assembleQuadData <- function(presences, backgroundpoints, sitecovariates,
                             wts, coord, parameters, control){

  ismulti <- checkMultispecies(presences)

  if(!ismulti) format <- "long"
  else format <- control$multispeciesFormat

  final_dat <- switch(format,
                      long=longdat(presences, backgroundpoints,
                                   sitecovariates,
                                   wts, coord),
                      wide=widedat(presences, backgroundpoints,
                                   sitecovariates,
                                   wts, coord),
                      list=listdat(presences, backgroundpoints,
                                   sitecovariates,
                                   wts, coord))

  return(list(modelmatrix=final_dat,parameters=parameters))

}


longdat <- function(presences, backgroundpoints, sitecovariates=NULL, wts, coord){

  if(!is.null(sitecovariates)){
      dat2 <- merge(wts[,-which(colnames(wts)%in%coord)],sitecovariates[!duplicated(sitecovariates$SiteID),],
                  by = "SiteID", sort=FALSE)
    } else {
      dat2 <- wts
    }
  if(is.null(presences)){ #presences true
      if(!is.null(sitecovariates)){ # covariates true
        dat2 <- cbind(backgroundpoints[,c(coord)],sitecovariates,pres=c(rep(0,nrow(backgroundpoints))))
      } else { # covariates false
        dat2 <- cbind(backgroundpoints[,c(coord)],pres=c(rep(0,nrow(backgroundpoints))))
      }
  }
  return(dat2)
}

listdat <- function(presence=NULL, backgroundpoints, sitecovariates, wts=NULL, coord){

  if(!is.null(presences)){ #presences true
    if(!is.null(sitecovariates)){ #covariates true
      dat2 <- list(presences=presences[,c(coord,"SpeciesID")],
                   background=backgroundpoints[,c(coord,"SpeciesID")],
                   covariates=sitecovariates,
                   pres=c(rep(1,nrow(presences)),rep(0,nrow(backgroundpoints))),
                   wts=wts)
    } else { #covariates false
      dat2 <- list(presences=presences[,c(coord,"SpeciesID")],
                   background=backgroundpoints[,c(coord,"SpeciesID")],
                   covariates=NULL,
                   pres=c(rep(1,nrow(presences)),rep(0,nrow(backgroundpoints))),
                   wts=wts)
    }
  } else { # presences false
    if(!is.null(sitecovariates)){
      dat2 <- list(presences=NULL,
                   background=backgroundpoints[,c(coord,"SpeciesID")],
                   covariates=sitecovariates,
                   pres=c(rep(0,nrow(backgroundpoints))),
                   wts=wts)
      } else { # covariates false
      dat2 <- list(presences=NULL,
                   background=backgroundpoints[,c(coord,"SpeciesID")],
                   covariates=NULL,
                   pres=c(rep(0,nrow(backgroundpoints))),
                   wts=wts)
    }
  }
  return(dat2)
}

widedat <- function(presence, backgroundpoints, sitecovariates, wts, coord){

  # Assemble a data.frame with all the bits we want.
  pamat <- fastwidemat(wts)
  presences_pamat <- pamat[-which(pamat[,"quad"]==0),-which(colnames(pamat)=='quad')]
  presences_pamat[presences_pamat==0]<-NA
  quad_pamat <- pamat[which(pamat[,"quad"]==0),-which(colnames(pamat)=='quad')]
  quad_pamat[is.na(quad_pamat)]<-0
  response_ppmmat <- as.data.frame(rbind(presences_pamat,quad_pamat))
  response_ppmmat$Const <- 1
  response_ppmmat$SiteID <- as.integer(rownames(response_ppmmat))

  df <- merge(response_ppmmat,sitecovariates[!duplicated(sitecovariates$SiteID),],
              by = "SiteID", sort=FALSE)
  wtsmat <- fastwidematwts(wts)
  ids <- wts[!duplicated(wts[,c('SpeciesID','DatasetID')]),c('SpeciesID','DatasetID')]
  ids <- ids[-which(ids$SpeciesID=='quad'),]
  idx_rows <- df$SiteID
  idx_cols <- match(colnames(quad_pamat),ids$SpeciesID)
  wtsmat <- wtsmat[idx_rows,idx_cols]
  colnames(wtsmat) <- colnames(quad_pamat)
  return(list(mm=df,wtsmat=wtsmat))
}

# gridMethod <- function(resolution=1, window, control){
# 
#   if(!inherits(window, c('RasterLayer','RasterStack','RasterBrick')))
#     stop("'grid' method currently only works a raster input as a 'window'")
# 
#   if(inherits(window, c('RasterLayer','RasterStack','RasterBrick'))){
# 
#     #set up the dissaggreation or aggregate
#     fct <- (res(window)/resolution)[1]
# 
#     #if fct is >= 1 dissaggregate, else aggregate
#     if(fct>=1) dd <- disaggregate(window, fct, na.rm=control$na.rm)
#     else dd <- aggregate(window, 1/fct, na.rm=control$na.rm)
# 
#     #create a dataframe of coordinates w/ area
#     grid <- as.data.frame(rasterToPoints(dd)[,-3])
#     colnames(grid) <- c('X','Y')
#   }
# 
#   newres <- res(dd)
# 
#   return( list( bkg_pts = as.data.frame( randpoints), newres = newres))
# 
# }

quasirandomMethod <- function(npoints, window, covariates=NULL, control, coord){

  #generate a set of potential sites for quasirandom generation
  if(!is.null(covariates)){
    rast_coord <- raster::xyFromCell(covariates,1:raster::ncell(covariates))
    rast_data <- raster::values(covariates)
    na_sites <- which(!complete.cases(rast_data))
    covariates_ext <- raster::extent(covariates)[1:4]
    potential_sites <- cbind(rast_coord,rast_data)
  } else {
    rast_coord <- raster::xyFromCell(window,1:raster::ncell(window))
    rast_data <- raster::values(window)
    na_sites <- which(!complete.cases(rast_data))
    covariates_ext <- raster::extent(window)[1:4]
    potential_sites <- rast_coord
  }
  
  ## setup the inclusion probs.
  N <- nrow(potential_sites)
  inclusion_probs <- rep(1/N, N)	
  inclusion_probs1 <- inclusion_probs/max(inclusion_probs)
  inclusion_probs1[na_sites] <- 0	

  samp <- MBHdesign::quasiSamp(n = npoints, dimension = control$quasiDim,
                               potential.sites = potential_sites[,1:control$quasiDim],
                               inclusion.probs = inclusion_probs1, nSampsToConsider = control$quasiSamps)

  colnames( samp[,1:2]) <- coord
  if(!is.null(covariates)){
    covars <- raster::extract(covariates,samp[,1:2])
    randpoints <- cbind(samp[,1:2],covars)
  } 

  return( list( bkg_pts = as.data.frame( randpoints)))
}

# randomMethod <- function(npoints, window, covariates = NULL){
# 
#   if(!inherits(window, c('RasterLayer','RasterStack','RasterBrick')))
#     stop("'grid' method currently only works a raster input as a 'window'")
# 
#   if(inherits(window, c('RasterLayer','RasterStack','RasterBrick'))){
# 
#     if(npoints>sum(!is.na(window[]))){
#       stop('Eeek! More background points than cells avaliable... maybe use a finer grid or less points')
#     }
# 
#     ## use dismo to get random points.
#     randpoints <- dismo::randomPoints(window,npoints)
# 
#     #create a dataframe of coordinates w/ area
#     colnames(randpoints) <- c('X','Y')
#   }
# 
#   if(!is.null(covariates)){
#     covars <- extract(covariates,randpoints)
#     randpoints <- cbind(randpoints,covars)
#   }
# 
#   return( list( bkg_pts = as.data.frame( randpoints), newres = raster::res( window)))
# }

coordMatchDim <- function (known.sites,dimension){

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
    colnames(df) <- c('X','Y','SpeciesID')
  } else {
    colnames(df) <- c('X','Y')
  }
  return (df)
}


## function to extract covariates for presence and background points.
getCovariates <- function(pbxy, covariates=NULL, interpolation, coord, control){
  if(is.null(covariates))return(NULL)
  covars <- raster::extract(covariates, pbxy[,coord], method=interpolation,
                            na.rm=control$na.rm, buffer=control$extractBuffer)
  covars <- cbind(SiteID=pbxy[,"SiteID"],pbxy[,coord],covars)
  return(covars)
}

# adjust the resolution to match desired number of background points
# guessResolution <- function(npoints,window){
#   message('Guessing resolution based on window resolution and approximately ',npoints,' background points')
#   reso <- raster::res(window)
#   ncello <-  sum(!is.na(window)[])
#   newres <- (((ncello*reso[1]))/npoints)
#   newres
# }

## Some checks, check yo self before you reck yo self. https://www.youtube.com/watch?v=bueFTrwHFEs
## check the covariates that go into building the quadrature scheme.
checkCovariates <- function(covariates){

  if(!is.null(covariates)){
    if(!inherits(covariates, c('RasterLayer','RasterStack','RasterBrick')))
      stop("Covariates must be a raster, rasterstack or rasterbrick of covariates which match the spatial window.")
    covars <- TRUE
  } else {
    covars <- FALSE
  }
  covars
}

## check to see if there are duplicated points per species.
## duplicated points are allowed across multiple species ala marked points.
checkDuplicates <- function(presences,coord){
  if(is.null(presences))return(NULL)
  dups <- duplicated(presences)
  if(sum(dups)>0){ message("There were ",sum(dups)," duplicated points unique to X, Y & SpeciesID, they have been removed.")
  dat <- presences[!dups,]
  } else {
  dat <- presences
  }
  dat <- as.data.frame(dat)
  colnames(dat) <- c(coord,"SpeciesID")
  dat
}

## check to see if the presences dataset is multispecies.
checkMultispecies <- function(presences){
  if(length(unique(presences[,"SpeciesID"]))>1) multi <- TRUE
  else multi <- FALSE
  multi
}

## check resolution and throw error if lots of quad points to be generated.
# checkResolution <- function(resolution,window,control,method){
#   reso <- raster::res(window)
#   ncello <-  sum(!is.na(window)[])
#   fct <- (reso/resolution)
#   fctprod <- prod(fct)
#   newncell <- ncello/(1/fctprod)
# if(method%in%"grid"){
#     if(newncell>control$maxpoints)stop(message("Hold up, the current resolution of ",resolution," will produce a a grid of approximately ",round(newncell),",\n choose a larger resolution. Limit is currently set to 500000 quadrature points"))
#   else message("Based on the provided resolution of ",resolution," a grid of approximately ",round(newncell)," quadrature points will be produced.")
#   }
# }


checkWindow <- function(presences,window,coord){

  if(is.null(presences)){
    presences <- data.frame(X=runif(100,0,100),Y=runif(100,0,100))
  }

  if (is.null(window)) {
    message("Window is NULL, a raster-based window will be generated based on the extent of 'presences'\n with a default resolution of 1 deg will be returned.")
    window <- defaultWindow(presences,coord)
  }
  window
}

defaultWindow <- function (presences,coord) {
  # get limits
  xlim <- range(presences[,coord[1]])
  ylim <- range(presences[,coord[1]])

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


fastwidemat <- function(dat){

  dat[,"SiteID"] <- factor(dat[,"SiteID"])
  dat[,"SpeciesID"] <- factor(dat[,"SpeciesID"])

  result <- with(dat, {
    out <- matrix(nrow=nlevels(SiteID), ncol=nlevels(SpeciesID),
                  dimnames=list(levels(SiteID), levels(SpeciesID)))
    out[cbind(SiteID, SpeciesID)] <- pres
    out
  })

  return(result)
}


fastwidematwts <- function(dat){

  dat[,"SiteID"] <- factor(dat[,"SiteID"])
  dat[,"DatasetID"] <- factor(dat[,"DatasetID"])
  # dat[,"SpeciesID"] <- factor(dat[,"SpeciesID"])

  wtsdat <- with(dat, {
    out <- matrix(nrow=nlevels(SiteID), ncol=nlevels(DatasetID),
                  dimnames=list(levels(SiteID), levels(DatasetID)))
    out[cbind(SiteID, DatasetID)] <- wts
    out
  })

  wtsdat
}

# a little function for making a finer scale window.
# disaggregateWindow <- function(window, npoints){
# 
#   #set up the dissaggreation or aggregate
#   nc <- length(window[!is.na(window[])])
# 
#   fct <- (nc/npoints)
#   if(fct<1){
#     dd <- disaggregate(window, 1/fct, na.rm=control$na.rm)
#     newres <- res(dd)
#   } else {
#     dd <- NULL
#     newres <- NULL
#   }
# 
#   return(list(ddwindow = dd, newres = newres))
# 
# }




