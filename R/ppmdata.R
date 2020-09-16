#' @name ppmData
#' @title Create a Point Process dataset for spatial modelling.
#' @description Creates a point process data frame for multiple species (marked) presences. 
#' Generates a quadrature scheme based on Berman & Turner 1992; Warton & Shepard 2010. 
#' The function can generate a quadrature scheme for a regular grid, quasi-random or random points.
#' @export
#' @param npoints The number of background points to generate.
#' @param presences a matrix, dataframe or SpatialPoints* object giving the coordinates of each species' presence in (should be a matrix of nsites * 3) 
#' with the three columns being c("X","Y","SpeciesID"), where X is longitude, Y is latitude and SpeciesID is a integer, character or factor which assoicates each point to a species.
#' If presences parameter is NULL then ppmDat will return the quadrature (background) points.
#' @param window a Raster* object giving the area over which to generate background points. NA cells are ignored and masked out of returned data.
#' If ignored, a rectangle defining the extent of \code{presences} will be used.
#' @param covariates A Raster* object containing covariates for modelling the point process (best use a Raster stack or Raster brick).
# #' @param resolution resolution setup grid for integration points (default is 1 deg) - but this need to be setup  with reference to original raster resolution.
## #' @param method the type of method that should be used to generate background points. The options are:
## #' 'quasirandom' generates quasirandom background points. See Bratley & Fox 1998 or Foster etal 2015 for details.
#' @param interpolation either 'simple' or 'bilinear' and this determines the interpolation method for interpolating data across different cell resolutions.
#' 'simple' is nearest neighbour, 'bilinear' is bilinear interpolation.
#' @param coord is the name of site coordinates. The default is c('X','Y').
#' @param mc.cores The number of cores to use in the processing to quasirandom points and weighting scheme. 
#' @param quasirandom.samples This set the total number of samples to consider in the BAS step (rejection sampling). 
#' The default is 10000, which means that 10000 halton numbers are drawn and then thinned according to the inclusion probabilities.
#' You will need to increase the number of samples if sampling across > 2 dimensions or selecting a large number of background points. 
#' The more quasirandomSample selected the slower the background point generation will be.
#' @param quasirandom.dimensions The the number of dimensions that the samples are located in. Equal to 2 for areal sampling. 
#' Care should be taken with large dimensions as :1) the number of potential sampling sites needed for effective coverage starts to explode (curse of dimensionality);
#' and 2) the well-spaced behaviour of the Halton sequence starts to deteriorate.

ppmData <- function(npoints = 10000,
                    presences = NULL,
                    window = NULL,
                    covariates = NULL,
                    interpolation='bilinear',
                    coord = c('X','Y'),
                    mc.cores = parallel::detectCores()-1,
                    quasirandom.samples = 10000,
                    quasirandom.dimensions = 2){

  ## if no resolution is provided guess the nearest resolution to return npoints for grid method.
  # if(is.null(resolution)) resolution <- guessResolution(npoints,window)
  # if(method=='quasirandom') control$quasiSamps <- ifelse(control$quasiSamps>npoints,control$quasiSamps,npoints*4)

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
   # backgroundpoints <- switch(method,
   #                            grid = gridMethod(resolution, window,control),
   #                            quasirandom = quasirandomMethod(npoints,  window, covariates, control, coord),
   #                            psuedorandom = randomMethod(npoints, window))
   backgroundpoints <- quasirandomMethod(npoints,  window, covariates, coord, mc.cores,
                                         quasirandom.samples,
                                         quasirandom.dimensions)
   
   sitecovariates <- getCovariates(backgroundpoints$bkg_pts[,coord],covariates,
                                   interpolation=interpolation,
                                   coord=coord,control=control)

  } else {

  presences <- coordMatchDim(presences,3)

  # create background points based on method.
  # backgroundpoints <- switch(method,
  #                            grid = gridMethod(resolution, window,control),
  #                            quasirandom = quasirandomMethod(npoints,  window, covariates, control, coord),
  #                            psuedorandom = randomMethod(npoints, window))
  backgroundpoints <- quasirandomMethod(npoints,  window, covariates, control, coord)
  
  ismulti <- checkMultispecies(presences)
  if(ismulti){
      message("Developing a quadrature scheme for multiple species (marked) dataset.")
      wts <- getMultispeciesWeights(presences, backgroundpoints$bkg_pts, coord, method, areaMethod=control$weightMethod, window, epsilon=control$epsilon)
    } else {
      message("Developing a quadrature scheme for a single species dataset.")
      ####  Will need to look at areaMethod -- hard-wired for now.
      wts <- getSinglespeciesWeights(presences, backgroundpoints$bkg_pts, coord, method, window, epsilon=control$epsilon)
    }

  sitecovariates <- getCovariates(wts,covariates,interpolation=interpolation,
                                  coord=coord,control=control)
  }

  parameters <- list(npoints=npoints,#resolution=resolution,
                     newresolution=backgroundpoints$newres,
                     # method=method,
                     interpolation=interpolation,
                     control=control)
  dat <- assembleQuadData(presences, backgroundpoints$bkg_pts, sitecovariates, wts,
                          coord, parameters, control=control)
  class(dat) <- "ppmdata"
  return(dat)
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

quasirandomMethod <- function(npoints, window, covariates=NULL, coord,
                              mc.cores,
                              quasirandom.samples,
                              quasirandom.dimensions){
  
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
  
  ####SDF CHANGE
  ids <- inclusion_probs1 > 0
  potential_sites <- potential_sites[ids,]
  inclusion_probs1 <- inclusion_probs1[ids]
  quasirandom.samples <- ifelse(npoints<1000,10000,10*npoints)
  ####
  
  samp <- suppressMessages(MBHdesign::quasiSamp(n = npoints, dimension = quasirandom.dimensions,
                               potential.sites = potential_sites[, seq_len(quasirandom.dimensions)],
                               inclusion.probs = inclusion_probs1, nSampsToConsider = quasirandom.samples))
  
  colnames( samp[,1:2]) <- coord
  if(!is.null(covariates)){
    covars <- raster::extract(covariates,samp[,1:2])
    randpoints <- cbind(samp[,1:2],covars)
  } 
  
  return( list( bkg_pts = as.data.frame( randpoints)))
}

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






