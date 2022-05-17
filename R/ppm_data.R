#' @name ppmData
#' @title Create a point process quadrature scheme for spatial modelling.
#' @description This package is a way to efficiently generate a quasirandom set
#' of background points for presence-only modelling of single or multiple
#' responses. The package was set up to model multiple species presence-only
#' datasets, but could be used for an point process spatial modelling.
#' Quasi-random points are a nice alternative to pseudo-random samples, this is
#' because we can generate a quasirandom sample across and areal region
#' (X and Y coordinates). This in turn should reduce autocorrelation in
#' quadrature scheme. The weight of each quadrature point is calculated using
#' Dirichlet (Voronoi) Tessellation written in c++ based on the duel-graph of a
#' Delaunay triangulation which is constructed based on the sweep algorithm.
#' @details The approach uses quasi-random sampling to generate a quadrature
#' scheme based (e.g Berman & Turner 1992; Warton & Shepard 2010;
#' Foster et al, 2017). The weights each quasi-random point in the quadrature
#' scheme is calculated using a Dirichlet tessellation (Turner 2020). To improve
#' computational efficiency we have rewritten the Delaunay triangulation and
#' Dirichlet tessellation in c++ using the sweep algorithm of .
#' @export
#' @param presences a matrix, dataframe or SpatialPoints object giving the
#' coordinates of each species' presence in (should be a matrix of nsites * 3)
#' with the three columns being c("X","Y","SpeciesID"), where X is longitude,
#' Y is latitude and SpeciesID is factor which associated each occurrence to a
#' species. If presences parameter is NULL then ppmDat will return the
#' quadrature scheme without presences.
#' @param window a raster object giving the area over which to generate
#' background points. NA cells are ignored and masked out of returned data.
#' If NULL, a rectangle bounding the extent of \code{presences} will be used as
#' the default window.
#' @param covariates A raster object containing covariates for modelling the
#' point process (best use a Raster stack or Raster brick). This should match
#' the resolution and extent of the window provided. If NULL, only the
#' coordinates of the presences and quadrature points are returned.
#' @param npoints The number of quadrature points to generate. If NULL, the
#' number of quadrature points is calculated based on a linear scaling. In
#' reality, the number of quadrature points needed to approximate the
#' log-likelihood function will depend on the data and likelihood function
#' being approximated. Typically, more quadrature the better the estimate,
#' but there is a trade off between computational efficiency and accuracy.
#' See Warton & Shepard (2010) or Renner et al., 2015 for useful discussions
#' on the location and number of quadrature points required to converge a
#' ppm likelihood.
#' @param coord is the name of site coordinates. The default is c('X','Y').
#' This should match the name of the coordinates in the presences dataset.
#' @param species.id is the name of the species ID in the presences dataset.
#' The default is "SpeciesID".
#' @param quad.method The quadrature generation method. Default is "quasi" for
#' quasi-random, "random" for pseudo-random (regular random) and "grid" for
#' a regular grid at a set resolution (with respect to the original window resolution).
#' @param mc.cores The number of cores to use in the processing. default is
#' parallel::detectCores()-1
#' @param quasirandom.samples This set the total number of samples to consider
#' in the BAS step (rejection sampling). The default is 10000, which means that
#' 10000 halton numbers are drawn and then thinned according to the inclusion
#' probabilities. You will need to increase the number of samples if selecting
#' a large number of quadrature points. The more quasirandomSamples selected the
#' slower the quasirandom quadrature scheme will be to generate.
#' @param interpolation The interpolation method to use when extracting
#' covariate data. Default is "bilinear", can also use "simple", this is based
#' on the raster package  \code{\link[raster]{extract}}.
#' @param unit The scale of the area weights, default is "geo" and returns the area based on the geographic coorinates
#' alternatively meters squared "m", kilometers squared "km", or hectars "ha" can be used (this is still experimental).
#' @param bufferNA If extract from \code{\link[raster]{extract}} returns NA for
#' point extract, do you want us to attempt to user buffer to calculate cells
#' which are NA.
#' @param bufferSize If you call 'bufferNA' what is the range of the buffer in
#' meters.
#' @param quiet If TRUE, do not print messages. Default is FALSE.
#' @importFrom graphics legend par points
#' @importFrom methods as
#' @importFrom stats complete.cases median runif
#' @importFrom utils txtProgressBar
#' @importFrom raster extract ncell raster extent extract res
#' @examples
#' \dontrun{
#' library(ppmData)
#' library(terra)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' preds <- rast(lst)
#' window <- preds[[1]]
#' presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#' bkgrid <- ppmData(npoints = 1000, presences=presences, window = window, covariates = preds)
#' }

ppmData <- function(presences,
                    window = NULL,
                    covariates = NULL,
                    npoints = NULL,
                    coord = c('X','Y'),
                    species.id = "SpeciesID",
                    quad.method = c("quasi","random","grid"),
                    mc.cores = 1,
                    quasirandom.samples = NULL,
                    interpolation = c("simple","bilinear"),
                    unit = c("geo","m","km","ha"),
                    bufferNA = FALSE,
                    bufferSize = NULL,
                    quiet=FALSE){

  # default methods
  quad.method <- match.arg(quad.method)
  interp.method <- match.arg(interpolation)
  unit <- match.arg(unit)

  ## Make sure the column ids are characters and check for missing/wrong named coord/species.id vars.
  if(!is.character(coord)) coord <- as.character(coord)
  if(all(!coord%in%colnames(presences))) stop(paste0('coord: "',coord[1],'" & "',coord[2],'" do not match any of the colnames in the presences data.frame'))
  if(!is.character(species.id)) species.id <- as.character(species.id)
  if(all(!species.id%in%colnames(presences))) stop(paste0('species.id: "',species.id,'" does not match any of the colnames in the presences data.frame'))

  # This should check the presences and make it returns the data in the correct format for the remaining function.
  pressies <- checkPresences(known.sites = presences,
                             window = window,
                             coord = coord,
                             species.id = species.id)

  # Check for duplicate presences - will remove duplicated points per species.
  pressies <- checkDuplicates(presences = pressies,
                              coord = coord,
                              species.id = species.id,
                              quiet = quiet)

  ## If npoints in NULL setup a default amount. This is taken from spatstat
  npoints <- ChechNumPoints(npoints = npoints,
                           presences = pressies,
                           species.id = species.id)

  ## If not window is provided provide a dummy window
  if(is.null(window)) default_window <- TRUE
  else default_window <- FALSE
  window <- checkWindow(presences = pressies,
                        window = window,
                        coord = coord,
                        quiet = quiet)

  ## Hold onto the species names from the species.id column
  sppNames <- getSppNames(presences = pressies,
                          species.id = species.id)

  ## Create some quadrature points
  bckpts <- quadMethod(quad.method = quad.method,
                       npoints = npoints,
                       window = window,
                       coord =  coord,
                       quasirandom.samples = quasirandom.samples)

  ## Sometimes the points are very close together, so let's jitter if needed.
  reswindow <- terra::res(window)[1]
  tmpPts <- jitterIfNeeded(pressiesJ = pressies, bckpts=bckpts$quasiPoints,
                           window=window,coord=coord,
                           species.id = species.id,
                           aBit=reswindow/2)

  pressies <- tmpPts$pressies
  bckptsQ <- tmpPts$bckpts
  bckptsD <- bckpts$quasiDummy

  # Check to see if the presences are for a single species or multiple.
  ismulti <- checkMultispecies(pressies, species.id)

  if(ismulti){
    if(!quiet)message("Developing a quadrature scheme for multiple species (marked) dataset.")
      wts <- getMultispeciesWeights(quad.method = quad.method,
                                    presences = pressies,
                                    quadrature = bckptsQ,
                                    quadDummy = bckptsD,
                                    window = window,
                                    coord = coord,
                                    speciesIdx = species.id,
                                    mc.cores = mc.cores,
                                    sppNames = sppNames,
                                    unit=unit)

      sitecovariates <- getCovariates(pbxy = wts,
                                      covariates,
                                      interpolation=interpolation,
                                      coord=coord,
                                      bufferNA = bufferNA,
                                      bufferSize = bufferSize)

    } else {
      if(!quiet)message("Developing a quadrature scheme for a single species dataset.")
      wts <- getSingleSpeciesWeights(quad.method= quad.method,
                                     presences = pressies,
                                     quadrature = bckptsQ,
                                     quadDummy = bckptsD,
                                     window = window,
                                     coord = coord,
                                     speciesIdx = species.id,
                                     unit=unit)
      # extract the covariate data
      sitecovariates <- getCovariates(pbxy = wts,
                                      covariates = covariates,
                                      interpolation=interpolation,
                                      coord=coord,
                                      bufferNA = bufferNA,
                                      bufferSize = bufferSize)
    }

  # Assemble data
  dat <- assembleQuadData(presences = pressies,
                          quadrature = bckptsQ,
                          sitecovariates = sitecovariates,
                          wts = wts,
                          coord = coord,
                          species.id = species.id,
                          sppNames = sppNames)

  if(!is.null(covariates)){
    covarNames <- names(covariates)
  } else {
    covarNames <- NULL
  }
  coordNames <- coord
  res <- list()
  if(ismulti){
     res$ppmData <- transposePPMdata(dat, sppNames, coordNames, covarNames)
     res$marked <- TRUE
  } else {
    res$ppmData <- dat
    res$marked <- FALSE
  }

  if(!is.null(presences)) res$presences <- presences
  res$window <- window
  res$params <- list(coord = coord,
                     species.id = species.id,
                     mc.cores = mc.cores,
                     quasirandom.samples = quasirandom.samples,
                     interpolation = interpolation,
                     dw = default_window)

  class(res) <- c("ppmData")

  if(!quiet)print(res)
  return(res)
}

jitterIfNeeded <- function( pressiesJ, bckpts, window, coord, species.id, aBit=1e-4){
  #the pressie bit first
  #are there any duplicates within a species?  If so, then jitter the duplicates
  for( jj in as.character( unique( pressiesJ[,species.id]))){ #I think that we have made the assumption that this is called SpeciesID...?
    sppJ <- which(  pressiesJ[,species.id]==jj)
    dupes <- which( duplicated( pressiesJ[sppJ,coord]))  #shouldn't need to round this as deldir reportedly uses duplicated

    if( length( dupes)>0){
      pressiesJ[dupes,coord[1]] <- jitter( pressiesJ[sppJ[dupes],coord[1]], amount=aBit)
      pressiesJ[dupes,coord[2]] <- jitter( pressiesJ[sppJ[dupes],coord[2]], amount=aBit)
      #fixing up those points that have been jittered outside of the window
      kount <- 1
      tmp <- terra::extract(window, pressiesJ[,coord])
      outOfWindow <- which( is.na( tmp))
      while( kount < 10 & length( outOfWindow)>0){
      	pressiesJ[outOfWindow,coord[1]] <- jitter( pressiesJ[sppJ[outOfWindow],coord[1]], amount=aBit)
      	pressiesJ[outOfWindow,coord[2]] <- jitter( pressiesJ[sppJ[outOfWindow],coord[2]], amount=aBit)
      	tmp <- terra::extract( window, pressiesJ[,coord])
      	outOfWindow <- which( is.na( tmp))
      	kount <- kount + 1
      }
    }
  }

  #background now
  #are there any points that duplicate a presence point?  If so, jitter.
  npres <- nrow( pressiesJ)
  dupes <- which( duplicated( rbind( pressiesJ[,coord], bckpts)))
  if( length( dupes)>0)
    dupes <- dupes[dupes>npres]
  if( length( dupes)>0){
    dupes <- dupes - npres
    bckpts[dupes,coord[1]] <- jitter( bckpts[dupes,coord[1]], amount=aBit)
    bckpts[dupes,coord[2]] <- jitter( bckpts[dupes,coord[2]], amount=aBit)
    #fixing up those points that have been jittered outside of the window
    kount <- 1
    tmp <- terra::extract( window, bckpts[,coord])
    outOfWindow <- which( is.na( tmp))
    while( kount < 10 & length( outOfWindow)>0){
      bckpts[outOfWindow,coord[1]] <- jitter( bckpts[,coord[1]], amount=aBit)
      bckpts[outOfWindow,coord[2]] <- jitter( bckpts[,coord[2]], amount=aBit)
      tmp <- terra::extract( window, bckpts[,coord])
      outOfWindow <- which( is.na( tmp))
      kount <- kount + 1
    }
  }

  return( list( pressies=pressiesJ, bckpts=bckpts))
}

assembleQuadData <- function(presences, quadrature, sitecovariates, wts, coord, species.id, sppNames){

  ismulti <- checkMultispecies(presences, species.id)

  if(!ismulti) type <- "long"
  else type <- "wide"

  final_dat <- switch(type,
                      long=longData(wts = wts,
                                   sitecovariates = sitecovariates,
                                   coord = coord),
                      wide=wideData(presence = presences,
                                   quadrature = quadrature,
                                   sitecovariates = sitecovariates, wts = wts,
                                   coord = coord,
                                   species.id = species.id,
                                   sppNames = sppNames))

  return(final_dat)

}


longData <- function(wts, sitecovariates=NULL, coord){

  if(!is.null(sitecovariates)) dat2 <- data.frame(wts[,coord],sitecovariates[,-1:-3],presence=wts$pres, weights=wts$wts)
  else dat2 <- data.frame(wts[,coord],presence=wts$pres,weights=wts$wts)

  return(dat2)
}

wideData <- function(presence, quadrature, sitecovariates, wts, coord, species.id, sppNames){

  # Assemble a data.frame with all the bits we want.
  pamat <- fastWideMatrix(wts, species.id)
  presences_pamat <- pamat[-which(pamat[,"quad"]==0),-which(colnames(pamat)%in%c('quad','dummy'))]
  quad_pamat <- pamat[which(pamat[,"quad"]==0),-which(colnames(pamat)%in%c('quad','dummy'))]
  quad_pamat[is.na(quad_pamat)]<-0
  response_ppmmat <- as.data.frame(rbind(presences_pamat,quad_pamat))
  response_ppmmat$Const <- 1
  response_ppmmat$OrigOrder <- as.integer(rownames(response_ppmmat))


  sitecovariates$OrigOrder <- wts$OrigOrder
  df <- merge(response_ppmmat,sitecovariates[!duplicated(sitecovariates$OrigOrder),], by = "OrigOrder", sort=TRUE)
  wtsmat <- fastWideMatrixWeights(wts, sppNames)

    if(!all.equal(colnames(wtsmat),colnames(quad_pamat))) message('names could be getting mixed up with factors in R')
  return(list(mm=df,wtsmat=wtsmat))
}


quadMethod <- function(quad.method, npoints, window, coord, quasirandom.samples=NULL){
  quad <- switch(quad.method,
                  quasi = quasiRandomQuad(npoints,
                                          window,
                                          coord,
                                          quasirandom.samples),
                  random = pseudoRandomQuad(npoints,
                                            window,
                                            coord),
                  grid = gridQuad(npoints,
                                  window,
                                  coord))
  return(quad)
}



ChechNumPoints <- function(npoints, presences, species.id){

  if(is.null(npoints)){

    ## Taken from spatstat.
    ## linear scaling of quad points compared to max n presences.
    npmx <- max(table(presences[,species.id]))
    nquad <- rep(pmax(32, 10 * ceiling(2 * sqrt(npmx)/10)),2)
    npoints <- prod(nquad)
  }
  return(npoints)
}

checkPresences <- function (known.sites, window, coord, species.id){

  # check for null sites
  if(is.null(known.sites))
    stop("This function requires a set of presences to run.\n")

  # check object classes
  expectClasses(known.sites, c('matrix','data.frame'),name = 'known.sites')

  if(!any(colnames(known.sites)%in%species.id))
    stop("'species.id': ",species.id," ,does not match any of the column names in your presences data.\n")
  if(!any(colnames(known.sites)%in%coord))
    stop("The coordinates names: ",paste(coord,collapse = ", ")," do not match any of the column names in your presences data.\n")

  #check to see if the points lie within the raster (not over an NA)
  # tmpCheck <- terra::extract( window, known.sites[,coord])
  # if( !all( !is.na( tmpCheck))){
    # warning("There are presence records outside the region window.  Removing them but please check.")
    # known.sites <- known.sites[!is.na( tmpCheck),]
  # }

  # try and sort out factors.
  known.sites[[species.id]] <- factor(known.sites[[species.id]], levels = unique(known.sites[[species.id]]))
  df <- known.sites[,c(coord,species.id)]

  return (df)
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
checkDuplicates <- function(presences, coord, species.id, quiet){
  if(is.null(presences))return(NULL)
  dups <- duplicated(presences)
  if(sum(dups)>0){ if(!quiet)message("There were ",sum(dups)," duplicated points unique to X, Y & SpeciesID, they have been removed.")
  dat <- presences[!dups,]
  } else {
  dat <- presences
  }
  dat <- as.data.frame(dat)
  colnames(dat) <- c(coord, species.id)
  dat
}

## check to see if the presences dataset is multispecies.
checkMultispecies <- function(presences, species.id){
  if(length(unique(presences[,species.id]))>1) multi <- TRUE
  else multi <- FALSE
  multi
}

checkWindow <- function(presences, window, coord, quiet){

  if(!is.null(window))
    expectClasses(window,"SpatRaster") #switching to terra (as it appears to be faster)

  if (is.null(window)) {
    if(!quiet) message("Window is NULL, a raster-based window will be generated based on the extent of 'presences'.\n")
    window <- defaultWindow(presences,coord)
  }
  window
}


## function to extract covariates for presence and background points.
getCovariates <- function(pbxy, covariates=NULL, interpolation, coord, bufferNA, bufferSize){
  if(is.null(covariates)){
    covars <- cbind(SiteID=pbxy[,"SiteID"],pbxy[,coord])
  } else {
    # expectClasses(covariates,c("RasterLayer","RasterStack","RasterBrick"),covariates)
    covars <- terra::extract(x = covariates,
                             y = data.frame(X=as.numeric(pbxy[,coord[1]]),
                                           Y=as.numeric(pbxy[,coord[2]])),
                             method=interpolation,
                             na.rm=TRUE)
    covars <- cbind(SiteID=pbxy[,"SiteID"],pbxy[,coord],covars)
  if(bufferNA){
    if(any(!complete.cases(covars))){
        message('NA cells generated during covariate extraction. Extracting values from nearest (1 step) neighbour -- might be prudent to check imputation (and why it was imputed).')
        missXY <- which(!complete.cases(covars))
        missCoord <- covars[missXY,coord]
        if(is.null(bufferSize)){
          if(terra::is.lonlat(covariates)) ltlnscale <- 100000
          else ltlnscale <- 1
          buff <- terra::global(terra::area(covariates),fun="mean")*ltlnscale
          # buff <- terra::cellS
        }else {
          buff <- bufferSize
        }
        buffCovars <- extract(x=covariates,y=missCoord,fun=mean,na.rm=TRUE,buffer=buff)
        covars[missXY,-1:-3] <- buffCovars
        }
      }
    }
  return(covars)
}



getSppNames <- function(presences, species.id){

  sppIdx <- list()
  sppIdx$sppNames <- unique(presences[,species.id])
  sppIdx$sppNames <- factor(sppIdx$sppNames, levels = unique(presences[,species.id]))
  sppIdx$sppNumber <- seq_along(sppIdx$sppNames)

  return(sppIdx)

}




defaultWindow <- function (presences, coord) {

  # get limits
  xlim <- range(presences[,coord[1]])
  ylim <- range(presences[,coord[2]])

  # buffer
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.1
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.1

  xlim[1] <- floor(xlim[1])
  xlim[2] <- ceiling(xlim[2])
  ylim[1] <- floor(ylim[1])
  ylim[2] <- ceiling(ylim[2])

  # reso <- round(diff(seq(ylim[1],ylim[2],length.out=50))[1],1)
  e <- extent(c(xlim,ylim))
  sa <- terra::rast(xmin=xlim[1],xmax=xlim[2],
             ymin=ylim[1],ymax=ylim[2],
             nrows=50,ncols=50,
             crs="+proj=longlat +datum=WGS84")
  # sa <- raster(e,res=reso, crs="+proj=longlat +datum=WGS84")
  terra::values(sa) <- 1#:terra::ncell(sa)
  return (sa)
}


fastWideMatrix <- function(dat, species.id){

  dat[,"OrigOrder"] <- factor(dat[,"OrigOrder"])
  # dat[,species.id] <- factor(dat[,"wts.dataset"])

  out <- matrix(nrow=nlevels(dat[,"OrigOrder"]),
                ncol=nlevels(dat[,"wts.dataset"]),
                dimnames=list(levels(dat[,"OrigOrder"]),levels(dat[,"wts.dataset"])))

  out[cbind(dat[,"OrigOrder"], dat[,"wts.dataset"])] <- dat[,"pres"]
  # out

  return(out)
}


fastWideMatrixWeights <- function(dat, sppNames){

  dat[,"OrigOrder"] <- factor(dat[,"OrigOrder"])
  dat[,"DatasetID"] <- factor(dat[,"DatasetID"])
  # [,"wts.dataset"]
  wtsdat <- with(dat, {
    out <- matrix(nrow=nlevels(OrigOrder),
                  ncol=nlevels(DatasetID),
                  dimnames=list(levels(OrigOrder),levels(DatasetID)));
    out[cbind(OrigOrder, DatasetID)] <- wts.area;
    out})

  colnames(wtsdat) <- sppNames$sppNames[which(colnames(wtsdat)%in%sppNames$sppNumber)]
  wtsdat
}

transposePPMdata <- function( dat, sppNames, coordNames, covarNames){

  dat1 <- list()
  dat1$wts <- dat$wtsmat
  ## look out fot any negative weights
  dat1$wts[!is.na(dat1$wts)& dat1$wts<=0] <- sqrt(.Machine$double.eps)
  my.ord <- match(sppNames$sppNames,colnames(dat1$wts))#gtools::mixedorder( colnames( dat1$wts))

  # responses
  dat1$y <- dat$mm[,colnames( dat$mm) %in% colnames( dat1$wts)]
  dat1$y <- dat1$y[,my.ord]
  dat1$y <- as.matrix( dat1$y)
  colnames(dat1$y) <- sppNames$sppNames[my.ord]

  #weights
  dat1$wts <- dat1$wts[,my.ord]
  colnames(dat1$wts) <- sppNames$sppNames[my.ord]

  # Model matrix add in coordinates in the model matrix
  if(!is.null(covarNames)){
    dat1$covars <- dat$mm[,covarNames]
    dat1$mm <- cbind(dat1$y,dat$mm[,coordNames],dat1$covars)
  } else {
    dat1$mm <- cbind(dat1$y,dat$mm[,coordNames])
  }

  # locations
  dat1$locations <- dat$mm[,coordNames] #passed to ppmdata as coord argument

  # expectation z
  dat1$z <- dat1$y / dat1$wts

  # colnames(dat1$mm) <- sppNames$sppNames[my.ord]

  dat1$bkg <- apply( dat1$y, 1, function(x) all( x==0))
  dat1$nspp <- ncol( dat1$wts)
  dat1$m <- nrow( dat1$wts)
  dat1$sppNames <- colnames( dat1$wts)
  dat1$nUniquePres <- sum( !dat1$bkg)
  dat1$nBkg <- sum( dat1$bkg)

  return( dat1)
}


plapply <- function (X, FUN, ..., .parallel = 1, .seed = NULL, .verbose = TRUE) {
  if (!(useCluster <- inherits(.parallel, "cluster"))) {
    stopifnot(length(.parallel) == 1L, is.vector(.parallel,
                                                 "numeric"), .parallel >= 1)
    .parallel <- as.vector(.parallel, mode = "integer")
    if (.Platform$OS.type == "windows" && .parallel > 1L) {
      useCluster <- TRUE
      .parallel <- parallel::makeCluster(.parallel)
      on.exit(parallel::stopCluster(.parallel))
    }
  }
  FUN <- match.fun(FUN)
  .FUN <- if (useCluster || is.primitive(FUN)) {
    FUN
  }
  else {
    verboseExpr <- if (isTRUE(.verbose)) {
      if (.parallel == 1L && interactive()) {
        env <- new.env(hash = FALSE, parent = environment(FUN))
        environment(FUN) <- env
        env$pb <- txtProgressBar(min = 0, max = length(X),
                                 initial = 0, style = 3)
        on.exit(close(env$pb), add = TRUE)
        quote(setTxtProgressBar(pb, pb$getVal() + 1L))
      }
      else {
        on.exit(cat("\n"), add = TRUE)
        quote(cat("."))
      }
    }
    else if (is.call(.verbose) || is.expression(.verbose)) {
      .verbose
    }
    else if (is.character(.verbose)) {
      on.exit(cat("\n"), add = TRUE)
      substitute(cat(.verbose))
    }
    do.call(add.on.exit, list(FUN, verboseExpr))
  }
  if (!is.null(.seed)) {
    if (useCluster) {
      parallel::clusterSetRNGStream(cl = .parallel, iseed = .seed)
    }
    else {
      if (!exists(".Random.seed", envir = .GlobalEnv,
                  inherits = FALSE)) {
        set.seed(NULL)
      }
      .orig.seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .orig.seed, envir = .GlobalEnv),
              add = TRUE)
      if (.parallel == 1L) {
        set.seed(seed = .seed)
      }
      else {
        stopifnot(requireNamespace("parallel", quietly = TRUE))
        set.seed(seed = .seed, kind = "L'Ecuyer-CMRG")
        parallel::mc.reset.stream()
      }
    }
  }
  if (useCluster) {
    parallel::parLapply(cl = .parallel, X = X, fun = .FUN,
                        ...)
  }
  else if (.parallel == 1L) {
    lapply(X = X, FUN = .FUN, ...)
  }
  else {
    parallel::mclapply(X = X, FUN = .FUN, ..., mc.preschedule = TRUE,
                       mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = .parallel)
  }
}

add.on.exit <- function (FUN, expr){
  FUN <- match.fun(FUN)
  if (is.null(expr <- substitute(expr))) {
    return(FUN)
  }
  if (is.primitive(FUN)) {
    stop("not implemented for primitive functions")
  }
  onexitexpr <- substitute(on.exit(expr))
  obody <- body(FUN)
  body(FUN) <- if (is.call(obody) && identical(as.name("{"),
                                               obody[[1L]])) {
    as.call(append(x = as.list(obody), values = onexitexpr,
                   after = 1L))
  }
  else {
    as.call(c(as.name("{"), onexitexpr, obody))
  }
  FUN
}

