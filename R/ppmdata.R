#' @name ppmdata
#' @title Create a point process quadrature scheme for spatial modelling.
#' @description This package is a way to efficiently generate a quasirandom set
#' of background points for presence-only modelling of single or multiple
#' responses. The package was set up to model multiple species presence-only
#' datasets, but could be used for an point process spatial modelling.
#' Quasi-random points are a nice alternative to pseudo-random samples, this is
#' because we can generate a quasirandom sample across and areal region
#' (X and Y coordinates), but we can also extend the dimensions of the
#' quasirandom sample to a N-dimensional hypervolume, which will allow users to
#' effectively sample the spatial and environmental space. This in turn should
#' reduce autocorrelation in quadrature scheme. The weight of each quadrature
#' point is calculated using Dirichlet (Voronoi) Tessellation as provided from
#' the \link[deldir]{deldir} function.
#' @details The approach uses quasi-random sampling to generate a quadrature
#' scheme based (e.g Berman & Turner 1992; Warton & Shepard 2010;
#' Foster et al, 2017). The weights each quasi-random point in the quadrature
#' scheme is calculated using a dirichlet tessellation (Turner 2020). To improve
#' computational efficiency \link[deldir]{deldir} function for a large number of
#' of quadrature points, we set up an approach which breaks up the problem into
#' manageable sub-windows. We do this by keeping each deldir call to less that
#' 5000 points (which appears to be the point where the algorithm slows
#' noticeably). To avoid edge effect (large areas on the edges of sub-areas),
#' we rotate the subregions three times, the first two use a the nearest largest
#' prime number closest to total number of points (presences+quadrature points)
#' divided by 5000, which allows us to rotate the window on the x and y axis
#' with an subset of the sub-windows. We then calculate a third set of sub-
#' windows using a even set of squares. We then take the median weight across
#' all weight calculated for each point. We then can calculate this in parallel
#' for each species to make it computationally more efficient.
#' @export
#' @param npoints The number of quadrature points to generate.
#' @param presences a matrix, dataframe or SpatialPoints object giving the
#' coordinates of each species' presence in (should be a matrix of nsites * 3)
#' with the three columns being c("X","Y","SpeciesID"), where X is longitude,
#' Y is latitude and SpeciesID is factor which assoicates each occurence to a
#' species. If presences parameter is NULL then ppmDat will return the
#' quadrature scheme without presences.
#' @param window a raster object giving the area over which to generate
#' background points. NA cells are ignored and masked out of returned data.
#' If ignored, a rectangle defining the extent of \code{presences} will be used.
#' @param covariates A raster object containing covariates for modelling the
#' point process (best use a Raster stack or Raster brick).
#' @param coord is the name of site coordinates. The default is c('X','Y').
#' @param speciesIdx is the name of the species ID in the presences dataset.
#' The default is "SpeciesID".
#' @param mc.cores The number of cores to use in the processing. default is
#' parallel::detectCores()-1
#' @param quasirandom.samples This set the total number of samples to consider
#' in the BAS step (rejection sampling). The default is 10000, which means that
#' 10000 halton numbers are drawn and then thinned according to the inclusion
#' probabilities. You will need to increase the number of samples if selecting
#' a large number of quadrature points. The more quasirandomSample selected the
#' slower the quasirandom quadrature scheme will be to generate.
#' @param interpolation The interpolation method to use when extracting covariate data. Default is "bilinear", can also use "simple".
#' @examples
#' \dontrun{
#' library(ppmdata)
#' library(raster)
#' path <- system.file("extdata", package = "ppmdata")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' preds <- stack(lst)
#' window <- preds[[1]]
#' presences <- snails
#' bkgrid <- ppmdata(npoints = 1000, presences=presences, window = window, covariates = preds)
#' }

ppmdata <- function(npoints = 10000,
                    presences = NULL,
                    window = NULL,
                    covariates = NULL,
                    coord = c('X','Y'),
                    speciesIdx = "SpeciesID",
                    mc.cores = parallel::detectCores()-1,
                    quasirandom.samples = NULL,
                    interpolation = "bilinear"){

  ####  If not window is provided provide a dummy window
  window <- checkWindow(presences,window,coord)

  if(is.null(presences)){
   message('Generating background points in the absence of species presences')
   bckpts <- quasirandomMethod(npoints = npoints,
                               window = window,
                               covariates =  covariates,
                               coord =  coord,
                               quasirandom.samples = quasirandom.samples)

   sitecovariates <- getCovariates(pbxy = bckpts,
                                   covariates = covariates,
                                   interpolation = interpolation,
                                   coord = coord)

  } else {

  # This should check the presences and make it returns the data in the correct format for the remaining function.
  pressies <- coordMatchDim(known.sites = presences,
                            coord = coord,
                            speciesIdx = speciesIdx)

  # Hold onto the species names from the speciesIdx column
  sppNames <- getSppNames(pressies, speciesIdx)

  # create some quasirandom background points
  bckpts <- quasirandomMethod(npoints = npoints,
                              window = window,
                              covariates = covariates,
                              coord = coord,
                              quasirandom.samples = quasirandom.samples)

  # Sometimes the points are very close together, so let's jitter if needed.
  reswindow <- raster::res(window)[1]
  tmpPts <- jitterIfNeeded( pressies=pressies, bckpts=bckpts$quasiPoints, coord=coord, aBit=reswindow/2)
  pressies <- tmpPts$pressies
  bckptsQ <- tmpPts$bckpts
  bckptsD <- bckpts$quasiDummy

  # Check to see if the presences are for a single species or multiple.
  ismulti <- checkMultispecies(pressies, speciesIdx)

  if(ismulti){
      message("Developing a quadrature scheme for multiple species (marked) dataset.")
      wts <- getMultispeciesWeights(presences = pressies,
                                    quadrature = bckpts$quasiPoints,
                                    quadDummy = bckpts$quasiDummy,
                                    window = window, coord = coord,
                                    mc.cores = mc.cores)

      sitecovariates <- getCovariates(pbxy = wts,
                                      covariates,
                                      interpolation=interpolation,
                                      coord=coord)

    } else {
      message("Developing a quadrature scheme for a single species dataset.")
      wts <- getSinglespeciesWeights(presences = pressies, quadrature = bckpts$quasiPoints,
                                     quadDummy = bckpts$quasiDummy, window = window,
                                     coord = coord, mc.cores = mc.cores)
      sitecovariates <- getCovariates(pbxy = wts, covariates = covariates, interpolation=interpolation, coord=coord)
    }


  }

  dat <- assembleQuadData(pressies, bckpts$quasiPoints, sitecovariates, wts, coord)

  if(!is.null(covariates)) covarNames <- names(covariates)
  coordNames <- coord
  if(ismulti){
     datOut <- transpose_ppmdata(dat, sppNames, coordNames, covarNames)
     class(datOut) <- c("ppmdata","multiple.species")
  } else {
    datOut <- list(dat)
    class(datOut) <- c("ppmdata","single.species")
  }

  return(datOut)
}

jitterIfNeeded <- function( pressies, bckpts, coord, aBit=1e-4){
  #the pressie bit first
  #are there any duplicates within a species?  If so, then jitter the duplicates
  for( jj in as.character( unique( pressies$SpeciesID))){ #I think that we have made the assumption that this is called SpeciesID...?
    sppJ <- which( pressies$SpeciesID==jj)
    dupes <- which( duplicated( pressies[sppJ,coord]))  #shouldn't need to round this as deldir reportedly uses duplicated
    if( length( dupes)>0){
      pressies[sppJ[dupes],coord[1]] <- jitter( pressies[sppJ[dupes],coord[1]], amount=aBit)
      pressies[sppJ[dupes],coord[2]] <- jitter( pressies[sppJ[dupes],coord[2]], amount=aBit)
    }
  }

  #background now
  #are there any points that duplicate a presence point?  If so, jitter.
  npres <- nrow( pressies)
  dupes <- which( duplicated( rbind( pressies[,coord], bckpts)))
  if( length( dupes)>0)
    dupes <- dupes[dupes>npres]
  if( length( dupes)>0){
    dupes <- dupes - npres
    bckpts[dupes,coord[1]] <- jitter( bckpts[dupes,coord[1]], amount=aBit)
    bckpts[dupes,coord[2]] <- jitter( bckpts[dupes,coord[2]], amount=aBit)
  }

  return( list( pressies=pressies, bckpts=bckpts))
}

assembleQuadData <- function(presences, quadrature, sitecovariates, wts, coord){

  ismulti <- checkMultispecies(presences, speciesIdx)

  if(!ismulti) type <- "long"
  else type <- "wide"

  final_dat <- switch(type,
                      long=longdat(wts, sitecovariates, coord),
                      wide=widedat(presences, quadrature,
                                   sitecovariates,
                                   wts, coord))

  return(final_dat)

}


longdat <- function(wts, sitecovariates=NULL, coord){

  if(!is.null(sitecovariates)) dat2 <- data.frame(wts[,coord],sitecovariates[,-1:-3],presence=wts$pres, weights=wts$wts)
  else dat2 <- data.frame(wts[,coord],presence=wts$pres,weights=wts$wts)

  if(length(nrow(dat2[!complete.cases(dat2), ]))>0) message("Your covariate data has ", nrow(dat2[!complete.cases(dat2), ])," rows with NAs in them - check before modelling.")

  return(dat2)
}

widedat <- function(presence, quadrature, sitecovariates, wts, coord){

  # Assemble a data.frame with all the bits we want.
  pamat <- fastwidemat(wts)
  presences_pamat <- pamat[-which(pamat[,"quad"]==0),-which(colnames(pamat)=='quad')]
  # presences_pamat[presences_pamat==0]<-NA
  quad_pamat <- pamat[which(pamat[,"quad"]==0),-which(colnames(pamat)=='quad')]
  quad_pamat[is.na(quad_pamat)]<-0
  response_ppmmat <- as.data.frame(rbind(presences_pamat,quad_pamat))
  response_ppmmat$Const <- 1
  response_ppmmat$OrigOrder <- as.integer(rownames(response_ppmmat))
  sitecovariates$OrigOrder <- wts$OrigOrder
  df <- merge(response_ppmmat,sitecovariates[!duplicated(sitecovariates$OrigOrder),], by = "OrigOrder", sort=FALSE)
  wtsmat <- fastwidematwts(wts)
  ids <- wts[!duplicated(wts[,c('SpeciesID','DatasetID')]),c('SpeciesID','DatasetID')]
  ids <- ids[-which(ids$SpeciesID=='quad'),]
  idx_rows <- df$OrigOrder
  idx_cols <- match(colnames(quad_pamat),ids$SpeciesID)
  wtsmat <- wtsmat[idx_rows,idx_cols]
  colnames(wtsmat) <- colnames(quad_pamat)
  return(list(mm=df,wtsmat=wtsmat))
}

quasirandomMethod <- function(npoints, window, covariates=NULL, coord, quasirandom.samples=NULL){

  #generate a set of potential sites for quasi-random generation
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

  if(is.null(quasirandom.samples)) quasirandom.samples <- 10*npoints

  exty <- raster::extent( window)
  study.area <- matrix( as.vector( exty)[c( 1,2,2,1, 3,3,4,4)], nrow=4, ncol=2)

  #this gives many many samples, unless the study region is an odd shape, which it shouldn't be (as it is just an extent at this stage)
  samp <- MBHdesign:::quasiSamp_fromhyperRect(nSampsToConsider=quasirandom.samples,
                                              randStartType=2,
                                              designParams=list(dimension=2,study.area=study.area))

  ## setup the inclusion probs.
  sampValues <- extract(window,samp[,1:2])
  NAsamps <- which(!complete.cases(sampValues))  ### these will give the

  ## generate a set of buffer points using the samp
  if(length(NAsamps) > 0){
    totalCells <- raster::ncell(window)
    NAcount <- raster::freq(window,value=NA)
    NonNAcount <- totalCells - NAcount
    ratio <- NAcount/NonNAcount
    samplesNAcells <- samp[NAsamps,]
    randpointsDummy <- samplesNAcells[1:(round(npoints*ratio)),1:2]
    colnames(randpointsDummy) <- coord
    randpointsDummy <- as.data.frame(randpointsDummy)
  } else {
    randpointsDummy <- NULL
  }

  ## generate the actuall quasi random points we want to generate for the model.
  Nsamps <- length(sampValues)
  inclusion_probs <- rep(1/Nsamps, Nsamps)
  inclusion_probs1 <- inclusion_probs/max(inclusion_probs)
  inclusion_probs1[NAsamps] <- 0

  #Spatially thin the sample which are NA sites in the raster window
  samp <- samp[samp[,3]<inclusion_probs1,1:2]
  if(nrow( samp) < npoints)
    stop("No set of background points found for this region.  Please increase the number of possible samples.")
  samp <- samp[1:npoints,]


  randpoints <- as.data.frame(samp[,1:2])
  colnames(randpoints ) <- coord


  return(list(quasiPoints = randpoints, quasiDummy = randpointsDummy))
}

coordMatchDim <- function (known.sites,coord,speciesIdx){

  # check object classes
  expectClasses(known.sites,
                c('matrix',
                  'data.frame'),
                name = 'known.sites')


  if(!any(colnames(known.sites)%in%speciesIdx))
    stop("'speciesIdx': ",speciesIdx," ,does not match any of the column names in your presences data.\n")
  if(!any(colnames(known.sites)%in%coord))
    stop("The coordinates names: ",paste(coord,collapse = ", ")," do not match any of the column names in your presences data.\n")

  df <- known.sites[,c(coord,speciesIdx)]

  return (df)
}


## function to extract covariates for presence and background points.
getCovariates <- function(pbxy, covariates=NULL, interpolation, coord){
  if(is.null(covariates))return(NULL)
  covars <- raster::extract(x = covariates,
                            y = data.frame(X=as.numeric(pbxy[,coord[1]]),
                                           Y=as.numeric(pbxy[,coord[2]])),
                            method=interpolation,
                            na.rm=TRUE)
  covars <- cbind(SiteID=pbxy[,"SiteID"],pbxy[,coord],covars)
  return(covars)
}



getSppNames <- function(presences, speciesIdx){

  sppIdx <- list()
  sppIdx$sppNames <- unique(presences[,speciesIdx])
  if(is.factor(sppIdx$sppNames)) sppIdx$sppNames <- droplevels(sppIdx$sppNames)
  sppIdx$sppNumber <- seq_along(sppIdx$sppNames)

  return(sppIdx)

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
# checkDuplicates <- function(presences, coord, speciesIdx){
#   if(is.null(presences))return(NULL)
#   dups <- duplicated(presences)
#   if(sum(dups)>0){ message("There were ",sum(dups)," duplicated points unique to X, Y & SpeciesID, they have been removed.")
#   dat <- presences[!dups,]
#   } else {
#   dat <- presences
#   }
#   dat <- as.data.frame(dat)
#   colnames(dat) <- c(coord, speciesIdx)
#   dat
# }

## check to see if the presences dataset is multispecies.
checkMultispecies <- function(presences, speciesIdx){
  if(length(unique(presences[,speciesIdx]))>1) multi <- TRUE
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
  ylim <- range(presences[,coord[2]])

  # buffer
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


fastwidemat <- function(dat, speciesIdx){

  dat[,"OrigOrder"] <- factor(dat[,"OrigOrder"])
  dat[,speciesIdx] <- factor(dat[,speciesIdx])

  result <- with(dat, {
    out <- matrix(nrow=nlevels(OrigOrder), ncol=nlevels(speciesIdx),
                  dimnames=list(levels(OrigOrder), levels(speciesIdx)))
    out[cbind(OrigOrder, SpeciesID)] <- pres
    out
  })

  return(result)
}


fastwidematwts <- function(dat){

  dat[,"OrigOrder"] <- factor(dat[,"OrigOrder"])
  dat[,"DatasetID"] <- factor(dat[,"DatasetID"])

  wtsdat <- with(dat, {
    out <- matrix(nrow=nlevels(OrigOrder), ncol=nlevels(DatasetID),
                  dimnames=list(levels(OrigOrder), levels(DatasetID)))
    out[cbind(OrigOrder, DatasetID)] <- wts.area
    out
  })

  wtsdat
}

transpose_ppmdata <- function( dat, sppNames, coordNames, covarNames){

  dat1 <- list()
  dat1$wts <- dat$wtsmat
  my.ord <- match(sppNames$sppNames,colnames(dat1$wts))#gtools::mixedorder( colnames( dat1$wts))
  idx <- complete.cases(dat$mm[,covarNames])
  dat1$y <- dat$mm[idx,colnames( dat$mm) %in% colnames( dat1$wts)]
  dat1$y <- dat1$y[,my.ord]
  dat1$y <- as.matrix( dat1$y)
  dat1$covars <- dat$mm[idx,covarNames]
  dat1$locations <- dat$mm[idx,coordNames] #passed to ppmdata as coord argument
  dat1$wts <- dat1$wts[idx,my.ord]
  dat1$z <- dat1$y / dat1$wts
  dat1$mm <- cbind(dat1$y,dat1$covars)
  # colnames(dat1$y) <- sppNames$sppNames[my.ord]
  # colnames(dat1$wts) <- sppNames$sppNames[my.ord]

  dat1$bkg <- apply( dat1$y, 1, function(x) all( x==0))
  dat1$nspp <- ncol( dat1$wts)
  dat1$m <- nrow( dat1$wts)
  dat1$sppNames <- colnames( dat1$wts)
  dat1$nUniquePres <- sum( !dat1$bkg)
  dat1$nBkg <- sum( dat1$bkg)

  return( dat1)
}


#'@rdname print.ppmdata
#'@name print.ppmdata
#'@title Print a summary of ppmdata object.
#'@param x A model object.
#'@param \\dots Ignored
#'@export

print.ppmdata <- function (x, ...){

    if(class(x)[2]%in%"single.species"){

      message("There are ", sum(x[[1]]$presence), " presence observations for this species")
      message("There are ", sum(x[[1]]$presence==0), " background quadrature (integration) points")
      message("There are a total of ", nrow(x[[1]]), " sites in the model.matrix")
      no_nans <- nrow(x[[1]][!complete.cases(x[[1]]),])
      if(no_nans>0)message("There are a total of ",no_nans, " NaNs in the covariates, check before modelling.")

    } else {

      message("There are ", x$nUniquePres, " presence observations for ", x$nspp," species")
      message("There are ", x$nBkg, " quadrature (integration) points for each of the ", x$nspp ," species")
      message("There are a total of ",x$m, " sites in the model.matrix")
      no_nans <- nrow(x$covars[!complete.cases(x$covars),])
      if(no_nans>0)message("There are a total of ",no_nans, " NaNs in the covariates, check before modelling.")

      }
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

