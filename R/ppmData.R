#' @name ppmData
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
#' scheme is calculated using a Dirichlet tessellation (Turner 2020). To improve
#' computational efficiency of the \link[deldir]{deldir} function for a large
#' number of quadrature points, we set up an approach which breaks up the problem into
#' manageable sub-windows. We do this by keeping each deldir call to less that
#' approximately 5000 points (which appears to be the point where the algorithm
#' slows noticeably). To avoid edge effect (large areas on the edges of sub-areas),
#' we rotate the sub-regions three times, the first two use a the nearest largest
#' prime number closest to total number of points (presences+quadrature points)
#' divided by 5000, which allows us to rotate the window on the x and y axis
#' with an subset of the sub-windows. We then calculate a third set of sub-
#' windows using a even set of squares. We then take the median weight across
#' all weight calculated for each point. We then can calculate this in parallel
#' for each species to make it computationally more efficient.
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
#' @param speciesIdx is the name of the species ID in the presences dataset.
#' The default is "SpeciesID".
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
#' library(raster)
#' path <- system.file("extdata", package = "ppmdata")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' preds <- stack(lst)
#' window <- preds[[1]]
#' presences <- snails
#' bkgrid <- ppmData(npoints = 1000, presences=presences, window = window, covariates = preds)
#' }

ppmData <- function(presences = NULL,
                    window = NULL,
                    covariates = NULL,
                    npoints = NULL,
                    coord = c('X','Y'),
                    speciesIdx = "SpeciesID",
                    mc.cores = 1,
                    quasirandom.samples = NULL,
                    interpolation = "bilinear",
                    bufferNA = FALSE,
                    bufferSize = NULL,
                    quiet=FALSE){

  # if npoints in NULL setup a default amount.
  if(is.null(npoints)){
    ## linear scaling of quad points compared to max n presences.
    ## count bump this up more is so desired.
    # xx <- seq(0,10000,100)
    # plot(xx,(10 * ceiling(2 * sqrt(xx)/10))^2)

    npmx <- max(table(presences[,speciesIdx]))
    nquad <- rep(pmax(32, 10 * ceiling(2 * sqrt(npmx)/10)),2)
    npoints <- prod(nquad)
  }

  ## Make sure the column ids are characters and check for missing/wrong named coord/speciesIdx vars.
  if(!is.character(coord)) coord <- as.character(coord)
  if(all(!coord%in%colnames(presences))) stop(paste0('coord: "',coord[1],'" & "',coord[2],'" do not match any of the colnames in the presences data.frame'))
  if(!is.character(speciesIdx)) speciesIdx <- as.character(speciesIdx)
  if(all(!speciesIdx%in%colnames(presences))) stop(paste0('speciesIdx: "',speciesIdx,'" does not match any of the colnames in the presences data.frame'))

  ##  If not window is provided provide a dummy window
  if(is.null(window)) default_window <- TRUE
  else default_window <- FALSE
  window <- checkWindow(presences,window,coord,quiet)

  if(is.null(presences)){
   if(!quiet)message('Generating background points in the absence of species presences')
   bckpts <- quasiRandomMethod(npoints = npoints,
                               window = window,
                               covariates =  covariates,
                               coord =  coord,
                               quasirandom.samples = quasirandom.samples)

   sitecovariates <- getCovariates(pbxy = bckpts,
                                   covariates = covariates,
                                   interpolation = interpolation,
                                   coord = coord,
                                   bufferNA = bufferNA,
                                   bufferSize = bufferSize)

  } else {

  # This should check the presences and make it returns the data in the correct format for the remaining function.
  pressies <- checkPresences(known.sites = presences,
                            coord = coord,
                            speciesIdx = speciesIdx)

  # Hold onto the species names from the speciesIdx column
  sppNames <- getSppNames(pressies, speciesIdx)

  # create some quasirandom background points
  bckpts <- quasiRandomMethod(npoints = npoints,
                              window = window,
                              covariates = covariates,
                              coord = coord,
                              quasirandom.samples = quasirandom.samples)

  # Sometimes the points are very close together, so let's jitter if needed.
  reswindow <- raster::res(window)[1]
  tmpPts <- jitterIfNeeded( pressies=pressies, bckpts=bckpts$quasiPoints, coord=coord, speciesIdx = speciesIdx, aBit=reswindow/2)
  pressies <- tmpPts$pressies
  bckptsQ <- tmpPts$bckpts
  bckptsD <- bckpts$quasiDummy

  # Check to see if the presences are for a single species or multiple.
  ismulti <- checkMultispecies(pressies, speciesIdx)

  if(ismulti){
    if(!quiet)message("Developing a quadrature scheme for multiple species (marked) dataset.")
      wts <- getMultispeciesWeights(presences = pressies,
                                    quadrature = bckptsQ,
                                    quadDummy = bckptsD,
                                    window = window,
                                    coord = coord,
                                    speciesIdx = speciesIdx,
                                    mc.cores = mc.cores,
                                    sppNames = sppNames)

      sitecovariates <- getCovariates(pbxy = wts,
                                      covariates,
                                      interpolation=interpolation,
                                      coord=coord,
                                      bufferNA = bufferNA,
                                      bufferSize = bufferSize)

    } else {
      if(!quiet)message("Developing a quadrature scheme for a single species dataset.")
      wts <- getSingleSpeciesWeights(presences = pressies,
                                     quadrature = bckptsQ,
                                     quadDummy = bckptsD,
                                     window = window,
                                     coord = coord,
                                     speciesIdx = speciesIdx)
      # extract the covariate data
      sitecovariates <- getCovariates(pbxy = wts,
                                      covariates = covariates,
                                      interpolation=interpolation,
                                      coord=coord,
                                      bufferNA = bufferNA,
                                      bufferSize = bufferSize)
    }


  }

  # Assemble data
  dat <- assembleQuadData(presences = pressies,
                          quadrature = bckptsQ,
                          sitecovariates = sitecovariates,
                          wts = wts,
                          coord = coord,
                          speciesIdx = speciesIdx,
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
                     speciesIdx = speciesIdx,
                     mc.cores = mc.cores,
                     quasirandom.samples = quasirandom.samples,
                     interpolation = interpolation,
                     dw = default_window)

  class(res) <- c("ppmData")

  if(!quiet)print(res)
  return(res)
}

jitterIfNeeded <- function( pressies, bckpts, coord, speciesIdx, aBit=1e-4){
  #the pressie bit first
  #are there any duplicates within a species?  If so, then jitter the duplicates
  for( jj in as.character( unique( pressies[,speciesIdx]))){ #I think that we have made the assumption that this is called SpeciesID...?
    sppJ <- which(  pressies[,speciesIdx]==jj)
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

assembleQuadData <- function(presences, quadrature, sitecovariates, wts, coord, speciesIdx, sppNames){

  ismulti <- checkMultispecies(presences, speciesIdx)

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
                                   speciesIdx = speciesIdx,
                                   sppNames = sppNames))

  return(final_dat)

}


longData <- function(wts, sitecovariates=NULL, coord){

  if(!is.null(sitecovariates)) dat2 <- data.frame(wts[,coord],sitecovariates[,-1:-3],presence=wts$pres, weights=wts$wts)
  else dat2 <- data.frame(wts[,coord],presence=wts$pres,weights=wts$wts)

  # if(length(nrow(dat2[!complete.cases(dat2), ]))>0) message("Your covariate data has ", nrow(dat2[!complete.cases(dat2), ])," rows with NAs in them - check before modelling.")

  return(dat2)
}

wideData <- function(presence, quadrature, sitecovariates, wts, coord, speciesIdx, sppNames){

  # Assemble a data.frame with all the bits we want.
  pamat <- fastWideMatrix(wts, speciesIdx)
  presences_pamat <- pamat[-which(pamat[,"quad"]==0),-which(colnames(pamat)%in%c('quad','dummy'))]
  quad_pamat <- pamat[which(pamat[,"quad"]==0),-which(colnames(pamat)%in%c('quad','dummy'))]
  quad_pamat[is.na(quad_pamat)]<-0
  response_ppmmat <- as.data.frame(rbind(presences_pamat,quad_pamat))
  response_ppmmat$Const <- 1
  response_ppmmat$OrigOrder <- as.integer(rownames(response_ppmmat))

  # if(!is.null(sitecovariates)){
    sitecovariates$OrigOrder <- wts$OrigOrder
    df <- merge(response_ppmmat,sitecovariates[!duplicated(sitecovariates$OrigOrder),], by = "OrigOrder", sort=TRUE)
  # } else {
    # df <- merge(response_ppmmat,wts[!duplicated(wts$OrigOrder),c("OrigOrder")], by="OrigOrder", sort=TRUE)

  # }
  wtsmat <- fastWideMatrixWeights(wts, sppNames)
  # ids <- wts[!duplicated(wts[,c(speciesIdx,'DatasetID')]),c(speciesIdx,'DatasetID')]
  # ids <- ids[-which( ids[,speciesIdx] == 'quad'),]
  # idx_rows <- df$OrigOrder
  # idx_cols <- match(colnames(quad_pamat), ids[,speciesIdx] )
  # wtsmat <- wtsmat[idx_rows,idx_cols]
  if(!all.equal(colnames(wtsmat),colnames(quad_pamat))) message('names could be getting mixed up with factors in R')
  return(list(mm=df,wtsmat=wtsmat))
}

quasiRandomMethod <- function(npoints, window, covariates=NULL, coord, quasirandom.samples=NULL){

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

checkPresences <- function (known.sites,coord,speciesIdx){


  # check object classes
  expectClasses(known.sites, c('matrix','data.frame'),name = 'known.sites')

  if(!any(colnames(known.sites)%in%speciesIdx))
    stop("'speciesIdx': ",speciesIdx," ,does not match any of the column names in your presences data.\n")
  if(!any(colnames(known.sites)%in%coord))
    stop("The coordinates names: ",paste(coord,collapse = ", ")," do not match any of the column names in your presences data.\n")

  # try and sort out factors.
  known.sites[[speciesIdx]] <- factor(known.sites[[speciesIdx]], levels = unique(known.sites[[speciesIdx]]))
  df <- known.sites[,c(coord,speciesIdx)]

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

checkWindow <- function(presences, window, coord, quiet){

  if(!is.null(window))
    expectClasses(window,"RasterLayer","window")

  if(is.null(presences)){
    presences <- data.frame(X=runif(100,0,100),Y=runif(100,0,100))
  }

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
    expectClasses(covariates,c("RasterLayer","RasterStack","RasterBrick"),covariates)
    covars <- raster::extract(x = covariates,
                              # y = pbxy[,"SiteID"])#,## let's try with cell numbers to see if they line up better
                              y = data.frame(X=as.numeric(pbxy[,coord[1]]),
                                           Y=as.numeric(pbxy[,coord[2]])),
                              method=interpolation,
                              na.rm=TRUE)
    covars <- cbind(SiteID=pbxy[,"SiteID"],pbxy[,coord],covars)
  if(bufferNA){
    if(any(!complete.cases(covars))){
        message('Buffering NA cells generated during covariate extraction.')
        missXY <- which(!complete.cases(covars))
        missCoord <- covars[missXY,coord]
        if(is.null(bufferSize)){
          if(raster::isLonLat(covariates)) ltlnscale <- 100000
          else ltlnscale <- 1
          buff <- raster::cellStats(raster::area(covariates),mean)*ltlnscale
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



getSppNames <- function(presences, speciesIdx){

  sppIdx <- list()
  sppIdx$sppNames <- unique(presences[,speciesIdx])
  sppIdx$sppNames <- factor(sppIdx$sppNames, levels = unique(presences[,speciesIdx]))
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

  reso <- round(diff(seq(ylim[1],ylim[2],length.out=20))[1],1)
  e <- extent(c(xlim,ylim))
  sa <- raster(e,res=reso, crs="+proj=longlat +datum=WGS84")
  sa[]<- 1:ncell(sa)
  return (sa)
}


fastWideMatrix <- function(dat, speciesIdx){

  dat[,"OrigOrder"] <- factor(dat[,"OrigOrder"])
  # dat[,speciesIdx] <- factor(dat[,"wts.dataset"])

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

