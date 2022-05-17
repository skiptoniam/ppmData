#' @importFrom terra xyFromCell values ext freq extract ncell
quasiRandomQuad <- function(npoints,
                            window,
                            coord,
                            quasirandom.samples=NULL){

  #generate a set of potential sites for quasi-random generation
  rast_coord <- terra::xyFromCell(window,ncell(window))
  rast_data <- terra::values(window)
  na_sites <- which(!complete.cases(rast_data))
  potential_sites <- rast_coord

  if(is.null(quasirandom.samples)) quasirandom.samples <- 10*npoints

  exty <- terra::ext(window)
  study.area <- matrix( as.vector( exty)[c( 1,2,2,1, 3,3,4,4)], nrow=4, ncol=2)

  #this gives many many samples, unless the study region is an odd shape, which it shouldn't be (as it is just an extent at this stage)
  samp <- quasiSampFromHyperRect(nSampsToConsider=quasirandom.samples,
                                 randStartType=2,
                                 designParams=list(dimension=2,study.area=study.area))

  ## setup the inclusion probs.
  sampValues <- terra::extract(window,samp[,1:2])
  NAsamps <- which(!complete.cases(sampValues))  ### these will give the

  ## generate a set of buffer points using the samp
  if(length(NAsamps) > 0){
    totalCells <- terra::ncell(window)
    NAcount <- terra::freq(window,value=NA)[3]
    NonNAcount <- totalCells - NAcount
    NAratio <- NAcount/NonNAcount
    samplesNAcells <- samp[NAsamps,]
    randpointsDummy <- samplesNAcells[1:(round(npoints*NAratio)),1:2]
    colnames(randpointsDummy) <- coord
    randpointsDummy <- as.data.frame(randpointsDummy)
  } else {
    randpointsDummy <- NULL
  }

  ## generate the actuall quasi random points we want to generate for the model.
  Nsamps <- nrow(sampValues)
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

#' quasiSampFromHyperRect
#'@author Scott Foster
#'@param nSampsToConsider Number of samples to consider
#'@param randStartType Random start type, default is 2.
#'@param designParams A list of design parameters for the sampler.
#'@importFrom randtoolbox halton
#'@importFrom mgcv in.out
#'@description Copied from MBHdesign as it wasn't an exported function.

quasiSampFromHyperRect <- function (nSampsToConsider, randStartType = 2, designParams){
  samp <- randtoolbox::halton(nSampsToConsider * 2, dim = designParams$dimension +
                                1, init = TRUE)
  if (randStartType == 1)
    skips <- rep(sample(1:nSampsToConsider, size = 1, replace = TRUE),
                 designParams$dimension + 1)
  if (randStartType == 2)
    skips <- sample(1:nSampsToConsider, size = designParams$dimension +
                      1, replace = TRUE)
  samp <- do.call("cbind", lapply(1:(designParams$dimension +
                                       1), function(x) samp[skips[x] + 0:(nSampsToConsider -
                                                                            1), x]))
  myRange <- apply(designParams$study.area, -1, range)
  for (ii in 1:designParams$dimension) samp[, ii] <- myRange[1,ii] + (myRange[2, ii] - myRange[1, ii]) * samp[, ii]
  if (designParams$dimension == 2) {
    tmp <- mgcv::in.out(designParams$study.area, samp[,1:designParams$dimension])
    samp <- samp[tmp, ]
  }
  return(samp)
}


quasiRandomWeights <- function(presences, quadrature, quadDummy, window, coord, speciesIdx){#, returnDirtess){

  if(!is.null(quadDummy)) {
    allpts.id <- rbind(presences, quadrature, quadDummy)
    allpts <- allpts.id[,coord]
  } else {
    allpts.id <- rbind(presences, quadrature)
    allpts <- rbind(presences[,coord], quadrature[,coord])
  }

  ## Get a bbox which has a small buffer
  bbox <- convert2pts(allpts)
  bbox[c(1,3)] <- bbox[c(1,3)] - 1e-10
  bbox[c(2,4)] <- bbox[c(2,4)] + 1e-10

  ## Tracking of site and species ids
  allpts$id <- 1:nrow(allpts)
  allpts$dataset <- allpts.id[,speciesIdx]

  ## Do the Dirichlet with me
  dirareas <- getDirichlet(allpts, bbox, coord)#, clip, polyclip)#, returnDirtess)

  ## merge with all pts
  res <- merge(allpts, dirareas, by='id', all=TRUE, sort=TRUE)

  ## remove dummy points if present
  res <- res[!res$dataset=="dummy",]

  return( res)
}

getDirichlet <- function(allpts, bbox, coord){#}, clip=FALSE, polyclip = NULL){#}, return_dirtess = TRUE ){

  res <- data.frame( id=1:nrow( allpts), area=NA)

  tess <- dirTess(as.matrix(allpts[,coord]), bbox = bbox)

  # if(clip){
  #   tess <- polygonise(tess,clip=TRUE, polyclip = polyclip, proj = proj)
  # }

  area_wts <- getAreasDirichlet(tess, nodummy=TRUE)

  res$area <- area_wts

  return(res)
}

getAreasDirichlet <- function(object, nodummy=TRUE){

  if(nodummy){
    area_wts <- sapply(object$polygons$poly,function(x)x$area)[1:nrow(object$coords)]
  } else {
    area_wts <- sapply(object$polygons$poly,function(x)x$area)
  }
  return(area_wts)
}



