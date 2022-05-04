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

#'@author Scott Foster
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
  for (ii in 1:designParams$dimension) samp[, ii] <- myRange[1,
                                                             ii] + (myRange[2, ii] - myRange[1, ii]) * samp[, ii]
  if (designParams$dimension == 2) {
    tmp <- mgcv::in.out(designParams$study.area, samp[,
                                                      1:designParams$dimension])
    samp <- samp[tmp, ]
  }
  return(samp)
}



quasiRandomWeights <- function( presences, quadrature, quadDummy, window, coord, speciesIdx){

  if(!is.null(quadDummy)) {
    allpts.id <- rbind(presences, quadrature, quadDummy)
    allpts <- allpts.id[,coord]
  } else {
    allpts.id <- rbind(presences, quadrature)
    allpts <- rbind(presences[,coord], quadrature[,coord])
  }
  window_ext <- convert2pts(allpts)
  window_ext[c(1,3)] <- window_ext[c(1,3)] - 1e-10
  window_ext[c(2,4)] <- window_ext[c(2,4)] + 1e-10
  allpts$id <- 1:nrow( allpts)
  allpts$dataset <- allpts.id[,speciesIdx]

  ##generate some dummy points around allpts - doesn't work very well for odd shaped things.
  ptsPerArea <- 5000
  primmy <- primest(ceiling( nrow( allpts)/ptsPerArea))
  prod.primmy <- outer( primmy, primmy)
  prod.primmy[upper.tri(prod.primmy)] <- Inf
  diag( prod.primmy) <- Inf

  tmp <- apply( prod.primmy, 1, function(xx) which.min( abs( xx*ptsPerArea - nrow( allpts))))
  tmp.id <- cbind( 1:nrow( prod.primmy), tmp)
  min.id <- which.min( abs( prod.primmy[tmp.id]*ptsPerArea-nrow( allpts)))
  ndivisions <- primmy[tmp.id[min.id,]]

  #cuts on unit square in arrangement 1
  cuts1 <- list( seq( from=window_ext[1], to=window_ext[2], length=ndivisions[1]+1), seq( from=window_ext[3], to=window_ext[4], length=ndivisions[2]+1), ndivisions)
  all.boxes <- getBoxes( cuts1)

  #voronoi areas for first rotation, set up a cluster for parallel
  tmp <- lapply(seq_len(nrow(all.boxes)), areasWithinBoxes, allpts=allpts, boxes=all.boxes)
  areas1 <- do.call( "rbind", tmp)

  #cuts on unit square in arrangement 2
  cuts2 <- list( seq( from=window_ext[1], to=window_ext[2], length=ndivisions[2]+1), seq( from=window_ext[3], to=window_ext[4], length=ndivisions[1]+1), ndivisions[2:1])
  all.boxes2 <- getBoxes( cuts2)
  tmp <- lapply(seq_len(nrow(all.boxes2)), areasWithinBoxes, allpts=allpts, boxes=all.boxes2)
  areas2 <- do.call( "rbind", tmp)

  #cuts on unit square in arrangement 3 (squares)
  ndivisions <- rep( ceiling( sqrt(nrow( allpts)/ptsPerArea)), 2)
  cuts3 <- list( seq( from=window_ext[1], to=window_ext[2], length=ndivisions[1]+1), seq( from=window_ext[3], to=window_ext[4], length=ndivisions[2]+1), ndivisions)
  all.boxes3 <- getBoxes( cuts3)
  tmp <- lapply(seq_len(nrow(all.boxes3)), areasWithinBoxes, allpts=allpts, boxes=all.boxes3)
  areas3 <- do.call( "rbind", tmp)

  #put them both together
  allAreas <- merge( areas1, areas2, all=TRUE, by="id", sort=TRUE)
  allAreas <- merge( allAreas, areas3, all=TRUE, by='id', sort=TRUE)
  # allAreas[allAreas==-99999] <- NA

  # remove really extreme areas - this is quite slow.
  res <- cbind( allAreas$id, apply(allAreas[,-1], 1, median, na.rm=TRUE))
  colnames( res) <- c("id","area")

  #form the return object
  res <- merge(allpts, res, by='id', all=TRUE, sort=TRUE)

  ## now drop the dummy sites
  res <- res[!res$dataset=="dummy",]

  return( res)
}

areasWithinBoxes <- function(kounter, allpts, boxes){
  res <- data.frame( id=1:nrow( allpts), area=NA)
  # for(kounter in 1:nrow(boxes)){
  subsetty <- (1:nrow( allpts))[allpts[,1]>=boxes[kounter,1] & allpts[,1]<boxes[kounter,2] & allpts[,2]>=boxes[kounter,3] & allpts[,2]<boxes[kounter,4]]
  # cat(length(subsetty)," ",kounter,"\n")
  # }
  #
  if( length( subsetty)>1){
    coords <- allpts[subsetty,-3]
    coordsDups <- coords[!duplicated(coords),]
    if(nrow(coordsDups)==1){
      res$area <- NA
    } else {
      res$area[subsetty] <- deldir::deldir( x=allpts[subsetty,1], y=allpts[subsetty,2], rw=boxes[kounter,],round=FALSE,digits=12, suppressMsge=TRUE)$summary$dir.area
    }
  } else {
    res$area <- NA
  }
  return( res[subsetty,])
}

getBoxes <- function(cutsPts){
  ids <- expand.grid( seq_len(cutsPts[[3]][1]), seq_len(cutsPts[[3]][2]))
  boxes <- matrix( NA, nrow=nrow( ids), ncol=4)
  for( kount in 1:nrow( ids)) boxes[kount,] <- c( cutsPts[[1]][ids[kount,1]+0:1], cutsPts[[2]][ids[kount,2]+0:1])
  return( boxes)
}

primest <- function(n){
  ### Stolen from https://stackoverflow.com/questions/3789968/generate-a-list-of-primes-up-to-a-certain-number
  p <- 2:n
  i <- 1
  while (p[i] <= sqrt(n)) {
    p <-  p[p %% p[i] != 0 | p==p[i]]
    i <- i+1
  }
  p
}
