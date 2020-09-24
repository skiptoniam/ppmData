
getWeights <- function( presences, backgroundpoints, window, coord, mc.cores){
  
  allpts <- rbind(presences[,coord], backgroundpoints[,coord])
  # window_ext <- convert2pts(window)
  window_ext <- convert2pts(allpts)
  window_ext[c(1,3)] <- window_ext[c(1,3)] - 1e-10
  window_ext[c(2,4)] <- window_ext[c(2,4)] + 1e-10
  allpts$id <- 1:nrow( allpts)
  
  ##generate some dummy points around allpts - doesn't work very well for odd shaped things.
  # ch <- allpts[chull(allpts[,coord]),coord]
  # dummyPts <- do.call(cbind,spatstat::spokes(ch[,1],ch[,2],nrad = 12))
  # dummyPts <- as.data.frame(dummyPts[chull(dummyPts[,1:2]),])

  ###Skip: this is hardwired to unit square
  #this is a bit of a dodge and exists due to < not being the same as <=.  Only a problem for the upper edge of of unit square
  # allpts[allpts[,1]==1,1] <- 1-1e-10
  # allpts[allpts[,2]==1,2] <- 1-1e-10
  ptsPerArea <- 5000
  primmy <- primest(ceiling( nrow( allpts)/ptsPerArea))
  prod.primmy <- outer( primmy, primmy)
  prod.primmy[upper.tri(prod.primmy)] <- Inf
  diag( prod.primmy) <- Inf
  
  tmp <- apply( prod.primmy, 1, function(xx) which.min( abs( xx - nrow( allpts))))
  tmp.id <- cbind( 1:nrow( prod.primmy), tmp)
  min.id <- which.min( abs( prod.primmy[tmp.id]-nrow( allpts)))
  ndivisions <- primmy[tmp.id[min.id,]]
  
  ###SKIP: this is hardwired to unit square
  #cuts on unit square in arrangement 1
  cuts1 <- list( seq( from=window_ext[1], to=window_ext[2], length=ndivisions[1]+1), seq( from=window_ext[3], to=window_ext[4], length=ndivisions[2]+1), ndivisions)
  all.boxes <- getBoxes( cuts1)
  
  ##############
  #voronoi areas for first rotation
  #set up a cluster for parallel
  cl <- parallel::makeCluster(mc.cores)
  # parallel::clusterExport( cl, "deldir", envir = .getNamespace("deldir"))
  parallel::clusterEvalQ(cl, library("deldir"))
  
  tmp <- parallel::parLapply( cl, seq_len(nrow(all.boxes)), areasWithinBoxes, allpts=allpts, boxes=all.boxes)
  areas1 <- do.call( "rbind", tmp)
  
  ###SKIP: this is hardwired to unit square
  #cuts on unit square in arrangement 2
  cuts2 <- list( seq( from=window_ext[1], to=window_ext[2], length=ndivisions[2]+1), seq( from=window_ext[3], to=window_ext[4], length=ndivisions[1]+1), ndivisions[2:1])
  all.boxes <- getBoxes( cuts2)
  tmp <- parallel::parLapply( cl, seq_len(nrow(all.boxes)), areasWithinBoxes, allpts=allpts, boxes=all.boxes)
  areas2 <- do.call( "rbind", tmp)
  
  ###SKIP: this is hardwired to unit square
  #cuts on unit square in arrangement 3 (squares)
  ndivisions <- rep( ceiling( sqrt(nrow( allpts)/ptsPerArea)), 2)
  cuts3 <- list( seq( from=window_ext[1], to=window_ext[2], length=ndivisions[1]+1), seq( from=window_ext[3], to=window_ext[4], length=ndivisions[2]+1), ndivisions)
  all.boxes <- getBoxes( cuts3)
  tmp <- parallel::parLapply( cl, seq_len(nrow(all.boxes)), areasWithinBoxes, allpts=allpts, boxes=all.boxes)
  areas3 <- do.call( "rbind", tmp)
  
  #tidying
  parallel::stopCluster(cl)
  
  #put them both together
  allAreas <- merge( areas1, areas2, all=T, by="id", sort=TRUE)
  allAreas <- merge( allAreas, areas3, all=TRUE, by='id', sort=TRUE)
  allAreas[allAreas==-99999] <- NA
  
  # remove really extreme areas - this is quite slow.
  # quant99 <- apply( allAreas[,-1], 1, quantile, .99, na.rm=TRUE)
  # res <- cbind( allAreas$id, apply(allAreas[,-1], 1, function(x) median(x[x>=quant99[parent.frame()$i[]]], na.rm=TRUE)))
  res <- cbind( allAreas$id, apply(allAreas[,-1], 1, median, na.rm=TRUE))
  colnames( res) <- c("id","area")
  
  #form the return object
  res <- merge( allpts, res, by='id', all=TRUE, sort=TRUE)
  
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
      res$area <- -99999
    } else {
      res$area[subsetty] <- deldir::deldir( x=allpts[subsetty,1], y=allpts[subsetty,2], rw=boxes[kounter,])$summary$dir.area
    }
  } else {
    res$area <- -9999
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

getSinglespeciesWeights <- function(presences, backgroundpoints, window, coord, mc.cores){

  backgroundpoints$SpeciesID <- "quad"
  wts <- getWeights(presences, backgroundpoints, window, coord, mc.cores)
  pbxy <- rbind(presences[,coord],backgroundpoints[,coord])
  pbxy$OrigOrder <- seq_len(nrow(pbxy))
  pbxy$DatasetID <- 1
  dat <- cbind(pbxy,
               pres=c(rep(1,nrow(presences)),
                     rep(0,nrow(backgroundpoints))),
               wts=wts)

  df <- getSiteID(dat,coord)
  return(df)
}


combineDF.fun <- function( ii, xxx, yyy, coords){
  ####  Assumes that the colnames of xxx are a subset of those from yyy
  ####  This is not a totally memory efficient implementation: data on all species is passed to all species...
  xxx <- xxx[[ii]]
  newdf <- as.data.frame( matrix( NA, nrow=nrow( xxx) + nrow( yyy), ncol=length( coords) + 2))
  colnames( newdf) <- c(coords,"SpeciesID","OrigOrder")
  newdf[1:nrow( xxx), ] <- xxx[,colnames( newdf)]
  newdf[nrow( xxx) + 1:nrow( yyy), ] <- yyy[,colnames( newdf)]

  return( newdf)
}

getMultispeciesWeights <- function(presences, backgroundpoints, window, coord, mc.cores){

  presences$OrigOrder <- seq_len(nrow(presences))
  nspp <- length(unique(presences[,"SpeciesID"]))
  spps <- unique(presences[,"SpeciesID"])
  backgroundpoints$SpeciesID <- "quad"
  backgroundpoints$OrigOrder <- seq_len(nrow(backgroundpoints))+max(presences$OrigOrder)
  sppdata <- lapply(seq_len(nspp), function(ii)presences[presences$SpeciesID==spps[ii],])

  # sppBckWtsList <- parallel::mclapply( seq_len(nspp), function(ii) {cat( ii, " "); getWeights( sppdata[[ii]], backgroundpoints, coord, mc.cores)})
  # sppBckDatList <- parallel::mclapply( seq_len(nspp), combineDF.fun, xxx=sppdata, yyy=backgroundpoints, coords=coord)
  # sppCounts <- parallel::mclapply(seq_len(nspp),function(ii)nrow(sppdata[[ii]]))
  # sppBckDatList <- parallel::mclapply(seq_len(nspp),function(ii){sppBckDatList[[ii]]$DatasetID <- ii;sppBckDatList[[ii]]})
  # sppWtsList <- parallel::mclapply(seq_len(nspp), function(ii)cbind(sppBckDatList[[ii]],
                                                                    # pres=c(rep(1,sppCounts[[ii]]),rep(0,nrow(backgroundpoints))),

  sppBckWtsList <- lapply( seq_len(nspp), function(ii) {cat( ii, " "); getWeights( sppdata[[ii]], backgroundpoints, window, coord, mc.cores)})
  sppBckDatList <- lapply( seq_len(nspp), combineDF.fun, xxx=sppdata, yyy=backgroundpoints, coords=coord)                                                                    # wts=sppBckWtsList[[ii]]))
  sppCounts <- lapply(seq_len(nspp),function(ii)nrow(sppdata[[ii]]))
  sppBckDatList <- lapply(seq_len(nspp),function(ii){sppBckDatList[[ii]]$DatasetID <- ii;sppBckDatList[[ii]]})
  sppWtsList <- lapply(seq_len(nspp), function(ii)cbind(sppBckDatList[[ii]],
                                                                    pres=c(rep(1,sppCounts[[ii]]),rep(0,nrow(backgroundpoints))),
                                                                    wts=sppBckWtsList[[ii]]))
 dat <- do.call(rbind,sppWtsList)
 df <- getSiteID(dat,coord)
 return(df)

}

getSiteID <- function(dat,coord){
  stidfn <- function(df, cols) {
    comb <- do.call(paste, c(as.list(df[cols]), sep = "."))
    df$SiteID <- match(comb, unique(comb))
    df
  }
  df <- stidfn(dat,c(coord,'pres'))
  return(df)
}

#estimate the total area of all cells in extent in km^2
estimateWindowArea <- function(window){
  if(raster::isLonLat(window)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- raster::area(window)
    area_study <- raster::mask(area_rast,window)
    area_of_region <- raster::cellStats(area_study,sum, na.rm=TRUE)*1e+6
  } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    area_of_region <- (raster::ncell(window[!is.na(window)])  * raster::xres(window) * raster::yres(window))
  }
  return(area_of_region)
}

convert2pts <- function( window) {
  #convert raster to extreme points (of bounding rectangle)
  #example only, *should* work for raster.  But what data type have we got?
  tmp <- sp::coordinates( window)
  tmp1 <- c( range( tmp[,1]), range( tmp[,2]))

  return( tmp1)
  #be careful with extents of rasters though, they tend to have slightly larger dimensions than they ought
}


