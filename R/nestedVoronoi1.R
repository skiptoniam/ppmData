
####ALL STUFF WITH FOUR "#" is just used for testing.
####The function getWeights() should be a direct replacement for the old getWeights() BUT there are some small differences in arguments
####This is not *beautiful* code, but it does seem to work, but only for the unit square.
####As discussed, can you please look at expanding it to arbitary shapes for boundaries.

####It does appear to be acceptably quick.
####It could be made quicker by, for example, by re-using a previously created compute cluster

#####already in qrbp
####convert2pts <- function( window){
####  #convert raster to extreme points (of bounding rectangle)
####  #example only, *should* work for raster.  But what data type have we got?
####  tmp <- sp::coordinates( window)
####  tmp1 <- c( range( tmp[,1]), range( tmp[,2]))
####
####  return( tmp1)
####  #be careful with extents of rasters though, they tend to have slightly larger dimensions than they ought
####}
####
####library( qrbp)
#####create some data
####bkpts <- ppmData( npoints = 10000,
####                  presences = sppDat[,c("x","y","spp")],
####                  window = Xbrick[[1]],
####                  covariates = Xbrick,
####                  resolution = 0.00001,
####                  method = "quasirandom",
####                  interpolation = "bilinear",
####                  coord = c("X","Y"),
####                  control = ppmData.control( cores=7))
####
#####call the function to do the bizzo

####nbkg <- 50000
####backy <- data.frame( x=runif( nbkg), y=runif( nbkg))
####backy[backy[,1]<0.75 & backy[,2] < 0.75,] <- NA
####backy <- backy[!is.na( backy[,1]),]
####presences <- sppDat[sppDat$spp==1,c("x","y")]
####presences[presences[,1]<0.75 & presences[,2] < 0.75,] <- NA
####presences <- presences[!is.na( presences[,1]),]
####
####par( pty='s'); plot( backy, pch='.', asp=1)
####points( presences, pch=20, cex=0.6, col='red')
####
####getWeights( presences[,c('x','y')], backy, coord = c("x","y"), window = Xbrick[[1]])


#the bizzo
getWeights <- function( presences, backgroundpoints, coord, window, epsilon=sqrt(.Machine$double.eps), mc.cores=parallel::detectCores()-1){
  window_ext <- convert2pts( window)
  allpts <- rbind(presences[,coord], backgroundpoints[,coord])
  allpts$id <- 1:nrow( allpts)
  ###Skip: this is hardwired to unit square
  #this is a bit of a dodge and exists due to < not being the same as <=.  Only a problem for the upper edge of of unit square
  allpts[allpts[,1]==1,1] <- 1-1e-10
  allpts[allpts[,2]==1,2] <- 1-1e-10
  ptsPerArea <- 5000
  #this doesn't scale well with very large numbers
  #but use it if there aren't too many points
  #actually don't use it at all...  It won't really play well with irregular shapes
#  if( nrow( presences) + nrow( backgroundpoints) < ptsPerArea){
#    tmptmp <- deldir::deldir( x=allpts[,1], y=allpts[,2], plot=FALSE, rw=window_ext)
#    return( tmptmp$summary$dir.area)#there may be some round-off error...
#  }
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
  cuts1 <- list( seq( from=0, to=1, length=ndivisions[1]+1), seq( from=0, to=1, length=ndivisions[2]+1), ndivisions)
  all.boxes <- getBoxes( cuts1)

  ##############
  #voronoi areas for first rotation

  #set up a cluster for parallel
  cl <- parallel::makeCluster(mc.cores)
#  parallel::clusterExport(cl, c("allpts", "window_ext", "boxes"), envir = environment())
  parallel::clusterExport( cl, "deldir", envir = as.environment( "package:deldir"))

  tmp <- parallel::parLapply( cl, seq_len(nrow(all.boxes)), areasWithinBoxes, allpts=allpts, boxes=all.boxes)
  areas1 <- do.call( "rbind", tmp)

  ###SKIP: this is hardwired to unit square
  #cuts on unit square in arrangement 2
  cuts2 <- list( seq( from=0, to=1, length=ndivisions[2]+1), seq( from=0, to=1, length=ndivisions[1]+1), ndivisions[2:1])
  all.boxes <- getBoxes( cuts2)
  tmp <- parallel::parLapply( cl, seq_len(nrow(all.boxes)), areasWithinBoxes, allpts=allpts, boxes=all.boxes)
  areas2 <- do.call( "rbind", tmp)

  ###SKIP: this is hardwired to unit square
  #cuts on unit square in arrangement 3 (squares)
  ndivisions <- rep( ceiling( sqrt(nrow( allpts)/ptsPerArea)), 2)
  cuts3 <- list( seq( from=0, to=1, length=ndivisions[1]+1), seq( from=0, to=1, length=ndivisions[2]+1), ndivisions)
  all.boxes <- getBoxes( cuts3)
  tmp <- parallel::parLapply( cl, seq_len(nrow(all.boxes)), areasWithinBoxes, allpts=allpts, boxes=all.boxes)
  areas3 <- do.call( "rbind", tmp)

  #tidying
  parallel::stopCluster(cl)

  #put them both together
  allAreas <- merge( areas1, areas2, all=T, by="id", sort=TRUE)
  allAreas <- merge( allAreas, areas3, all=TRUE, by='id', sort=TRUE)
  res <- cbind( allAreas$id, apply( allAreas[,-1], 1, median, na.rm=TRUE))
  colnames( res) <- c("id","area")

  #form the return object
  res <- merge( allpts, res, by='id', all=TRUE, sort=TRUE)

  return( res)

}

areasWithinBoxes <- function(kounter, allpts, boxes){
  res <- data.frame( id=1:nrow( allpts), area=NA)
  subsetty <- (1:nrow( allpts))[allpts[,1]>=boxes[kounter,1] & allpts[,1]<boxes[kounter,2] & allpts[,2]>=boxes[kounter,3] & allpts[,2]<boxes[kounter,4]]
  if( length( subsetty)>0)
    res$area[subsetty] <- deldir::deldir( x=allpts[subsetty,1], y=allpts[subsetty,2], rw=boxes[kounter,], eps=1e-5)$summary$dir.area
  else
    res$area <- -9999
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


