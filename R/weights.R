# this is code out of ppm lasso package.
getTileWeights <- function (presences, backgroundpoints, coord = c("X", "Y")){

  sp.col = c(which(names(presences) == coord[1]), which(names(presences) == coord[2]))
  back.col = c(which(names(backgroundpoints) == coord[1]), which(names(backgroundpoints) == coord[2]))
  X.inc = sort(unique(backgroundpoints[, back.col[1]]))[2] - sort(unique(backgroundpoints[, back.col[1]]))[1]
  Y.inc = sort(unique(backgroundpoints[, back.col[2]]))[2] - sort(unique(backgroundpoints[, back.col[2]]))[1]
  back.0X = min(backgroundpoints[, back.col[1]]) - floor(min(backgroundpoints[, back.col[1]])/X.inc) * X.inc
  back.0Y = min(backgroundpoints[, back.col[2]]) - floor(min(backgroundpoints[, back.col[2]])/Y.inc) * Y.inc
  X = c(presences[, back.col[1]], backgroundpoints[, back.col[1]])
  Y = c(presences[, back.col[2]], backgroundpoints[, back.col[2]])
  round.X = round((X - back.0X)/X.inc) * X.inc
  round.Y = round((Y - back.0Y)/Y.inc) * Y.inc
  round.id = paste(round.X, round.Y)
  round.table = table(round.id)
  wts = X.inc * Y.inc/as.numeric(round.table[match(round.id, names(round.table))])
  wts
}

getWeights <- function( presences, backgroundpoints, window, epsilon=sqrt(.Machine$double.eps), method='voronoi'){
  window_ext <- convert2pts( window)  #Don't really know what window is, nor what it contains
  if( method=="voronoi"){
    #this doesn't scale well with very large numbers
    allpts <- rbind(presences, backgroundpoints)
    tmptmp <- deldir::deldir( x=allpts[,1], y=allpts[,2], plot=FALSE, rw=window_ext)
    return( tmptmp$summary$dir.area)#there may be some round-off error...
  }
  # if(method == "Scott2"){
  # ###  #using Scott's made-up method based on expected
  # ###  #based on David's made-up method based on expected (but assumes lots of points)
  #  total_area <- qrbp:::estimateWindowArea(window)
  #  npres <- nrow( presences)
  #  nbkg <- nrow( backgroundpoints)
  #  bkwts.single <- total_area/nbkg#assumes that all bkg points occupy the same space
  #  tmptmp <- deldir::deldir( x=presences[,1], y=presences[,2], plot=FALSE, rw=window_ext)
  #  pres.areas_nobkg <- tmptmp$summary$del.area
  #  pres.props_nobkg <- pres.areas_nobkg / sum( pres.areas_nobkg)
  # }
}

# getRandomWeights <- function(presences, backgroundpoints, window, epsilon=sqrt(.Machine$double.eps)){
#
#   # xy <- rbind(presences[,coord],backgroundpoints[,coord])
#   preswts <- rep(epsilon,nrow(presences))
#   total_area <- estimateWindowArea(window)
#   nbkg <- nrow(backgroundpoints)
#   bkwts <- rep(total_area/nbkg,nbkg)
#   wts <- c(preswts,bkwts)
#
#   return(wts)
# }

getSinglespeciesWeights <- function(presences, backgroundpoints, coord, method, window, epsilon){

  backgroundpoints$SpeciesID <- "quad"
  if(method%in%"grid")  wts <- getTileWeights(presences,backgroundpoints,coord)
  if(method%in%c("quasirandom","psuedorandom")) wts <- qrbp:::getWeights(presences[,coord], backgroundpoints, window, epsilon)
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

getMultispeciesWeights <- function(presences, backgroundpoints, coord, method, window, epsilon){

  presences$OrigOrder <- seq_len(nrow(presences))
  nspp <- length(unique(presences[,"SpeciesID"]))
  spps <- unique(presences[,"SpeciesID"])
  backgroundpoints$SpeciesID <- "quad"
  backgroundpoints$OrigOrder <- seq_len(nrow(backgroundpoints))+max(presences$OrigOrder)
  sppdata <- lapply(seq_len(nspp), function(ii)presences[presences$SpeciesID==spps[ii],])

  if(method%in%"grid")  sppBckWtsList <- parallel::mclapply(seq_len(nspp), function(ii)getTileWeights(sppdata[[ii]],backgroundpoints,coord))
  if(method%in%c("quasirandom","psuedorandom"))  sppBckWtsList <- parallel::mclapply(seq_len(nspp), function(ii)getWeights(sppdata[[ii]],backgroundpoints, window, epsilon))

  sppCounts <- parallel::mclapply(seq_len(nspp),function(ii)nrow(sppdata[[ii]]))
  sppBckDatList <- parallel::mclapply(seq_len(nspp),function(ii)rbind(sppdata[[ii]],backgroundpoints))
  sppBckDatList <- parallel::mclapply(seq_len(nspp),function(ii){sppBckDatList[[ii]]$DatasetID <- ii;sppBckDatList[[ii]]})
  sppWtsList <- parallel::mclapply(seq_len(nspp), function(ii)cbind(sppBckDatList[[ii]],
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


# estimateSiteArea <- function(window,site_coords){
#   if(raster::isLonLat(window)){
#     #calculate area based on area function
#     #convert kms to ms
#     area_rast <- raster::area(window)
#     area_study <- raster::mask(area_rast,window)
#     total_area <- cellStats(area_study,sum,na.rm=TRUE)
#     wts <- total_area/raster::extract(area_study,site_coords,na.rm=TRUE)
#   } else {
#     # calculate area based on equal area cell resolution
#     # mode equal area should be in meters
#     cell_area <- raster::res(window)[1]*raster::res(window)[2]
#     n_cell <- length(window[!is.na(window[])])
#     wts <- rep((n_cell*cell_area)/nrow(site_coords),nrow(site_coords))/1000
#   }
#   return(wts)
# }
#


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
  tmp <- raster::coordinates( window)
  tmp1 <- c( range( tmp[,1]), range( tmp[,2]))

  return( tmp1)
  #be careful with extents of rasters though, they tend to have slightly larger dimensions than they ought
}


