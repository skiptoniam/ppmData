# this is code out of ppm lasso package.
getWeights <- function (presences, backgroundsites, coord = c("X", "Y")){
  sp.col = c(which(names(presences) == coord[1]), which(names(presences) == coord[2]))
  back.col = c(which(names(backgroundsites) == coord[1]), which(names(backgroundsites) == coord[2]))
  X.inc = sort(unique(backgroundsites[, back.col[1]]))[2] - sort(unique(backgroundsites[, back.col[1]]))[1]
  Y.inc = sort(unique(backgroundsites[, back.col[2]]))[2] - sort(unique(backgroundsites[, back.col[2]]))[1]
  back.0X = min(backgroundsites[, back.col[1]]) - floor(min(backgroundsites[, back.col[1]])/X.inc) * X.inc
  back.0Y = min(backgroundsites[, back.col[2]]) - floor(min(backgroundsites[, back.col[2]])/Y.inc) * Y.inc
  X = c(presences[, back.col[1]], backgroundsites[, back.col[1]])
  Y = c(presences[, back.col[2]], backgroundsites[, back.col[2]])
  round.X = round((X - back.0X)/X.inc) * X.inc
  round.Y = round((Y - back.0Y)/Y.inc) * Y.inc
  round.id = paste(round.X, round.Y)
  round.table = table(round.id)
  wts = X.inc * Y.inc/as.numeric(round.table[match(round.id, names(round.table))])
  wts
}


getMultispeciesWeights <- function(presences, backgroundsites, coord = c("X", "Y")){

  nspp <- length(unique(presences[,"SpeciesID"]))
  sppBckWtsList <- lapply(seq_len(nspp), function(ii)getWeights(presences[presences$SpeciesID==ii,],
                                                         backgroundsites,coord))
  sppCounts <- parallel::mclapply(seq_len(nspp),function(ii)nrow(presences[presences$SpeciesID==ii,]))
  sppCountsVec <- do.call(c,sppCounts)
  sppWtsList <- parallel::mclapply(seq_len(nspp), function(ii)sppBckWtsList[[ii]][seq_len(sppCountsVec[ii])])
  sppWtsVec <- do.call(c,sppWtsList)
  bckWts <- sppBckWtsList[[1]][-1:-sppCountsVec[1]]

  dat <- data.frame(wts=c(sppWtsVec,bckWts),
             pres=c(rep(1,length(sppWtsVec)),rep(0,length(bckWts))),
             SpeciesID = c(rep(unique(presences$SpeciesID),sppCountsVec),rep('quad',length(bckWts))))

 return(dat)

}

## old function.
# getWeights <- function(presences, backgroundsites, window, coord=c("X","Y")){
#
#   xy <- rbind(presences[,coord],backgroundsites[,coord])
#   areas <- estimateSiteArea(window = window, site_coords = xy)
#   cell_id <- raster::cellFromXY(window, xy)
#   count_pts <- table(cell_id)
#   wts <- areas/as.numeric(count_pts[match(cell_id, names(count_pts))])
#   # weights <- data.frame(x=xy[,1],y=xy[,2],weights=weights)
#
#   return(wts)
# }

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
# #estimate the total area of all cells in extent in km^2
# estimateWindowArea <- function(window){
#   if(raster::isLonLat(window)){
#     #calculate area based on area function
#     #convert kms to ms
#     area_rast <- area(window)
#     area_study <- mask(area_rast,window)
#     area_of_region <- cellStats(area_study,sum, na.rm=TRUE)
#   } else {
#     # calculate area based on equal area cell resolution
#     # mode equal area should be in meters
#     area_of_region <- (ncell(window[!is.na(window)])  * xres(window) * yres(window))/1000
#   }
#   return(area_of_region)
# }


