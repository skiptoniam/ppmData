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

getSinglespeciesWeights <- function(presences,backgroundsites,coord){

  backgroundsites$SpeciesID <- "quad"
  wts <- qrbp:::getWeights(presences,backgroundsites,coord)
  dat <- cbind(rbind(presences,backgroundsites),
               pres=c(rep(1,sppCounts[[ii]]),
                     rep(0,nrow(backgroundsites))),
               wts=wts)
  return(dat)
}

getMultispeciesWeights <- function(presences, backgroundsites, coord = c("X", "Y")){

  presences$oo <- seq_len(nrow(presences))
  nspp <- length(unique(presences[,"SpeciesID"]))
  backgroundsites$SpeciesID <- "quad"
  backgroundsites$oo <- seq_len(nrow(backgroundsites))+max(presences$oo)
  sppdata <- lapply(seq_len(nspp), function(ii)presences[presences$SpeciesID==ii,])
  sppBckWtsList <- parallel::mclapply(seq_len(nspp), function(ii)qrbp:::getWeights(sppdata[[ii]],backgroundsites,coord))
  sppCounts <- parallel::mclapply(seq_len(nspp),function(ii)nrow(sppdata[[ii]]))
  sppBckDatList <- parallel::mclapply(seq_len(nspp),function(ii)rbind(sppdata[[ii]],backgroundsites))
  sppBckDatList <- parallel::mclapply(seq_len(nspp),function(ii){sppBckDatList[[ii]]$DatasetID <- ii;sppBckDatList[[ii]]})
  sppWtsList <- parallel::mclapply(seq_len(nspp), function(ii)cbind(sppBckDatList[[ii]],
                                                                    pres=c(rep(1,sppCounts[[ii]]),rep(0,nrow(backgroundsites))),
                                                                    wts=sppBckWtsList[[ii]]))
 dat <- do.call(rbind,sppWtsList)
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


