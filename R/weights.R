# this is code out of ppm lasso package.
ppmWeights <- function (presences, backgroundsites, coord = c("X", "Y")){
  sp.col = c(which(names(sp.xy) == coord[1]), which(names(sp.xy) == coord[2]))
  quad.col = c(which(names(quad.xy) == coord[1]), which(names(quad.xy) == coord[2]))
  X.inc = sort(unique(quad.xy[, quad.col[1]]))[2] - sort(unique(quad.xy[, quad.col[1]]))[1]
  Y.inc = sort(unique(quad.xy[, quad.col[2]]))[2] - sort(unique(quad.xy[, quad.col[2]]))[1]
  quad.0X = min(quad.xy[, quad.col[1]]) - floor(min(quad.xy[, quad.col[1]])/X.inc) * X.inc
  quad.0Y = min(quad.xy[, quad.col[2]]) - floor(min(quad.xy[, quad.col[2]])/Y.inc) * Y.inc
  X = c(sp.xy[, quad.col[1]], quad.xy[, quad.col[1]])
  Y = c(sp.xy[, quad.col[2]], quad.xy[, quad.col[2]])
  round.X = round((X - quad.0X)/X.inc) * X.inc
  round.Y = round((Y - quad.0Y)/Y.inc) * Y.inc
  round.id = paste(round.X, round.Y)
  round.table = table(round.id)
  wt = X.inc * Y.inc/as.numeric(round.table[match(round.id, names(round.table))])
  wt
}

## old function.
getWeights <- function(presences, backgroundsites, window, coord=c("X","Y")){

  xy <- rbind(presences,background_sites)
  areas <- estimateSiteArea(window = window, site_coords = xy)
  cell_id <- raster::cellFromXY(window, xy)
  count_pts <- table(cell_id)
  weights <- areas/as.numeric(count_pts[match(cell_id, names(count_pts))])
  weights <- data.frame(x=xy[,1],y=xy[,2],weights=weights)

  return(weights)
}

estimateSiteArea <- function(window,site_coords){
  if(raster::isLonLat(window)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- raster::area(window)
    area_study <- raster::mask(area_rast,window)
    total_area <- cellStats(area_study,sum,na.rm=TRUE)
    wts <- total_area/raster::extract(area_study,site_coords,na.rm=TRUE)
  } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    cell_area <- raster::res(window)[1]*raster::res(window)[2]
    n_cell <- length(window[!is.na(window[])])
    wts <- rep((n_cell*cell_area)/nrow(site_coords),nrow(site_coords))/1000
  }
  return(wts)
}

#estimate the total area of all cells in extent in km^2
estimateWindowArea <- function(window){
  if(raster::isLonLat(window)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- area(window)
    area_study <- mask(area_rast,window)
    area_of_region <- cellStats(area_study,sum, na.rm=TRUE)
  } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    area_of_region <- (ncell(window[!is.na(window)])  * xres(window) * yres(window))/1000
  }
  return(area_of_region)
}


