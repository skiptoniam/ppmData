#' @title random_background_weights
#' @description Quick way to generate psuedo-random points using the terria
#' package. Weights are taken as the the area(window)/npoints. The default unit
#' is km^2, but other units can be used such as meters squared "m" or hectars
#' "ha". There is an internal check to see if the input data is in lon/lat and
#' if so, it returns the area of the cell in which the point lies.
#' @param window A SpatRaster from terra package which will represent the extent
#'  and resolution of the point process model.
#' @param npoints The number of background points to use. The default is 10000,
#' but it is generally recommended to use more when fitting a ppm.
#' @param unit The scale of the area weights, default is kilometers squared "km"
#' but meters squared "m" or hectars "ha" can be used.
#' @examples
#' library(ppmData)
#' library(terra)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' window <- terra::rast(lst[1])
#' res <- random_background_points(window,10000)
#' plot(window)
#' points(res[,1:2],pch=".")


random_background_points <- function(window = NULL,
                                     npoints = 10000,
                                     unit = "km"){

  # testing
  # r.list <- list.files("/home/woo457/unimelb_local/data/global/elevation/",full.names = TRUE)
  # window <- terra::rast(r.list)
  # require(terra)

  if(is.null(window)) stop("This function requires a window (terra raster) to work.")
  if(class(window)[1]!="SpatRaster") stop("'window' needs to be a 'SpatRaster' from the 'terra' package.")

  background_sites <- terra::spatSample(x = window,
                                        size = npoints,
                                        na.rm = TRUE,
                                        as.df = TRUE,
                                        xy = TRUE)

  ## expanse broke for large raster
  areas <- terra::cellSize(window, unit="km")

  window_area <- terra::global(areas, "sum", na.rm=TRUE)

  if(is.lonlat(window)){
    bck_wts <- terra::extract(areas,background_sites[,-3])
    res <- data.frame(background_sites[,-3],weights=bck_wts[,-1])
  }else{
    bck_wts <- rep(as.numeric(window_area)/npoints,npoints)
    res <- data.frame(background_sites[,-3],weights=bck_wts)
  }

  return(res)

}
