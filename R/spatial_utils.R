#' Convert a terra SpatRaster* object to an sf object

#' Convert a terra SpatRaster* object to an sf object.
#' This creates a binary mask for the raster NA vs data and then converts it to
#' POLYGONS. For use internally to turn a rast object in to a mask for polygon
#' clipping.
#' @param x rast* to be converted into an object class \code{sf}
#' @param ... further specifications, see \link{st_as_sf}
#' @name rast_to_sf
#' @importFrom terra ifel as.polygons
#' @importFrom sf st_as_sf st_make_valid st_sf
#' @export
rast_to_sf = function(x, ...) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("package sf required, please install it first")

  ## converts it to binary layer
  window <- terra::ifel(is.na(x),NA,1)

  ## create the sf mulitpolygon
  pols <- terra::as.polygons(window)
  sf.out <- sf::st_as_sf(pols,...)
  sf.valid <- sf::st_make_valid(sf.out)
  x.out <- sf::st_sf(sf.valid)

  return(x.out)
}


#' convert at terra raster to a spatstat image taken from maptools and converted to terra
#' @name terra2im
#' @param rast rast A terra raster to convert to a spatstat image
#' @param factor.col.name string Names of factors in the raster to convert to window
#' @export
#' @importFrom terra res ext hasValues
#' @importFrom spatstat.geom transmat im
#' @examples
#' library(ppmData)
#' library(terra)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' covariates <- rast(lst[1])
#' im.r <- terra2im(covariates)
#'
#' plot(im.r)
terra2im <- function(rast, factor.col.name = NULL){
  if (!requireNamespace("spatstat.geom", quietly = TRUE))
    stop("package spatstat.geom required for coercion")
  if (!requireNamespace("terra", quietly = TRUE))
    stop("package raster required for coercion")
  if (!terra::hasValues(rast))
    stop("values required in RasterLayer object")
  # if (raster::rotated(rast)) {
  #   stop("\n Cannot coerce because the object is rotated.\n Either coerce to SpatialPoints* rast\n or first use the \"rectify\" function")
  # }
  rs <- terra::res(rast)
  orig <- terra::ext(rast)[c(1,3)] + 0.5 * rs
  dm <- dim(rast)[2:1]
  xx <- unname(orig[1] + cumsum(c(0, rep(rs[1], dm[1] - 1))))
  yy <- unname(orig[2] + cumsum(c(0, rep(rs[2], dm[2] - 1))))
  val <- terra::values(rast)
  if (is.factor(rast)) {
    lev <- levels(rast)[[1]]
    if (!is.null(factor.col.name)) {
      if (factor.col.name %in% colnames(lev)) {
        factor.col <- which(colnames(lev) == factor.col.name)
      }
      else {
        stop("'factor.col.name' is not a column name of the raster 'rast'")
      }
    }
    else {
      factor.col <- length(lev)
    }
    val <- factor(val, levels = lev$ID, labels = lev[[factor.col]])
  }
  dim(val) <- dm
  val <- spatstat.geom::transmat(val, from = list(x = "-i",y = "j"), to = "spatstat")
  im <- spatstat.geom::im(val, xcol = xx, yrow = yy)
  return(im)
}

