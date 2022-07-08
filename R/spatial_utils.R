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
