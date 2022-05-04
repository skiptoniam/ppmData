##### Terra untilities #####
#' @name xyFromWindow
#' @param window
#' @param coord.name
#' @param mask.na
#' @export
#' @importFrom terra xyFromCell ncell mask
#' @examples
#' library(ppmData)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' covariates <- rast(lst)
#' xy.r <- xyFromWindow(covariates[[1]],mask.na=TRUE)

xyFromWindow <- function(window, coord.name = c("X","Y"), mask.na=FALSE){

  grid.locations <- terra::xyFromCell(window,1:ncell(window))
  X <- Y <- window
  X[] <- grid.locations[,1]
  Y[] <- grid.locations[,2]
  names(X) <- coord.name[1]
  names(Y) <- coord.name[2]

  st.lonlat <- c(X,Y)

  if(mask.na){
    st.lonlat <- mask(st.lonlat,window)
  }

  return(st.lonlat)

}

## convert at terra raster to a spatstat image taken from maptools and converted to terra
#' @name terra2im
#' @param rast
#' @param factor.col.name
#' @export
#' @importFrom terra res ext hasValues
#' @importFrom spatstat.geom transmat im
#' @examples
#' library(ppmData)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' covariates <- rast(lst[1])
#' im.r <- terra2im(covariates)
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
