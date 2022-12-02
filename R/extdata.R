#' Raster data set of environmental and social covariates clipped to Tasmania
#'
#' @name raster_datasets
#' @details The a subset of the ANUCLIM v6 layers are included in the package, but have been clipped to the Tasmanian extent. A distance from roads layer is also included which is calculated as the distance from main roads, based on the road network dataset from OpenStreetMap (OSM) across the geographic area of Tasmania
#'
#' \itemize{
#' \item{ precip_of_driest_month.tif}{ Precipitation of the driest month}
#' \item{ min_temp_coldest_month.tif}{ Minimum temperature of the coldest month}
#' \item{ max_temp_coldest_month.tif}{ Maximum temperature of the warmest month}
#' \item{ annual_mean_temp.tif}{ Annual mean temperature}
#' \item{ annual_mean_precip.tif}{ Annual mean precipitation}
#' \item{ distance_from_main_roads.tif}{ Distance from main roads}
#' }
#' @references Xu, Tingbao, and Michael Hutchinson. "ANUCLIM version 6.1 user
#' guide." The Australian National University, Fenner School of Environment and
#' Society, Canberra 90 (2011).
#'
#' OpenStreetMap, (2020): OpenStreetMap - Road Network (Australia) 2020;
#' accessed from AURIN.
#' @examples
#' \dontrun{
#' ## This is how you would load these data
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' covariates <- rast(lst)
#' }
NULL
