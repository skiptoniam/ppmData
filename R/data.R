#' Presence only observations of snails from the family Rhytididae located in Tasmania, Australia.
#'
#' @source Atlas of Living Australia. \url{http://www.ala.org/}
#' @format A data frame with columns:
#' \describe{
#'  \item{X}{Longitude.}
#'  \item{Y}{Latitude.}
#'  \item{SpeciesID}{Species taxonomic name.}
#' }
#' @examples
#' \dontrun{
#'  snails
#' }
"snails"

#' Raster stack of environmental and bias covariates for Tasmania, Australia.
#'
#' @source Environmental covariates come from [Worldclim](http://www.worldclim.org), human population density for Tasmania comes from [NASA](https://neo.sci.gsfc.nasa.gov/view.php?datasetId=SEDAC_POP), distance from roads was estimated from main roads based on this dataset from the tasmanian goverment [Listmap](http://maps.thelist.tas.gov.au/listmap/app/list/map).

#' @format A raster stack with the layers:
#' \describe{
#'  \item{"annual_mean_temperature"}{Annual Mean Temperature from worldclim.}
#'  \item{"annual_precipitation"}{Annual Precipitation from worldclim.}
#'  \item{"max_temperature_of_warmest_month"}{Max temperature of warmest month from worldclim.}
#'  \item{"min_temperature_of_coldest_month"}{Min temperature of coldest month from worldclim.}
#'  \item{"precipitation_of_driest_month"}{Precipition of driest month from worldclim}
#'  \item{"human_population_density"}{Human population density from NASA.}
#'  \item{"distance_from_main_roads"}{Distance from main roads calulated from Tasmanian road network.}
#'  }
#'  @name tassie_preds
"tassie_preds"
