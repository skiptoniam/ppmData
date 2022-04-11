#' @title psuedoRandomQuad
#' @description A relatively quick way to generate psuedo-random quadrature set
#' for ppm models.
#' @param npoints The number of background points to use. The default is 1000,
#' but it is generally recommended to use more when fitting a ppm.
#' @param window A SpatRaster from terra package which will represent the extent
#'  and resolution of the point process model.
#' @param coord The names of the coordinates. Default is c("X","Y").
#' @export
#' @author Skipton Woolley
#' @examples
#' library(ppmData)
#' library(terra)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' window <- terra::rast(lst[1])
#' res <- psuedoRandomQuad(window,10000)
#' plot(window)
#' points(res[,1:2],pch=".")

## test on big data
# r.list <- list.files("/home/woo457/unimelb_local/data/global/elevation/",full.names = TRUE)
# window <- terra::rast(r.list)
# coord = c('X','Y')
# covariates <- terra::rast(rep(r.list,3))
# npoints = 1000
# pseudoRandomQuad(100000,window,coord)


pseudoRandomQuad <- function(npoints,
                             window,
                             coord = c("X","Y")){

  if(is.null(window)) stop("This function requires a window (terra raster) to work.")
  if(class(window)[1]!="SpatRaster") stop("'window' needs to be a 'SpatRaster' from the 'terra' package.")

  background_sites <- terra::spatSample(x = window,
                                        size = npoints,
                                        na.rm = TRUE,
                                        as.df = TRUE,
                                        values = FALSE,
                                        xy = TRUE)

   bk_sites <- background_sites
   colnames(bk_sites) <- coord

  return(list(quasiPoints = bk_sites, quasiDummy = NULL))

}

#' @title random_background_weights
#' @description Quick way to generate psuedo-random points using the terria
#' package. Weights are taken as the the area(window)/npoints. The default unit
#' is km^2, but other units can be used such as meters squared "m" or hectars
#' "ha". There is an internal check to see if the input data is in lon/lat and
#' if so, it returns the area of the cell in which the point lies.
#' @param presence The observed presence locations
#' @param quadrature The quadrature locations
#' @param window A SpatRaster from terra package which will represent the extent
#'  and resolution of the point process model.
#' @param unit The scale of the area weights, default is kilometers squared "km"
#' but meters squared "m" or hectars "ha" can be used.
#' @author Skipton Woolley
#' @importFrom terra init res extract cellSize mask global
#' @export
#' @examples
#' library(ppmData)
#' library(terra)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' window <- terra::rast(lst[1])
#' presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")[,1:2]
#' quadrature <- pseudoRandomQuad(1000,window)
#' res <- pseudoRandomWeights(presences,quad,window,unit="geo")
#' plot(window)
#' points(res[res$presence==1,1:2],pch=16,cex=0.5,col='blue')
#' points(res[res$presence==0,1:2],pch='.')

pseudoRandomWeights <- function(presences,
                                quadrature,
                                window,
                                coord,
                                unit = c("geo","m","km","ha")){

  # work out what unit you want the areas in.
  unit <- match.arg(unit)

  if(is.null(presences)) stop("This function needs a set of presence sites to work")
  if(is.null(quadrature)) stop("This function needs a set of quadrature sites to work")

  if(is.null(window)) stop("This function requires a window (terra raster) to work.")
  if(class(window)[1]!="SpatRaster") stop("'window' needs to be a 'SpatRaster' from the 'terra' package.")

  ## expanse broke for large raster
  npres <- nrow(presences)
  nquad <- nrow(quadrature)

  # if in geographic coordinates return area in unit of interest.
  if(is.lonlat(window)){
    if(unit=="geo"){ # useful if you just want areas on the geographic scale (this is what deldir gives you)
      areas <- terra::init(window,prod(terra::res(window)))
      bck_wts <- terra::extract(areas,quadrature)
      bck_wts_num <- list2numeric(bck_wts)
    }
    areas <- terra::cellSize(window, unit=unit)
    bck_wts <- terra::extract(areas,quadrature)
    bck_wts_num <- list2numeric(bck_wts)
  }else{
    if(unit=="geo"){
      areas <- terra::cellSize(window, unit=unit)
      areas.m <- terra::mask(areas,window)
      window_area <- terra::global(areas.m, "sum", na.rm=TRUE)
      bck_wts <- rep(as.numeric(window_area)/npoints,npoints)
    }
    areas <- terra::init(window,prod(terra::res(window)))
    window_area <- terra::global(areas, "sum", na.rm=TRUE)
    bck_wts <- rep(as.numeric(window_area)/npoints,npoints)
  }

  pres_wts <- rep(sqrt(.Machine$double.eps),npres)
  all_sites <- rbind(presences,quadrature)
  all_weights <- c(pres_wts,bck_wts_num)
  pres_bck <- c(rep(1,npres),rep(0,nquad))

  res <- data.frame(all_sites,presence=pres_bck,weights=all_weights)

  return(res)

}


