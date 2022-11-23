#' Function to covert a spatstat quadrature
#' @name as.spatstat
#' @param object A ppmData object
#' @param \\dots Ignored
#' @export
#' @importFrom spatstat.geom ppp quad
#' @examples
#' \dontrun{
#' library(ppmData)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' preds <- rast(lst)
#' window <- preds[[1]]
#' presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#' object <- ppmData(npoints = 10000, presences=presences, window = window, covariates = preds)
#' ss.object <- as.spatstat(object)
#' ppmfit <- ppm(ss.object$Q~1+poly(min_temp_coldest_month,2)+poly(annual_mean_precip,2),data=ss.object$covariates,Poisson())
#' pred <- predict(ppmfit, covariates = ss.object$covariates,type="trend")
#' }

as.spatstat <- function(object, ...){

  if(!isa(object,"ppmData"))
    stop("This function requires a ppmData object")

  ## create the owin
  quad.win <- as.owin(object)

  ## if covariates transform to spatstat image
  if(!is.null(object$covariates)){
    covariates <- lapply(as.list(object$covariates), function(ii) terra2im(ii))
    names(covariates) <- names(object$covariates)
  } else {
    covariates <- NULL
  }

  ## coordinate names
  xname <- object$params$coord[1]
  yname <- object$params$coord[2]

  ## presences
  px <- object$presences.cleaned[,xname]
  py <- object$presences.cleaned[,yname]

  ## quadrature
  qx <- object$ppmData[object$ppmData$presence==0,xname]
  qy <- object$ppmData[object$ppmData$presence==0,yname]

  ## weights
  wt <- object$ppmData$weights

  ## create a spatstat quad
  pres.dat <- spatstat.geom::ppp(px, py, window = quad.win, check = FALSE)
  quad.dat <- spatstat.geom::ppp(qx, qy, window = quad.win, check = FALSE)

  if(object$params$quad.method%in%c("quasi.random","psuedo.random")){
    quad.out <- spatstat.geom::quad(data = pres.dat, dummy = quad.dat, w = wt,
                                    param = list(weight = list(method = "dirichlet")))
  }

  if(object$params$quad.method%in%c("grid")){
    quad.out <- spatstat.geom::quad(data = pres.dat, dummy = quad.dat, w = wt,
                                    param = list(weight = list(method = "grid")))
  }

  return(list(Q = quad.out, covariates = covariates))

}

#' Function to covert a spatstat quadrature
#' @name as.ppp
#' @param X A ppmData object
#' @param \\dots Ignored
#' @param fatal	Logical value specifying what to do if the data cannot be converted. See Details.
#' @importFrom spatstat.geom as.ppp
#' @export
#' @examples
#' \dontrun{
#' library(ppmData)
#' path <- system.file("extdata", package = "ppmData")
#' lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#' preds <- rast(lst)
#' window <- preds[[1]]
#' presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#' object <- ppmData(npoints = 10000, presences=presences, window = window, covariates = preds)
#' X <- as.ppp(object)
#' }

as.ppp.ppmData <- function(X, ..., fatal=TRUE){

  if(isa(X$window,"SpatRaster")){
    quad.win <- spatstat.geom::as.owin(terra2im(X$window))
  } else {
    bb <- get_bbox(X$presences.original)
    quad.win <- spatstat.geom::owin(bb[1:2],bb[3:4])
  }

  ## coordinate names
  xname <- object$params$coord[1]
  yname <- object$params$coord[2]

  ## presences
  px <- object$presences.cleaned[,xname]
  py <- object$presences.cleaned[,yname]

  ## create a spatstat quad
  ppp.out <- spatstat.geom::ppp(px, py, window = quad.win, check = FALSE)

  return(ppp.out)

}

#' Convert a ppmData object to an owin from the SpatStat package.
#' @param W Window.
#' @param ... Unused.
#' @param fatal Unused.
#' @importFrom spatstat.geom as.owin
#' @export
as.owin.ppmData <- function(W,...,fatal=TRUE){

  if(isa(W$window,"SpatRaster")){
  win <- spatstat.geom::as.owin(terra2im(object$window))
  } else {
  bb <- get_bbox(W$presences.original)
  win <- spatstat.geom::owin(bb[1:2],bb[3:4])
  }

 return(win)
}

