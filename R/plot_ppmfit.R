#' plotting function for a fitted weighted poisson
#' @name plot.ppmfit
#' @title Plotting a fitted spatial point process mode.
#' @description This package is a way to efficiently generate a quasirandom set
#' of background points for presence-only modelling of single or multiple
#' responses. The package was set up to model multiple species presence-only
#' datasets, but could be used for an point process spatial modelling.
#' Quasi-random points are a nice alternative to pseudo-random samples, this is
#' because we can generate a quasirandom sample across and areal region
#' (X and Y coordinates), but we can also extend the dimensions of the
#' quasirandom sample to a N-dimensional hypervolume, which will allow users to
#' effectively sample the spatial and environmental space. This in turn should
#' reduce autocorrelation in quadrature scheme. The weight of each quadrature
#' point is calculated using Dirichlet (Voronoi) Tessellation as provided from
#' the \link[deldir]{deldir} function.
#' @details The approach uses quasi-random sampling to generate a quadrature
#' scheme based (e.g Berman & Turner 1992; Warton & Shepard 2010;
#' Foster et al, 2017). The weights each quasi-random point in the quadrature
#' scheme is calculated using a Dirichlet tessellation (Turner 2020). To improve
#' computational efficiency of the \link[deldir]{deldir} function for a large
#' number of quadrature points, we set up an approach which breaks up the problem into
#' manageable sub-windows. We do this by keeping each deldir call to less that
#' approximately 5000 points (which appears to be the point where the algorithm
#' slows noticeably). To avoid edge effect (large areas on the edges of sub-areas),
#' we rotate the sub-regions three times, the first two use a the nearest largest
#' prime number closest to total number of points (presences+quadrature points)
#' divided by 5000, which allows us to rotate the window on the x and y axis
#' with an subset of the sub-windows. We then calculate a third set of sub-
#' windows using a even set of squares. We then take the median weight across
#' all weight calculated for each point. We then can calculate this in parallel
#' for each species to make it computationally more efficient.
#' @param model a fitted ppm model object.
#' @param window a raster object giving the area over which to generate
#' background points. NA cells are ignored and masked out of returned data.
#' If NULL, a rectangle bounding the extent of \code{presences} will be used as
#' the default window.
#' @param covariates A raster object containing covariates for modelling the
#' point process (best use a Raster stack or Raster brick). This should match
#' the resolution and extent of the window provided. If NULL, only the
#' coordinates of the presences and quadrature points are returned.
#' @export
#' @examples

# library(ppmData)
# path <- system.file("extdata", package = "ppmData")
# lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
# preds <- raster::stack(lst)
# presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#
# ## Plot study area/first raster and snails.
# plot(preds[[1]])
# points(presences,pch=16,cex=0.5,col=gray(0.4))
#
# ## Set up the raster stack (here I'm include coordinates as predictors)
# xy.st <- lonlat_from_window(preds[[1]],mask.na = TRUE)
# preds2 <- stack(xy.st,preds)
# ppmdata <- ppmData(npoints = 1000,presences=presences, window = preds[[1]], covariates = preds2)
#
# predict.ppmfit <- function(model,
#                         window,
#                         covarites,
#                         type = "response"){
#
#
#
# }


