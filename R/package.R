#' @name ppmData
#' @description ppmData is a package for setting up quadrature to implement
#' spatial Poisson Point process models and extensions. The approach uses quasi-
#' random sampling (Grafston & Tille, 2013, Foster et al., 2018) to generate a
#' quadrature scheme for numerical approximation of a Poisson point process
#' model (Berman & Turner 1992; Warton & Shepard 2010). Quasi-random sampling
#' quadrature are form of spatially-balanced survey design or point
#' stratification that aims to reduce the frequency of placing samples close to
#' each other (relative to pseudo-random or grid designs). A quasi-random
#' quadrature design improves efficiency of background point sampling (and
#' subsequent modelling) by reducing the amount of spatial auto-correlation
#' between data implying that each sample is providing as much unique
#' information as possible (Grafston & Tille, 2013, Foster et al., 2018) and
#' thus reducing low errors for geostatistical prediction (Diggle & Ribeiro,
#' 2007). Because the quasi-random design is not on a regular grid we use
#' Dirichlet tessellation to generate polygons for each point in the quadrature
#' scheme. Areal weights are then derived from these polygons.
#'
#' @author Skipton Woolley <skip.woolley@csiro.au> & Scott Foster <scott.foster@data61.csiro.au>
#'
#' @references Diggle, P. J., P. J. Ribeiro, Model-based Geostatistics. Springer Series in Statistics. Springer, 2007.
#'
#' Foster, S.D., Monk, J., Lawrence, E., Hayes, K.R., Hosack, G.R. and Przeslawski, R., 2018. Statistical considerations
#' for monitoring and sampling. Field manuals for marine sampling to monitor Australian waters, pp.23-41.
#'
#' Grafstrom, Anton, and Yves Tille. Doubly balanced spatial sampling with spreading and restitution of auxiliary totals.
#' Environmetrics 24.2 (2013): 120-131.
#'
#' Warton, D. I., and L. C. Shepherd. Poisson point process models solve the pseudo-absence problem for presence-only data #'in ecology. The Annals of Applied Statistics 4.3 (2010): 1383-1402.
#'
#' @docType package
#' @useDynLib "ppmData", .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

