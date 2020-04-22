context('Test qrbp for null species')


testthat::test_that('Test background generation for null species - i.e. no species presences', {

  library(qrbp)
  library(raster)
  path <- system.file("extdata", package = "qrbp")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- stack(lst)
  presences <- NULL
  window <- preds[[1]]
  covariates <- NULL#preds
  method <- 'grid'
  interpolation <- 'bilinear'
  npoints <- 1000
  resolution <- 0.1
  covariates <- NULL

  bkgrid <- ppmData(npoints = npoints,
                   resolution = resolution,
                   presences = presences,
                   window = window,
                   covariates = covariates,
                   method = method,
                   interpolation = interpolation)

  method <- "quasirandom"
  bkquasi <- ppmData(npoints = npoints,
                     resolution = resolution,
                     presences = presences,
                     window = window,
                     covariates = covariates,
                     method = method,
                     interpolation = interpolation)

  method <- "random"
  bkrandom <- ppmData(npoints = npoints,
                      resolution = resolution,
                      presences = presences,
                      window = window,
                      covariates = covariates,
                      method = method,
                      interpolation = interpolation)

  testthat::expect_message()
  testthat::expect_message()
})
