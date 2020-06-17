context('Test ppmData for multiple species')

testthat::test_that('Test background generation for multiple species - i.e. for joint/mixture models', {

  library(raster)
  path <- system.file("extdata", package = "qrbp")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- raster::stack(lst)

  presences <- snails
  window <- preds[[1]]
  covariates <- NULL#preds
  method <- 'grid'
  interpolation <- 'bilinear'
  npoints <- 1000
  resolution <- res(preds)[1]
  covariates <- preds
  control <- ppmData.control()
  coord <- c("X","Y")


  ## test built in functions.
  presences$SpeciesID <- factor(presences$SpeciesID)
  presences <- qrbp:::checkDuplicates(presences,coord)

  backgroundsites <- switch(method,
                            grid = qrbp:::gridMethod(resolution, window,control),
                            quasirandom = qrbp:::quasirandomMethod(npoints,  window, covariates,control,coord),
                            random = qrbp:::randomMethod(npoints,  window, covariates))
  testthat::expect_is(backgroundsites,'list')

  wts <- qrbp:::getMultispeciesWeights(presences,backgroundsites$grid,coord)
  testthat::expect_is(wts,'data.frame')

  sitecovariates <-  qrbp:::getCovariates(wts, covariates = preds,
                                          interpolation, coord, control)
  testthat::expect_is(sitecovariates,'data.frame')

  parameters <- list(npoints=npoints,resolution=resolution,
                     newresolution=backgroundsites$newres,method=method,
                     interpolation=interpolation,control=control)
  testthat::expect_is(parameters,'list')

  dat <- qrbp:::assembleQuadData(presences, backgroundsites$grid, sitecovariates, wts,
                          coord, parameters, control=control)
  testthat::expect_is(dat,'list')

  rm(list=ls())

  ## test function alone
  path <- system.file("extdata", package = "qrbp")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- raster::stack(lst)

  presences <- snails
  window <- preds[[1]]
  covariates <- NULL#preds
  method <- 'grid'
  interpolation <- 'bilinear'
  npoints <- 1000
  resolution <- res(preds)[1]
  covariates <- preds
  control <- ppmData.control()
  coord <- c("X","Y")

  bkgrid<- ppmData(npoints = npoints,
                   resolution = resolution,
                   presences = presences,
                   window = window,
                   covariates = covariates,
                   method = method,
                   interpolation = interpolation,
                   control=ppmData.control())
  testthat::expect_type(bkgrid,"list")
  testthat::expect_that(bkgrid, testthat::is_a("ppmdata"))

  ## quasirandom
  method <- "quasirandom"
  backgroundsites <- switch(method,
                            grid = qrbp:::gridMethod(resolution, window),
                            quasirandom = qrbp:::quasirandomMethod(npoints,  window, covariates,control,coord),
                            random = qrbp:::randomMethod(npoints,  window, covariates))
  testthat::expect_is(backgroundsites,'list')

  wts <- qrbp:::getMultispeciesWeights(presences,backgroundsites$grid)
  testthat::expect_is(wts,'data.frame')

  sitecovariates <-  qrbp:::getCovariates(wts, covariates = preds,
                                          interpolation, coord, control)
  testthat::expect_is(sitecovariates,'data.frame')

  parameters <- list(npoints=npoints,resolution=resolution,
                     newresolution=backgroundsites$newres,method=method,
                     interpolation=interpolation,control=control)
  testthat::expect_is(parameters,'list')

  dat <- qrbp:::assembleQuadData(presences, backgroundsites$grid, sitecovariates, wts,
                                 coord, parameters, control=control)
  testthat::expect_is(dat,'list')

  bkquasi <- ppmData(npoints = npoints,
                   resolution = resolution,
                   presences = presences,
                   window = window,
                   covariates = covariates,
                   method = method,
                   interpolation = interpolation,
                   control=ppmData.control())
  testthat::expect_type(bkquasi,"list")
  testthat::expect_that(bkquasi, testthat::is_a("ppmdata"))

})
