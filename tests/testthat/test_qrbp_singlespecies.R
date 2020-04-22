context('Test ppmData for single species')

testthat::test_that('Test ppm data generation for a single species - i.e. for joint/mixture models', {

  set.seed(42)
  library(raster)
  library(qrbp)

  species <- subset(snails,SpeciesID=='Victaphanta lampra')
  path <- system.file("extdata", package = "qrbp")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- stack(lst)

  presences <- species
  window <- preds[[1]]
  covariates <- NULL
  interpolation <- 'bilinear'
  npoints <- 1000
  resolution <- 0.1
  control <- ppmData.control()
  coord <- c("X","Y")

  ## grid
  method <- 'grid'
  backgroundsites <- switch(method,
                            grid = qrbp:::gridMethod(resolution, window),
                            quasirandom = qrbp:::quasirandomMethod(npoints,  window, covariates),
                            random = qrbp:::randomMethod(npoints,  window, covariates))
  testthat::expect_is(backgroundsites,'list')

  wts <- qrbp:::getSinglespeciesWeights(presences,backgroundsites$grid,coord)
  testthat::expect_is(wts,'data.frame')

  pbxy <- wts

  sitecovariates <-  qrbp:::getCovariates(pbxy,covariates = preds, interpolation, coord, control)
  testthat::expect_is(sitecovariates,'data.frame')

  parameters <- list(npoints=npoints,resolution=resolution,
                     newresolution=backgroundsites$newres,method=method,
                     interpolation=interpolation,control=control)
  testthat::expect_is(parameters,'list')

  dat <- qrbp:::assembleQuadData(presences, backgroundsites$grid, sitecovariates, wts,
                                 coord, parameters, control=control)
  testthat::expect_is(dat,'list')

  bkgrid<- ppmData(npoints = npoints,
                   resolution = resolution,
                   presences = presences,
                   window = window,
                   covariates = covariates,
                   method = method,
                   interpolation = interpolation)

  testthat::expect_type(bkgrid,"list")
  testthat::expect_that(bkgrid, testthat::is_a("ppmdata"))

  ## random
  method <- "random"
  backgroundsites <- switch(method,
                            grid = qrbp:::gridMethod(resolution, window),
                            quasirandom = qrbp:::quasirandomMethod(npoints,  window, covariates),
                            random = qrbp:::randomMethod(npoints,  window, covariates))
  testthat::expect_is(backgroundsites,'list')

  wts <- qrbp:::getSinglespeciesWeights(presences,backgroundsites$grid,coord)
  testthat::expect_is(wts,'data.frame')

  sitecovariates <-  qrbp:::getCovariates(wts,covariates = preds, interpolation, coord, control)
  testthat::expect_is(sitecovariates,'data.frame')

  parameters <- list(npoints=npoints,resolution=resolution,
                     newresolution=backgroundsites$newres,method=method,
                     interpolation=interpolation,control=control)
  testthat::expect_is(parameters,'list')

  dat <- qrbp:::assembleQuadData(presences, backgroundsites$grid, sitecovariates, wts,
                                 coord, parameters, control=control)
  testthat::expect_is(dat,'list')

  bkrandom <- ppmData(npoints = npoints,
                   resolution = resolution,
                   presences = presences,
                   window = window,
                   covariates = covariates,
                   method = method,
                   interpolation = interpolation)

  testthat::expect_type(bkrandom,"list")
  testthat::expect_that(bkrandom, testthat::is_a("ppmdata"))

  ## quasirandom
  method <- "quasirandom"
  backgroundsites <- switch(method,
                            grid = qrbp:::gridMethod(resolution, window),
                            quasirandom = qrbp:::quasirandomMethod(npoints,  window, covariates,
                                                                   control=control,coord=coord),
                            random = qrbp:::randomMethod(npoints,  window, covariates))
  testthat::expect_is(backgroundsites,'list')

  wts <- qrbp:::getSinglespeciesWeights(presences,backgroundsites$grid,coord)
  testthat::expect_is(wts,'data.frame')

  sitecovariates <-  qrbp:::getCovariates(wts,covariates = preds, interpolation, coord, control)
  testthat::expect_is(sitecovariates,'data.frame')

  parameters <- list(npoints=npoints,resolution=resolution,
                     newresolution=backgroundsites$newres,method=method,
                     interpolation=interpolation,control=control)
  testthat::expect_is(parameters,'list')

  dat <- qrbp:::assembleQuadData(presences, backgroundsites$grid, sitecovariates, wts,
                                 coord, parameters, control=control)
  testthat::expect_is(dat,'list')

  bkqrand<- ppmData(npoints = npoints,
                   resolution = resolution,
                   presences = presences,
                   window = window,
                   covariates = covariates,
                   method = method,
                   interpolation = interpolation)

  testthat::expect_type(bkqrand,"list")
  testthat::expect_that(bkqrand, testthat::is_a("ppmdata"))



})
