context('Test ppmData for single species')

testthat::test_that('Test ppm data generation for a single species - i.e. for a single species model', {

  library(raster)
  library(ppmData)
  species <- subset(snails,SpeciesID=='Victaphanta lampra')
  path <- system.file("extdata", package = "ppmData")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- raster::stack(lst)

  presences <- species
  window <- preds[[1]]
  covariates <- preds
  interpolation <- 'bilinear'
  npoints <- 1000
  resolution <- 0.1
  control <- ppmData.control()
  coord <- c("X","Y")

  presences <- ppmData:::checkDuplicates(presences,coord)
  ppmData:::checkResolution(resolution,window,control,method)
  window <- ppmData:::checkWindow(presences,window)

  ## grid
  method <- 'grid'
  backgroundpoints <- switch(method,
                            grid = ppmData:::gridMethod(resolution, window,control),
                            quasirandom = ppmData:::quasirandomMethod(npoints,  window, covariates),
                            psuedorandom = ppmData:::randomMethod(npoints,  window, covariates))
  testthat::expect_is(backgroundpoints,'list')


  wts <- ppmData:::getSinglespeciesWeights(presences,backgroundpoints$grid,coord,method)
  testthat::expect_is(wts,'data.frame')

  pbxy <- wts

  sitecovariates <-  ppmData:::getCovariates(pbxy,covariates = preds, interpolation, coord, control)
  testthat::expect_is(sitecovariates,'data.frame')

  parameters <- list(npoints=npoints,resolution=resolution,
                     newresolution=backgroundpoints$newres,method=method,
                     interpolation=interpolation,control=control)
  testthat::expect_is(parameters,'list')

  dat <- ppmData:::assembleQuadData(presences, backgroundpoints$grid, sitecovariates, wts,
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

  ## psuedorandom
  method <- "psuedorandom"
  backgroundpoints <- switch(method,
                            grid = ppmData:::gridMethod(resolution, window),
                            quasirandom = ppmData:::quasirandomMethod(npoints,  window, covariates),
                            psuedorandom = ppmData:::randomMethod(npoints,  window, covariates))
  testthat::expect_is(backgroundpoints,'list')

  wts <- ppmData:::getSinglespeciesWeights(presences,backgroundpoints$grid,coord,method,window)
  testthat::expect_is(wts,'data.frame')

  sitecovariates <-  ppmData:::getCovariates(wts,covariates = preds, interpolation, coord, control)
  testthat::expect_is(sitecovariates,'data.frame')

  parameters <- list(npoints=npoints,resolution=resolution,
                     newresolution=backgroundpoints$newres,method=method,
                     interpolation=interpolation,control=control)
  testthat::expect_is(parameters,'list')

  dat <- ppmData:::assembleQuadData(presences, backgroundpoints$grid, sitecovariates, wts,
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
  backgroundpoints <- switch(method,
                            grid = ppmData:::gridMethod(resolution, window),
                            quasirandom = ppmData:::quasirandomMethod(npoints,  window, covariates,
                                                                   control=control,coord=coord),
                            psuedorandom = ppmData:::randomMethod(npoints,  window, covariates))
  testthat::expect_is(backgroundpoints,'list')

  wts <- ppmData:::getSinglespeciesWeights(presences,backgroundpoints$grid,coord)
  testthat::expect_is(wts,'data.frame')

  sitecovariates <-  ppmData:::getCovariates(wts,covariates = preds, interpolation, coord, control)
  testthat::expect_is(sitecovariates,'data.frame')

  parameters <- list(npoints=npoints,resolution=resolution,
                     newresolution=backgroundpoints$newres,method=method,
                     interpolation=interpolation,control=control)
  testthat::expect_is(parameters,'list')

  dat <- ppmData:::assembleQuadData(presences, backgroundpoints$grid, sitecovariates, wts,
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
