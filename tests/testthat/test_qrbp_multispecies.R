context('Test ppmData for multiple species')
library(qrbp)

testthat::test_that('Test background generation for multiple species - i.e. for joint/mixture models', {

  set.seed(42)
  library(sdm)
  library(raster)

  file <- system.file("external/", package="qrbp") #
  species <- shapefile(file)
  path <- system.file("external", package="sdm") # path to the folder contains the data
  lst <- list.files(path=path,pattern='asc$',full.names = T)
  preds <- stack(lst)
  projection(preds) <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +units=m +no_defs"

  presences <- cbind(coordinates(species[species$Occurrence == 1,]),SpeciesID=1)
  presences2 <- cbind(coordinates(species[species$Occurrence == 1,]),SpeciesID=2)
  colnames(presences)[1:2] <- c("X","Y")
  presences <- as.data.frame(rbind(presences,presences2))
  presences$SpeciesID <- sample(10,nrow(presences),replace = TRUE)
  # presences <- NULL
  window <- preds[[1]]
  covariates <- NULL#preds
  method <- 'grid'
  interpolation <- 'bilinear'
  npoints <- 1000
  resolution <- ceiling(res(preds)[1])
  covariates <- preds
  control <- ppmData.control()

  backgroundsites <- switch(method,
                            grid = qrbp:::gridMethod(resolution, window),
                            quasirandom = qrbp:::quasirandomMethod(npoints,  window, covariates),
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

})
