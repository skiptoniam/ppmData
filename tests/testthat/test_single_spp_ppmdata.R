# context('ppmData single species')
# library(ppmData)


testthat::test_that('ppmdata single species quasi', {

library(ppmData)
library(terra)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
preds <- rast(lst)
presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
npoints <- 1000

## some test a-roos!
# test with just presences
ppmdata0 <- ppmData(npoints = npoints, presences = presences)
plot(ppmdata0)


# test with just pres & window
ppmdata1 <- ppmData(presences=presences, window = preds[[1]], npoints = npoints)
plot(ppmdata1)

# test with just pres, window & covars
ppmdata2 <- ppmData(presences=presences, window = preds[[1]], covariates = preds, npoints = npoints)
plot(ppmdata2)

# test with consistant masks.
x <- sum(preds)
preds2 <- mask(preds,x)
ppmdata2a <- ppmData(npoints = npoints, presences=presences, window = preds2[[1]], covariates = preds2)

# test with no npoints
ppmdata3 <- ppmData(presences=presences, window = preds[[1]])
plot(ppmdata3)


})

testthat::test_that('ppmdata single species random', {

  library(ppmData)
  library(terra)
  path <- system.file("extdata", package = "ppmData")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- rast(lst)
  presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
  npoints <- 1000

  ## some test a-roos!
  # test with just presences
  ppmdata0 <- ppmData(npoints = npoints, presences = presences, quad.method = 'pseudo.random')
  plot(ppmdata0)


  # test with just pres & window
  ppmdata1 <- ppmData(presences=presences, window = preds[[1]],
                      npoints = 10000, quad.method = 'pseudo.random')
  plot(ppmdata1)

  # test with just pres, window & covars
  ppmdata2 <- ppmData(presences=presences, window = preds[[1]],
                      covariates = preds, npoints = npoints,
                      quad.method = 'pseudo.random')
  plot(ppmdata2)

  # test with consistant masks.
  x <- sum(preds)
  preds2 <- mask(preds,x)
  ppmdata2a <- ppmData(npoints = npoints, presences=presences,
                       window = preds2[[1]], covariates = preds2,
                       quad.method = 'pseudo.random')

  plot(ppmdata2a)

  # test with no npoints
  ppmdata3 <- ppmData(presences=presences, window = preds[[1]],
                      quad.method = 'pseudo.random')
  plot(ppmdata3)


})

testthat::test_that('ppmdata single species grid', {

  library(ppmData)
  library(terra)
  path <- system.file("extdata", package = "ppmData")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- rast(lst)
  presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
  npoints <- 1000

  ## some test a-roos!
  # test with just presences
  ppmdata0 <- ppmData(npoints = npoints, presences = presences, quad.method = 'grid')
  plot(ppmdata0)


  # test with just pres & window
  ppmdata1 <- ppmData(presences=presences, window = preds[[1]],
                      npoints = npoints, quad.method = 'grid')
  plot(ppmdata1)
  sum(ppmdata1$ppmData$weights)


  # test with just pres, window & covars
  ppmdata2 <- ppmData(presences=presences, window = preds[[1]],
                      covariates = preds, npoints = npoints, quad.method = 'grid')
  plot(ppmdata2)

  # test with consistant masks.
  x <- sum(preds)
  preds2 <- mask(preds,x)
  ppmdata2a <- ppmData(npoints = 50000, presences=presences,
                       window = preds2[[1]], covariates = preds2,
                       quad.method = 'grid')

    # test with no npoints
  ppmdata3 <- ppmData(presences=presences, window = preds[[1]], quad.method = 'grid')
  plot(ppmdata3)


})
