context('ppmData multiple species')

testthat::test_that('ppmData multiple species', {

library(ppmData)
library(raster)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
preds <- raster::stack(lst)
presences <- snails
npoints <- 1000

# test with just presences
ppmdata0 <- ppmData(npoints = npoints, presences=presences)

# test with just pres & window
# window <- raster(extent(preds))
# window[] <- 0
ppmdata1 <- ppmData(npoints = npoints, presences=presences, window = preds[[1]])

# test with just pres, window & covars
ppmdata2 <- ppmData(npoints = npoints, presences=presences, window = preds[[1]], covariates = preds)

# test with no npoints
ppmdata3 <- ppmData(presences=presences, window = preds[[1]])

})
