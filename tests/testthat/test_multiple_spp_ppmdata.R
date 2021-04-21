context('ppmData multiple species')

testthat::test_that('ppmData multiple species', {

library(ppmData)
library(raster)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
preds <- raster::stack(lst)
presences <- snails
npoints <- 1000
ppmdata <- ppmData(npoints = npoints, presences=presences, window = preds[[1]], covariates = preds)

})
