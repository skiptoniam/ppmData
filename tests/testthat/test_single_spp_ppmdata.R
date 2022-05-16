context('ppmData single species')
library(ppmData)


testthat::test_that('species mix ippm', {

library(ppmData)
library(terra)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
preds <- rast(lst)
presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
npoints <- 1000

## some test a-roos!
# test with just presences
ppmdata0 <- ppmData(presences = presences, window = preds[[1]])

# test with just pres & window
ppmdata1 <- ppmData(npoints = npoints, presences=presences, window = preds[[1]])

# test with just pres, window & covars
ppmdata2 <- ppmData(npoints = npoints, presences=presences, window = preds[[1]], covariates = preds)

# test with consistant masks.
x <- sum(preds)
preds2 <- mask(preds,x)
ppmdata2a <- ppmData(npoints = npoints, presences=presences, window = preds2[[1]], covariates = preds2)

# test with no npoints
ppmdata3 <- ppmData(presences=presences, window = preds[[1]])

})
