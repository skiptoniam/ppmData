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

true_count <- table( presences$SpeciesID)
ppm_count <- apply( ppmdata0$ppmData$y, 2, function(x) sum( x!=0, na.rm=TRUE))
testthat::expect_identical(c(true_count),c(ppm_count))
apply( ppmdata0$ppmData$wts[!ppmdata0$ppmData$bkg,], 2, function(x) sum( !is.na( x)))

# test with just pres & window
ppmdata1 <- ppmData(npoints = npoints, presences=presences, window = preds[[1]])
ppm_count <- apply( ppmdata1$ppmData$y, 2, function(x) sum( x!=0, na.rm=TRUE))
testthat::expect_identical(c(true_count),c(ppm_count))

# test with just pres, window & covars
ppmdata2 <- ppmData(npoints = npoints, presences=presences, window = preds[[1]], covariates = preds)
ppm_count <- apply( ppmdata2$ppmData$y, 2, function(x) sum( x!=0, na.rm=TRUE))
testthat::expect_identical(c(true_count),c(ppm_count))

# test with no npoints
ppmdata3 <- ppmData(presences=presences, window = preds[[1]])

})
