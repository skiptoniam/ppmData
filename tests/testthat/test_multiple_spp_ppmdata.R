testthat::test_that('ppmData multiple species', {

library(ppmData)
library(terra)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
preds <- rast(lst)
presences <- snails
npoints <- 1000

# test with just presences
ppmdata0 <- ppmData(npoints = npoints, presences=presences)
testthat::expect_s3_class(ppmdata0,"ppmData")

true_count <- table( presences$SpeciesID)
ppm_count <- table( ppmdata0$presences.original$SpeciesID)
testthat::expect_identical(c(true_count),c(ppm_count))

# test with just pres & window
ppmdata1 <- ppmData(npoints = npoints, presences=presences, window = preds[[1]])
testthat::expect_s3_class(ppmdata1,"ppmData")

# test with just pres, window & covars
ppmdata2 <- ppmData(npoints = npoints, presences=presences, window = preds[[1]], covariates = preds)
testthat::expect_s3_class(ppmdata2,"ppmData")

# test with no npoints
ppmdata3 <- ppmData(presences=presences, window = preds[[1]])
testthat::expect_s3_class(ppmdata3,"ppmData")

## run in parallel
ppmdata4 <- ppmData(presences=presences, window = preds[[1]], control=list(mc.cores=2))
testthat::expect_s3_class(ppmdata4,"ppmData")


})

testthat::test_that('ppmData multiple species plot', {

  library(ppmData)
  library(terra)
  path <- system.file("extdata", package = "ppmData")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- rast(lst)
  presences <- snails
  npoints <- 1000

  # test with just pres & window
  ppmdata <- ppmData(npoints = npoints, presences=presences, window = preds[[1]])
  p <- plot(ppmdata)
  expect_type(p,"list")

})

