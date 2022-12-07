testthat::test_that('test spatstat utils', {

  library(ppmData)
  library(terra)
  path <- system.file("extdata", package = "ppmData")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- rast(lst)
  presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
  npoints <- 1000

  pd <- ppmData(presences=presences, window = preds[[1]], covariates = preds, npoints = npoints)

  ## check to see the convertion to spatstat works
  ss <- as.spatstat(pd)
  expect_s3_class(ss$Q,"quad")
  expect_s3_class(ss$covariates$annual_mean_precip,"im")

  ## check to see ppmdata to owin works
  expect_s3_class(as.owin(pd),"owin")

  ## check to see that ppmdata to ppp works
  expect_s3_class(class(as.ppp(pd)),"ppp")

  })
