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
testthat::expect_s3_class(ppmdata0,"ppmData")

## error if no presences
testthat::expect_error(ppmData(npoints = npoints))

## error if presences in the wrong format message to tell correct format.
testthat::expect_error(ppmData(presences = "a",npoints = npoints))

# test with just pres & window
ppmdata1 <- ppmData(presences=presences, window = preds[[1]], npoints = npoints)
testthat::expect_s3_class(ppmdata1,"ppmData")

## test with a dodgy window
testthat::expect_error(ppmData(presences=presences, window = "window", npoints = npoints))
## show throw an error if a bounding box is given rather than a spatial raster
win <- ppmData:::get_bbox(presences[,1:2])
testthat::expect_error(ppmData(presences=presences, window = win, npoints = npoints))

# test with just pres, window & covars
ppmdata2 <- ppmData(presences=presences, window = preds[[1]], covariates = preds, npoints = npoints)
testthat::expect_s3_class(ppmdata2,"ppmData")
testthat::expect_error(ppmData(presences=presences, window = preds[[1]], covariates = as.data.frame(preds), npoints = npoints))


# test with consistant masks.
x <- sum(preds)
preds2 <- mask(preds,x)
ppmdata2a <- ppmData(npoints = npoints, presences=presences, window = preds2[[1]], covariates = preds2)
testthat::expect_s3_class(ppmdata2a,"ppmData")

# test with no npoints
ppmdata3 <- ppmData(presences=presences, window = preds[[1]])
testthat::expect_s3_class(ppmdata3,"ppmData")


ppmdata4 <- ppmData(npoints=1000,
                    presences=presences,
                    window = preds[[1]],
                    unit="m")
testthat::expect_s3_class(ppmdata4,"ppmData")

ppmdata5 <- ppmData(npoints=1000,
                    presences=presences,
                    window = preds[[1]],
                    unit="km")
testthat::expect_s3_class(ppmdata5,"ppmData")

ppmdata6 <- ppmData(npoints=1000,
                    presences=presences,
                    window = preds[[1]],
                    unit="ha")
testthat::expect_s3_class(ppmdata6,"ppmData")

## check areal sums
testthat::expect_gt(sum(ppmdata4$ppmData$weights),6e10)
testthat::expect_gt(sum(ppmdata5$ppmData$weights),6e4)
testthat::expect_gt(sum(ppmdata6$ppmData$weights),6e6)

## throw some errors if wrong unit
testthat::expect_error(ppmData(npoints=1000,
                               presences=presences,
                               window = preds[[1]],
                               unit="meters"))



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
  ppmdata0 <- ppmData(npoints = npoints,
                      presences = presences,
                      quad.method="pseudo.random")
  testthat::expect_s3_class(ppmdata0,"ppmData")

  ## error if no presences
  testthat::expect_error(ppmData(npoints = npoints,
                                 quad.method="pseudo.random"))

  ## error if presences in the wrong format message to tell correct format.
  testthat::expect_error(ppmData(presences = "a",npoints = npoints,
                                 quad.method="pseudo.random"))

  # test with just pres & window
  ppmdata1 <- ppmData(presences=presences,
                      window = preds[[1]],
                      npoints = npoints,
                      quad.method="pseudo.random")
  testthat::expect_s3_class(ppmdata1,"ppmData")

  ## test with a dodgy window
  testthat::expect_error(ppmData(presences=presences,
                                 window = "window",
                                 npoints = npoints,
                                 quad.method="pseudo.random"))
  ## show throw an error if a bounding box is given rather than a spatial raster
  win <- ppmData:::get_bbox(presences[,1:2])
  testthat::expect_error(ppmData(presences=presences,
                                 window = win,
                                 npoints = npoints,
                                 quad.method="pseudo.random"))

  # test with just pres, window & covars
  ppmdata2 <- ppmData(presences=presences,
                      window = preds[[1]],
                      covariates = preds,
                      npoints = npoints,
                      quad.method="pseudo.random")
  testthat::expect_s3_class(ppmdata2,"ppmData")
  testthat::expect_error(ppmData(presences=presences,
                                 window = preds[[1]],
                                 covariates = as.data.frame(preds),
                                 npoints = npoints,
                                 quad.method="pseudo.random"))


  # test with consistant masks.
  x <- sum(preds)
  preds2 <- mask(preds,x)
  ppmdata2a <- ppmData(npoints = npoints,
                       presences=presences,
                       window = preds2[[1]],
                       covariates = preds2,
                       quad.method="pseudo.random")
  testthat::expect_s3_class(ppmdata2a,"ppmData")

  # test with no npoints
  ppmdata3 <- ppmData(presences=presences,
                      window = preds[[1]],
                      quad.method="pseudo.random")
  testthat::expect_s3_class(ppmdata3,"ppmData")


  ppmdata4 <- ppmData(npoints=1000,
                      presences=presences,
                      window = preds[[1]],
                      unit="m",,
                      quad.method="pseudo.random")
  testthat::expect_s3_class(ppmdata4,"ppmData")

  ppmdata5 <- ppmData(npoints=1000,
                      presences=presences,
                      window = preds[[1]],
                      unit="km",
                      quad.method="pseudo.random")
  testthat::expect_s3_class(ppmdata5,"ppmData")

  ppmdata6 <- ppmData(npoints=1000,
                      presences=presences,
                      window = preds[[1]],
                      unit="ha",
                      quad.method="pseudo.random")
  testthat::expect_s3_class(ppmdata6,"ppmData")


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
  ppmdata0 <- ppmData(npoints = npoints,
                      presences = presences,
                      quad.method="grid")
  testthat::expect_s3_class(ppmdata0,"ppmData")

  ## error if no presences
  testthat::expect_error(ppmData(npoints = npoints,
                                 quad.method="grid"))

  ## error if presences in the wrong format message to tell correct format.
  testthat::expect_error(ppmData(presences = "a",npoints = npoints,
                                 quad.method="grid"))

  # test with just pres & window
  ppmdata1 <- ppmData(presences=presences,
                      window = preds[[1]],
                      npoints = npoints,
                      quad.method="grid")
  testthat::expect_s3_class(ppmdata1,"ppmData")

  ## test with a dodgy window
  testthat::expect_error(ppmData(presences=presences,
                                 window = "window",
                                 npoints = npoints,
                                 quad.method="grid"))
  ## show throw an error if a bounding box is given rather than a spatial raster
  win <- ppmData:::get_bbox(presences[,1:2])
  testthat::expect_error(ppmData(presences=presences,
                                 window = win,
                                 npoints = npoints,
                                 quad.method="grid"))

  # test with just pres, window & covars
  ppmdata2 <- ppmData(presences=presences,
                      window = preds[[1]],
                      covariates = preds,
                      npoints = npoints,
                      quad.method="grid")
  testthat::expect_s3_class(ppmdata2,"ppmData")
  testthat::expect_error(ppmData(presences=presences,
                                 window = preds[[1]],
                                 covariates = as.data.frame(preds),
                                 npoints = npoints,
                                 quad.method="grid"))


  # test with consistant masks.
  x <- sum(preds)
  preds2 <- mask(preds,x)
  ppmdata2a <- ppmData(npoints = npoints,
                       presences=presences,
                       window = preds2[[1]],
                       covariates = preds2,
                       quad.method="grid")
  testthat::expect_s3_class(ppmdata2a,"ppmData")

  # test with no npoints
  ppmdata3 <- ppmData(presences=presences,
                      window = preds[[1]],
                      quad.method="grid")
  testthat::expect_s3_class(ppmdata3,"ppmData")


  ppmdata4 <- ppmData(npoints=1000,
                      presences=presences,
                      window = preds[[1]],
                      unit="m",,
                      quad.method="grid")
  testthat::expect_s3_class(ppmdata4,"ppmData")

  ppmdata5 <- ppmData(npoints=1000,
                      presences=presences,
                      window = preds[[1]],
                      unit="km",
                      quad.method="grid")
  testthat::expect_s3_class(ppmdata5,"ppmData")

  ppmdata6 <- ppmData(npoints=1000,
                      presences=presences,
                      window = preds[[1]],
                      unit="ha",
                      quad.method="grid")
  testthat::expect_s3_class(ppmdata6,"ppmData")


})
