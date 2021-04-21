context('ppmData single species')
library(ppmData)


testthat::test_that('species mix ippm', {

library(ppmData)
library(raster)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
preds <- raster::stack(lst)
presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
coords <- xyFromCell(preds,1:ncell(preds))
X <- preds[[1]]
Y <- preds[[1]]
X[] <- coords[,1]
Y[] <- coords[,2]
names(X) <- "X"
names(Y) <- "Y"
preds2 <- stack(X,Y,preds)

## some test a-roos!

#
ppmdata <- ppmData(npoints = 10000,
                   presences=presences)#,
                   # window = preds2[[1]])

})
