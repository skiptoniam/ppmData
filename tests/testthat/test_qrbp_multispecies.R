context('Test ppmData for multiple species')
library(qrbp)

testthat::test_that('Test background generation for null species - i.e. no species presences', {

  set.seed(42)
  library(sdm)
  library(raster)

  file <- system.file("external/species.shp", package="sdm") #
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
  resolution <- 16000
  covariates <- NULL

  bkgrid<- ppmData(npoints = npoints,
                   resolution = resolution,
                   presences = presences,
                   window = window,
                   covariates = covariates,
                   method = method,
                   interpolation = interpolation,
                   control=list(control$multispeciesFormat))

})
