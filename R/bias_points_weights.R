# bias points weights for a single species while taking into account bias points (actual observations)
# weight <-  'total area of study extent'/'number of bias points (actual observations)'
# bias points should be all the observational data - eg. Fithian 2014.

# we need to know the area of the study region - if equal area cells,
# this is easy - we just sum areas !is.na*area of cell - other wise we use cellStats(area(study_areas)) - they should give the same answer.

# if using bias option - provide the the exisiting survey sites as bg pointspoints'
# library(raster)
# N <- 100
# known.sites <- as.data.frame(cbind(x1=runif(N, min=-10, max=10),x2=runif(N, min=-10, max=10)))
# study_area <- raster(nrows=100, ncols=100, xmn=-10, xmx=10,ymn=-10,ymx=10)
# study_area[]<-rnorm(10000)
# projection(study_area) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

library(dismo)
bradypus <- read.csv(paste0(system.file(package="dismo"),
                            "/ex/bradypus.csv"))[, -1]
covariates <- stack(list.files(paste0(system.file(package="dismo"), '/ex'),
                               pattern='grd', full.names=TRUE))

dat <- qrbp2(number_of_background_points = 10000,
             known_sites = bradypus,
             study_area = covariates[[1]],
             model_covariates = covariates)

plot(covariates[[1]])
points(bradypus, pch = 16, cex = 0.4)
points(dat[dat$presence==0,2:3], col='red',pch = 16, cex = 0.4)

library(gbm)
m <- gbm(presence ~ offset(log(weights)) +
           bio1 +
           bio5 +
           bio6 +
           bio7 +
           bio8,
         data = dat,
         n.trees = 10000,
         cv.folds = 5,
         distribution = 'poisson')

# optimal number of trees
trees <- gbm.perf(m, plot.it = FALSE, method = 'cv')


# predict to a raster
p <- predict(covariates, m,
             type = 'response',
             const = data.frame(weights = 1),
             n.trees = trees)

p_cell <- area(p)*1000 * p # in m^2
plot(p_cell)
points(bradypus, pch = 16, cex = 0.4)
cellStats(p_cell, sum)

# need to then create weight for observations and background/bias integration points
qrbp2 <- function(number_of_background_points = 2000, # number of back ground points to create.
                  known_sites = NULL,  #a set of coordinates,
                  study_area = NULL,#raster
                  model_covariates = NULL, # a set of covariates.
                  bias_points = NULL){
  eps <- sqrt(.Machine$double.eps)
  if(is.null(known_sites))stop('come on mate, you need some occurrence points.')

  #bias points will attempt to use other information from other species to inform sampling bias.
  #this could be a set of multi-species observations or a layer of observation bias/points.
  #hopefully in the future I can setup a fitian style offset based on collection method. e.g PO, PA, Abundance, ect.
  if(!is.null(bias_points)){
    message("Note - this isn't set up correctly...\nUsing 'bias_points', which should be a set of occurrences from other species observations.")
    eps <- sqrt(.Machine$double.eps)
    ppm_weights <-  c(rep(eps,known_sites),bias_points)
  }

  #if there is no study_area generate potential sites for halton function and a study area.
  if (is.null(study_area)) {
      message("no study area provided, generating a regular grid based on extent of 'known_sites'.")
      N <- 100
      potential_sites <- as.matrix(expand.grid(as.data.frame(matrix(rep(1:N,
                                                                        times = 2), ncol = 2)))/N -
                                     1/(2 * N))
      study_area <- as.matrix(expand.grid(as.data.frame(matrix(c(rep(0,
                                                                     2), rep(1, 2)), nrow = 2, byrow = TRUE))))
      colnames(potential_sites) <- colnames(study_area) <- paste0("2",
                                                                  1:2)
      study_area <- study_area[c(1, 3, 4, 2), ]
      study_area_ext <- c(study_area[,1])
  } else {
  potential_sites <- raster::rasterToPoints(study_area)[,-3]
  study_area_ext <- extent(study_area)[1:4]
  }
  N <- nrow(potential_sites)
  dimension <- dim(known_sites)[2]
  mult <- 20
  samp <- randtoolbox::halton(sample(1:10000, 1), dim = dimension +
                                1, init = TRUE)
  njump <- number_of_background_points * mult
  samp <- randtoolbox::halton(max(njump, 5000), dim = dimension +
                                1, init = FALSE)
  myRange <- t(matrix(study_area_ext,2,2,byrow = TRUE))
  for (ii in 1:2) samp[, ii] <- myRange[1, ii] + (myRange[2,ii] - myRange[1, ii]) * samp[, ii]
  sampIDs <- rep(NA, nrow(samp))
  sampIDs <- class::knn1(potential_sites,samp[, -(dimension + 1), drop = FALSE], 1:nrow(potential_sites))
  sampIDs <- sampIDs[!duplicated(sampIDs)]
  sampIDs <- sampIDs[1:number_of_background_points]
  samp <- as.data.frame(cbind(potential_sites[sampIDs, , drop = FALSE],
                              sampIDs))
  colnames(samp) <- c(colnames(potential_sites), "ID")
  # if using bias layer - we generate background points based on known sites
  # this requires information on other species where is the sampling bias in the landscape is...
  area_rast <- area(study_area)
  area_study <- mask(area_rast,study_area)
  region.size <- cellStats(area_study,sum, na.rm=TRUE)*1000 # size of the study extent in meters^2

if (!is.null(covariates)){
  bk_cvr <- extract(covariates,samp[,1:2],method='simple',na.rm=TRUE)
  bk_pts <- cbind(samp,bk_cvr)
  bk_pts2 <- na.omit(bk_pts)
  n <- nrow(bk_pts2)
  weight <- region.size/n
  bk_wts <- rep(weight,n)
  cat(nrow(bk_pts)-nrow(bk_pts2),"background points had NA data,\n now only",
    nrow(bk_pts2),"background points from the original\n",
    nrow(bk_pts),"background points")
  po_wts <- rep(eps,nrow(known_sites))
  po_cvr <- extract(covariates,known_sites,method='simple',na.rm=TRUE)
  po_pts <- cbind(known_sites,po_cvr)
}

names(po_pts)[1:2]<-c('x','y')
covar_data <- rbind(po_pts,bk_pts2[,-3])
dat <- data.frame(presence=c(rep(1,nrow(known_sites)),rep(0,nrow(bk_pts2))),
                  covar_data,
                  weights=c(po_wts,bk_wts))
return(dat)
}
