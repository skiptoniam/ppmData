#' @importFrom terra xyFromCell values ext freq extract ncell
quasiRandomQuad <- function(npoints,
                            window,
                            coord,
                            control){

  #generate a set of potential sites for quasi-random generation
  rast_coord <- terra::xyFromCell(window,ncell(window))
  rast_data <- terra::values(window)
  na_sites <- which(!complete.cases(rast_data))
  potential_sites <- rast_coord

  if(is.null(control$quasirandom.samples)) control$quasirandom.samples <- 10*npoints

  exty <- terra::ext(window)
  study.area <- matrix( as.vector( exty)[c( 1,2,2,1, 3,3,4,4)], nrow=4, ncol=2)

  #this gives many many samples, unless the study region is an odd shape, which it shouldn't be (as it is just an extent at this stage)
  samp <- quasiSampFromHyperRect(nSampsToConsider=control$quasirandom.samples,
                                 randStartType=2,
                                 designParams=list(dimension=2,study.area=study.area))

  ## setup the inclusion probs.
  sampValues <- terra::extract(window,samp[,1:2])
  NAsamps <- which(!complete.cases(sampValues))  ### these will give the

  ## generate the actuall quasi random points we want to generate for the model.
  Nsamps <- nrow(sampValues)
  inclusion_probs <- rep(1/Nsamps, Nsamps)
  inclusion_probs1 <- inclusion_probs/max(inclusion_probs)
  inclusion_probs1[NAsamps] <- 0

  #Spatially thin the sample which are NA sites in the raster window
  samp <- samp[samp[,3]<inclusion_probs1,1:2]
  if(nrow( samp) < npoints)
    stop("No set of background points found for this region.  Please increase the number of possible samples.")
  samp <- samp[1:npoints,]

  randpoints <- as.data.frame(samp[,1:2])
  colnames(randpoints ) <- coord


  return(randpoints)#, quasiDummy = randpointsDummy))
}

#' quasiSampFromHyperRect
#'@author Scott Foster
#'@param nSampsToConsider Number of samples to consider
#'@param randStartType Random start type, default is 2.
#'@param designParams A list of design parameters for the sampler.
#'@importFrom randtoolbox halton
#'@importFrom mgcv in.out
#'@description Copied from MBHdesign as it wasn't an exported function.

quasiSampFromHyperRect <- function (nSampsToConsider, randStartType = 2, designParams){
  samp <- randtoolbox::halton(nSampsToConsider * 2, dim = designParams$dimension + 1, init = TRUE)
  if (randStartType == 1)
    skips <- rep(sample(1:nSampsToConsider, size = 1, replace = TRUE), designParams$dimension + 1)
  if (randStartType == 2)
    skips <- sample(1:nSampsToConsider, size = designParams$dimension + 1, replace = TRUE)
  samp <- do.call("cbind", lapply(1:(designParams$dimension + 1), function(x) samp[skips[x] + 0:(nSampsToConsider - 1), x]))
  myRange <- apply(designParams$study.area, -1, range)
  for (ii in 1:designParams$dimension) samp[, ii] <- myRange[1,ii] + (myRange[2, ii] - myRange[1, ii]) * samp[, ii]
  if (designParams$dimension == 2) {
    tmp <- mgcv::in.out(designParams$study.area, samp[,1:designParams$dimension])
    samp <- samp[tmp, ]
  }
  return(samp)
}


quasiRandomWeights <- function(presences,
                               quadrature,
                               window,
                               coord,
                               mark.id,
                               unit,
                               crs = sf::st_crs("EPSG:4326"),
                               control){#, returnDirtess){

  ## Merge the presences and absences
  allpts.id <- rbind(presences, quadrature)
  allpts <- rbind(presences[,coord], quadrature[,coord])


  ## Get a bbox which has a small buffer
  bbox <- convert2pts(allpts)
  bbox[c(1,3)] <- bbox[c(1,3)] - 1e-10
  bbox[c(2,4)] <- bbox[c(2,4)] + 1e-10

  ## Tracking of site and species ids
  allpts$id <- seq_len(nrow(allpts))
  allpts$dataset <- allpts.id[,mark.id]

  ## Do the Dirichlet with me
  dirareas <- getDirichlet(allpts,
                           coord,
                           unit,
                           clippy = TRUE,
                           window,
                           crs = crs,
                           control = control)

  ## merge with all pts
  res <- merge(allpts, dirareas$data, by='id', all=TRUE, sort=TRUE)

  return( res)
}

#'@importFrom stats quantile
getDirichlet <- function(allpts,
                         coord,
                         unit,
                         clippy = TRUE,
                         window,
                         crs = sf::st_crs("EPSG:4326"),
                         control){#}, return_dirtess = TRUE ){

  ## set up the data.frame to catch the results.
  df <- data.frame(id = seq_len(nrow(allpts)), area = NA)

  ## Run the tessellation
  tess <- dirTess(as.matrix(allpts[,coord]))#, bbox = bbox)


  if(control$approx){
  quick_area <- sapply(tess$polygons$poly,function(x)x$area)[seq_len(nrow(allpts))]

  ## remove polygons with very large areas and make them approx equal edge polys
  approx_area <- quantile(quick_area,0.975)

  quick_area <- ifelse(quick_area>approx_area,approx_area,quick_area)

  df$area <- quick_area

  } else {
  # Take my dodgy polys and make them into sf ones.
  tess.out <- polygonise(x = tess,
                         window = window,
                         clippy = clippy,
                         crs = crs,
                         unit = unit)
  # tess.out$polygons.areas[seq_len(nrow(tess$coords))]

  df$area <- tess.out$polygons.areas[seq_len(nrow(tess$coords))]
  }


  res <- list()
  # res$polygons <- tess.out$polygons
  res$data <- df

  return(res)
}

