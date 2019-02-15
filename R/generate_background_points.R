#' @name generate_background_points
#' @title create background points for study area or point data
#' @export
#' @param number_of_background_points number of back ground points to create - only applicable to quasirandom.
#' @param known_sites a matrix, dataframe or SpatialPoints* object giving the
#'   coordinates of the points to use in sampling (of size Nxdimension).
#'   Note: SpatialPoints* will only be suitiable for \code{dimension=2}.
#' @param study_area a Raster* object giving the area over which to generate background points.
#'  If ignored, a rectangle defining the extent of \code{known_sites} will be used.
#' @param model_covariates an optional Raster* object containing covariates for
#'   modelling the point process (best use a Raster* stack or Raster* brick)
#' @param resolution resolution setup grid for integration points (default is 1 deg).
#' @param multispecies_presences a community matrix of species x sites - which has the columns "x","y","sp1,...,"spN". Presences as denoted with a 1, while non-presences are denoted with an NA, this is to distingish between non-species specific observations and generated background points.
#' @param method the type of method that should be used to generate background points. The options are:
#' 'single_species_grid' generates a regular grid and weights for a single species. See Warton 2010 for details.
#' 'single_species_quasi' generates quasirandom background points based on input data for a single species. See Foster 2015 for details.
#' 'multispecies_grid' this method produces species-specific bacground points and weights multispecies presence-only using a quadrature grid #' method. It should return a list, that contains a model matrix of (n_unique_presence_only_sites + n_background_points) x (n_sps + const=1 #' + n_covariates), along with a  (n_unique_presence_only_sites + n_background_points) x (n_sps) matrix of weights.
#' 'multispecies_quasi' this method produces species-specific bacground points and weights multispecies presence-only using a quasirandom design method. It should return a list, that contains a model matrix of (n_unique_presence_only_sites + n_background_points) x (n_sps + const=1 #' + n_covariates), along with a  (n_unique_presence_only_sites + n_background_points) x (n_sps) matrix of weights.
#' 'multispecies_quasi'
#' @param interpolation_method either 'simple' or 'bilinear' and this determines the interpolation method for selecting covariate data. 'simple' is nearest neighbour, 'bilinear' is bilinear interpolation.
#' @param coords is the name of site coordinates. The default is c('x','y').
#'
#' @importFrom mgcv in.out
#' @importFrom raster extract

generate_background_points <- function(number_of_background_points = 10000,
                                       known_sites = NULL,  #a set of coordinates,
                                       study_area = NULL,  #raster
                                       model_covariates = NULL, # a set of covariates.
                                       resolution = 1, #resolution to setup quadriature points in lat/lon
                                       multispecies_presences = NULL, #a community matrix of speciesXsites.
                                       method = c('single_species_grid','single_species_quasi',
                                                  'multispecies_grid','multispecies_quasi'),
                                       interpolation_method='bilinear',
                                       coords = c("x","y")){

  #this function was built to match the dimensions of coordinates and other dimensions
  known_sites <- coords_match_dim(known_sites,dim(known_sites)[2])
  # eps <- sqrt(.Machine$double.eps)

  if(is.null(known_sites))stop('come on mate, you need some occurrence points.')

  #multispecies points will attempt to use other information from other species to inform sampling bias.
  #this could be a set of multi-species observations or a layer of observation bias/points.
  #hopefully in the future I can setup a fitian style offset based on collection method. e.g PO, PA, Abundance, ect.
  if(method == 'multispecies'){
    if(is.null(multispecies_presences))
    stop("Using 'multispecies_presences' or method='multispecies' requires a matrix of sites x presences. \n
         Note NA in this matrix indicate a non-presence record for that species.")
    if(!is.null(multispecies_presences))
      message("'known_sites' must be all the unique coordinates for known species' presences")
    if(nrow(known_sites)!=nrow(multispecies_presences))
      stop("nrows of 'known_sites' must equal nrows of 'multispecies_presences'")
  }

  #if there is no study_area generate potential sites for halton function and a study area.
  if (is.null(study_area)) {
      message("no study area provided, a raster based on the extent of 'known_sites'\n with a default resolution of 1 deg will be returned.")
      study_area <- default_study_area(known_sites)
  }

  # create background points based on method.
  background_sites <- switch(method,
                single_species_grid = grid_method(resolution, study_area),
                single_species_quasi = quasirandom_covariates_method(number_of_background_points,  study_area, model_covariates),
                multispecies_grid = grid_method(resolution,  study_area),
                multispecies_quasi = quasirandom_covariates_method(number_of_background_points, study_area, model_covariates))

  site_weights <- switch(method,
                         single_species_grid = get_weights(known_sites,background_sites[,1:2],study_area,coords),
                         single_species_quasi = get_weights(known_sites,background_sites[,1:2],study_area,coords),
                         multispecies_grid = get_weights(multispecies_presences,background_sites[,1:2],study_area,coords),
                         multispecies_quasi = get_weights(multispecies_presences,background_sites[,1:2],study_area,coords))

  #id have provided covaraites use this to extract out environmental data from rasterstack
   if (!is.null(model_covariates)){
     if(!inherits(study_area, c('RasterLayer','RasterStack','RasterBrick')))
       stop("model_covariates must be a raster, rasterstack or rasterbrick of model_covariates desired for modelling")

       # extract covariate data for background points
       if(method=='multispecies_grid'|method=='multispecies_quasi') {
         covars <- extract(model_covariates,site_weights$multispecies_presence[,coords],method=interpolation_method,na.rm=TRUE)
         NAsites <- which(!complete.cases(covars))
         if(length(NAsites)>0){
           print(paste0('A total of ',length(NAsites),' sites where removed from the background data, because they contained NAs, check environmental data and species sites data overlap.'))
            covars <- covars[-NAsites,,drop=FALSE]
            dat <- list()
            dat$model_matrix <- data.frame(site_weights$multispecies_presence[-NAsites,-1:-2],const=1,site_weights[-NAsites,coords],covars)
            dat$species_weights <- site_weights$multispecies_weights[-NAsites,-1]
         } else {
           dat <- list()
           dat$model_matrix <- data.frame(site_weights$multispecies_presence[,-1:-2],const=1,site_weights,covars)
           dat$species_weights <- site_weights$multispecies_weights[,-1]
         }
       } else {
         # print(head(site_weights[,coords]))
         covars <- extract(model_covariates,site_weights[,coords],method=interpolation_method,na.rm=TRUE)
         NAsites <- which(!complete.cases(covars))
         if(length(NAsites)>0){
           print(paste0('A total of ',length(NAsites),' sites where removed from the background data, because they contained NAs, check environmental data and species sites data overlap.'))
           covars <- covars[-NAsites,,drop=FALSE]
           dat <- data.frame(presence=c(rep(1,nrow(known_sites)),rep(0,nrow(background_sites)))[-NAsites],
                             site_weights[-NAsites,coords],
                             covars,
                             weights=site_weights$weights[-NAsites])#replace this with a function 'get_weights'
         } else {
           dat <- data.frame(presence=c(rep(1,nrow(known_sites)),rep(0,nrow(background_sites))),
                             site_weights,
                             covars,
                             weights=site_weights$weights)
        }
      }
    } else {
    if(method=='multispecies_grid'|method=='multispecies_quasi'){
       dat <- list()
       dat$model_matrix <- data.frame(site_weights$multispecies_presence[,-1:-2],
                                       x=c(known_sites$x,background_sites$x),
                                       y=c(known_sites$y,background_sites$y))
       dat$species_weights <- site_weights$multispecies_weights[,-1]
     } else {
     # create an entire dataset
     dat <- data.frame(presence=c(rep(1,nrow(known_sites)),rep(0,nrow(background_sites))),
                       x=site_weights$x,y=site_weights$y,
                       weights=site_weights$weights)#replace this with a function 'get_weights'
   }
  }

  if(method!='multispecies_grid'|method!='multispecies_quasi') dat <- rm_na_pts(dat)
  return(dat)
}

grid_method <- function(resolution=1, study_area){

  if(!inherits(study_area, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'grid' method currently only works a raster input as a 'study_area'")

  if(inherits(study_area, c('RasterLayer','RasterStack','RasterBrick'))){

    #set up the dissaggreation or aggregate
    fct <- (res(study_area)/resolution)[1]

    #if fct is >= 1 dissaggregate, else aggregate
    if(fct>=1) dd <- disaggregate(study_area, fct, na.rm=TRUE)
    else dd <- aggregate(study_area, 1/fct, na.rm=TRUE)

    #create a dataframe of coordinates w/ area
    grid <- as.data.frame(rasterToPoints(dd)[,-3])
    colnames(grid) <- c('x','y')
  }
  return(grid)
}


# still working on this method.
quasirandom_covariates_method <- function(number_of_background_points, study_area, covariates=NULL){

  if(!inherits(covariates, c('RasterLayer','RasterStack','RasterBrick')))
    stop("'quasirandom_covariates' method currently only works a raster input as a 'covariates'")

  #generate a set of potential sites for quasirandom generation
  if(!is.null(covariates)){
    rast_coords <- raster::xyFromCell(covariates,1:ncell(covariates))
    rast_data <- raster::values(covariates)
    na_sites <- which(is.na(rast_data[,1]))
    covariates_ext <- extent(covariates)[1:4]
  } else {
    rast_coords <- raster::xyFromCell(study_area,1:ncell(study_area))
    rast_data <- raster::values(study_area)
    na_sites <- which(is.na(rast_data))
    covariates_ext <- extent(study_area)[1:4]
    }

  potential_sites <- cbind(rast_coords,rast_data)
  potential_sites[na_sites,-1:-2] <- 0

  #if number_of_background_points>nrow(potential_sites) recursively create a finer grid.
  if(number_of_background_points>nrow(potential_sites)){
    stop('more background points than cells avaliable')
  }

  # covariates_ext <- extent(covariates)[1:4]

  #setup the dimensions need to halton random numbers - let's keep it at 2 for now, could expand to alternative dimension in the future
  dimension <- dim(potential_sites)[2]

  N <- nrow(potential_sites)
  inclusion_probs <- rep(1/N, N)
  inclusion_probs[na_sites] <- 0
  inclusion_probs1 <- inclusion_probs/max(inclusion_probs)
  mult <- 10
  samp <- randtoolbox::halton(sample(1:10000, 1), dim = dimension +
                                1, init = TRUE)
  njump <- number_of_background_points * mult
  samp <- randtoolbox::halton(max(njump, 5000), dim = dimension +
                                1, init = FALSE)
  myRange <- apply(potential_sites[-na_sites,], -1, range)
  for (ii in 1:dimension) samp[, ii] <- myRange[1, ii] + (myRange[2,
                                                                  ii] - myRange[1, ii]) * samp[, ii]
  if (dimension == 2) {
    tmp <- mgcv::in.out(study_area, samp[, 1:dimension])
    samp <- samp[tmp, ]
  }
  sampIDs <- rep(NA, nrow(samp))
  kount <- 0
  flag <- TRUE
  while (flag & (kount < nrow(samp))) {
    if (kount == 0)
      message("Number of samples considered (number of samples found): ",
              njump, "(0) ", sep = "")
    else message(kount + njump, "(", length(sampIDs.2),
                 ") ", sep = "")
    sampIDs[kount + 1:min(njump, nrow(samp) - kount)] <- class::knn1(potential_sites,
                                                                     samp[kount + 1:min(njump, nrow(samp) - kount), -(dimension + 1), drop = FALSE],
                                                                     1:nrow(potential_sites))
    sampIDs.2 <- which(samp[1:(kount + min(njump, nrow(samp) - kount)), dimension + 1] < inclusion_probs1[sampIDs[1:(kount +
                                                                                                     min(njump, nrow(samp) - kount))]])
    if (length(sampIDs.2) >= number_of_background_points) {
      sampIDs <- sampIDs[sampIDs.2][1:number_of_background_points]
      flag <- FALSE
    }
    kount <- kount + njump
  }
  # message("Finished\n")
  if (kount > nrow(samp))
    stop("Failed to find a design. It is likely that the inclusion probabilities are very low and uneven. Please try again OR make inclusion probabilities more even")
  samp <- as.data.frame(cbind(potential_sites[sampIDs, , drop = FALSE],
                              inclusion_probs[sampIDs], sampIDs))
  colnames(samp) <- c(colnames(potential_sites), "inclusion.probabilities",
                      "ID")

  return(samp)
}

coords_match_dim <- function (known.sites,dimension){

  # check object classes
  expectClasses(known.sites,
                c('matrix',
                  'data.frame',
                  'SpatialPoints',
                  'SpatialPointsDataFrame'),
                name = 'known.sites')

  if (is.matrix(known.sites) | is.data.frame(known.sites)) {

    # for matrices/dataframes, make sure there are only two columns
    if (ncol(known.sites) != dimension ) {
      stop (sprintf('known.sites should match the number of dimensions used in quasi-random sampling.
                    The object passed had %i columns, while the sampling dimensions are %i',NCOL(known.sites),dimension))
    }

    # otherwise, coerce into a data.frame and rename the columns
    df <- data.frame(known.sites)

    } else {
      # otherwise, for SpatialPoints* objects, just grab the coordinates
      df <- data.frame(known.sites@coords)
    }

  # set column names
  if (ncol(known.sites)>2) {
    if(!is.null(colnames(known.sites))) colnames(df) <- c("x","y",colnames(known.sites)[-1:-2])
    colnames(df) <- c("x","y",paste0("var",seq_len(ncol(known.sites)-2)))
  } else {
    colnames(df) <- c("x","y")
  }
  return (df)
}

default_study_area <- function (known_sites) {
  # get limits
  xlim <- range(known_sites$x)
  ylim <- range(known_sites$y)

  # add on 5%
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.1
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.1

  xlim[1] <- floor(xlim[1])
  xlim[2] <- ceiling(xlim[2])
  ylim[1] <- floor(ylim[1])
  ylim[2] <- ceiling(ylim[2])

  e <- extent(c(xlim,ylim))
  sa <- raster(e,res=1, crs="+proj=longlat +datum=WGS84")
  sa[]<- 1:ncell(sa)
  return (sa)
}

# remove NA values from a pts object
rm_na_pts <- function (pts) {

  # remove any points that are NA and issue a warning
  if (any(is.na(pts))) {

    # copy the original data
    pts_old <- pts
    # remove NAs
    pts <- na.omit(pts)

    # report which ones were removed
    rm <- as.vector(attributes(pts)$na.action)
    which_rm <- pts_old$points[rm]

    if (any(which_rm == 0)) {
      warning(sprintf('Removed %i integration points for which covariate values could not be assigned.
                      This may affect the results of the model, please check alignment between model_covariates and area.',sum(which_rm == 0)))
    }

    if (any(which_rm == 1)) {
      warning(sprintf('Removed %i observed points for which covariate values could not be assigned.
                      This may affect the results of the model, please check alignment between model_covariates and coords.',sum(which_rm == 0)))
    }

  }

  return (pts)

}

estimate_area <- function(study_area,site_coords){
  if(raster::isLonLat(study_area)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- raster::area(study_area)
    area_study <- raster::mask(area_rast,study_area)
    total_area <- cellStats(area_study,sum,na.rm=TRUE)
    wts <- total_area/raster::extract(area_study,site_coords,na.rm=TRUE)
    } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    cell_area <- raster::res(study_area)[1]*raster::res(study_area)[2]
    n_cell <- length(study_area[!is.na(study_area[])])
    wts <- rep((n_cell*cell_area)/nrow(site_coords),nrow(site_coords))/1000
    }
   return(wts)
}

#estimate the total area of all cells in extent in km^2
estimate_area_region <- function(study_area){
  if(raster::isLonLat(study_area)){
    #calculate area based on area function
    #convert kms to ms
    area_rast <- area(study_area)
    area_study <- mask(area_rast,study_area)
    area_of_region <- cellStats(area_study,sum, na.rm=TRUE)
  } else {
    # calculate area based on equal area cell resolution
    # mode equal area should be in meters
    area_of_region <- (ncell(study_area[!is.na(study_area)])  * xres(study_area) * yres(study_area))/1000
  }
  return(area_of_region)
}

get_weights <- function(known_sites,background_sites,study_area,coords){

  # if single species known sites will be just coords
  # if multispeices matrix will be x, y, sp1, sp2 matrix (min four columns)
  # weights are the realtive area that point represents within the domain.


  if(ncol(known_sites)>2) {

    known_sites <- as.matrix(known_sites)
    if(!check_pres_bk_names(known_sites[,1:2],background_sites))colnames(background_sites) <- colnames(known_sites[,1:2])


    # how many species are there?
    n_sp <- ncol(known_sites) # this is actually this number - 2
    # get the species specific presences
    sp_presence_coords <- lapply(3:n_sp,function(x)known_sites[!is.na(known_sites[,x]),coords])
    # what
    sp_xy <- lapply(sp_presence_coords,function(x)rbind(x,background_sites))
    # sp specific areas
    sp_areas <- lapply(sp_xy,function(x)qrbp:::estimate_area(study_area,x))
    # sp cell ids.
    sp_cell_ids <- lapply(sp_xy,function(x)raster::cellFromXY(study_area,x))
    # sp specific counts of points
    sp_count_pts <- lapply(sp_cell_ids,table)
    # sp specific weights
    sp_weights <- lapply(1:c(n_sp-2),function(x)sp_areas[[x]]/as.numeric(sp_count_pts[[x]][match(sp_cell_ids[[x]],
                                                                                                 names(sp_count_pts[[x]]))]))

    ## now I need to put this all back in to a matix that matches the xy.
    sp_weights_mat <- matrix(NA,nrow(known_sites)+nrow(background_sites),n_sp-1)
    sp_pres_mat <- matrix(NA,nrow(known_sites)+nrow(background_sites),n_sp)
    all_sites_ids <- raster::cellFromXY(study_area,rbind(known_sites[,coords],background_sites[,coords]))
    sp_weights_mat[,1]<- all_sites_ids

    for (ii in 1:c(n_sp-2)){
      sp_pres_mat[,c(ii+2)] <- c(known_sites[,c(ii+2)],rep(0,nrow(background_sites)))
      sp_weights_mat[which(!is.na(sp_pres_mat[,ii+2])),c(ii+1)]<-sp_weights[[ii]]
    }
    sp_pres_mat[,c(1,2)] <- as.matrix(rbind(known_sites[,coords],background_sites[,coords]))
    colnames(sp_pres_mat) <- colnames(known_sites)
    colnames(sp_weights_mat) <- c("cell_id",colnames(known_sites)[-1:-2])
    weights <- list()
    weights$multispecies_weights <- sp_weights_mat
    weights$multispecies_presence <- sp_pres_mat

  } else {
    xy <- rbind(known_sites,background_sites)
    areas <- estimate_area(study_area = study_area, site_coords = xy)
    cell_id <- raster::cellFromXY(study_area, xy)
    count_pts <- table(cell_id)
    weights <- areas/as.numeric(count_pts[match(cell_id,
                                                names(count_pts))])
    weights <- data.frame(x=xy[,1],y=xy[,2],weights=weights)
  }
  return(weights)
}


check_pres_bk_names <- function(pres_sites,back_sites){

  return(all(colnames(pres_sites)==colnames(back_sites)))



}
