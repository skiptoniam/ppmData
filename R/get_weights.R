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
