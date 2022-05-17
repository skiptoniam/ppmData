#' Helper function to convert df to polygons
#' @param df data.frame of polygon coordinates and ids.
#' @param keys id of each polygon
#' @param coords names of the coordinates
#' @param proj The CRS to project the polygons in
#' @importFrom sp Polygon Polygons SpatialPolygons CRS
#' @importFrom methods is
df2polygons <- function(df,keys,coords,proj) {

  ## Basic checks
  if(!is(df,"data.frame")) stop("df needs to be a data frame")
  if(!is(keys,"character")) stop("keys needs to be of class character")
  if(!is(coords,"character")) stop("coords needs to be of class character")
  if(!all(keys %in% names(df))) stop("All keys needs to be labels in data frame")
  if(!all(coords %in% names(df))) stop("All coordinate labels needs to be labels in data frame")
  if(!methods::is(proj,"CRS")) stop("proj needs to be of class CRS")

  ## dfun takes a data frame with coordinates for 1 polygon, and makes one POLYGON object from it
  ## with a UID from the polygon key
  dfun <- function(d) {
    sp::Polygons(list(sp::Polygon(d[coords])),
             as.character(d[keys]))
  }

  ## Now apply dfun to all polygons in data frame
  df_poly <- plyr::dlply(df,keys,dfun)

  ## Frorm a SpatialPolygons object from all the returned Polygons
  Sr <- sp::SpatialPolygons(df_poly,             # Polygons
                        1:length(df_poly),   # plotting order
                        proj4string=proj)    # CRS
}

cpp_to_df <- function(x){
  df <- as.data.frame(matrix(x,ncol=2,byrow=TRUE))
  colnames(df) <- c("x","y")
  return(df)
}
