#' @title delTri
#' @description  Computes the Delaunay triangulation for large data sets. The
#' main use of this function is to quickly calculate the Delaunay triangulation
#' so we can calculate the duel-graph and Dirichlet tessellation to use in point
#' process modelling, where the area of each Dirichlet tessellation is used as a
#' areal weight in a statistical model.
#' @param coords The coordinates of the points to triangulate. Needs to be a two column matrix or data.frame with longitude and latitude as columsn 1 and 2 respectively.
#' @export
#' @examples
#' coordin <- matrix(runif(2000),ncol=2)
#' tri <- delTri(coordin)
# plot(tri)

"delTri" <- function(coords){

  ## Do some stuff with the coordinates
  if(is.data.frame(coords)) coords.in <- as.numeric(coords)
  if(is.matrix(coords)) coords.in <- as.numeric(coords)

  dt.out <- deltri_cpp(coords.in)

  res <- list()
  res$coords_cpp <- dt.out$coords
  res$coordinates <- matrix(dt.out$coords, ncol=2)
  res$triangles <- dt.out$triangle  # make the indexing R friendly.
  res$tri.mat <- matrix(dt.out$triangle, ncol=3)
  res$halfedges <- dt.out$halfedges
  res$convex.hull <- dt.out$convexhull
  res$convex.hull.area <- dt.out$convexhull.area

  class(res) <- "delTri"
  return(res)

}

#' Printing function for delTri
#' @param x A delTri object
#' @param \\dots Additional print calls.
#' @export
"print.delTri" <- function(x, ...){
  message(paste0("Delaunay triangulation of ", nrow(x$coordinates)," sites."))
  message(paste0("There are a total of ", nrow(x$tri.mat)," unique triangles for this solution."))
}




