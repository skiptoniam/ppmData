#' @title deltri
#' @description  Computes the Delaunay triangulation for large data sets. The
#' main use of this function is to quickly calculate the Delaunay triangulation
#' so we can calculate the duel-graph and Dirichlet tessellation to use in point
#' process modelling, where the area of each Dirichlet tessellation is used as a
#' areal weight in a statistical model.
#' @param coords The coordinates of the points to triangulate.
#' @export
#' @examples
#' coordin <- matrix(runif(2000),ncol=2)
#' tri <- deltri(coordin)
# plot(tri)

deltri <- function(coords){

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

  class(res) <- "deltri"
  return(res)

}

#' @export
print.deltri <- function(x, ...){
  message(paste0("Delaunay triangulation of ", nrow(x$coordinates)," sites."))
  message(paste0("There are a total of ", nrow(x$tri.mat)," unique triangles for this solution."))
}

# x <- tri
# next_half_edge <- function(e) { #edge index starting at zero.
#   tmp <- ifelse(e %% 3 == 2, e - 2, e + 1)
#   return(tmp)
# }
#
#
# tri_edges <- function(x, ...) {
#
#   for (e in 1:length(x$triangles)) {
#     if (e > x$halfedges[e]) {
#       p <- x$coordinates[x$triangles[e],]
#       q <- x$coordinates[x$triangles[next_half_edge(e)]];
#
#      }
#   }
# }

#' @export
plot.deltri <- function(x, add = FALSE, axis = FALSE, boxed = FALSE, ...){

  p <- x$coordinates
  # if (!is.matrix(p)) {
  #   p = cbind(p, p2)
  # }
  xlim = range(p[, 1])
  ylim = range(p[, 2])
  if (!add) {
    plot.new()
    plot.window(xlim, ylim, ...)
  }
  if (boxed) {
    box()
  }
  if (axis) {
    axis(1)
    axis(2)
  }
  m = rbind(tri[, -1], tri[, -2], tri[, -3])
  segments(p[m[, 1], 1], p[m[, 1], 2], p[m[, 2], 1], p[m[,
                                                         2], 2], ...)
  return(invisible(list(tri = tri, p = p)))
}




