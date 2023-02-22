#' @title dirTess
#' @description Computes the Dirichlet tessellation using the dual-graph of the
#' Delaunay triangulation. The Delaunay triangulation is calculated internally
#' for this function using a sweep algorithm.
#' @param coords The coordinates of the points to triangulate.
#' @param bbox Optional - doesn't actually do anything anymore :)
#' @export
#' @examples
#' coords <- matrix(runif(2000),ncol=2)
#' tess <- dirTess(coords)

"dirTess" <- function(coords, bbox=NULL){

  if(!is.matrix(coords)) stop('"coords" must be a matrix of spatial points for tessellation.')
  if(ncol(coords)!=2) stop('"coords" must be a two column matrix with longitude "x" and latitude "y".')

  # ## Add in a function which checks for bounding box.
  if(is.null(bbox)){
    bbox <- get_bbox(coords,buffer=0)
  }

  ## check bbox is right, vectorise.
  bbox <- check_bbox(bbox, coords)

  ## Create dummy points on the edges of bounding box.
  # bbox <- get_bbox(coords,buffer=0.05)
  bbox_dummy <- get_bbox(coords, buffer = 0.05)
  dummycoords <- get_bbox_dummy_coords(coords,bbox_dummy,n=9)
  coords.in <- rbind(coords,dummycoords)

  ## run dirichlet tessellation
  res <- dirtess_cpp(coords = as.numeric(t(coords.in)))#, as.integer(ncoords))

  # res <- list()
  res$coords <- coords
  res$ncoords <- nrow(coords)
  res$dummy.coords <- dummycoords
  res$bbox <- bbox

  class(res) <- "dirTess"

  return(res)

}

"check_bbox" <- function(bbox, coords, ...){

  if(length(bbox)!=4)
    stop("bbox should be of length four in xmin, xmax, ymin, ymax format")
  if(bbox[1]>=bbox[2])
    stop("xmin >= xmax check inputs")
  if(bbox[3]>=bbox[4])
    stop("ymin >= ymax check inputs")

  xmn <- min(coords[,1])
  xmx <- max(coords[,1])
  ymn <- min(coords[,2])
  ymx <- max(coords[,2])

  if(bbox[1]>xmn)
    stop("bbox xmin > than coordinate xmin  ")
  if(bbox[2]<xmx)
    stop("bbox xmax < that coordinate xmax  ")
  if(bbox[3]>ymn)
    stop("bbox ymin > than coordinate ymin  ")
  if(bbox[4]<ymx)
    stop("bbox ymax < than coordinate ymax  ")

  return(bbox)

}

"get_bbox" <- function(coords,buffer=0.05,...){

  ## coord ranges
  xmn <- min(coords[,1])
  xmx <- max(coords[,1])
  ymn <- min(coords[,2])
  ymx <- max(coords[,2])

  ## Add buffer to bbox
  xr <- diff(range(coords[,1]))
  yr <- diff(range(coords[,2]))
  xb <- xr*buffer
  yb <- yr*buffer
  xmn <- xmn - xb
  xmx <- xmx + xb
  ymn <- ymn - yb
  ymx <- ymx + yb

  bbox <- c('xmin'=xmn,'xmax'=xmx,'ymin'=ymn,'ymax'=ymx)
  return(bbox)

}

"get_bbox_dummy_coords" <- function( coords, bbox, n = NULL,...){

  if(is.null(n))
    n <- ceiling(sqrt(nrow(coords))) ## approx spacing of points along edge

  ## get seq of points along bbox edges
  xpts <- seq(bbox[1],bbox[2],len=n)
  ypts <- seq(bbox[3],bbox[4],len=n)

  dpx <- expand.grid(xpts,bbox[3:4])
  dpy <- expand.grid(bbox[1:2],ypts)

  dp <- rbind(dpx,dpy)

  dp <- dp[!duplicated(dp),]
  dp <- dp[order(dp[,1]),]

  dupdp <- which(as.matrix(dp)%in%coords)
  if(length(dupdp)>0){
    dp[-dupdp,]
  }

  return(as.matrix(dp))

}


#'plot a Dirichlet tessellation
#'@param x dirTess object
#'@param \\dots Additional plotting arguments
#'@importFrom grDevices gray
#'@importFrom graphics polygon points
#'@export
"plot.dirTess" <- function(x, ...){

  xlim <- range(x$coords[,1]) + abs(diff(range(x$coords[,1])))*0.1*c(-1,1)
  ylim <- range(x$coords[,2]) + abs(diff(range(x$coords[,2])))*0.1*c(-1,1)
  plot(x$coords, type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=ylim,...)
  for(i in 1:x$ncoords){
    graphics::polygon(x$polygons$poly[[i]]$poly.coords,border = grDevices::gray(0.3,0.5),
            col = sample(grDevices::hcl.colors(x$ncoords,alpha = 0.5),1))
  }
  graphics::points(x$coords, pch=16, cex=.8)
}


#' Convert the Dirichlet tessellation to polygons and clip if so desired.
#' @rdname polygonise
#' @export polygonise
#' @param x Dirichlet tessellation object
#' @param window polygon A polygon to clip the dirichlet tessellation - intersection between boundary polygons and polygon
#' @param clippy Logical Clip to the bounding box or polgon if polyclip = TRUE.
#' @param crs CRS Projection of coordinates default is "EPSG:4326" (lon/lat wgs84)
#' @param unit Character The type of area to return. The default is "geo" and
#' returns the area based on the euclidean distance between geographic
#' coordinates. This will default to the values of the raster and presence
#' coordinate system. Alternatively, meters squared "m", kilometers squared "km", or hectares "ha" can be used.
#' @param \\dots Additional arguments for a polygonise function
#' @importFrom sf st_crs st_as_sfc st_bbox st_intersection st_area
#' @examples
#' library(sf)
#' coords <- matrix(runif(2000),ncol=2)
#' tess <- dirTess(coords)
#' tess.polys <- polygonise(tess)
#' plot(st_geometry(tess.polys$polygons), col=hcl.colors(100))
#' points(coords,col='tomato',pch=16, cex=.3)

"polygonise" <- function (x,  window=NULL, clippy=TRUE, crs = sf::st_crs("EPSG:4326"), unit, ...){
  UseMethod("polygonise", x)
}

#' @export
"polygonise.dirTess" <- function(x, window=NULL, clippy=TRUE, crs = sf::st_crs("EPSG:4326"), unit="geo", ...){

  # convert the vector of coordinates per polygon into polygons
  bbox <- x$bbox
  res <- list()
  tes.polys <- dirichletPolygons(x = x, crs = crs)
  # res$polygons <- tes.polys
  if(clippy){
    clip.poly <- sf::st_as_sfc(st_bbox(c(bbox[1], bbox[2], bbox[4], bbox[3]), crs = crs))
    # clip.poly <- rgeos::bbox2SP(bbox[4],bbox[3],bbox[1],bbox[2],proj4string = proj)
    if(!is.null(window)){
      if(classTrue(window,"SpatRaster"))
        clip.poly <- rast_to_sf(window)
      else if(classTrue(window,"sf"))
        clip.poly <- window
      else
        stop("window needs to be a sf polygons or a terra SpatRaster object.")
    }

    ## get the intersection between the two polygons
    tes.clip <- suppressWarnings(sf::st_intersection(tes.polys,clip.poly))

    polygon.clipped.areas <-  switch(unit,
                                     geo = as.numeric(geoArea(tes.clip)),
                                     m = as.numeric(sf::st_area(tes.clip)),
                                     km = as.numeric(sf::st_area(tes.clip))/1e6,
                                     ha = as.numeric(sf::st_area(tes.clip))/1e4)
    # polygon.clipped.areas <- suppressWarnings(st_area(tes.clip))
    res$polygons <- tes.clip
    res$polygons.areas <- polygon.clipped.areas
  } else {
    polygon.clipped.areas <-  switch(unit,
                                     geo = as.numeric(geoArea(tes.polys)),
                                     m = as.numeric(sf::st_area(tes.polys)),
                                     km = as.numeric(sf::st_area(tes.polys))/1e6,
                                     ha = as.numeric(sf::st_area(tes.polys))/1e4)
    res$polygons.areas <- polygon.clipped.areas
  }

  return(res)

}

"geoArea" <- function(x){

    ## data.frame of coordinates
    coords <- suppressWarnings(sf::st_coordinates(sf::st_cast(x,"POLYGON")))

    ## ids of polygons
    ids <- seq_along(unique(coords[,"L2"]))
    ea <- sapply(ids,function(ii)dirtess_poly_area(coords[coords[,"L2"]==ii,"X"],coords[coords[,"L2"]==ii,"Y"]))

    return(ea)

}

"dirichletPolygons" <- function(x, crs = sf::st_crs("EPSG:4326"),...){

  if(!isa(x,"dirTess"))
    stop("Must be a dirichlet tessellation object from 'dirTess'")

  ## extract the polygon lists from cpp object
  poly.list <- lapply(x$polygons$poly,function(i)i$poly.coords)

  ## close the polygon for sf
  poly.list2 <- lapply(poly.list,function(x)rbind(x,x[1,]))

  ## keep track of bogus polygons
  dodgy_polys <- which(sapply(poly.list2,nrow)<4)
  if(length(dodgy_polys)==0)
    dodgy_polys <- -seq_len(length(poly.list2))

  ## Generate a list of polygons
  sfpoly <- lapply(poly.list2[-dodgy_polys],function(x)sf::st_polygon(list(x)))

  ## Convert them to a single sf object
  sfc_poly <- sf::st_sfc(sfpoly)
  sf_polys <- sf::st_sf(sfc_poly)

  ## add an id
  poly.ids <- sapply(x$polygons$poly[-dodgy_polys],function(i)i$poly.id)
  sf_polys$id <- poly.ids
  sf_polys <- sf::st_set_crs(sf_polys,crs)

  return(sf_polys)

}
