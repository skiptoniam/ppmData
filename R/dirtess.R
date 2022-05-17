#' @title dirtess
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

  ## Add in a function which checks for bounding box.
  if(is.null(bbox)){
    bbox <- get_bbox(coords,buffer=0)
  }

  ## check bbox is right, vectorise.
  bbox <- check_bbox(bbox, coords)

  ## Create dummy points on the edges of bounding box.
  # bbox <- get_bbox(coords,buffer=0.05)
  bbox_dummy <- get_bbox(coords,buffer = 0.05)
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
#'@export
"plot.dirTess" <- function(x, ...){

  xlim <- range(x$coords[,1]) + abs(diff(range(x$coords[,1])))*0.1*c(-1,1)
  ylim <- range(x$coords[,2]) + abs(diff(range(x$coords[,2])))*0.1*c(-1,1)
  plot(x$coords, type="n",axes=FALSE,xlab="",ylab="",xlim=xlim,ylim=ylim,...)
  for(i in 1:x$ncoords){
    polygon(x$polygons$poly[[i]]$poly.coords,border = gray(0.3,0.5),
            col = sample(hcl.colors(x$ncoords,alpha = 0.5),1))
  }
  points(x$coords, pch=16, cex=.8)
}


#' Convert the Dirichlet tessellation to polygons and clip if so desired.
#' @rdname polygonise
#' @export polygonise
#' @param x Dirichlet tessellation object
#' @param clip Logical Clip to the bounding box or polgon if polyclip = TRUE.
#' @param polyclip polygon A polygon to clip the dirichlet tessellation - intersection between boundary polygons and polygon
#' @param proj CRS Projection of coordinates default is "EPSG:4326" (lon/lat wgs84)
#' @param \\dots Additional arguments for a polygonise function
#' @importFrom sp CRS
#' @importFrom rgeos gArea gIntersection
#' @examples
#' library(sp)
#' coords <- matrix(runif(2000),ncol=2)
#' tess <- dirTess(coords)
#' tess.polys <- polygonise(tess)
#' plot(tess.polys$polygons, col=hcl.colors(100))

"polygonise" <- function (x, clip=TRUE, polyclip = NULL, proj = sp::CRS("EPSG:4326"), ...){
  UseMethod("polygonise", x)
}

#' @export
"polygonise.dirTess" <- function(x, clip=TRUE, polyclip = NULL, proj = sp::CRS("EPSG:4326"), ...){

  # convert the vector of coordinates per polygon into polygons
  bbox <- x$bbox
  res <- list()
  tes.polys <- dirichlet.polygon(x = x, proj = proj)
  res$polygons <- tes.polys
  if(clip){
    clip.poly <- rgeos::bbox2SP(bbox[4],bbox[3],bbox[1],bbox[2],proj4string = proj)
    if(!is.null(polyclip)){
      clip.poly <- polyclip
    }
    tes.clip <- suppressWarnings(rgeos::gIntersection(tes.polys,clip.poly,byid = TRUE))
    polygon.clipped.areas <- suppressWarnings(rgeos::gArea(tes.clip, byid=TRUE))
    res$polygons <- tes.clip
    res$polygons.areas <- polygon.clipped.areas
  } else {
    res$polygons.areas <-suppressWarnings(rgeos::gArea(tes.polys, byid=TRUE))
  }

  return(res)

}


"dirichlet.polygon" <- function(x, proj = sp::CRS("EPSG:4326"),...){

  if(class(x)!="dirTess")
    stop("Must be a dirichlet tessellation object from 'dirTess'")

  poly.list <- lapply(x$polygons$poly,function(i)i$poly.coords)
  poly.size <- sapply(poly.list,nrow)

  id <- sapply(x$polygons$poly,function(i)i$poly.id)
  ids <- rep(id,poly.size)

  poly.df <- do.call(rbind,poly.list)
  poly.df <- data.frame(x=poly.df[,1],y=poly.df[,2],id=ids)

  polys <- suppressWarnings(df2polygons(poly.df,"id",c("x","y"),proj)) ## warnings for three/four point polygons

  return(polys)
}


