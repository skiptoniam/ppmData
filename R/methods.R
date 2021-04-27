#' @rdname plot.ppmData
#' @name plot.ppmData
#' @title Plot a ppmData object
#' @param x A model object.
#' @param \\dots Ignored
#' @export
#' @importFrom raster extent

plot.ppmData <- function(x, ...){


  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))

  cols <- c("#1B9E77","#D95F02","#7570B3",
            "#E7298A","#66A61E","#E6AB02",
            "#A6761D","#666666")
  # cols <- rep(cols,100)

  pchs <- 15:20
  colshp <-  expand.grid(cols,pchs)

  if(x$marked){

    marks <- x$presences
    quad <- x$ppmData$locations[x$ppmData$bkg,]

    par(mar=c(7,7,7,7))
    if(x$params$dw){
      e <- extent(x$window)
      p <- as(e, 'SpatialPolygons')
      raster::plot(p,main='quadrature scheme')
    } else {
      p <- x$window
      raster::plot(p,axes=FALSE, box=FALSE,legend=FALSE,main='quadrature scheme')
    }
    points(quad,pch='.')
    points(marks,col=as.character(colshp[as.numeric(as.factor(marks[,x$params$speciesIdx])),1]),
           pch=colshp[as.numeric(as.factor(marks[,x$params$speciesIdx])),2],cex=0.5)
    legend(x="left",
           legend=unique(marks[,x$params$speciesIdx]),
           col=as.character(colshp[1:length(unique(marks[,x$params$speciesIdx])),1]),
           pch=colshp[1:length(unique(marks[,x$params$speciesIdx])),2],
           cex=0.75,
           xpd = TRUE, horiz = FALSE, inset = c(-.2, 0),
           bty="n")


  } else {

    pressies <- x$presences
    quad <- x$ppmData[x$ppmData$presence%in%0,x$params$coord]

    par(mar=c(7,7,7,7))
    if(x$params$dw){
      e <- extent(x$window)
      p <- as(e, 'SpatialPolygons')
      raster::plot(p,main='quadrature scheme')
    } else {
      p <- x$window
      raster::plot(p,axes=FALSE, box=FALSE,legend=FALSE,main='quadrature scheme')
    }
    points(quad,pch='.')
    points(pressies,col='dodgerblue',pch=16,cex=0.5)
    legend(x="bottom",
           legend=unique(pressies[,x$params$speciesIdx]),
           col='dodgerblue',
           pch=16,
           cex=0.75,
           xpd = TRUE, horiz = TRUE, inset = c(0, -0.2),
           bty="n")
  }

}

#'@rdname print.ppmData
#'@name print.ppmData
#'@title Print a summary of ppmData object.
#'@param x A model object.
#'@param \\dots Ignored
#'@export

print.ppmData <- function (x, ...){

  if(!x$marked){

    message("There are ", sum(x$ppmData$presence), " presence observations for this species")
    message("There are ", sum(x$ppmData$presence==0), " background quadrature (integration) points")
    message("There are a total of ", nrow(x$ppmData), " sites in the model.matrix")
    no_nans <- nrow(x$ppmData[!complete.cases(x$ppmData),])
    if(no_nans>0)
      message("There are a total of ",no_nans, " NaNs in the covariates, check before modelling.")

  } else {

    message("There are ", x$ppmData$nUniquePres, " presence observations for ", x$ppmData$nspp," species")
    message("There are ", x$ppmData$nBkg, " quadrature (integration) points for each of the ", x$ppmData$nspp ," species")
    message("There are a total of ",x$ppmData$m, " sites in the model.matrix")
    if(!is.null(x$ppmData$covars)){
      num_nans <- nrow(x$ppmData$covars[!complete.cases(x$ppmData$covars),])
      if(num_nans>0)message("There are a total of ",num_nans, " NaNs in the covariates, check before modelling.")
    }
  }
}



