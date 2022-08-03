#' @rdname plot.ppmData
#' @name plot.ppmData
#' @title Plot a ppmData object
#' @param x A ppmData object.
#' @param main Title of the plot, default is 'Quadrature Scheme'
#' @param rcols Vector of colours or colour codes for raster
#' @param pch.cols Vector of point shapes
#' @param pchs Vector of pch codes
#' @param plt.legend boolean plot the legend or not, default is TRUE.
#' @param \\dots Ignored
#' @export


plot.ppmData <- function(x,
                         main='Quadrature Scheme',
                         rcols = gray(0.5,0.75),
                         pch.cols =  c("#1B9E77","#D95F02","#7570B3",
                                   "#E7298A","#E6AB02",
                                   "#A6761D","#666666"),
                         pchs = 15:20,
                         plt.legend=TRUE, ...){


  # op <- graphics::par(no.readonly = TRUE)
  # on.exit(graphics::par(op))

  colshp <-  expand.grid(pch.cols,pchs)

  if(x$marked){

    marks <- x$presences.cleaned
    quad <- x$ppmData$locations[x$ppmData$bkg,]

    # par(mar=c(7,7,7,7))
    if(x$params$dw){
      e <- terra::ext(x$window)
      p <- terra::as.polygons(e)
      terra::plot(p,col=rcols,main=main,axes=FALSE, legend=FALSE)
    } else {
      p <- x$window
      terra::plot(p,col=rcols,axes=FALSE, legend=FALSE, main=main)
    }
    points(quad,pch='.')
    points(marks,col=as.character(colshp[as.numeric(as.factor(marks[,x$params$species.id])),1]),
           pch=colshp[as.numeric(as.factor(marks[,x$params$species.id])),2], ...)
    if(plt.legend){
    legend(x="bottom",
           legend=unique(marks[,x$params$species.id]),
           col=as.character(colshp[1:length(unique(marks[,x$params$species.id])),1]),
           pch=colshp[1:length(unique(marks[,x$params$species.id])),2],
           cex=0.75,
           xpd = TRUE, horiz = FALSE, inset = c(-.2, 0),
           bty="n")
    }

  } else {

    pressies <- x$presences.cleaned
    quad <- x$ppmData[x$ppmData$presence%in%0,x$params$coord]

    # par(mar=c(7,7,7,7))
    if(x$params$dw){
    e <- terra::ext(x$window)
    p <- terra::as.polygons(e)
    terra::plot(p,col=rcols,main=main,axes=FALSE, legend=FALSE)
    } else {
      p <- x$window
      terra::plot(p,col=rcols,axes=FALSE, legend=FALSE,main=main)
    }
    points(quad,pch='.')
    points(pressies,col=pch.cols[1],pch=16, ...)
    if(plt.legend){
    legend(x="bottom",
           legend=unique(pressies[,x$params$species.id]),
           col=pch.cols[1],
           pch=16,
           cex=0.75,
           xpd = TRUE, horiz = TRUE, inset = c(0, -0.2),
           bty="n")
    }
  }

}

#'@rdname print.ppmData
#'@name print.ppmData
#'@title Print a summary of ppmData object.
#'@param x A ppmData object.
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



