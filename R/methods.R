#' @rdname plot.ppmData
#' @name plot.ppmData
#' @title Plot a ppmData object
#' @param x A model object.
#' @param \\dots Ignored
#' @export

# x <- ppmdata1

plot.ppmData <- function(x, ...){

  # if(x$marked){
  #
  #   marks <- x$presences
  #   # mark.sites <- x$p,]
  #   # cbind(mark.sites, marks)
  #
  #
  # }

  raster::plot(x$window,legend=FALSE)
  points(x$ppmData$locations[x$ppmData$bkg,],pch='.')
  points(x$presences,col=as.numeric(as.factor(x$presences$SpeciesID)),pch=16,cex=0.8)



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
    if(no_nans>0)message("There are a total of ",no_nans, " NaNs in the covariates, check before modelling.")

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
