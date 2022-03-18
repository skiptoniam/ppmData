#' @name lonlat_from_window
#' @param window
#' @param coord.name
#' @param mask.na
#' @export
#' @importFrom raster xyFromCell ncell mask stack
#' @examples
#' xy.st <- lonlat_from_window(preds[[1]],mask.na=TRUE)

lonlat_from_window <- function(window, coord.name = c("X","Y"),mask.na=FALSE){

          grid.locations <- raster::xyFromCell(window,1:ncell(window))
          X <- Y <- window
          X[] <- grid.locations[,1]
          Y[] <- grid.locations[,2]
          names(X) <- coord.name[1]
          names(Y) <- coord.name[2]

          st.lonlat <- stack(X,Y)

          if(mask.na){
            st.lonlat <- mask(st.lonlat,window)
          }

          return(st.lonlat)

}

#'@export
is.formula <- function(x) inherits(x, "formula")

#'@export
merge.formula <- function(x, y, ...){
  if(!is.formula(x) || length(x) != 3)
    stop("First argument is invalid")
  if(!is.formula(y)) stop("Second argument is invalid")
  if(length(list(...))) warning("extraneous arguments discarded")
  is.gEnv <- function(e) identical(e, .GlobalEnv)

  str <- paste(c(deparse(x[[2]]), "~",
                 deparse(x[[3]]), "+",
                 deparse(y[[length(y)]])), collapse = "")
  f <- as.formula(str)
  ## MM: try to keep environment (where reasonable)
  ex <- environment(x)
  ey <- environment(y)
  if(!is.gEnv(ex)) {
    environment(f) <- ex
    if(!is.gEnv(ey) && !identical(ex,ey)) {
      warning("`x' and `y' have different environments; x's is used")
    }
  } else if(!is.gEnv(ey))
    environment(f) <- ey
  f
}

## Helper functions
#'@export
list2numeric <- function(x){
  as.numeric(as.matrix(x))
}
