## need a generic function which can take weights presences and quad points and generate weights.
#'@export
getWeights <- function(method,
                      presences,
                      quadrature,
                      quadDummy, #a hack to reduce large weights on edge of dirichlet.
                      window,
                      coord,
                      speciesIdx,
                      unit){

  wts <- switch(method,
                grid=gridWeights(presences,quadrature,window,coord,unit),
                pseudo=pseudoRandomWeights(presences,quadrature,window,coord,unit),
                quasi=quasiRandomWeights(presences, quadrature, quadDummy, window, coord, speciesIdx))

  return(wts)
}



getSingleSpeciesWeights <- function(method,
                                    presences,
                                    quadrature,
                                    quadDummy,
                                    window,
                                    coord,
                                    speciesIdx,
                                    unit){

  quadrature[[speciesIdx]] <- "quad"
  if(!is.null(quadDummy))quadDummy[[speciesIdx]] <- "dummy"
  wts <- getWeights(method, presences, quadrature, quadDummy, window, coord, speciesIdx, unit)
  wts$OrigOrder <- wts$id
  wts$DatasetID <- 1
  wts$pres <- ifelse(wts$dataset=="quad",0,1)
  colnames(wts)[which(colnames(wts)=="area")] <- "wts"
  df <- getSiteID(wts,coord)
  return(df)
}


combineDF.fun <- function( ii, xxx, yyy, coords, speciesIdx){
  ####  Assumes that the colnames of xxx are a subset of those from yyy
  ####  This is not a totally memory efficient implementation: data on all species is passed to all species...
  xxx <- xxx[[ii]]
  newdf <- as.data.frame( matrix( NA, nrow=nrow( xxx) + nrow( yyy), ncol=length( coords) + 2))
  colnames( newdf) <- c(coords,speciesIdx,"OrigOrder")
  newdf[1:nrow( xxx), ] <- xxx[,colnames( newdf)]
  newdf[nrow( xxx) + 1:nrow( yyy), ] <- yyy[,colnames( newdf)]

  return( newdf)
}

getMultispeciesWeights <- function(method, presences, quadrature, quadDummy, window,
                                   coord, speciesIdx, mc.cores, sppNames, unit){

  presences$OrigOrder <- seq_len(nrow(presences))
  nspp <- length(unique(presences[,speciesIdx]))
  spps <- sppNames$sppNames
  # print(spps)
  quadrature[[speciesIdx]] <- "quad"
  quadrature[["OrigOrder"]] <- seq_len(nrow(quadrature))+max(presences[["OrigOrder"]])
  if(!is.null(quadDummy)){
    quadDummy[[speciesIdx]] <- "dummy"
    quadDummy[["OrigOrder"]] <- seq_len(nrow(quadDummy))+max(quadrature[["OrigOrder"]])
  }
  sppdata <- lapply(seq_len(nspp), function(ii)presences[presences[,speciesIdx]==spps[ii],])

  sppBckWtsList <- plapply(seq_len(nspp), function(ii) {getWeights( method, sppdata[[ii]], quadrature, quadDummy, window, coord, speciesIdx, unit)}, .parallel = mc.cores, .verbose = FALSE)
  sppBckDatList <- plapply(seq_len(nspp), combineDF.fun, xxx=sppdata, yyy=quadrature, coords=coord, speciesIdx = speciesIdx, .parallel = mc.cores, .verbose = FALSE)
  sppCounts <-  plapply(seq_len(nspp),function(ii)nrow(sppdata[[ii]]),.parallel = mc.cores, .verbose = FALSE)
  sppBckDatList <- plapply(seq_len(nspp),function(ii){sppBckDatList[[ii]]$DatasetID <- ii;sppBckDatList[[ii]]},.parallel = mc.cores, .verbose = FALSE)
  sppWtsList <- plapply(seq_len(nspp), function(ii)cbind(sppBckDatList[[ii]],
                                                         pres=c(rep(1,sppCounts[[ii]]),rep(0,nrow(quadrature))),
                                                         wts=sppBckWtsList[[ii]]),.parallel = mc.cores, .verbose = FALSE)

  dat <- do.call(rbind,sppWtsList)
  df <- getSiteID(dat,coord)
  return(df)

}

getSiteID <- function(dat,coord){
  stidfn <- function(df, cols) {
    comb <- do.call(paste, c(as.list(df[cols]), sep = "."))
    df$SiteID <- match(comb, unique(comb))
    df
  }
  df <- stidfn(dat,c(coord,'pres'))
  return(df)
}

convert2pts <- function(window) {
  #convert raster to extreme points (of bounding rectangle)
  #example only, *should* work for raster.  But what data type have we got?
  tmp <- sp::coordinates( window)
  tmp1 <- c( range( tmp[,1]), range( tmp[,2]))

  return( tmp1)
  #be careful with extents of rasters though, they tend to have slightly larger dimensions than they ought
}


