#' wrapper around ppm.fit using a mccv

mccv <- function(species_formula,
                 bias_formula,
                 ppmdata,
                 method,
                 control,
                 nsim=10,
                 ncores=1,...){


    #need update the quadrature per-fit

    apply()


    ppm.fit()

}



#' @name p_thin
#' @rdname ppm_cv
#' @description We need to run a monte-carlo cross validation on the observed
#' point process. To do so requires the thinning of the observed pp. We can do
#' this via independent p-thinning. If the thinning is independent then the
#' resulting thinned point process can be used as another point process!
#' The current realization we have is the observed species locations. So if we
#' thin this based on a p-thinning and use the retained points as the validation
#' set and the thinned (removed) as the testing set we should be right.
#' @param sites The recorded presences (realisation) of a species point process
#' @param p The thinning probability (Bernoulli variable).
#' @author Skipton Woolley
#' @examples
#' library(ppmData)
#' presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#' p_thin(presences[,-3],p=0.25)

p_thin <- function(sites, p=0.25, seed = NULL){

  if(!ncol(sites)==2)stop("'sites' should just be the coordinates of the species presences.")

  # The number of
  n <- nrow(sites)

  # Bernoulli variable of known points being thinned
  if(!is.null(seed))set.seed(seed)

  boole_thin <- runif(n)<p;
  boole_retain <- !boole_thin

  #locations of thinned points
  thinned=sites[boole_thin,]
  retained=sites[boole_retain,]

  # return as list
  res <- list(thinned_sites=thinned,
              retained_sites=retained)

  return(res)

}


#' @name block_sample
#' @rdname ppm_cv
#' @description Do a block sample of sites
#' @param p The thinning probability (Bernoulli variable).
#' @author Skipton Woolley
#' @examples
#' library(ppmData)
#' presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#' block_sample(presences[,-3],p.blocks=0.25,block.res=0.1)

block_sample <- function(sites, coords=c("X","Y"), p.blocks = 0.25, block.res = 1, seed = NULL){


  cell.group <- rep(0, length(sites[,coords[1]]))
  n.groups   <- ceiling(c(max((sites[,coords[1]] - min(sites[,coords[1]]))/block.res), max((sites[,coords[2]] - min(sites[,coords[2]]))/block.res)))

  xq <- block.res*(1:n.groups[1]) + min(sites[,coords[1]])
  for (i.group in 1:n.groups[1]){
    cell.group <- cell.group + as.numeric(sites[,coords[1]] > xq[i.group])
  }

  yq <- block.res*(1:n.groups[2]) + min(sites[,coords[2]])
  for (i.group in 1:n.groups[2]){
    cell.group <- cell.group + n.groups[1] * as.numeric(sites[,coords[2]] > yq[i.group])
  }

  block.group <- factor(cell.group)

  cat("There are a total of",length(levels(block.group)), "blocks at a",block.res,"resolution")

  if(!is.null(seed))set.seed(seed)

  levs  <- levels(block.group)
  nlevs <- length(levs)
  p.blocks.in <- ceiling(p.blocks*nlevs)
  lev.sample <- sample(levs,p.blocks.in)

  block.sites <- which(block.group%in%lev.sample)

  #locations of thinned points
  block_train=sites[-block.sites,]
  block_test=sites[block.sites,]

  # return as list
  res <- list(block_train_sites=block_train,
              block_test_sites=block_test)

  return(res)
}
