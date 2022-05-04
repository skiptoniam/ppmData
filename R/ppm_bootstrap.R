#' @title Non-parametric bootstrap to assess uncertainty in ppm parameters
#' @param object A ppm model object
#' @param nboot Default is 10; The number of bootstraps to run.
#' @param quiet Default is TRUE; Run the bootstrap with out any trace.
#' @param mc.cores Default is 1; the number of cores if run in parallel
#' @param \dots Ignored
#' @description This function uses a non-parametric bootstrap to estimate
#' uncertainty in a ppm. Conditional on X_j = {X_i,...,X_n}, let N*
#' have a Poisson distribution with Lambda_j. Lambda_j can be estimated as the
#' integral for each species intensities (lambda_ij), conditional on the
#' archetypes. For each species we can then draw \eqn{X*_{1},...,X*_{N}} by
#' sampling randomly with replacement N* times from X. We refit the ppm with
#' this new realisation of the ppm data. Estimated coefficients are returned
#' from each of B bootstrap model fits. This bootstrap object can then be used
#' in the to predict \link[stats]{predict} ppm and generate confidence intervals
#' function. If the model is very large (many po observations + quadrature sites)
#' vthis will be slow, and will take nboot times the time it takes to fit a
#' single ppm with the same data. So be warned for large datasets. We implement
#' method two from Cowling et al., (1996).
#' @references Cowling, A., Hall, P. and Phillips, M.J., 1996. Bootstrap
#' confidence regions for the intensity of a Poisson point process. Journal of
#' the American Statistical Association, 91(436), pp.1516-1524.
#' @references Woolley, S.N.C., Dunstan, P.K., Warton, D., and Foster, S.,
#' InPrep. Multiple species modelling of presence-only data-sets with
#' mixture-of-regressions models.

#' @export

#' @rdname bootstrap
#' @export bootstrap
"bootstrap" <- function (object, nboot=10, quiet=TRUE, mc.cores=1, ...){
  UseMethod("bootstrap", object)
}

#' @export
"bootstrap.ppm" <- function (object, nboot=10, quiet=TRUE, mc.cores=1, ...){

  if( is.null(object$titbits$X))
    stop( "Argument object has been created with the control option titbits=FALSE.\n
         These are now needed. Please re-fit with titbits=TRUE.")

  # fit0 <-
  fit0 <- predict(object = object,)
  # wtd_fitty <- fit0*object$ppmdata$ppmData$wts[-object$na.sites,]
  Lambda_j <- colSums(wtd_fitty,na.rm=TRUE)

  coef_star <- plapply(X=1:nboot,
                       FUN=ppmsam_bootFun,
                       Lambda_j=Lambda_j,
                       object=object,
                       quiet=quiet,
                       .parallel = mc.cores,
                       .verbose = TRUE)

  boot.estis <- do.call("rbind", coef_star)
  class( boot.estis) <- "species_mix_ppm.bootstrap"
  return(boot.estis)

}


"ppmsam_bootFun" <- function(X, Lambda_j, object, quiet = TRUE){
  #simulate N_j from Lambda_j
  N_star <- rpois(n=object$S, lambda=Lambda_j)
  #find the scaled intensity for each of the observed spp locations
  dat_star <- object$ppmdata
  for( jj in 1:object$S){
    ids <- sort( which( !is.na( object$ppmdata$ppmData$y[1:object$ppmdata$ppmData$nUniquePres,jj])))
    ids_star <- sample( ids, size=N_star[jj], replace=TRUE)
    bootWts <- table( factor( ids_star, levels=ids))
    #the bootstrap outcome
    dat_star$ppmData$y[ids,jj] <- bootWts

    #the bootstrap working variate (note that this extends the Berman and Turner as there are repeat presences in the bootstrap sample)
    dat_star$ppmData$z[ids,jj] <- dat_star$ppmData$y[ids,jj] / dat_star$ppmData$wts[ids,jj]

  }

  controlly <- object$titbits$control
  if(quiet){
    controlly$quiet <- TRUE
    controlly$trace <- 0
  }

  etta_boots_gerald <- additive_logistic(object$pi,inv = TRUE)[seq_len(object$G-1)]
  names(etta_boots_gerald) <- names(object$eta)
  object$eta <- etta_boots_gerald
  tmp_star <- species_mix_ppm.fit(y = dat_star$ppmData$y[-object$na.sites,],
                                  X = object$titbits$X,
                                  W = object$titbits$W,
                                  U = object$titbits$U,
                                  G = object$G, S = object$S,
                                  site_spp_weights = dat_star$ppmData$wts[-object$na.sites,],
                                  offset = object$titbits$offset,
                                  y_is_na = is.na(dat_star$ppmData$y[-object$na.sites,]),
                                  control =  controlly,
                                  inits = object$coefs)

  boot.covars <- c(alpha=tmp_star$alpha,beta=tmp_star$beta,tmp_star$eta,
                   gamma=tmp_star$gamma,delta=tmp_star$delta)
  gc()
  return(boot.covars)
}
