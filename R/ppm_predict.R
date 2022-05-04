#'@title ppmFit
#'@description This function should predict an intensity surface based on the
#'on the model fitting in the ppmFit.
#'@param object A fitted ppmFit object.
#'@param newdata SpatRaster. A terra raster stack of covariates for the model, or it can be a data.frame.
#'@param type Character. Either "response","link", "unit" & "tloglog". The type of response variable to return. The default is 'response' which is on the intensity scale or 'link' which is one the linear predictor scale (log).
#'@param offset Numeric vector or raster. If an offset is used in the model. Either an observed offset at prediction sites. If an offset is used and this is not known at prediction sites something like the mean offset used to fit the model can be used.
#'@param slambda Character Either 'lambda.min' or 'lambda.1se'. Value(s) of the penalty parameter lambda at which predictions are required. Default is "lambda.min".
#'@param quad.only Logical. If TRUE prediction is only done at the quadrature locations - useful for some of the diagnostic tools.
#'@param cores Integer. The number of cores to use in the prediction, useful for large rasters.
#'@export
#'@examples
#'library(ppmData)
#'library(terra)
#'path <- system.file("extdata", package = "ppmData")
#'lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#'covariates <- rast(lst)
#'bias <- covariates[[1]]
#'names(bias) <- "bias"
#'covariates <- c(covariates,bias)
#'presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#'ppmdata <- ppmData(npoints = 1000,presences=presences, window = covariates[[1]], covariates = covariates)
#'sp_form <- presence/weights ~ poly(annual_mean_precip,2) + poly(annual_mean_temp,2) + poly(distance_from_main_roads,2)# + offset(log(bias))
#'## Fit a ppmlasso
#' ft.ppm <- ppmFit(species_formula = sp_form, ppmdata=ppmdata)
#' pred <- ppmData:::predict.ppmFit(ft.ppm, covariates)
#'## Fit a ppmlasso with offset


predict.ppmFit <- function(object,
                           # bootobject=NULL,
                           newdata = NULL,
                           type = c("response","link","unit","cloglog"),
                           offset = NULL,
                           slambda= c("lambda.min","lambda.1se"),
                           quad.only = TRUE,
                           cores = 1,
                           ...){

  # newdata <- covariates
  type <- match.arg(type)
  model <- class(object[[1]])[1]
  slambda <- match.arg(slambda)
  object.mod <- object[[1]]

  ## check if data is supplied
  if (is.null(newdata)) {
      newdata <- getPredQuad(object.mod,quad.only)
      newdata <- newdata$newdata
      wts <- newdata$wts
  }

  ## set up the offset, for use as a thinned process.
  if(is.null(offset)){
    offy <- getPredOffset(object.mod, newdata, quad.only)
  } else {
    offy <- offset
  }

  if(class(newdata)=="SpatRaster"){
    ## if you pass a spatRaster stack
    ## let's predict as a raster
    ## Do predictions directly on the rasters with terra
    if(model=="ppmlasso")
      pred <- terra::predict(object = newdata, model=object.mod,
                             const = data.frame(wt = 1),
                             fun=pred.fun.ppmlasso, na.rm=TRUE, cores=cores)
    if(model=="glm")
      pred <- terra::predict(object = newdata, model=object.mod,
                             # const = data.frame(weight = 1),
                             na.rm=TRUE, cores=cores)
    if(model=="gam")
      pred <- terra::predict(object = newdata, model=object.mod,
                             # const = data.frame(weight = 1),
                             na.rm=TRUE, cores=cores)
    if(model=="lasso")
      pred <- terra::predict(object = newdata, model=object.mod,
                             # const = data.frame(weight = 1),
                             fun=pred.fun.glmnet, na.rm=TRUE, cores=cores)
    if(model=="ridge")
      pred <- terra::predict(object = newdata, model=object.mod,
                             # const = data.frame(weight = 1),
                             fun=pred.fun.glmnet, na.rm=TRUE, cores=cores)

    ## change the type
    if(type=="link")
      pred <- log(pred);
    if(type=="unit")
      pred <- pred*prod(res(pred))
    if(type=="cloglog"){
      cell.pred <- pred*prod(res(pred))
      Lambda <- terra::global(cell.pred,"sum",na.rm=TRUE)
      pred <- 1-exp(-pred/as.numeric(Lambda))
    }

  } else {
    ## Do prediction on a data.frame
    if(model=="ppmlasso")
      pred <- ppmlasso::predict.ppmlasso(object = object.mod, newdata = newdata)
    if(model=="glm")
      pred <- predict(object = object.mod, newdata = newdata, type = type)
    if(model=="gam")
      pred <- predict(object = object.mod, newdata = newdata, type = type)
    if(model=="lasso")
      pred <- glmnet::predict.glmnet(object = object.mod, newx = newdata, type=type, s=slambda)
    if(model=="ridge")
      pred <- glmnet::predict.glmnet(object = object.mod, newx = newdata, type = type, s=slambda)

    ## change the type
    if(type=="link")
      pred <- log(pred);
    if(type=="unit")
      pred <- pred*wts
    if(type=="cloglog"){
      tmp.pred <- pred*wts
      Lambda <- sum(tmp.pred)
      pred <- 1-exp(-pred/Lambda)
    }
  }
  pred
}

## This is a wrapper for fitting a ppmlasso models using terra rasters - this will be important for big data.
pred.fun.ppmlasso <- function(object, dat, type, offy, ...) {
  ppmlasso::predict.ppmlasso(object = object, newdata = dat, ...)
}

## Wrapper for predicting glmnet to terra rast
pred.fun.glmnet <- function(object, dat, ...) {
  glmnet::predict.glmnet(object = object, newx = dat, type = type, s = slambda, newoffset = offy, ...)
}

predict.ppmlasso2 <- function(object, newdata=NULL, type = c("response","link"),
                              offset=NULL, interactions = NA, ...){

  if (any(lapply(newdata, class) == "factor")) {
    unpacknewdata = CatConvert(newdata)
    newdata = unpacknewdata$X
    cat.names = setdiff(unique(unpacknewdata$cat.names),
                        NA)
    use.form = as.character(object$formula)[2]
    for (i in 1:length(cat.names)) {
      use.form = gsub(cat.names[i], paste(names(newdata)[which(unpacknewdata$cat.names ==
                                                                 cat.names[i])], collapse = " + "), use.form)
    }
    object$formula = as.formula(paste("~", use.form))
  }
  var.0 = which(apply(newdata, 2, var) == 0)
  if (length(var.0) > 0) {
    newdata[, var.0] = newdata[, var.0] + rnorm(dim(newdata)[1],
                                                0, 1e-08)
  }
  mf = model.frame(object$formula, data = newdata)
  mt = attr(mf, "terms")
  X.des = if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(0, length(object$mu), 0L)
  X.var = X.des[, -1]
  if (is.null(object$s.means) == FALSE) {
    X.var = scale(X.var, center = object$s.means, scale = object$s.sds)
    X.des = cbind(1, X.var)
  }
  if (object$family == "area.inter") {
    if (is.na(interactions) == TRUE) {
      if (is.null(object$s.means) == FALSE) {
        X.des = cbind(X.des, min(scale(object$pt.interactions)))
      }
      if (is.null(object$s.means) == TRUE) {
        X.des = cbind(X.des, 0)
      }
    }
    if (is.na(interactions) == FALSE) {
      if (is.null(object$s.means) == FALSE) {
        X.des = cbind(X.des, scale(interactions, center = mean(object$pt.interactions),
                                   scale = sd(object$pt.interactions)))
      }
      if (is.null(object$s.means) == TRUE) {
        X.des = cbind(X.des, interactions)
      }
    }
  }

  link <- make.link('log')
  eta <- as.matrix(X.des) %*% object$beta + offy

  pred.int <- switch(type,
                     reponse = link$linkinv(eta),
                     link = eta)
  return(pred.int)

}

getPredQuad <- function(object.mod, quad.only){

  if(class(object.mod)[1]=="ppmlasso"){
    if(quad.only){
      X <- as.data.frame(object.mod$data[object.mod$pres==0,])
      wts <- object.mod$wt[object.mod$pres==0]
    } else {
      X <- as.data.frame(object.mod$data)
      wts <- object.mod$wt
    }
  }

  return(list(newdata=X,wts=wts))
}

getPredOffset <- function(object.mod, newdata, quad.only){

                      if(!is.null(newdata)){
                        if(class(newdata)=="SpatRaster"){
                          offy <- newdata[[1]]
                          id <- !is.na(offy[])
                          roffy <- ifelse(id,0,NA)
                          offy <- setValues(offy,roffy)
                          } else {
                          offy <- rep(0,nrow(newdata))
                        }
                      } else {
                        if(quad.only){
                          if(class(object.mod)[1]=="ppmlasso"){
                            offy <- rep(0,nrow(object.mod$data[object.mod$pres==0,]))
                          }
                        } else {
                          if(class(object.mod)[1]=="ppmlasso"){
                            offy <- rep(0,nrow(object.mod$data))
                          }
                        }
                      }
  return(offy)
}





