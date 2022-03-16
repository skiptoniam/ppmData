#'@param species_formula The ppm model formula.
#'@param bias_formula Default is NULL. The idea will be to implement this as part of an integrated model.
#'Currently its just a nice way to keep track of which covariates are used in species or bias variables.
#'@param ppmdata A ppmData data object
#'@param method A method to fit the a ppm. Default is 'glm'. Others options are: 'gam','lasso','ridge'
#'@param control Options to pass to fitting functions.
#'@author Skipton Woolley
#'@references Warton, D.I. and Shepherd, L.C., 2010. Poisson point process models solve the" pseudo-absence problem" for presence-only data in ecology. The Annals of Applied Statistics, pp.1383-1402. \url{https://doi.org/10.1214/10-AOAS331}
#'@example
#'library(ppmData)
#'path <- system.file("extdata", package = "ppmData")
#'lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#'preds <- raster::stack(lst)
#'presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#'ppmdata <- ppmData(npoints = 1000,presences=presences, window = preds[[1]], covariates = preds)
#'sp_form <- presence/weights ~ poly(X,2) + poly(Y,2) + poly(max_temp_hottest_month,2) + poly(annual_mean_precip,2) + poly(annual_mean_temp,2) + poly(distance_from_main_roads,2)
#'ft.ppm <- ppm.fit(species_formula = sp_form, ppmdata=ppmdata)
#'ft.ppm <- ppm.fit(species_formula = sp_form, ppmdata=ppmdata, method='gam')
#'ft.ppm <- ppm.fit(species_formula = sp_form, ppmdata=ppmdata, method='lasso')
#'ft.ppm <- ppm.fit(species_formula = sp_form, ppmdata=ppmdata, method='ridge')

ppm.fit <- function(species_formula = presence/weights ~ 1,
                    bias_formula = NULL,
                    ppmdata,
                    method="glm",
                    control=list()){

  # lambda will be a vector of
  if(!class(ppmdata)=="ppmData")
    stop("'ppm.fit' requires a 'ppmData' object to run.")
  if(!any(c("glm","gam","lasso","ridge","gibbs")%in%method))
    stop("'ppm.fit' must used one of the follow methods\n 'glm','gam','lasso','ridge' to run.")

  ## merge formulas
  if(!is.null(bias_formula)){
   form <- merge.formula(species_formula,bias_formula)
  } else {
   form <- species_formula
  }
  ## check if response is ok
  if(form[[2]]!="presence/weights"){
    form <- update(form,presence/weights ~ .)
  }

  ## setup the model matrix/frames
  if(method!="gam"){
    ## set up the data for ppm fit
    ppp <- ppmdata$ppmData
    mf <- model.frame(form,ppp,weights=ppmdata$ppmData$weights)
    mt <- terms(mf)
    x <- model.matrix(mt,mf)
    y <- model.response(mf)
    wts <- model.weights(mf)
    offy <- model.offset(mf)
  }

  if(is.null(offy))
    offy <- rep(0,length(y))

  if(method=="glm"){
    ft <- glm.fit(x = x, y = y/wts, weights = wts, offset = offy, family = poisson())
  }
  if(method=="gam"){
    ft <- mgcv::gam(formula = form, data = ppmdata$ppmData, weights = x$weights, offset = offy, family = poisson())
  }
  if(method=="lasso"){
    # ft <- ppmlasso::ppmlasso() ## maybe could sub in ppmlasso
    ft <- glmnet::glmnet(x=x, y=y/wts, weights = wts, offset = offy, family = "poisson", alpha = 1) #lasso
  }
  if(method=="ridge"){
    ft <- glmnet::glmnet(x=x, y=y/wts, weights = wts, offset = offy, family = "poisson", alpha = 0) #lasso
  }


  return(ft)

}


