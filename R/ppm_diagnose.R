# convert ppm.fit to a spatstat obje
# take a glm/glmnet and try to do diagnostics.
ppm.ss <- function (fit) {
  is.ai = is.numeric(fit$pt.interactions)
  pres.x = fit$x[fit$pres > 0]
  pres.y = fit$y[fit$pres > 0]
  quad.x = fit$x[fit$pres == 0]
  quad.y = fit$y[fit$pres == 0]
  ux = unique(quad.x)
  uy = unique(quad.y)
  ux = sort(ux)
  uy = sort(uy)
  nx = length(ux)
  ny = length(uy)
  quad.mask = matrix(NA, ny, nx, dimnames = list(uy, ux))
  col.ref = match(quad.x, ux)
  row.ref = match(quad.y, uy)
  all.vec = rep(NA, max(row.ref) * max(col.ref))
  vec.ref = (col.ref - 1) * max(row.ref) + row.ref
  all.vec[vec.ref] = 1
  num.vec = all.vec
  num.vec[is.na(all.vec)] = 0
  quad.mask = matrix(all.vec, max(row.ref), max(col.ref),
                     dimnames = list(uy, ux))
  quad.im = im(quad.mask, xcol = ux, yrow = uy)
  quad.win = as.owin(quad.im)
  pres.dat = ppp(pres.x, pres.y, window = quad.win, check = FALSE)
  quad.dat = ppp(quad.x, quad.y, window = quad.win, check = FALSE)
  Q = quad(data = pres.dat, dummy = quad.dat, w = fit$wt,
           param = list(weight = list(method = "grid", ntile = c(length(ux),
                                                                 length(uy)), npix = NULL, areas = num.vec)))
  num.var = dim(fit$data)[2] - 1 - is.ai
  cov.list = vector("list", num.var)
  for (var in 1:num.var) {
    x.dat = fit$data[fit$pres == 0, (var + 1)]
    v.vec = rep(NA, max(row.ref) * max(col.ref))
    v.vec[vec.ref] = x.dat
    x.mat = matrix(v.vec, max(row.ref), max(col.ref), dimnames = list(uy,
                                                                      ux))
    assign(paste("var.im.", var, sep = ""), im(x.mat, xcol = ux,
                                               yrow = uy))
    cov.list[[var]] = im(x.mat, xcol = ux, yrow = uy)
  }
  names(cov.list) = paste("V", 1:num.var, sep = "")
  trend = "~V1"
  for (var in 2:num.var) {
    trend = paste(trend, " + V", var, sep = "")
  }
  glmfit.form = as.formula(paste(".mpl.Y", trend))
  if (is.ai) {
    glmfit.form = as.formula(paste(".mpl.Y", trend, " + Interaction"))
  }
  trend = as.formula(trend)
  call.1 = quote(ppm(Q = Q, trend = trend, covariates = cov.list,
                     interaction = Poisson(), correction = "none"))
  if (is.ai) {
    call.1 = bquote(ppm(Q = Q, trend = trend, covariates = cov.list,
                        interaction = AreaInter(.(fit$r)), correction = "none"))
  }
  class(fit) = "ppm"
  fit$fitter = "glm"
  fit$coef = fit$beta
  fit$method = "mpl"
  fit$projected = FALSE
  fit$trend = trend
  fit$interaction = NULL
  if (is.ai) {
    fit$interaction = AreaInter(fit$r)
  }
  fit.int = list(name = "Poisson process", creator = "Poisson",
                 family = NULL, pot = NULL, par = NULL, parnames = NULL)
  class(fit.int) = "interact"
  fit$fitin = list(interaction = Poisson(), coefs = fit$beta,
                   Vnames = character(0), IsOffset = logical(0))
  if (is.ai) {
    fit$fitin = list(interaction = AreaInter(fit$r), coefs = fit$beta,
                     Vnames = "Interaction", IsOffset = FALSE)
  }
  class(fit$fitin) = c("fii", "list")
  fit$Q = Q
  fit$maxlogpl = fit$likelihood
  fit$covariates = cov.list
  fit$covfunargs = list()
  fit$correction = "none"
  fit$rbord = 0
  glmdata = data.frame(fit$wt, fit$pres/fit$wt, fit$data[,
                                                         -1], TRUE)
  if (is.ai == FALSE) {
    names(glmdata) = c(".mpl.W", ".mpl.Y", names(cov.list),
                       ".mpl.SUBSET")
  }
  if (is.ai) {
    names(glmdata) = c(".mpl.W", ".mpl.Y", names(cov.list),
                       "Interaction", ".mpl.SUBSET")
  }
  fit$version = list(major = 1, minor = 31, release = 0)
  fit$problems = list()
  fit$call = call.1
  fit$callstring = "character"
  fit$callframe = environment()
  terms.int = terms(glmfit.form)
  terms.ppm = terms(trend)
  fit$terms = terms.ppm
  fit$internal$glmfit$terms = terms.int
  glm.1 = glm(fit$pres/fit$wt ~ as.matrix(fit$data[, -1]),
              weights = fit$wt, family = poisson())
  glm.1$coefficients = fit$beta
  glm.1$fitted.values = fit$mu
  glm.1$residuals = (fit$pres/fit$wt - fit$mu)/fit$mu
  glm.1$linear.predictors = log(fit$mu)
  fam = poisson()
  vari = fam$variance(glm.1$fitted.values)
  deriv = 1/vari
  wt = fit$wt
  weii = wt * 1/(deriv^2 * vari)
  w.1 = weii
  Xw = t(as.vector(sqrt(weii)) * t(t(fit$data)))
  q.1 = qr(t(Xw))
  glm.1$qr = q.1
  glm.1$weights = w.1
  glm.1$family = quasi(log)
  glm.1$data = glmdata
  glm.1$formula = glmfit.form
  fit$internal = list(glmfit = glm.1, glmdata = glmdata)
  if (is.ai) {
    fit$internal = list(glmfit = glm.1, glmdata = glmdata,
                        Vnames = "Interaction", IsOffset = FALSE)
  }
  fit
}
