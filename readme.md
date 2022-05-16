
[![Travis build
status](https://travis-ci.org/skiptoniam/ppmData.svg?branch=master)](https://travis-ci.org/skiptoniam/ppmData)

## ppmData is an R package that can be used to set up a quadrature scheme for multiple species point processes modelling and could be used in [ecomix](https://github.com/skiptoniam/ecomix) or for other multiple species Point Processes.

The approach uses quasi-random sampling (Grafston & Tille, 2013, Foster
et al., 2018) to generate a quadrature scheme for numerical
approximation of a Poisson point process model (Berman & Turner 1992;
Warton & Shepard 2010). Quasi-random sampling quadrature are form of
spatially-balanced survey design or point stratification that aims to
reduce the frequency of placing samples close to each other (relative to
pseudo-random or grid designs). A quasi-random quadrature design
improves efficiency of background point sampling (and subsequent
modelling) by reducing the amount of spatial auto-correlation between
data implying that each sample is providing as much unique information
as possible (Grafston & Tille, 2013, Foster et al., 2018) and thus
reducing low errors for geostatistical prediction (Diggle & Ribeiro,
2007). Because the quasi-random design is not on a regular grid we use a
Dirichlet tessellation to generate polygons for each point. Areal
weights are then derived from these polygons. Under this approach we
**do not** use a targeted background scheme for generating
‘pseudo-absences’, but rather can account for bias in the sighting
process of presences using covariates or a offset when modelling (e.g
Warton et al., 2013; Fithian et al, 2015).

### References

Diggle, P. J., P. J. Ribeiro, Model-based Geostatistics. Springer Series
in Statistics. Springer, 2007.

Foster, S.D., Monk, J., Lawrence, E., Hayes, K.R., Hosack, G.R. and
Przeslawski, R., 2018. Statistical Considerations for Monitoring and
Sampling.

Grafstrom, Anton, and Yves Tille. “Doubly balanced spatial sampling with
spreading and restitution of auxiliary totals.” Environmetrics 24.2
(2013): 120-131.

Renner, I.W., Elith, J., Baddeley, A., Fithian, W., Hastie, T.,
Phillips, S.J., Popovic, G. and Warton, D.I., 2015. Point process models
for presence only analysis. Methods in Ecology and Evolution, 6(4),
p.366-379.

Warton, D. I., and L. C. Shepherd. “Poisson point process models solve
the pseudo-absence problem for presence-only data in ecology.” The
Annals of Applied Statistics 4.3 (2010): 1383-1402.

Warton, D.I., Renner, I.W. and Ramp, D., 2013. Model-based control of
observer bias for the analysis of presence-only data in ecology. PloS
one, 8(11), p.e79168.
