[![Travis build
status](https://travis-ci.org/skiptoniam/ppmData.svg?branch=master)](https://travis-ci.org/skiptoniam/ppmData)

## ppmData is an R package that can be used to set up a quadrature scheme for multiple species point processes modelling and could be used in [ecomix](https://github.com/skiptoniam/ecomix) or for other multiple species Point Processes.

The approach uses quasi-random sampling to generate a quadrature scheme
based (e.g Berman & Turner 1992; Foster et al, 2018). Quasi-random
sampling quadrature are form of spatially-balanced survey design or
point stratification that aims to reduce the frequency of placing
samples close to each other (relative to pseudo-random or grid designs).
A quasi-random quadrature design improves efficiency of background point
sampling (and subsequent modelling) by reducing the amount of spatial
auto-correlation between data implying that each sample is providing as
much unique information as possible (Grafston & Tille, 2013, Foster et
al., 2018) and thus reducing low errors for geostatistical prediction
(Diggle & Ribeiro, 2007). Because the quasi-random design is not on a
regular grid we use a Dirichlet tessellation to generate polygons for
each point. Areal weights are then derived from these polygons.

### References

Diggle, P. J., P. J. Ribeiro, Model-based Geostatistics. Springer Series
in Statistics. Springer, 2007.

Foster, S.D., Monk, J., Lawrence, E., Hayes, K.R., Hosack, G.R. and
Przeslawski, R., 2018. Statistical Considerations for Monitoring and
Sampling.

Grafstrom, Anton, and Yves Tille. “Doubly balanced spatial sampling with
spreading and restitution of auxiliary totals.” Environmetrics 24.2
(2013): 120-131.

Warton, D. I., and L. C. Shepherd. “Poisson point process models solve
the pseudo-absence problem for presence-only data in ecology.” The
Annals of Applied Statistics 4.3 (2010): 1383-1402.
