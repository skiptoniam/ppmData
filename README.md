qrbp is an r package for generating quasi random background points for Poisson point process models.
----------------------------------------------------------------------------------------------------

The package aims to generate background points that...

``` r
devtools::install_github('skiptoniam/qrbp')
```

Estimating the probabilty of presence from a series of spatial points. The probability of *absence in an area of size A* according to the poisson distribution is: *p**r*(*y* = 0)=*e**x**p*(−*λ*(*u*)\**A*)

The prob of *presence* is then *p**r*(*y* = 1)=1 − *p**r*(*y* = 0) =1 − *e**x**p*(−*λ*(*u*)\**A*) where *λ*(*u*) = the intensity value at point *u* and *A* is the area of the sampling unit (cell size). This is estimated using
