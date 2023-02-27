---
title: "The ppmData R-package for setting up spatial point process models"
tags:
  - Point Process
  - Quasi-random Sampling
  - Dirichlet Tesselation
  - Inhomogeneous Poisson Point Process
  - R package
authors:
  - name: Skipton N.C. Woolley
    orcid: 0000-0001-7022-0509
    affiliation: "1, 2"
  - name: Scott D. Foster
    orcid: 0000-0002-6719-8002
    affiliation: "3"
affiliations:
 - name: Environment, CSIRO
   index: 1   
 - name: School of Ecosystems and Forestry Science, The University of Melbourne
   index: 2
 - name: Data61, CSIRO
   index: 3    
citation_author: Woolley & Foster    
date: 14 December 2022
year: 2022
bibliography: paper.bib
header-includes: 
 - \usepackage{amssymb}
 - \newcommand{\R}{\mathbb{R}}
 - \newcommand\numberthis{\addtocounter{equation}{1}\tag{\theequation}}
 - \usepackage{color}
 - \usepackage{bm}
output: rticles::joss_article
csl: apa.csl
journal: JOSS
editor_options: 
  chunk_output_type: console
---

# Summary

Spatial data, the locations of objects in a spatial region, is a common type of data recorded in ecology, environmental sciences and epidemiology amongst many others. An appropriate statistical approach for data which can be described as points, such as plants or individuals of an animal population, are spatial point process models [@cressie_statistics_1993; @baddeley_practical_2000]. Spatial point processes can be used to understand the location of objects to inform decision making and improve ecological, health and other socio-economic outcomes. For example, the locations of individuals from a species' population are typically not distributed completely at random. Locations might be associated with favorable habitat conditions for the niche of that species. We might expect these species locations and the underlying spatial intensity are distributed according to relevant covariates, such as environmental conditions at those locations. 

Inhomogeneous Poisson Process Models allow for the spatial intensity of a point process to vary across space, usually as a function of spatially referenced covariates. One challenge when fitting an Inhomogeneous Poisson Process Model is that the integral used to described the intensity surface in the log-likelihood (as described below) has no known closed analytic solution. For Inhomogeneous Poisson Process Models, the Berman-Turner device [@berman_approximating_1992] is commonly used to approximate the model using quadrature within a weighted Poisson generalized linear model. The type of quadrature scheme can have important ramifications for overall model fitting and inference [@warton_poisson_2010]. Approximating the integral accurately often comes at an increased computation cost, where a large number of quadrature points are required to robustly estimate the integral [@renner_point_2015].

Our `ppmData` package is setup to use a quasi-random quadrature approach as an alternative to pseudo-random sampling [@phillips_sample_2009] and regular grid [@warton_poisson_2010] quadrature approaches. Quasi-random sampling, a base-case of Balanced Acceptance Sampling [@robertson_bas_2013], can be used to efficiently perform numerical integration [@halton_efficiency_1960] and create spatially balanced survey designs. Quasi-random sampling is an efficient form of spatial sampling as it approximately balances over all spatially smooth covariates [@grafstrom_doubly_2013], even if they are not included or considered in the quadrature generation. When used in the fitting of a Inhomogeneous Poisson Process Models, quasi-random (spatially balanced) quadrature also reduces the influence of spatial autocorrelation between quadrature points and is likely to improve numerical approximation of the intensity [@foster_spatially_2017; @liu_bayesian_2020].

Here we present the `ppmData` R package that is designed to setup a quadrature scheme using Dirichlet tessellation for fitting spatial point process models. `ppmData` can setup a quadrature scheme of point process (e.g the locations of species) or a marked point process (e.g the locations of multiple species, where each location is associated with a specific species). In this paper, we demonstrate how to set up quadrature for inhomogeneous Poisson Process Models and provide a simple example of how we can fit an inhomogeneous Poisson Process Model using a `ppmData` object.

# Inhomogeneous Poisson Process Models

An Inhomogeneous Poisson Process Model is a statistical model used to describe a random point pattern for $n$ observed locations, ($\bm{y}=(y(\bm{s}_1),y(\bm{s}_2),\ldots,y(\bm{s}_n))$, in a spatial window $\mathcal{A}\subset R^2$. An Inhomogeneous Poisson Process Model assumes that the number of points in $\mathcal{A}$, $n$, is a Poisson random variable with total intensity $\Lambda$ defined as $\Lambda = \int_{\bm{s}\in\mathcal{A}}\lambda(\bm{s})d\bm{s}$. The spatial function $\lambda(\bm{s})$ is usually a non-constant function and defines the spatial surface that the spatial locations $\bm{y}$ are independently drawn from [@cressie_statistics_1993]. For convenience, we use subscript notation for particular spatial locations, for example, $\lambda_i=\lambda(\bm{s}_i)$ is the intensity surface of the point process at the $i^{th}$ location with coordinates $\bm{s}_i$. The intensity surface can be defined as $\log\{\lambda_i\}=\bm{x}_i^{\top}\bm{{\beta}}$ where $\bm{x}_i$ is a vector of spatially located covariates and $\bm{{\beta}}$ is a corresponding vector of coefficients.

The approximate log-likelihood $\ell(\bm{{\beta}};\{\bm{x}_i\})$ can be written as a weighted Poisson likelihood [@berman_approximating_1992]

```{=tex}
\begin{align}
\ell(\bm{{\beta}}|\bm{y}) & = \sum_{i=1}^n\log(\lambda_i) - \int_{\bm{s} \in \mathcal{A}}\lambda(\bm{s})d\bm{s} - \log(n!) \nonumber\\
& \approx \sum_{i=1}^n\log(\lambda_i) - \sum_{i=1}^m w_i\lambda_i \nonumber\\
& = \sum_{i=1}^m w_i\{z_i \log(\lambda_i) - \lambda_i\} \numberthis \label{eq:two}
\end{align}
```
where $z_i$ is an indicator variable identifying if site $i$ is one of the $n$ observed points from the random point pattern or one of the $m$ quadrature points. $w_i$ stores the quadrature weights which sum to the total area of the region $|\mathcal{A}|$. In \eqref{eq:two}, we can see the integral in the log-likelihood $\ell(\beta|\mathbf{y})$ is being estimated using quadrature [@baddeley_practical_2000; @diggle_geostatistical_2010]. @warton_poisson_2010 demonstrated for presence-only species distribution models, that using numerical quadrature with a point process framework makes the models scale invariant and treats the quadrature purely as a tool for approximating the integral.

\begin{figure}

{\centering \includegraphics{paper_files/figure-latex/fig1-1} 

}

\caption{Quadrature schemes generated using the \texttt{ppmData} package for the species \textit{Tasmaphena sinclairi} located within Tasmania, Australia. The red points represent the known locations of \textit{Tasmaphena sinclairi}. The black points represent the quadrature locations. A) Quasi-random quadrature where the integration points are generated using a quasi-random areal sample. Weights for each integration point are calculated as the dual of the subsequent Dirichlet tessellation of both the presence and integration points. B) Pseudo-random quadrature generated using random points from within the window. C) Grid quadrature where quadrature points are generated on a regular grid.}\label{fig:fig1}
\end{figure}

# Statement of need

There are a number of ways to setup a quadrature scheme to approximate the integral $\Lambda = \int_{\bm{s}\in\mathcal{A}}\lambda(\bm{s})d\bm{s}$. One common approach is to use a regular grid [@warton_poisson_2010], which creates a regular grid at a known resolution within a spatial window ($\mathcal{A} \subset R^2$). Despite the simplicity of the grid approach, it is not without detraction though. In particular, it is unlikely to handle edge effects well (when there is a consistent and potentially large pattern near the boundaries of $\mathcal{A}$). Then the grid will either be defined away from the edge, missing important trends, or defined on the edge, over-representing this trend. Also, another potential problem stems from the unlikely event that the spatial surface $\lambda(\bm{s}_i)$ has periodicity in its pattern aligned to the direction of the grid, such as repeated locations of mountain ridges. In that case, the grid-based quadrature will over- or under-estimate the integral depending on whether the grid coincides with peaks or troughs of the surface.

Another common approach is to use a pseudo-random spatial samples and generate approximate areal weights per quadrature point [@phillips_sample_2009; @renner_point_2015]. Typically, this is done by generating pseudo-random points within the study window $\mathcal{A} \subset R^2$ and assuming that all quadrature points have equal weights, generally calculated as $w_i=\frac{m}{|\mathcal{A}|}$ [@renner_point_2015]. This approach is unlikely to suffer either of the potential problems from the grid approach, however it may be inefficient -- requiring more points to achieve the same level of accuracy [e.g. @robertson_bas_2013].

Here we use the quasi-random quadrature as a method to trade-off the robustness of pseudo-random sampling and the efficiency of the grid. Quasi-random quadrature spreads the quadrature points in space but simultaneously retains some of the important properties of random sampling. Quasi-random quadrature schemes are available in other R packages, like the excellent `spatstat` [@baddeley_spatial_2015]. The `rQuasi` function can be used for generating a quasi-random spatial sample within a window. However, if the window is not a rectangle, the number of points will be reduced to include only the ones remaining inside the window boundary, so specifying the number of active points is guess work. Our package `ppmData` allows the exact number of quadrature points to be specified within a spatial window.

The package `spatstat` also allows a user to develop a quasi-random quadrature scheme, but this appears to be quite slow when used with Dirichlet weights for larger numbers of dummy (quadrature) points (\> 10000). For example, a regular window with 100000 quadrature points takes approximately 2.5 minutes to run in `spatstat` using `quadscheme` with `method=Dirichlet` with quasi-random dummy points. Our package computes the same Dirichlet tessellation in about 7.5 seconds (Comparisons were from a Linux operating system using a Intel i7 8 core processor, with 16 GB of RAM, see package readme for a quick comparison). Typically, calculating a Dirichlet tessellation increases super-linearly and can become very slow for a large numbers of points. Part of the challenge when developing `ppmData` was to make this step efficient by a C++ implementation of a Delaunay triangulation radial sweep algorithm [@sinclair_s-hull_2016]. The Dirichlet tessellations are then calculated as the dual graph of the Delaunay triangulation. Areal weights for the quadrature points are calculated as the area of each Dirichlet tessellation polygon.

## Generating a quasi-random quadrature scheme

In the `ppmData` package we have tried to make the generation of a quadrature scheme of a single or marked point process as easy as possible. We provide a simple interface for generating a Poisson point process data object and base all inputs and outputs using `terra` [@hijmans_terra_2021] and `sf` [@pebesma_simple_2018] packages. Users need to pass a set of coordinates from an observed point pattern, a window to define the point pattern region and a set of covariates observed across this region. The window is defined using a `terra` SpatRaster object, and it is used to identify the point pattern region and where to generate the quadrature locations. Covariates are included as a multiple layered `terra` SpatRaster object that share the same extent and resolution as the window. The values of the covariates will be extracted at the locations of the point pattern and quadrature scheme. Whilst not unique to this package, the `ppmData` package can also generate grid and pseudo-random quadrature schemes as presented in Fig. 1.

\begin{figure}

{\centering \includegraphics{paper_files/figure-latex/multiple species quad-1} 

}

\caption{Multiple species quadrature scheme. Here we can see the presences of five different snail species across Tasmania. For the multiple species quadrature there is a common set of background quadrature sites, but the object returns species specific weights.}\label{fig:multiple species quad}
\end{figure}

## Multiple species (marked) point process quadrature

One of the main reasons for creating the `ppmData` package was to develop quadrature schemes for multiple species (marked) point process. The multiple species quadrature scheme can be used in multiple species models available for use in the `ecomix` package [@woolley_ecomix_2022], or models like joint species distribution models [@ovaskainen_making_2011]. Here is a quick example of how we would set up a quasi-random quadrature scheme for multiple species (Fig. 2). The main difference is the column `speciesID` contains multiple species, `ppmData` recognises this internally and subsequently sets up a marked point process quadrature scheme. The quadrature scheme contains a common set of quadrature locations, but the quadrature weights are calculated on a species-specific basis $w_{ij}$, where $i$ is a site index and $j$ is a species index.

## Fitting a model using quadratures generated from `ppmData`

Here we demonstrate how to fit an Inhomogeneous Poisson Process Model using the `glm` function and a `ppmData` object. This approach is very similar to the method presented in @warton_poisson_2010, except we are using a quasi-random quadrature. We use `glm` to fit an Inhomogeneous Poisson Process Model for simplicity, but if one wanted to fit an Inhomogeneous Poisson Process Model using different statistical machinery, then one could -- examples are: penalized regressions using `glmnet` [@friedman_regularization_2010] to achieve a `maxent` analysis [@renner_equivalence_2013], generalized additive models [ @wood_generalized_2017] or statistical learning approaches [@hastie_elements_2009].


```r
#load relevant libraries
library(ppmData)
#read the data into R
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
preds <- rast(lst)
#Format the data (subset and scale) 
presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
preds <- scale(preds)
#Obtain quadrature points, variables and weights
ppmdata <- ppmData(npoints = 20000, presences = presences,
                   window = preds[[1]], covariates = preds)
```

We present a snippet of example code for fitting an Inhomogeneous Poisson Process Model using `glm` from the R `stats` package [@r_core_team_r_2013]. We just need to extract the data.frame from the `ppmData` object and to specify a formula for the model. The only extra step is making the formula response equal to `presence/weights` (which gives us $z$ in the @berman_approximating_1992 approximation). The right side of the formula will represent the covariates and their functional forms. In this example we specify independent 2$^{nd}$ degree polynomials for each covariate, except for the spatial coordinates X & Y which we include in the model as interacting 2$^{nd}$ degree polynomials. Once this is done, we define the weights in the `glm` function call and we now can fit an IPPM using `ppmData` and `glm`.


```r
ppp <- ppmdata$ppmData
form  <- presence/weights ~ poly(X,Y, degree = 2) + 
                            poly(max_temp_hottest_month, degree = 2) + 
                            poly(annual_mean_precip, degree = 2) + 
                            poly(annual_mean_temp, degree = 2) + 
                            poly(distance_from_main_roads, degree = 2)

ft.ppm <- glm(formula = form, data = ppp,
              weights = as.numeric(ppp$weights),
              family = poisson())
```

Finally, we show how to predict the model using `terra` and a `glm` object. This will return the expected intensity when all raster cell areas equal to one. To get the intensity of *Tasmaphena sinclairi* per cell, we need to re-scale the intensity based on the area represented by each raster cells. The re-scaled intensity returns the expected count of points (species presences) per raster cell. See Fig. 3 for an example output.

\begin{figure}

{\centering \includegraphics{paper_files/figure-latex/unnamed-chunk-2-1} 

}

\caption{Predicted intensity per for $\approx $ 0.6 km$^2$ sized raster cell for the snail species \textit{Tasmaphena sinclairi}. The points are occurrence records for \textit{Tasmaphena sinclairi} in the region.}\label{fig:unnamed-chunk-2}
\end{figure}

<!-- # Future work -->

<!-- One common problem when developing a quadrature is identifying the right number of quadrature points and identify locations which contribute the most information to log-likelihood. A common approach is to choose an arbitrary large number of quadrature locations [@renner_point_2015] or fit multiple models at different grid resolutions to assess log-likelihood convergence [@warton_poisson_2010]. One way to improve computation efficiency of quadrature schemes would be to use inclusion probabilities [@foster_spatially_2017]. There will be regions which do not contribute very much to the log-likelihood [@warton_advancing_2013]. These regions could be give sparser quasi-random samples, while areas which are likely to contribute larger amounts of information could be targeted. For our future work, we envisage an adaptive scheme, where we start with a relatively small number of quadrature points and update the scheme in regions of the spatial domain $\mathcal{A} \in R^2$ where the gradient of the integral is highly variable. Inclusion probabilities can be thus updated to increase the probability of a quasi-random sample in regions where the integral gradient (or other appropriate metrics) are complex. An adaptive approach would hone a quadrature scheme to improve computation efficient and target areas that will contribute more information to the log-likelihood. -->

# Acknowledgments

We thank Piers Dunstan for comments on an earlier draft of this manuscript.

# References
