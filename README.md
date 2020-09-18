# Nonstat_Cov

# Estimation of Spatial Deformation for Nonstationary Processes via Variogram Alignment
<img align="right" src="https://github.com/ghulamabdul/Nonstat_Cov/blob/master/p1cover.png" alt="drawing" width="400"/>


This repository provides reproducible code for the manuscript titled *Estimation of Spatial Deformation for Nonstationary Processes via Variogram Alignment* by Ghulam A. Qadir, Sebastian Kurtek, and Ying Sun. The manuscript describes a new approach to model nonstationary covariance function of a univariate spatial process.

## Abstract

In modeling spatial processes, a second-order stationarity assumption is often made. However, for spatial data observed on a vast domain, the covariance function often varies over space, leading to a heterogeneous spatial dependence structure, therefore requiring nonstationary modeling. Spatial deformation is one of the main methods for modeling nonstationary processes, assuming the nonstationary process has a stationary counterpart in the deformed space. The estimation of the deformation function poses severe challenges. Here, we introduce a novel approach for nonstationary geostatistical modeling, using space deformation, when a single realization of the spatial process is observed. Our method is based on aligning regional variograms, where warping variability of the distance from each subregion explains the spatial nonstationarity. We propose to use multi-dimensional scaling to map the warped distances to spatial locations. We assess the performance of our new method using multiple simulation studies. Additionally, we illustrate our methodology on precipitation data to estimate the heterogeneous spatial dependence and to perform spatial predictions.
## Requirements

The codes are written in R, and reproducing would require installing and loading the following R-packages: `fields`,`fdasrvf`,`geoR`,`plot3D`,`rgl`,`scatterplot3d`,`scoringRules`,`doParallel`,`mvtnorm`,`ggplot2`,`ggmap`, `viridis`,`rgdal` and `viridis`. 

##
