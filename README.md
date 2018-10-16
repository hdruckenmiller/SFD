# SFD

Code for computing the spatial first differences (SFD) estimator. SFD is a research design that exploits the spatial structure of the data to address unobservable heterogeneity in cross section. Please see Druckenmiller & Hsiang, 2018. 

# Installation 
Currently, this package exists in a development version on GitHub. To use the package, you need to install it directly from GitHub using the `install_github` function from devtools.
```
library(devtools)
install_github("hdruckenmiller/SFD")
library(SFD)
```

# Usage 
The package can be used to compute spatial first differences (SFD) or spatial double differences (SDD). 
For example, you can compute SFD in your variables of interest using the `SFD_vars` function.
```
SFD_data <- SFD_vars(spatial_df = us.map, 
    n_channel = 50, 
    obs_var = "fips", 
    dependent_var = "logYield", 
    independent_vars = c("temp", "prec"))
```
The above code calculates SFD for the variables `logYield`, `temp`, and `prec`. 
The inputs are a SpatialPolygonsDataFrame `us.map`, a unique identifer for the observational unit `fips`, 
and a specification for the number of sampling channels `50`. 
The number sampling channels should be chosen such that the height of your sampling channels 
is approximately equal to the height of your observational units. 
The `n_channels` function inputs your SpatialPolygonsDataFrame and provides a recommendation for the number of sampling channels. 
```
n_channels <- SFD_vars(spatial_df = us.map)
```
Alternatively, the user can directly compute the SFD estimator using the function `SFD_lm`. This function performs a linear regression on the data after taking SFD. 
```
sfdmodel <- SFD_lm(spatial_df = us.map, 
    n_channel = 50, 
    obs_var = "fips", 
    dependent_var = "logYield", 
    independent_vars = c("temp", "prec"), 
    plot = TRUE)
```
The `plot=TRUE` option produces a plot the adjacent observations used to calculate SFD in your sample. 
We recommend checking this plot to make sure the algorithm identified reasonable neighbors. 
