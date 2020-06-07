# logdensity
Local polynomial estimates of the logarithm of the density and its derivatives. For details, please see Pinkse, J. and Schurter, K. (2020) "[Estimates of derivatives of (log) densities and related objects.](http://arxiv.org/abs/2006.01328)" Code is currently provided for R and Julia. Stata, Matlab, and Python are coming soon.

## Installation
You can manually download and install the code from this repository. The MATLAB function is defined in `logdensity.m`. Documentation for the MATLAB function can be accessed within MATLAB using `help logdensity` if the function file is in the current directory or in MATLAB's search path. A Julia function can be found in `logdensity.jl`. The R functions can be installed as a package from source by downloading all other files in the repository or by using the following R code.

    install.package('devtools') #If not already installed
    devtools::install_github("kschurter/logdensity")

## System requirements
For the convenience of Linux and Mac OS users, the R function `logdensity` parallelizes the computation if its argument `mc.cores` is greater than 1 by calling the R function `mcmapply`. This will not work with Windows, so make sure the argument `mc.cores` is equal to 1 (the default) if you are on a Windows machine. Refer to the package `parallel` for more information.
