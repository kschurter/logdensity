# logdensity
Local polynomial estimates of the logarithm of the density and its derivatives. For details, please see Pinkse, J. and Schurter, K. (2020) "Estimates of derivatives of (log) densities and related objects."

## Installation
You can manually download and install the code from this repository or use

    install.package('devtools') #If not already installed
    devtools::install_github("kschurter/logdensity")

## System requirements
For Linux and Mac OS users, parallelization via forking is implemented using the R function `mcmapply`. This will not work with Windows, so make sure the argument `mc.cores` is equal to 1 (the default) if you are on a Windows machine.
