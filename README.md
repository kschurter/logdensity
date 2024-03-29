# logdensity
Local polynomial estimates of the logarithm of the density and its derivatives. For details, please see Pinkse, J. and Schurter, K. (2022) "[Estimates of derivatives of (log) densities and related objects.](https://doi.org/10.1017/S0266466621000529)" Code is provided for MATLAB, Python, R, and Stata. The [Julia package](https://github.com/NittanyLion/LogDensity.jl) has its own repository, and can be installed using Julia's package manager (e.g. type the character `]` in the REPL and then run the command `add LogDensity`.)

## Installation
You can manually download and install the code from this repository. The MATLAB function is defined in `logdensity.m`. Documentation for the MATLAB function can be accessed within MATLAB using `help logdensity` if the function file is in the current directory or in MATLAB's search path.

The Stata command is in `logdensity.ado`, with accompanying documentation in `logdensity.sthlp` and `logdensity_sthlp.pdf`.   `logdensity.ado` and `logdensity.sthlp` should be saved in the working directory or Stata's personal ado directory. The location of your personal ado directory depends on the version of Stata and your operating system. Run the command `sysdir` in Stata if you are not sure where your personal ado directory is. The pdf version of the Stata help file in `logdensity_sthlp.pdf` is just for the convenience of viewing outside of Stata and does not need to be downloaded.

A Python module can be found in `logdensity.py`. This module depends on NumPy and SciPy, which do not come with the standard Python installation. If you already have the NumPy and SciPy packages, you can simply save `logdensity.py` anywhere along Python's search path and get started with `import logdensity` and `help(logdensity)`.

The R functions can be installed as a package from source by downloading all other files in the repository or by using the following R code.

    install.packages('devtools') #If not already installed
    devtools::install_github("kschurter/logdensity")

## System requirements
For the convenience of Linux and Mac OS users, the R function `logdensity` parallelizes the computation if its argument `mc.cores` is greater than 1 by calling the R function `mcmapply`. This will not work with Windows, so make sure the argument `mc.cores` is equal to 1 (the default) if you are on a Windows machine. Refer to the R package `parallel` for more information.
