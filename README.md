## BayesTwin
BayesTwin is an R package that can be used to perform bayesian analysis of twin data. 

Caution! This is work under progress! The master function (master.R) that calls the subroutines is still being developed. 

You can use subroutines (e.g. ace.R), but although most functions are based on publications of the package authors, we cannot guarantee that all scripts are bug free as the software is not yet fully tested. 


###Installation

install.packages("devtools")

library(devtools)

install_github("ingaschwabe/BayesTwin")

Note: An error might occur when (an earlier version) of the package is already installed on your system (especially when you're using a Mac).
In this case, deleting the package folder often solves the problem. The library path can be retrieved by calling the function .libPaths()