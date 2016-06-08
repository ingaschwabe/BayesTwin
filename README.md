## BayesTwin
BayesTwin is an R package that can be used to perform bayesian analysis of twin data. The most common univariate and multivariate biometric models are included, as well as a method to study concordance rates. In addition, it includes functions that help
plot relevant information in figures and compute posterior statistics such as 95% HPD intervals. Some models are directly implemented in R, but for most subroutines you will need the MCMC software package JAGS, which is freely obtainable at http://mcmc-jags.sourceforge.net.

Caution! This is work under progress! 

You can use subroutines (e.g. ACE decomposition with GxE and IRT model: ge_irt.R), but although most functions are based on publications of the package authors, we cannot guarantee that all scripts are bug free as the software is not yet fully tested. All subroutines can be found in the folder R. 


###Installation
Install JAGS (http://mcmc-jags.sourceforge.net) and Rtools (https://cran.r-project.org/bin/windows/Rtools/index.html). 
Then run in R:

install.packages("devtools")
library(devtools)
install_github("ingaschwabe/BayesTwin")

Note: An error might occur when (an earlier version) of the package is already installed on your system (especially when you're using a Mac). In this case, deleting the package folder often solves the problem. The library path can be retrieved by calling the function .libPaths() 
