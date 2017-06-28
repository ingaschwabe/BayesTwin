## BayesTwin
BayesTwin is an R package that can be used for Bayesian inference of item-level twin data. Simultaneously with the biometric model, an item response theory (IRT) measurement model is estimated. For dichotomous item data, the 1 parameter model (1PL) or the 2 parameter model (2PL) can be used and for polytomous item data, the partial credit model (PCM) or the generalized partial credit model (GPCM). Functions are included that help plot relevant information in figures and compute posterior statistics such as HPD intervals. 
 
Caution! To use this package, you need to install the MCMC program JAGS which can be freely obtained at http://mcmc-jags.sourceforge.net.

For an example analysis, see http://www.ingaschwabe.com/Exampleanalysis.html

### Installation

#### Stable CRAN version
Install JAGS (http://mcmc-jags.sourceforge.net) and R (https://www.r-project.org/).

Then run in R:

install.packages("BayesTwin")

library(BayesTwin)

#### Developmental Github version
Install JAGS (http://mcmc-jags.sourceforge.net) and R (https://www.r-project.org/).

Then run in R:

install.packages("devtools")

library(devtools)

install_github("ingaschwabe/BayesTwin")

Note: An error might occur when (an earlier version) of the package is already installed on your system (especially when you're using a Mac). In this case, deleting the package folder often solves the problem. The library path can be retrieved by calling the function .libPaths() 
