## BayesTwin
BayesTwin is an R package that can be used for Bayesian inference of item-level twin data. Simultaneously with the biometric model, an item response theory (IRT) measurement model is estimated. For dichotomous item data, the 1 parameter model (1PL) or the 2 parameter model (2PL) can be used and for polytomous item data, the partial credit model (PCM) or the generalized partial credit model (GPCM). Functions are included that help plot relevant information in figures and compute posterior statistics such as HPD intervals. 
 
_Caution! To use this package, you need to install the MCMC program JAGS which can be freely obtained at http://mcmc-jags.sourceforge.net._

You can find an example analysis using simulated data in the vignette map. 
Depending on the sample-size and number of items included in the analysis, the analysis is computationally intensive and can take several hours to complete. In this case, it is recommended to run the analysis on a computer cluster (see vignette/computer cluster analysis for an example analysis).

### Installation

#### Stable CRAN version
Install JAGS (http://mcmc-jags.sourceforge.net) and R (https://www.r-project.org/).

Then run in R:

```r
install.packages("BayesTwin")

library(BayesTwin)
```

#### Developmental Github version
Install JAGS (http://mcmc-jags.sourceforge.net) and R (https://www.r-project.org/).

Then run in R:

```r
install.packages("devtools")

library(devtools)

install_github("ingaschwabe/BayesTwin")
```

Note: An error might occur when (an earlier version) of the package is already installed on your system (especially when you're using a Mac). In this case, deleting the package folder often solves the problem. The library path can be retrieved by calling the function .libPaths() 
