\name{BayesTwin}
\alias{BayesTwin}
\docType{package}
\title{
    Bayesian analysis of item-level twin data
}
\description{
This package can be used to perform bayesian analysis of item-level twin data. Functions are included that help plot relevant information in figures and compute posterior statistics such as HPD intervals. 

Caution! The subroutines of this package rely on the program JAGS, which can be freely obtained from http://mcmc-jags.sourceforge.net.
}
\details{
    \tabular{ll}{
        Package: \tab BayesTwin\cr
        Type: \tab Package\cr
        Version: \tab 0.1.0\cr
        Date: \tab 2017-01-06\cr
        License: \tab GPL-2\cr
    }
    
    
The main function \code{IRT_twin} can be used to analyse item-level twin data under the ACE, ADE or AE model. Simultaneously with the biometric model, an item response theory (IRT) measurement model is estimated. For dichotomous item data, the 1 parameter model (1PL) or the 2 parameter model (2PL) can be used and for polytomous item data, the partial credit model (PCM) or the generalized partial credit model (GPCM). Optionally, genotype by environment interaction (GxE) can be estimated and/or covariates can be included in the analysis. GxE is assessed in the case that the unique environment features as latent (i.e., unmeasured) variable, using the parametrization as described in Schwabe & van den Berg (2014, 2016) and Schwabe, Boomsma & van den Berg (2017). The function returns MCMC samples as well as posterior means, standard deviations and 95\% HPD intervals for all variance components. Objects returned from the function are assigned the class "Bayestwin" which has its own \code{plot} method that can be used to create sampling plots or plot posterior distributions. 

The function \code{simulate_twin_data} can be used to simulate item-level twin data under all common biometric models (ACE/ADE or AE). For the simulation of the item data, a 1PL, 2PL, PCM or GPCM can be used. Optionally, the data can be simulated with GxE, using the parametrization as described in Schwabe & van den Berg (2014, 2016) and Schwabe, Boomsma & van den Berg (2017).
}
\author{
   Inga Schwabe <bayestwin@gmail.com>
}
\references{
Schwabe, I. & van den Berg, S.M. (2014). Assessing Genotype by Environment Interaction in Case of Heterogeneous Measurement Error, Behavior Genetics, 44 (4), 394-406.

Schwabe, I., Jonker, W. & van den Berg, S.M. (2016). Genes, Culture and Conservatism - A Psychometric-Genetic Approach, Behavior Genetics, 46 (4), 516-52.

Schwabe, I., Boomsma, D.I. & van den Berg, S.M. (2017). Increased Environmental Sensitivity in High Mathematics Performance. Learning and Individual Differences, 54, 196-201.
}
\note{
To work, this packages requires the program JAGS to be in the PATH variable. JAGS can be freely obtained from http://mcmc-jags.sourceforge.net.
}
\examples{
\dontrun{
#Simulate twin data under an ACE model
simulate_twin_data(nmz = 300, ndz = 200)
}
}