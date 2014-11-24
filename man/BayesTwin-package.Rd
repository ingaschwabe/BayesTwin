\name{BayesTwin}
\alias{BayesTwin}
\docType{package}
\title{
    Bayesian analysis of twin data
}
\description{
This package can be used to perform bayesian analysis of twin data. Using Gibbs sampling (Geman & Geman 1984;       Gelfand & Smith 1990; Gelman et al. 2004), twin data can be analysed on basis of an ACE or ADE model (Jinks & Fulker, 1970). When item data is available, it is also possible to analyse the data on the latent niveau. In that case, a measurement model (1 PL Rasch model, Rasch 1960) is incorporated into the biometric (ACE/ADE) model. 
}
\details{
    \tabular{ll}{
        Package: \tab BayeTwin\cr
        Type: \tab Package\cr
        Version: \tab 0.1.0\cr
        Date: \tab 2014-11-09\cr
        License: \tab GPL-2\cr
    }
The main functions can be used to .. With \code{simulate_twin_data} you can simulate twin data. 
}
\author{
   Inga Schwabe <i.schwabe@utwente.nl>, Stephanie M. van den Berg
}
\references{
Jinks, J. L. & Fulker, D.W. (1970). Comparison of the biometrical genetical, mava, and classical approaches to the analysis of human behavior. Psychol Bull, 73(5), 410--423.
}
\note{
    Analysis of genotype by environment interaction requires the program "JAGS" (Just Another Gibbs Sampler) to be in the PATH variable. This program can be obtained from http://mcmc-jags.sourceforge.net.
}
\examples{
\dontrun{
    #Simulate twin data
    simulate_twin_data(nmz = 300, ndz = 200)
        }
}
