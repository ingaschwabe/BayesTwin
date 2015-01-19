\name{BayesTwin}
\alias{BayesTwin}
\docType{package}
\title{
    Bayesian analysis of twin data
}
\description{
This package can be used to perform bayesian analysis of twin data. he most common univariate and multivariate biometric models are included, as well as a method to study concordance rates. In addition, it includes functions that help plot relevant information in figures and compute posterior statistics such as HPD intervals.
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
}
\note{
    Most subroutines require the program "JAGS" (Just Another Gibbs Sampler) to be in the PATH variable. This program can be obtained from http://mcmc-jags.sourceforge.net.
}
\examples{
\dontrun{
    #Simulate twin data
    simulate_twin_data(nmz = 300, ndz = 200)
        }
}
