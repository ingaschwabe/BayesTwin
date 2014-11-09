\name{BayesTwin}
\alias{BayesTwin}
\docType{package}
\title{
    Bayesian twin analaysis.
}
\description{
    This package can be used to perform bayesian analysis of twin models (ACE, ADE). It can be used for sumscore data as well as item data.
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
   Inga Schwabe <i.schwabe@utwente.nl>, Stephanie M. van den Berg.
}
\references{
}
\note{
    This function requires the program "JAGS" (Just Another Gibbs Sampler) to be in the PATH variable. This program can be obtained from http://mcmc-jags.sourceforge.net.
}
\examples{
    \dontrun{
        #simulate data
        simulate_twin_data()
    }
}
