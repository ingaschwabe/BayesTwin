<<<<<<< HEAD
\name{BayesTwin}
\alias{BayesTwin}
\docType{package}
\title{
    Bayesian analysis of item-level twin data
}
\description{
This package can be used to perform bayesian analysis of item-level twin data. Simultaneously to the biometric model, an item response theory (IRT) measurement model is estimated. For dichotomous item data, the 1 parameter model (1PL) or the 2 parameter model (2PL) can be used and for polytomous item data, the generalized partial credit model (PCM) or the partial credit model (GPCM). Functions are included that help plot relevant information in figures and compute posterior statistics such as HPD intervals. 

Caution! The subroutines of this package rely on the program JAGS (Just Another Gibbs Sampler). JAGS can be freely obtained from http://mcmc-jags.sourceforge.net.
}
\details{
    \tabular{ll}{
        Package: \tab BayeTwin\cr
        Type: \tab Package\cr
        Version: \tab 0.1.0\cr
        Date: \tab 2017-20-06\cr
        License: \tab GPL-2\cr
    }
    
    
The main function \code{IRT_twin} can be used to analyse item-level twin data under the ACE, ADE or AE model. For dichotomous data, the 1PL or 2PL can be used and for polytomous item data the GPCM or PCM. Optionally, a genotype by environment interaction (GxE) can be estimated and/or covariates can be included. The function returns MCMC samples as well as posterior means, standard deviations and 95\% HPD intervals for all variance components and item parameters. Objects returned from the function are assigned the class "Bayestwin" for which the \code{plot} function can be used to create sampling plots or plot posterior distributions. 

The function \code{simulate_twin_data} can be used to simulate item-level twin data under all common biometric models with or without GxE, using an 1PL/2PL/PCM or GPCM. 
}
\author{
   Inga Schwabe <bayestwin@gmail.com>
}
\references{
Schwabe, I. & van den Berg, S.M. (2014). Assessing Genotype by Environment Interaction in Case of Heterogeneous Measurement Error, Behavior Genetics, 44 (4), 394-406.

Schwabe, I., Jonker, W. & van den Berg, S.M. (2016). Genes, Culture and Conservatism - A Psychometric-Genetic Approach, Behavior Genetics, 46 (4), 516-52.

Schwabe, I., Boomsma, D.I. & van den Berg, S.M. (2017). Increased Environmental Sensitivity in High Mathematics Performance, Learning and Individual Differences, 54, 196-201.
}
\note{
To work, this packages requires the program "JAGS" (Just Another Gibbs Sampler) to be in the PATH variable. JAGS can be freely obtained from http://mcmc-jags.sourceforge.net.
}
\examples{
\dontrun{
    #Simulate twin data
    simulate_twin_data(nmz = 300, ndz = 200)
        }
}
=======
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
>>>>>>> 580e32ad986f43fb99925137796f29266474c179
