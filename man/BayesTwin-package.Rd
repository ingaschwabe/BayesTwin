\name{BayesTwin}
\alias{BayesTwin}
\docType{package}
\title{
    Bayesian analysis of twin data
}
\description{
    This package can be used to perform bayesian analysis of twin data. Using Gibbs sampling (Geman & Geman 1984;       Gelfand & Smith 1990; Gelman et al. 2004), twin data can be analysed on basis of an ACE or ADE model (Jinks & Fulker, 1970). It is also possible to estimate genotype by environment interaction. 
    
    When item data is available, it is also possible to analyse the data on the latent niveau. In that case, a measurement model (1 PL Rasch model, Rasch 1960) is incorporated into the biometric (ACE/ADE) model. 
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
   Inga Schwabe <i.schwabe@utwente.nl>, Stéphanie M. van den Berg
}
\references{
Jinks, J. L. & Fulker, D.W. (1970). Comparison of the biometrical genetical, mava, and classical approaches to the analysis of human behavior. Psychol Bull, 73(5), 311–349.

Gelfand A.E & Smith, A.F.M. (1990). Sampling-based approaches to calculating marginal densities. J Am Stat Assoc, 85(410), 398–409.

Gelman A., Carlin, J.B., Stern, H.S. & Rubin, D.B. (2004). Bayesian data analysis, 2nd edn. Chapman and Hall, London

Gelman A. & Rubin, D.B. (1992). Inference from iterative simulation using multiple sequences. Stat Sci 7(4), 457–511.

Geman S, Geman D. (1984) Stochastic relaxation, gibbs distributions and the bayesian restoration of images. IEEE Trans Pattern Anal Mach Intell 6(6), 721–741.

Rasch, G. (1960). Probabilistic models for some intelligence and attainment tests. Danish Institute or Educational Research, Copenhagen

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
