\name{BayesTwin}
\alias{BayesTwin}
\docType{package}
\title{Bayesian analysis of item-level twin data}
\description{
This package can be used to perform bayesian analysis of item-level twin data. Simulatenously with the biometric model (ACE/ADE or AE), an item response theory (IRT) model is estimated to take into account properties of the measurement scale. Functions are included in the package that help plot relevant information in figures and compute posterior statistics such as HPD intervals. 

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
    
    
The main function \code{IRT_twin} can be used to analyse item-level twin data under the ACE, ADE or AE model. Simultaneously with the biometric model, an item response theory (IRT) measurement model is estimated. For dichotomous item data, the 1 parameter model (1PL) or the 2 parameter model (2PL) can be used and for polytomous item data, the partial credit model (PCM) or the generalized partial credit model (GPCM). Optionally, genotype by environment interaction (GxE) can be estimated and/or covariates can be included in the analysis. GxE is assessed in the case that the unique environment features as latent (i.e., unmeasured) variable, using the parametrization as described in Schwabe & van den Berg (2014, 2016) and Schwabe, Boomsma & van den Berg (2017). The function returns MCMC samples as well as posterior means, standard deviations and 95\% HPD intervals for all variance components. Objects returned from the function are assigned the class "Bayestwin" which has its own \code{summary} method that can be used to summarize the output. The function \code{plotbayestwin} can be used to create sampling plots or plot posterior distributions and the function \code{geplot} to plot the 95\% credibility region of the GxE interaction effect.

The function \code{simulate_twin_data} can be used to simulate item-level twin data under all common biometric models (ACE/ADE or AE). For the simulation of the item data, a 1PL, 2PL, PCM or GPCM can be used. Optionally, the data can be simulated with GxE, using the parametrization as described in Schwabe & van den Berg (2014, 2016) and Schwabe, Boomsma & van den Berg (2017).
}
\author{
   Inga Schwabe <bayestwin@gmail.com>
}
\references{
Schwabe, I. & van den Berg, S.M. (2014). Assessing Genotype by Environment Interaction in Case of Heterogeneous Measurement Error, Behavior Genetics, 44 (4), 394-406.

Molenaar, D. & Dolan, C.V. (2014). Testing Systematic Genotype by Environment Interactions Using
Item Level Data, Behavior Genetics, 44(3), 212-231.

Schwabe, I., Jonker, W. & van den Berg, S.M. (2016). Genes, Culture and Conservatism - A Psychometric-Genetic Approach, Behavior Genetics, 46 (4), 516-52.

Schwabe, I., Boomsma, D.I. & van den Berg, S.M. (2017). Increased Environmental Sensitivity in High Mathematics Performance. Learning and Individual Differences, 54, 196-201.
}
\note{
To work, this packages requires the program JAGS to be in the PATH variable. JAGS can be freely obtained from http://mcmc-jags.sourceforge.net.
}
\examples{
data(results)
summary(results)

#Using the output to obtain the 95\% HPD for additive genetic variance 
HPD(results$samples_var_a)

#Using the output to obtain the 90\% HPD for all item difficulty parameters
apply(results$samples_item_b, 1, function (x) HPD(x, 0.90))

\dontrun{
##Simulate Item-level twin data under the 1PL Rasch model 
data = simulatetwin(irt_model = "1PL", var_a = 0.5, var_c = 0.3, ge_beta0 = log(0.2), 
                    ge = TRUE)
                
data_mz = data$y_mz
data_dz = data$y_dz

##Analyse the simulated data under an 1PL model with GxE
results = IRTtwin(data_mz, data_dz, 1:20, 21:40, ge = TRUE)

##Summarize results: 
summary(results)

#Using the output to obtain the 95\% HPD for additive genetic variance 
HPD(results$samples_var_a)

#Using the output to obtain the 90\% HPD for all item difficulty parameters
apply(results$samples_item_b, 1, function (x) HPD(x, 0.90))

##Plot trace lines for var(A)
plotbayestwin(results$samples_var_a, type = "trace")

##Plot posterior distribution of var(A)
plotbayestwin(results$samples_var_a)

##Plot 95\% credibility region of GxE interaction effect
geplot(results$var_a, results$samples_beta0, results$samples_beta1)
}
}