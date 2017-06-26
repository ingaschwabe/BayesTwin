\name{IRTtwin}
\alias{IRTtwin}
\title{Bayesian analysis of item-level twin data}
\description{
    Bayesian analysis of item-level twin data under the AE, ACE or ADE model. Covariates can be included an a GxE interaction effect can be estimated.}
\usage{IRTtwin(data_mz, data_dz, 
         twin1_datacols_p, twin2_datacols_p, 
         twin1_datacols_cov = NA, twin2_datacols_cov = NA, 
         decomp_model = "ACE", irt_model = "1PL", ge = FALSE, 
         n_iter = 5000, n_burnin = 7000, n_chains = 1, 
         fit_stats = FALSE, var_prior = "INV_GAMMA", 
         N_cov = 0, inits = NA, Nk = 0)
}
\arguments{
\item{data_mz}{Phenotypic as well as covariate data of MZ twins in matrix form.}
\item{data_dz}{Phenotypic as well as covariate data of DZ twins in matrix form.}
\item{twin1_datacols_p}{Columns of data_mz or respectively data_dz in which phenotypic data is stored for the     first twin of all MZ and DZ twins, e.g. c(1,2,3) if they are stored in the first three columns or c(20:40) if stored in column 20 through 40.}
\item{twin2_datacols_p}{Columns of data_mz or respectively data_dz in which phenotypic data is stored for the second twin of all MZ and DZ twins.}
\item{twin1_datacols_cov}{Columns of data_mz or respectively data_dz in which covariate data is stored for the first twin of all MZ and DZ twins.}
\item{twin2_datacols_cov}{Columns of data_mz or respectively data_dz in which covariate data is stored for the second twin of all MZ and DZ twins.}
\item{decomp_model}{Choice of the biometric model under which the data is analyzed, either "ACE", "ADE" or "AE". Default model is "ACE" model.}
\item{irt_model}{Choice of the measurement (IRT) model under which the item data is simulated. "1PL" = 1 PL model, "2PL" = 2PL model, "PCM" = Partial credit model, "GPCM" = Generalized partial credit model..}
\item{ge}{Indication whether data is analyzed with (ge = TRUE) or without (ge = FALSE) GxE. Default value is FALSE.}
\item{n_iter}{Total number of MCMC iterations (without burn-in period). Default number is 5.000}
\item{n_burnin}{Total number of burn-in iterations for the MCMC analysis. Default number is 7.000}
\item{n_chains}{Number of MCMC chains. Default number is 1.}
\item{fit_stats}{Should be set to TRUE when fit statistics are needed (e.g., DIC, BIC). In this case, multiple chains need to be used (the parameter n_chains should be set to at least 2)}
\item{var_prior}{Choice of prior distributions used for the variance compoenents, either "INV_GAMMA" for an inverse gamma distribution or "UNIF" for uniform distribution. For a low total phenotypic variance, uniform prior distributions are recommended. Default value is "INV_GAMMA".}
\item{N_cov}{Number of covariates. Default value is 0.}
\item{inits}{Initial values for the MCMC analysis. Default is NULL.}
\item{Nk}{Number of categories in case of polytomous item data. As the default IRT model is a 1PL, the default value is set 0.}


}
\details{Bayesian analysis of item-level twin data under all common biometric models (i.e., A(C)E, ADE) with or without genotype by environment interaction (GxE). For the simulation of dichotomous item data, the one parameter logistic model (1PL) or the 2 parameter logistisc model (2PL) can be used. For polytomous item data, the partial credit model (PCM) or the generalized partial credit model (GPCM) can be used. Optionally, covariates can be included in the analysis.
    
GxE is defined modelled as heterogeneous unique-environmental variance and is defined as \eqn{\sigma^2Eij = exp (\beta0 + \beta1 Aij)} under the ACE model and as \eqn{\sigma^2Eij = exp (\beta0 + \beta1 Gij)} under the ADE model where \eqn{\beta0} refers to an intercept (unique-environmental variance when \eqn{Aij} or respectively \eqn{Gij = 0}) and \eqn{\beta1} represents the GxE interaction effect. \eqn{Aij} represents the genotypic value for twin \eqn{j} from family \eqn{i}, likewise \eqn{G} represents the total genotypic value (both additive and non-additive genetic effects) for twin \eqn{j} from family \eqn{i}. For more detail see Schwabe & van den Berg (2014), Behavior Genetics, 44(4), 394-406 and Schwabe, Boomsma & van den Berg (2017), Learning and Individual Differences, 54, 196-201. 
}

\value{
Output of the class "bayestwin" including samples, posterior point estimates, standard deviations and the 95\% HPD interval for variance components and, if applicable, regression parameters. Call str(results) if the output of \code{IRT_twin} is stored in an object called results to see a list of all stored variables. 

Results are printed on the fly for variance components and, if applicable, regression parameters. 
}

\author{Inga Schwabe}

\examples{
data(results)

#Summarize results
summary(results)

##Plot posterior distribution of var(A)
plotbayestwin(results$samples_var_a)

##Plot 95\% credibility region of GxE interaction effect
geplot(results$var_a, results$samples_beta0, results$samples_beta1)

\dontrun{
##Simulate data
data = simulate(irt_model = "1PL", var_a = 0.5, var_c = 0.3, ge_beta0 = log(0.2), 
                ge = TRUE)
data_mz = data$y_mz
data_dz = data$y_dz

##Run analysis
results = IRTtwin(data_mz, data_dz, 1:20, 21:40, ge = TRUE)

##Summarize results: 
summary(results)

##Plot trace lines for var(A)
plotbayestwin(results$samples_var_a, type = "trace")

##Plot posterior distribution of var(A)
plotbayestwin(results$samples_var_a)

##Plot 95\% credibility region of GxE interaction effect
geplot(results$var_a, results$samples_beta0, results$samples_beta1)
}
}
\keyword{IRTtwin}