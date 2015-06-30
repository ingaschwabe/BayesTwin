\name{simulate_twin_data}
\title{Simulate twin data}
\description{
  Simulation of item twin data under the AE, ACE or ADE model.
  Optional with genotype by environment interaction (GxE) and/or environmental covariates.
}
\usage{
simulate_twin_data(n_mz, n_dz, var_a = 0.5, var_c = 0.3, 
var_e = 0.2, var_d = 0, model = "ACE", n_items = 0, n_var = 0, 
ge = FALSE, ge_beta0 = log(0.5), ge_beta1 = 1.5)
}
\arguments{
  \item{n_mz}{Number of monozygotic (MZ) twin pairs}
  \item{n_dz}{Number of dizygotic (DZ) twin pairs}
  \item{var_a}{Variance explained by additive genetic influences. Default value is 0.5.}
  \item{var_c}{Variance explained by shared environmental influences. Default value is 0.3.}
  \item{var_e}{Variance explained by unique environmental influences. Default value is 0.2.}
  \item{var_d}{Variance explained by dominance genetic influences. As the default model is an "ACE""
  model, the default value is set to 0.}
  \item{model}{Biometric model used for the simulation. Default model is "ACE" model.}
  \item{n_items}{Number of test items. When set above 0, 

}
\details{
  For each parameter the interval is constructed from the empirical cdf
  of the sample as the shortest interval for which the difference in
  the ecdf values of the endpoints is the nominal probability.  Assuming
  that the distribution is not severely multimodal, this is the HPD interval.
}
\value{
  For an \code{"mcmc"} object, a matrix with columns \code{"lower"} and
  \code{"upper"} and rows corresponding to the parameters.  The
  attribute \code{"Probability"} is the nominal probability content of
  the intervals.  A list of such matrices is returned for an
  \code{"mcmc.list"} object.
}
\author{Douglas Bates}
\examples{
data(line)
HPDinterval(line)
}
\keyword{univar}
\keyword{htest}