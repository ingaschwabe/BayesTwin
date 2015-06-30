\name{simulate_twin_data}
\title{Simulate twin data}
\description{
  Simulation of item twin data under the AE, ACE or ADE model.
  Optional with genotype by environment interaction (GxE) and/or environmental covariates.
}
\usage{
simulate_twin_data(n_mz, n_dz, var_a = 0.5, var_c = 0.3, 
var_e = 0.2, var_d = 0, model = "ACE", n_items = 20, n_var = 0, 
ge = FALSE, ge_beta0 = log(0.5), ge_beta1 = 1.5)
}
\arguments{
  \item{n_mz}{Number of monozygotic (MZ) twin pairs}
  \item{n_dz}{Number of dizygotic (DZ) twin pairs}
  \item{var_a}{Variance explained by additive genetic influences. Default value is 0.5.}
  \item{var_c}{Variance explained by shared environmental influences. Default value is 0.3.}
  \item{var_e}{Variance explained by unique environmental influences. Default value is 0.2.}
  \item{var_d}{Variance explained by dominance genetic influences. As the default model is an "ACE"
  model, the default value is set to 0.}
  \item{model}{Biometric model used for the simulation. Default model is "ACE" model.}
  \item{n_items}{Number of test items. Default value is 20.}
  \item{n_var}{Number of environmental covariates. Default number is 0.}
  \item{ge}{Indication whether data is simulated with (ge = TRUE) or without (ge = FALSE) GxE. Default value is FALSE}
  \item{ge_beta0}{When data is simulated with ge = TRUE, average environmental variance (when A = 0) has to be defined. Default value is log(0.5)}
  \item{ge_beta1}{When data is simulated with ge = TRUE, the GxE interaction parameter has to be defined. Default value is 1.5}
}
\details{
 Item twin data is simulated under different biometric models ("AE", "ACE" or "ADE") with or without genotype environment interaction. Optionally, environmental covariates can be simulated. Data is simulated with an average phenotypic value of 0. The 1PL IRT model is used as a measurement model to estimate item data. 
}
\value{Two matrices are returned, \code{y_mz} for MZ twin pairs and \code{y_dz} for DZ twin pairs. 
For a total of \code{n_items}, the result is a matrix of \code{n_mz} rows for the i-th MZ family and 2*\code{n_items} columns with the item answers of the first (columns 1:\code{n_items}) and second twin (columns \code{n_items} + 1:2*\code{n_items}) of a family. For example, \code{y_mz[1,22]} is the response of the first twin from family 1 to item 22 and \code{y_mz[1,23]} is the response of the seocnd twin2 from family 1 to item 1 if \code{n_items} = 22. The same logic applies to the data of DZ families.}
\author{Inga Schwabe}
\examples{
simulate_twin_data(n_mz = 500, n_dz = 800)
}
\keyword{univar}