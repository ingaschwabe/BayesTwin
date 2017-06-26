\name{simulatetwin}
\alias{simulatetwin}
\alias{simulatetwindata}
\alias{Simulatedata}
\title{Simulate item-level twin data}
\description{
  Simulation of item-level twin data under the AE, ACE or ADE model with or without genotype by environment interaction (GxE).}
\usage{
simulatetwin(n_mz = 140, n_dz = 360, var_a = 0.5, var_c = 0.3, 
var_e = 0.2, var_d = 0, model = "ACE", n_items = 20, n_cat = 0,
ge = FALSE, ge_beta0 = 0, ge_beta1 = 0, irt_model = "1PL")
}
\arguments{
  \item{n_mz}{Number of monozygotic (MZ) twin pairs. Default value is 140.}
  \item{n_dz}{Number of dizygotic (DZ) twin pairs. Default value is 360.}
  \item{var_a}{Variance explained by additive genetic influences. Default value is 0.5.}
  \item{var_c}{Variance explained by shared environmental influences. Default value is 0.3.}
  \item{var_e}{Variance explained by unique environmental influences. Default value is 0.2.}
  \item{var_d}{Variance explained by dominance genetic influences. As the default model is an "ACE"
  model, the default value is set to 0.}
  \item{model}{Biometric model used for the simulation, either "ACE", "ADE" or "AE". Default model is "ACE" model.}
  \item{n_items}{Number of test items. Default value is 20.}
  \item{n_cat}{Number of categories for the simulation of polytomous item data. As the default IRT model is a 1PL, the default is set 0.}
  \item{ge}{Indication whether data is simulated with (ge = TRUE) or without (ge = FALSE) GxE. Default value is FALSE}
  \item{ge_beta0}{When data is simulated with ge = TRUE, average environmental variance (when A = 0) has to be defined. Default value is 0.}
  \item{ge_beta1}{When data is simulated with ge = TRUE, the GxE interaction parameter has to be defined. Default value is 0.}
  \item{irt_model}{Choice of the measurement (IRT) model under which the item data is simulated. "1PL" = 1 PL model, "2PL" = 2PL model, "PCM" = Partial credit model, "GPCM" = Generalized partial credit model.}
}
\details{
 Item twin data is simulated under all common biometric models (i.e., A(C)E, ADE) with or without genotype by environment interaction (GxE). Data is simulated with an average phenotypic value of 0. For the simulation of dichotomous item data, the 1 PL or 2 PL can be used. For polytomous item data with three or more categories, the partial credit model (PCM) or the generalized partial credit model (GPCM) can be used. 
 
Genotype by environment interaction is defined modelled as heterogeneous unique-environmental variance and is defined as \eqn{\sigma^2Eij = exp (\beta0 + \beta1 Aij)} under the ACE model and as \eqn{\sigma^2Eij = exp (\beta0 + \beta1 Gij)} under the ADE model where \eqn{\beta0} refers to an intercept (unique-environmental variance when \eqn{Aij} or respectively \eqn{Gij = 0}) and \eqn{\beta1} represents the GxE interaction effect. \eqn{Aij} represents the genotypic value for twin \eqn{j} from family \eqn{i}, likewise \eqn{G} represents the total genotypic value (both additive and non-additive genetic effects) for twin \eqn{j} from family \eqn{i}. For more detail see Schwabe & van den Berg (2014, 2016) or Schwabe, Boomsma & van den Berg (2017).
}

\value{
Two matrices are returned, \code{y_mz} for MZ twin pairs and \code{y_dz} for DZ twin pairs. 
For a total of \code{n_items}, the result is a matrix of \code{n_mz} rows for the i-th MZ family and 2*\code{n_items} columns with the item answers of the first (columns 1:\code{n_items}) and second twin (columns \code{n_items} + 1:2*\code{n_items}) of a family. For example, \code{y_mz[1,22]} is the response of the first twin from family 1 to item 22 and \code{y_mz[1,23]} is the response of the seocnd twin2 from family 1 to item 1 if \code{n_items} = 22. The same logic applies to the data of DZ families.
}

\author{Inga Schwabe}
\references{
Schwabe, I. & van den Berg, S.M. (2014). Assessing Genotype by Environment Interaction in Case of Heterogeneous Measurement Error, Behavior Genetics, 44 (4), 394-406.

Schwabe, I., Jonker, W. & van den Berg, S.M. (2016). Genes, Culture and Conservatism - A Psychometric-Genetic Approach, Behavior Genetics, 46 (4), 516-52.

Schwabe, I., Boomsma, D.I. & van den Berg, S.M. (2017). Increased Environmental Sensitivity in High Mathematics Performance. Learning and Individual Differences, 54, 196-201.
}
\examples{
#100 MZ twins, 200 DZ twins, ACE model, no GxE, 1 PL: 
simulatetwin(n_mz = 100, n_dz = 200) 

\dontrun{
#500 MZ twins, 800 DZ twins, ACE model, no GxE, 1 PL: 
simulatetwin(n_mz = 500, n_dz = 800) 

#140 MZ twins, 360 DZ twins, ADE model, GPCM: 
simulatetwin(var_a = 0.4, var_d = 0.2, var_e = 0.4, 
             model = "ADE", irt_model = "GPCM")
}
}
\keyword{utilities}