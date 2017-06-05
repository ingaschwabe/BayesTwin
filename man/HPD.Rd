\name{HPD}
\alias{hpd}
\alias{HPD}
\title{Calculate highest posterior density interval}
\description{This function calculates the Bayesian highest posterior density interval (HPD) based on a parameters' posterior sample.}
\usage{
HPD(sample, cred_int = 0.95)
}
\arguments{
  \item{sample}{A vector representing draws from the target distribution of the paramter of interest, as produced by the main function \code{IRT_twin} of this package.}
  \item{cred_int}{The desired accuracy of the HPD. Default value is 0.95 for 95\%.}
}
\details{
The highest posterior density interval (HPD, see e.g. Box & Tia, 1992) contains the required mass such that all points within the interval have a higher probability density than points outside of the interval. 

The function expects as input a vector representing draws from the target distribution of the paramter of interest, such as produced by the main function \code{IRT_twin} of this package.

The result is a vector consisiting of two values, the first value representing the lower bound of the HPD and the second value representing the upper bound. 
}

\value{
A vector of length 2 with the lower (first value) and upper (second value) bound of the HPD.
}

\author{Inga Schwabe}
\references{
Box, G., & Tiao, G. (1992). Bayesian inference in statistical analysis. New York: John Wiley & Sons.
}
\examples{\dontrun{
#Simulate twin data under an ACE and 1PL model: 
data = simulate_twin_data(irt_model = "ACE"", model = "1PL"", n_items = 20)
data_mz = data$y_mz
data_dz = data$y_dz

#Analyse the data using the IRT_twin function: 
results = IRT_twin(data_mz = data$y_mz, data_dz = data$y_dz,
                  twin1_datacols_p = c(1:20), twin2_datacols_p = c(21:40))

#Using the output to obtain the 95\% HPD for additive genetic variance 
HPD(results$samples_var_a)

#Using the output to obtain the 90\% HPD for all item difficulty parameters
apply(results$samples_item_b, 1, function (x) HPD(x, 0.90))
}
}
\keyword{IRT_twin}