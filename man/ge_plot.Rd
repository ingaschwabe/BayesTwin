\name{ge_plot}
\alias{ge_plot}
\alias{GE_plot}
\alias{geplot}
\title{Plots 95\% credibility interval of GxE interaction}
\description{This function plots the 95\% credibility region of a GxE interaction effect.}
\usage{
ge_plot(var_a, samples_beta0, samples_beta1, main, xlab, ylab, col, ...)}
\arguments{
  \item{var_a}{A single value representing the posterior point estimate for genetic variance}
 \item{samples_beta0}{A vector representing draws from the posterior distribution of beta0, as produced by the main function \code{IRT_twin} of this package}
\item{samples_beta1}{A vector representing draws from the posterior distribution of beta1, as produced by the main function \code{IRT_twin} of this package} 
\item{main}{An overall title for the plot}
\item{xlab}{A title for the x axis}
\item{ylab}{A title for the y axis}
\item{col}{Color of the lines}
\item{...}{Further arguments for the default S3 plot method}
}
\details{The 95 \% credibility region for the GxE interaction effect is displayed for the entire
range of estimated genotypic values in the interval of -2 and 2 standard deviations of the posterior point estimate of genetic variance. 

The function expects as input the posterior point estimate for genetic variance and two vectors representing draws from the posterior distributions of beta0 and beta1 respectively, as produced by the main function \code{IRT_twin} of this package.
}
\value{
Plot of the 95\% credibility region of the GxE interaction.
}
\author{Inga Schwabe}
\examples{\dontrun{
#Simulate twin data under an ACE and 1PL model: 
data = simulate_twin_data(irt_model = "ACE"", model = "1PL"", n_items = 20, ge = TRUE)
data_mz = data$y_mz
data_dz = data$y_dz

#Analyse the data using the IRT_twin function: 
results = IRT_twin(data_mz = data$y_mz, data_dz = data$y_dz,
                  twin1_datacols_p = c(1:20), twin2_datacols_p = c(21:40))

#Using the output to plot the 95\% credibility region 
ge_plot(results$var_a, results$samples_beta0, results$samples_beta1)
}
}
\keyword{IRT_twin}