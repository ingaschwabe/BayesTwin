\name{geplot}
\alias{geplot}
\alias{GEplot}
\alias{ge_plot}
\title{Plots 95\% credibility interval of GxE interaction}
\description{This function plots the 95\% credibility region of a GxE interaction effect.}
\usage{
geplot(var_a, samples_beta0, samples_beta1, main, xlab, ylab, col, ...)}
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
\examples{
#Using the output to plot the 95\% credibility region 
data(results)
geplot(results$var_a, results$samples_beta0, results$samples_beta1)
}
\keyword{utilities}