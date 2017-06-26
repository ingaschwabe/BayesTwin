\name{plotbayestwin}
\alias{plotbayestwin}
\alias{plotbayestwin}
\alias{plotBayestwin}
\title{Summary plots of BayesTwin MCMC-objects}
\description{\code{plotbayestwin} returns either a density or trace plot for the target parameter.}
\usage{
plotbayestwin(sample, t = "density", main, xlab, ylab, legend = TRUE, lines = TRUE, ...)}
\arguments{
  \item{sample}{A vector representing draws from the posterior distribution of the parameter of interest, as produced by the main function \code{IRT_twin} of this package}
 \item{t}{Type of plot that is produced, either "density" to get a histogram of the posterior distribution or "trace"to get a trace plot. Default value is "density"}
\item{main}{An overall title for the plot}
\item{xlab}{A title for the x axis}
\item{ylab}{A title for the y axis}
\item{legend}{Indicating if a legend should be added. Default value is TRUE.}
\item{lines}{Indicating if lines for mean, median and lower and upper limit of the 95\% HPD should be added. Default value is TRUE} 
\item{...}{Further arguments for the generic S3 method histogram (for the density plot) or plot (for the trace plot) method}
}
\details{Creates either a density or trace plot for the target parameter. 

The function expects as input a vector representing draws from the posterior distributions of the parameter of interest, as produced by the main function \code{IRT_twin} of this package.
}
\author{Inga Schwabe}
\examples{
data(results)

##Plot posterior distribution of var(A)
plotbayestwin(results$samples_var_a)

}
\keyword{utilities}