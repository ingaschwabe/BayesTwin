\name{plot_bayestwin}
\alias{plot_bayestwin}
\alias{plotbayestwin}
\alias{plotBayestwin}
\title{Summary plots of BayesTwin MCMC-objects}
\description{\code{plot_bayestwin} returns either a density or trace plot for the target parameter.}
\usage{
plot_bayestwin(sample, type = "density", main, xlab, ylab, legend = TRUE, ...)}
\arguments{
  \item{sample}{A vector representing draws from the posterior distribution of the parameter of interest, as produced by the main function \code{IRT_twin} of this package}
 \item{type}{Type of plot that is produced, either "density" to get a histogram of the posterior distribution or "trace"to get a trace plot. Default value is "density"}
\item{main}{An overall title for the plot}
\item{xlab}{A title for the x axis}
\item{ylab}{A title for the y axis}
\item{legend}{Indicating if a legend should be added. Default value is TRUE.}
\item{...}{Further arguments for the generic S3 method histogram (for the density plot) or plot (for the trace plot) method}
}
\details{Creates either a density or trace plot for the target parameter. 

The function expects as input a vector representing draws from the posterior distributions of the parameter of interest, as produced by the main function \code{IRT_twin} of this package.
}
\author{Inga Schwabe}
\keyword{IRT_twin}