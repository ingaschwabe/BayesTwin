\name{summary.bayestwin}
\alias{summary.bayestwin}
\alias{summarybayestwin}
\alias{summaryBayestwin}
\title{Summary statistics for BayesTwin analysis}
\description{\code{summary.bayestwin} produces summary statistics for variance components and, if applicable, regression parameters: The posterior point estimate and its standard deviation and the 95\% HPD.}
\usage{\method{summary}{bayestwin}(object, ...)}
\arguments{
    \item{object}{An object storing output from the main function \code{IRT_twin}}
    \item{...}{Further arguments for the default S3 summary method}
}
\author{Inga Schwabe}
\examples{
data(results)
summary(results)}
\keyword{IRTtwin}