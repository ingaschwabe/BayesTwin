\name{simulate_twin_data}
\alias{simulate_twin_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of twin data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function simulates twin data on basis of the ACE or ADE (coming in the future) model (Jinks & Fulker, 1970). It can be use to simulate single responses (sumscores), but also to simuulate item data. Item data is simulated on basis of the 1PL Rasch model (Rasch, 1960). Answers to covariates can be simulated as well. 
}
\usage{
simulate_twin_data(nmz, ndz, 
                   var_a = 0.5, var_c = 0.3,  var_e = 0.2,
                   model = "ACE", n_items = 0, n_var = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nmz}{a numeric value, total number of MZ twins}
  \item{ndz}{a numeric value, total number of DZ twins}
  \item{var_a}{a numeric value, variance due to genetic effects. 
               Defaults to 0.5}
  \item{var_c}{a numeric value, variance due to shared-environmental effects. 
               Defaults to 0.3}
  \item{var_e}{a numeric value, variance due to unique-environmental effects. 
               Defaults to 0.2}
  \item{model}{specify the twin model that is used. 
               At this moment only "ACE" available. "ADE" coming in the future}
  \item{n_items}{numeric value, number of items when item data is simulated. Defaults to 0}
  \item{n_var}{numeric value, number of covariates that are simulated. Defaults to 0}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{y_mz}{
A matrix of nmz rows and 2 columns with the answers of all MZ twins. Answers of twin 1 in the first column, answers of twin 2 in the second column. 

When item data is simulated for a total of k items, the result is a matrix of nmz rows for the i-th MZ family and 2 columns of k items for twin 1 [1:k] and twin 2 [k+1:2*k] For example, y_mz[1,22] is the response of 
MZ twin 1 from family 1 to item 22 and y_mz[1,23] is the response of MZ twin 2 from family 1 to item 1 if n_items = 22
}
\item{y_dz}{
A matrix with the answers of all DZ twins. Results are in the same structure as for y_mz
}

\item{cov_mz}{
Responses of all MZ twins to covariates. A matrix of nmz rows and n_var columns. 
}

\item{cov_dz}{
Responses of all MZ twins to covariates. A matrix of nmz rows and n_var columns. 
}

\references{
%% ~put references to the literature/web site here ~
Jinks, J. L. & Fulker, D.W. (1970). Comparison of the biometrical genetical, mava, and classical approaches to the analysis of human behavior. Psychol Bull, 73(5), 311-€“349.

Rasch G (1960) Probabilistic models for some intelligence and attainment tests. Danish Institute or Educational Research, Copenhagen
}

\author{
Inga Schwabe <i.schwabe@utwente.nl>, StÃ©phanie M. van den Berg
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ace_decomp}}
}

\examples{

\dontrun{

#Generate data of 300 MZ twins and 500 DZ twins to 5 covariates
simulate_twin_data(nmz = 300, ndz = 500, n_var = 5)

#Generate data of 300 MZ twins and 500 DZ twins to 50 items
simulate_twin_data(nmz = 300, ndz = 500, n_items = 50)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
