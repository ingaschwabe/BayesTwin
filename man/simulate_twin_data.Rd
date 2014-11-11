\name{simulate_twin_data}
\alias{simulate_twin_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of twin data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function simulates twin data. It can be use to simulate single responses (sumscores), but also to 
simuulate item data. Item data is simulated on basis of the 1PL Rasch model (Rasch, 1960). 
}
\usage{
simulate_twin_data(nmz = 300, ndz = 500, var_a = 0.5, 
                   var_c = 0.3,  var_e = 0.2,
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
  \item{n_items}{numeric value, number of items when item data is simulated}
  \item{n_var}{numeric value, number of covariates that are simulated}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Rasch G (1960) Probabilistic models for some intelligence and attainment tests. Danish Institute or Educational Research, Copenhagen
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
