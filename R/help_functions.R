#==========================================================
# help_functions.R
#
# Help functions for BayesTwin pacakge
# Subroutines to calculate probabilites for every response
# categories under the GPCM and PCM, check if response is binary
#
# BayesTwin package
#==========================================================

#For warning in IRT_twin (check if dummy variables need to be used)
is.categorical = function(x){
    length(unique(x)) <= 7 #allow for values 0,1 and NA. 
}


#Help functions to simulate data under the GPCM/PCM
numerator <- function(a, b, theta, k, D=1.7) exp(sum(D*a*(theta-b[1:k])))

denominator <- function(a, b, theta, D=1.7) {
    sumval <- 0
    for (i in 1:length(b)) sumval <- sumval + numerator(a,b,theta,i,D)
    return(sumval)
}

#Category characteristic curve
prm <- function(a, b, theta, k, D=1.7) numerator(a, b, theta, k, D)/denominator(a, b, theta, D)

#cumulative probability function
cprm <- function(a, b, theta, D=1.7) {
    ncat = length(b)
    p=1
    for (i in 1:ncat) p[i]<- prm(a,b,theta,i,D)
    pmf = p
    for (i in 1:ncat) pmf[1:ncat>i]<-pmf[1:ncat>i]+p[i]
    return(pmf)
}