#==========================================================
# help_simulate_twin_data.R
# Help functions for simulate_twin_data functions
# Subroutines to calculate probabilites for every response
# categories under the GPCM and PCM. 
# BayesTwin package
#==========================================================

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