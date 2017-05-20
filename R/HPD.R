#==========================================================
# HPD.R
# Calculate highest posterior density (HPD)
#
# BayesTwin package
#==========================================================

HPD <- function(sample, cred_int = 0.95) {
    cred_int <- (1 - cred_int)/2 #calculate range outside of credibility region (both sides 2.5 in this case)
    lower <- round(length(sample) * cred_int, 0) 
    upper <- round(length(sample) * (1 - cred_int), 0)
    diff.int <- upper - lower
    HPDo <- sample[order(sample)][1:lower]
    HPDb <- sample[order(sample)][(diff.int + 1):upper]
    HPDI <- round(c(HPDo[order(HPDb - HPDo)[1]], HPDb[order(HPDb - HPDo)[
        1]]), 5)
    return(HPDI)
}