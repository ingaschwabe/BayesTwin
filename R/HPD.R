#==========================================================
# Calculate highest posterior density (HPD)
# BayesTwin pacakge
#==========================================================

HPD <- function(sample1, rel.int) {
    rel.int <- (1 - rel.int)/2 #calculate range outside of credibility region (both sides 2.5 in this case)
    lower <- round(length(sample1) * rel.int, 0) 
    upper <- round(length(sample1) * (1 - rel.int), 0)
    diff.int <- upper - lower
    HPDo <- sample1[order(sample1)][1:lower]
    HPDb <- sample1[order(sample1)][(diff.int + 1):upper]
    HPDI <- round(c(HPDo[order(HPDb - HPDo)[1]], HPDb[order(HPDb - HPDo)[
        1]]), 5)
    #CI <- round(c(sample1[order(sample1)][lower], sample1[order(sample1)][
    #  upper]), 3)
    return(HPDI)
}
