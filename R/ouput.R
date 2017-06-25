#==========================================================
# output.R
# All output and auxiliarly functions
#
# BayesTwin package
#==========================================================

#==========================================================
# summary.bayestwin
# S3 Method for the generic summary function
# to plot results of the IRT_twin.R function. 
#==========================================================
summary.bayestwin <- function(object, ...){
    if(exists("object$results_b") == TRUE){
        return(list(object$results, object$results_b))
    }else{
        return(object$results)
    }
}

#==========================================================
#plotbayestwin
# Plot sampling plots and histograms of posterior dist. 
#==========================================================
plotbayestwin <- function(sample, t = "density", 
                          main, xlab, ylab, legend = TRUE, lines = TRUE, ...){
    
    if(t == "density"){
        #Calculate HPD: 
        hpd = HPD(sample, 0.95)
    
        if(missing(main)){
            main=paste("Posterior distribution of ",deparse(substitute(sample)),sep=" ")
        }
        
        if(missing(xlab)){xlab=deparse(substitute(sample))}
        if(missing(ylab)){ylab = "Frequency"}
        
        #Plot histogram of posterior samples: 
        hist(sample, main=main, xlab=xlab, ylab=ylab, ...)
    
        #Lines to indicate lowest + highest region
        if (lines == TRUE){
            abline(v = hpd[1], lwd = 2, col = "black")
            abline(v = hpd[2], lwd = 2, col = "black")
            abline(v= mean(sample), lwd = 2, col = "blue")
            abline(v= median(sample), lwd = 2, col = "orange")
        }
        
        #Add legend
        if(legend == TRUE){
            legend("topright",
                    legend=c("95% HPD","Mean","Median"),
                    col=c("black","blue","orange"),lty=1,lwd=2,
                    cex=0.9, pt.cex = 21)
        }
        
    } else if (t == "trace"){
        #Plot iteration-history:
        if(missing(main)){
            main=paste("Iteration history of ",deparse(substitute(sample)),sep=" ")
        }
        
        if (missing(ylab)){ylab=deparse(substitute(sample))}
        if(missing(xlab)){xlab = "Iteration"}
        
        plot.default(sample, type = "l", ylab=ylab, xlab=xlab, main=main, ...)
        
    } else {
        plot.default(sample, ...)
    }
}

#==========================================================
# geplot
# Plot 95% CI of the GxE interaction
#==========================================================
geplot = function(var_a, samples_beta0, samples_beta1, main, xlab, ylab, col, ...){
    max_a <- 2 * sqrt(var_a); min_a <- - 2 * sqrt(var_a)
    var_e = exp(log(samples_beta0) + samples_beta1 %*% t(seq(min_a, max_a, by = .1)))
    quantiles <- apply(var_e, 2, function(x) quantile(x,probs = c(0.05,0.50,0.95)))
    
    if(missing(main)){main = "95% Credibility region GE interaction"}
    if(missing(xlab)){xlab = "Genetic value"}
    if(missing(ylab)){ylab = "Environmental variance"}
    if(missing(col)){col = "black"}
    
    matplot(main = main, seq(min_a,max_a,by=.1), t(quantiles),
            xlab = xlab, ylab = ylab, type = "l", col = col, 
            lty=c(1,2,1), ...)
}

#==========================================================
# HPD
# Calculate highest posterior density (HPD)
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
