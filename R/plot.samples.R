#==========================================================
# plot.samples.R
# S3 Method for the generic plot function
# to plot sampling plots and histograms of posterior dist. 
# BayesTwin package
#==========================================================

plot.samples <- function(x, type = "Histogram", ...){
    
    if(type == "Histogram"){
        #Calculate HPD: 
        hpd = HPD(x, 0.95)
        
        #Plot histogram of posterior samples: 
        hist(x,
             main=paste("Posterior distribution of ",deparse(substitute(x)),sep=" "),
             xlab=deparse(substitute(x)), ...)
        
        #Lines to indicate lowest + highest region
        abline(v = hpd[1], lwd = 2, col = "red")
        abline(v = hpd[2], lwd = 2, col = "red")
        abline(v= mean(x), lwd = 2, col = "green")
        abline(v= median(x), lwd = 2, col = "yellow")
        
        
        #Add legend
        legend("topright",
               legend=c("95% HPD","Mean","Median"),
               col=c("red","green","yellow"),lty=1,lwd=2,
               cex=0.9, pt.cex = 21)
    }
    
    
    if(type == "Sampling plot"){
        
        #Plot iteration-history: 
        plot.default(x, type = "l", ylab=deparse(substitute(x)), 
                     xlab = "Iteration")
    }
}
