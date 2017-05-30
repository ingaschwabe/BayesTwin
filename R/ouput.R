#==========================================================
# output.R
# All output functions
#
# BayesTwin package
#==========================================================

#==========================================================
# Coerce general R object into bayestwin object
#==========================================================
as.bayestwin <- function(x, ...){
    
    class(x) = "bayestwin"
}

#==========================================================
# S3 Method for the generic summary function
# to plot results of the IRT_twin.R function. 
#==========================================================
summary.bayestwin <- function(x, ...){
    x = as.bayestwin(x)
    print(x$results)
}


#==========================================================
# traceplot.R
#
# S3 Method for the generic plot function
# to plot sampling plots and histograms of posterior dist. 
# BayesTwin package
#==========================================================

plot.bayestwin <- function(x, type = "dist", ...){
    if(type == "dist"){
        #Calculate HPD: 
        hpd = HPD(x, 0.95)
        
        if(missing(main)){
            main=paste("Posterior distribution of ",deparse(substitute(x)),sep=" ")
        }
        
        if(missing(xlab)){
            xlab=deparse(substitute(x))
        }
        
        #Plot histogram of posterior samples: 
        hist(x,
             main=main, xlab = xlab, ...)
        
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
    
    
    if(type == "density"){
        
        #Plot iteration-history: 
        plot.default(x, type = "l", ylab=deparse(substitute(x)), 
                     xlab = "Iteration")
    }
}
