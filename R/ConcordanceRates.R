# Functions for a Bayesian analysis of concordance rates, based on 
# Van den Berg SM & Hjelmborg J (2012). Genetic analysis of rare disorders: 
# Bayesian estimation of twin concordance rates. Behavior Genetics, 42, 857-65. 
# doi:10.1007/s10519-012-9547-9


#library(xtable)

#if(!suppressMessages(require(LearnBayes))){ install.packages('LearnBayes'); suppressMessages(require(LearnBayes))}
#if(!suppressMessages(require(coda))){ install.packages('coda'); suppressMessages(require(coda))}
#if(!suppressMessages(require(KernSmooth))){ install.packages('KernSmooth'); suppressMessages(require(KernSmooth))}
#if(!suppressMessages(require(xtable))){ install.packages('xtable'); suppressMessages(require(xtable))}


## Function to compute highest posterior density interval for a parameter 
## sample1 : sample parameter of interest
## rel.int  : reliability interval, e.g. 0.95
HPD <- function(sample1, rel.int)
{
	rel.int <- (1 - rel.int)/2
	lower <- round(length(sample1) * rel.int, 0)
	upper <- round(length(sample1) * (1 - rel.int), 0)
	diff.int <- upper - lower
	HPDo <- sample1[order(sample1)][1:lower]
	HPDb <- sample1[order(sample1)][(diff.int + 1):upper]
	HPDI <- round(c(HPDo[order(HPDb - HPDo)[1]], HPDb[order(HPDb - HPDo)[
		1]]), 5)
        return(HPDI)
}
 
 
 
# reparameterization functions
# delta defined as p(aff|cotwin aff) - p(aff | cotwin not aff)
# q defined as conditional probability of being affected given that co-twin is affected
# pi defined as population prevalence of affected status
delta <- function(pi, q){  q-  pi*(1-q)/(1-pi)   }
q <- function(pi, delta){ return( delta*(1-pi)+ pi )   }

# A function that calculates the logaritm of the joint posterior distribution 
logpost<- function(theta, par)
{
        ymz <- par$data[,1]; ydz <- par$data[,2];
        lambda<- theta[1]; mu.mz<- theta[2]; mu.dz<- theta[3]
        alpha1 <- par$alpha1;alpha2<- par$alpha2
        beta1 <- par$beta1;beta2<- par$beta2; gamma1 <- par$gamma1; gamma2<- par$gamma2
        delta.mz <- (exp(mu.mz)-1)/(1+exp(mu.mz));delta.dz <- (exp(mu.dz)-1)/(1+exp(mu.dz))
        pi <- exp(lambda)/(1+exp(lambda))
        logl.mz1 <- ymz[1]*(log(pi)+log(delta.mz+ pi*(1-delta.mz)))
        logl.mz2 <-ymz[2]*(log(pi)+log(1-delta.mz*(1-pi)-pi))
        logl.mz3 <- ymz[3]*log(1+pi*(delta.mz*(1-pi)+pi-2))
        logl.dz1 <- ydz[1]*(log(pi)+ log(delta.dz + pi*(1-delta.dz)))
        logl.dz2 <-ydz[2]*(log(pi)+log(1-delta.dz*(1-pi)-pi))
        logl.dz3 <- ydz[3]*log(1+pi*(delta.dz*(1-pi)+pi-2))
        logpriorlambda <- (alpha1-1)*log(pi) + (alpha2-1)*log(1- pi)
        logpriormumz <- (beta1-1)* log((exp(delta.mz)-1)/ (1+exp(delta.mz))+1) + (beta2-1)*log(1- (exp(delta.mz)-1)/ (1+ exp(delta.mz)))
        logpriormudz <- (gamma2-1)* log((exp(delta.dz)-1)/ (1+ exp(delta.dz)) +1) + (gamma2-1)*log(1- (exp(delta.dz)-1)/ (1+ exp(delta.dz)))
        logJacobian<- lambda-2*log(1+ exp(lambda))+mu.mz-2*log(1+exp(mu.mz))+mu.dz-2*log(1+ exp(mu.dz))
        log <-logl.mz1 + logl.mz2 + logl.mz3 + logl.dz1 + logl.dz2 + logl.dz3
        log<- log + logpriorlambda + logpriormumz + logpriormudz + logJacobian
        return(log)
} # end function logpost

# Gather data:
data<- matrix(c(dataMZ, dataDZ), 3,2)

# Get parameters values for joint posterior:
par <- list(alpha1=1, alpha2=1, beta1=1, beta2=1, gamma1=1,gamma2=1, data=data)

# Compute a laplace approximation of posterior mode and variance:
outlaplace <- laplace(logpost, c(-1,0,0), par)
proposal = list(var=outlaplace$var, scale=1)

# sample parameters from unconstrained posterior
out <-rwmetrop(logpost, proposal, start=outlaplace$mode, NIterations, par)
pi <- 1/ (1+ exp(-1*out$par[,1])) # transforming mu into pi
delta.mz <- (exp(out$par[,2])-1) / (1+ exp(out$par[,2]))
delta.dz <- (exp(out$par[,3])-1) / (1+ exp(out$par[,3]))

# transforming pi and deltas into concordance rates
qmz <- q(pi,delta.mz) 
qdz <- q(pi,delta.dz)

 
# select posterior samples that satisfy constraint
satisfied<- which( (max(-1, (pi-1)/pi) < delta.mz) & (delta.mz< 1) & (max(-1, (pi-1)/pi) < delta.dz) & (delta.dz< 1))
posterior <- cbind(pi,qmz, qdz)[satisfied,]
posterior <- posterior[(NBurnin+1):NIterations, ] # first 1000 is burn-in


# Make a PDF with trace plots
pdf('PosteriorSamples.pdf')
par(mfrow=c(3,1))
for (i in 1:3) {plot((NBurnin+1):NIterations,posterior[,i], main = paste("Posterior samples of" , c('pi', 'q.mz','q.dz')[i]),ylab=c('pi', 'q.mz','q.dz')[i], type='l', xlab="iteration")}
dev.off() 

# Make a PDF with marginal posterior density plots
pdf('PosteriorDensities.pdf')
par(mfrow=c(1,3))
for (i in 1:3) 
        {
        if (i==1) 
                {
                prut<-density(posterior[,i])
                plot( function(x) dbeta(x, NAffectedPreviousStudies+1, NNotAffectedPreviousStudies+1) , main = '',xlab=c('pi', 'q.mz','q.dz')[i], ylab='Density',col=2,ylim=c(0,max(prut$y)) )
                points(density(posterior[,i]), type ='l')
                } else 
                        plot(density(posterior[,i]), main = '',xlab=c('pi', 'q.mz','q.dz')[i], ylab='Density')
                        lines(x=c(0,1), y=c(0,0), col=2, lty=1)
        }
dev.off()


# Determine 95% Highest Posterior Density intervals
HPDintervals <-t( matrix( c( HPD(posterior[,1],0.95), HPD(posterior[,2],0.95) , HPD(posterior[,3],0.95)), 2, 3) )
 