#==========================================================
# Analysis of genotype by environment interaction 
# in twin data. With IRT measurement model
# function ge_irt.R for BayesTwin package
#==========================================================

ge_irt <- function(y_mz, y_dz, n_mz, n_dz, n_items, burnin, n_iterations){
      
    #Install packages when necessary (later niet meer nodig!!)
    if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
    if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
    if(!require(MCMCpack)){ install.packages('MASS'); require(MASS)}    
    if(!require(MCMCpack)){ install.packages('rjags'); require(rjags)}    
    
    #==========================================================
    # I. Write JAGS model file
    #==========================================================
    jags_model <- "model{
        ##MZ twins
        for (fam in 1:n_mz){
            c_mz[fam] ~ dnorm(mu, tau_c)
            f_mz[fam] ~ dnorm(c_mz[fam], tau_a) 

            for (twin in 1:2){
                pheno_mz[fam,twin]~ dnorm(f_mz[fam],tau_e)   
            }

                #1pl model twin1
                for (k in 1:n_items){
                    logit(p[fam,k]) <- pheno_mz[fam,1] - b[k]
                    y_mz[fam,k] ~ dbern(p[fam,k])
                }   		
    
                #1pl model twin2
                for (k in (n_items+1):(2*n_items)){
                    logit(p[fam,k]) <- pheno_mz[fam,2] - b[k-n_items]
                    y_mz[fam,k] ~ dbern(p[fam,k])
                }
        }

        ##DZ twins
        for (fam in 1:n_dz){
            c_dz[fam] ~ dnorm(mu, tau_c)
            f0_dz[fam] ~ dnorm(c_dz[fam], doubletau_a)

            for (twin in 1:2){										
                f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
                pheno_dz[fam,twin] ~ dnorm(f_dz[fam,twin], tau_e)
            }

                #1pl model twin1 (DZ)
                for (k in 1:n_items){
                    logit(p2[fam,k]) <- pheno_dz[fam,1] - b[k]
                    y_dz[fam,k] ~ dbern(p2[fam,k])
                }

                #1pl model twin2 (DZ)
                for (k in (n_items+1):(2*n_items)){
                    logit(p2[fam,k]) <- pheno_dz[fam,2] - b[k-n_items]
                    y_dz[fam,k] ~ dbern(p2[fam,k])
                }

        }


    #Priors
    mu <- 0 #to identify scale

    doubletau_a <- 2*tau_a
    tau_a ~ dgamma(1,1)   
    tau_c ~ dgamma(1,.5)
    tau_e ~ dgamma(1,1)

    for (i in 1:n_items){
        b[i] ~ dnorm(0,.1)
    }

    }"


    jags_file_irt <- tempfile(fileext=".txt")
    write(jags_model,jags_file_irt)
    
    #==========================================================
    # II. Run JAGS analysis
    #==========================================================
    inits = list(tau_a = 2, tau_c = 5, tau_e = 2)
    jags_data <- list(y_mz, y_dz, n_mz, n_dz, n_items)
    names(jags_data)<- c("y_mz", "y_dz", "n_mz", "n_dz", "n_items") 
    jags <- jags.model(jags_file_irt, jags_data, inits, n.chains = 1, quiet=FALSE)
    update(jags, burnin)
    out <- jags.samples(jags, c("tau_a", "tau_c", "tau_e", "b"), n_iterations)
    
    return(out)
}


#Test function: 
'library(rjags)
data = simulate_twin_data(50, 20, n_items = 3)
y_mz = data$y_mz
y_dz = data$y_dz
n_mz = 50; n_dz = 20; n_items = 3; burnin = 200; n_iterations = 200
xx = ge_irt(n_mz = n_mz, n_dz = n_dz, n_items = n_items, burnin = burnin, n_iterations = n_iterations,
       y_mz = y_mz, y_dz = y_dz)

x = xx$samples_var_c
plot(x, type = "Sampling plot")'