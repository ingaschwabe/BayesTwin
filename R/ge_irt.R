#==========================================================
# Analysis of genotype by environment interaction 
# in twin data. With IRT measurement model
# function ge_irt.R for BayesTwin package
#==========================================================


# This function writes a JAGS script for IRT model with GE interaction
# and draws from the posterior
# and organizes output
ge_irt <- function(data_mz=data_mz, data_dz=data_dz, n_burnin=n_burnin, n_iter=n_iter){
    
        
        
    #Install packages when necessary (later niet meer nodig!!)
    if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
    if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
    if(!require(MCMCpack)){ install.packages('MASS'); require(MASS)}    
    if(!require(MCMCpack)){ install.packages('rjags'); require(rjags)}    
    
    # determine number of twin pairs
    n_mz <- nrow(data_mz) ; n_dz<- nrow(data_dz)
    
    # determine number of phenotypic variables
    n_var <- ncol(data_mz)/2

    
    #==========================================================
    # I. Write JAGS model file
    #==========================================================
    jags_model_ge <- "model{
        ##MZ twins
        for (fam in 1:n_mz){
            c_mz[fam] ~ dnorm(mu, tau_c)
            f_mz[fam] ~ dnorm(c_mz[fam], tau_a) 
            a_mz[fam] <- f_mz[fam] - c_mz[fam] 
            tau_e[fam] <- 1/(exp(beta0 + (beta1*a_mz[fam])))

            for (twin in 1:2){
                pheno_mz[fam,twin]~ dnorm(f_mz[fam],tau_e[fam])   
            }

                #1pl model twin1
                for (k in 1:n_var){
                    logit(p[fam,k]) <- pheno_mz[fam,1] - b[k]
                    data_mz[fam,k] ~ dbern(p[fam,k])
                }   		
    
                #1pl model twin2
                for (k in (n_var+1):(2*n_var)){
                    logit(p[fam,k]) <- pheno_mz[fam,2] - b[k-n_var]
                    data_mz[fam,k] ~ dbern(p[fam,k])
                }
        }

        ##DZ twins
        for (fam in 1:n_dz){
            c_dz[fam] ~ dnorm(mu, tau_c)
            f0_dz[fam] ~ dnorm(c_dz[fam], doubletau_a)

            for (twin in 1:2){										
                f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
                a_dz[fam,twin] <- f_dz[fam,twin] - c_dz[fam]
                tau_e_dz[fam,twin] <- 1/(exp(beta0 + (beta1*a_dz[fam,twin])))
                pheno_dz[fam,twin] ~ dnorm(f_dz[fam,twin], tau_e_dz[fam,twin])
            }

                #1pl model twin1 (DZ)
                for (k in 1:n_var){
                    logit(p2[fam,k]) <- pheno_dz[fam,1] - b[k]
                    data_dz[fam,k] ~ dbern(p2[fam,k])
                }

                #1pl model twin2 (DZ)
                for (k in (n_var+1):(2*n_var)){
                    logit(p2[fam,k]) <- pheno_dz[fam,2] - b[k-n_var]
                    data_dz[fam,k] ~ dbern(p2[fam,k])
                }

        }


    #Priors
    mu <- 0 #to identify scale 
    beta0 ~ dnorm(-1,0.5)
    beta1 ~ dnorm(0,.1)

    doubletau_a <- 2*tau_a
    tau_a ~ dgamma(1,1)   
    tau_c ~ dgamma(1,.5)

    for (i in 1:n_var){
        b[i] ~ dnorm(0,.1)
    }

    }"


    jags_file_ge_irt <- tempfile(fileext=".txt")
    write(jags_model_ge,jags_file_ge_irt)
    
    #==========================================================
    # II. Run JAGS analysis
    #==========================================================
    inits = list(tau_a = 2, tau_c = 5)
    jags_data <- list(data_mz, data_dz, n_mz, n_dz, n_var)
    names(jags_data)<- c("data_mz", "data_dz", "n_mz", "n_dz", "n_var") 
    jags <- jags.model(jags_file_ge_irt, jags_data, inits, n.chains = 1, quiet=FALSE)
    update(jags, n_burnin)
    out <- jags.samples(jags, c("tau_a", "tau_c", "beta0", "beta1", "b"), n_iter)
    
    #==========================================================
    # III. Organize results  
    #==========================================================
    
    samples_var_a = 1/out$tau_a[,,1]
    samples_var_c = 1/out$tau_c[,,1]
    samples_beta0 = out$beta0[,,1]
    samples_beta1 = out$beta1[,,1]
    samples_b = out$b[,,1]
    
    #Calculate mean values: 
    var_a = mean(samples_var_a); var_c = mean(samples_var_c)
    beta0 = mean(samples_beta0); beta1 = mean(samples_beta1)
    b = apply(samples_b,1,mean)
    
    #Calculate SDs: 
    sd_var_a = sd(samples_var_a); sd_var_c = sd(samples_var_c)
    sd_beta0 = sd(samples_beta0); sd_beta1 = sd(samples_beta1)
    sd_b = apply(samples_b,1, sd)
    
    #Calculate HPD: 
    hpd_var_a = HPD(samples_var_a, 0.95); hpd_var_c = HPD(samples_var_c, 0.95)
    hpd_beta0 = HPD(samples_beta0, 0.95); hpd_beta1 = HPD(samples_beta1, 0.95)
    hpd_b = apply(samples_b, 1,function(x) HPD(x,0.95))
    
    #Put results in a table 
    results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, beta0, sd_beta0, beta1, sd_beta1), 2, 4)
    colnames(results) <- c("varA","varC","beta0", "beta1")
    rownames(results) <- c("Posterior mean","Posterior standard deviation")
    results = as.table(results)
    
    #Print results on the fly: 
    cat("\n") 
    print(results)
    cat("\n newline")
    
    #Remind the user of the issue of convergence: 
    
    
    #And save output in a list: 
    output = list(results = results, samples_var_a = samples_var_a, samples_var_c = samples_var_a, 
                  samples_beta0 = samples_beta0, samples_beta1 = samples_beta1, samples_b = samples_b,
                  hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1)
    
    #Change class of objects in order to use right plot method 
    class(output) <- c("ge_irt","list")
    class(output$samples_var_a) <- "samples"
    class(output$samples_var_c) <- "samples"
    class(output$samples_beta0) <- "samples"
    class(output$samples_beta1) <- "samples" 
    
    
    return(output)
}


# #Test function: 
# library(rjags)
# data = simulate_twin_data(50, 20, n_items = 3)
# y_mz = data$y_mz
# y_dz = data$y_dz
# n_mz = 50; n_dz = 20; n_items = 3; n_burnin = 200; n_iter = 200
# xx = ge_irt( n_burnin = n_burnin, n_iter = n_iter,
#        data_mz = y_mz, data_dz = y_dz)
# 
# x = xx$samples_var_c
# plot(x, type = "Sampling plot")'







