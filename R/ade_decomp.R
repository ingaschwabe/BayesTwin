#==========================================================
# Analysis of twin data, using sum scores, ADE model. 
# function ade_decomp.R for BayesTwin package
#
# This function writes a JAGS script for the model,
# draws from the posterior distribution and then 
# organizes the output
#==========================================================

sumscores_ade <- function(data_mz, data_dz, n_burnin, n_iter){
    
    #Install packages when necessary (later niet meer nodig!!)
    if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
    if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
    if(!require(MCMCpack)){ install.packages('MASS'); require(MASS)}    
    if(!require(MCMCpack)){ install.packages('rjags'); require(rjags)}    
    
    #Load packages:
    require(rjags, MCMCpack, R.utils)
    
    # Determine number of twin pairs
    n_mz <- nrow(data_mz); n_dz <- nrow(data_dz)
    
    #==========================================================
    # I. Write JAGS model file
    #==========================================================
    jags_model_sumscores_ade <- "model{
        #First compute sum and difference scores from observed data
        #for MZ twins:
        for (i in n_mz){
            dif_mz[i] <- data_mz[i,1] - data_mz[i,2]
            sum_mz[i] <- data_mz[i,1] + data_mz[i,2]
        }
        
        #and for DZ twins: 
        for (i in n_dz){
            dif_dz[i] <- data_dz[i,1] - data_dz[i,2]
            sum_dz[i] <- data_dz[i,1] + data_dz[i,2]
        }

        #model sum and difference scores: 
        for (i in n_mz){
            dif_mz[i] ~ dnorm(0, tau_mz_dif)
            sum_mz[i] ~ dnorm(double_mu, tau_mz_sum)
        }

        for (i in 1:n_dz){
            dif_dz[i] ~ dnorm(0, tau_dz_dif)
            sum_dz[i] ~ dnorm(double_mu, tau_dz_sum)
        }

        #model for singleton: nu nog niet in dit script
        #for (i in 1:n_singleton){
        #    singleton[i] ~ dnorm(mu, tau_singleton)
        #}
        
        tau_mz_dif <- 1/var_mz_dif
        tau_dz_dif <- 1/var_dz_dif
        tau_mz_sum <- 1/var_mz_sum
        tau_dz_sum <- 1/var_dz_sum
        #tau_singleton <- 1/var_singleton
        
        var_mz_dif <- 2 * var_e
        var_dz_dif <- var_a + 1.5 * var_d + 2 * var_e
        var_mz_sum <- 4 * var_a + 4 * var_d + 2 * var_e
        var_dz_sum <- 3 * var_a + 2.5 * var_d + 2 * var_e
        #var_singleton <- var_a + var_d + var_e
        
        var_a <- pow(std_a,2)
        std_a ~ dunif(0,100)
        var_d <- pow(std_d,2)
        std_d ~ dunif(0,100)
        var_e <- pow(std_e,2)
        std_e ~ dunif(0,100)
        mu ~ dnorm(0, .01)
        doubl_emu <- 2 * mu 
    }"
        
    jags_file_sumscores_ade <- tempfile(fileext=".txt")
    write(jags_model_sumscores_ade, jags_file_sumscores_ade)
    
    #==========================================================
    # II. Run JAGS analysis
    #==========================================================
    inits = list(tau_a = 2, tau_d = 5)
    jags_data <- list(data_mz, data_dz, n_mz, n_dz)
    names(jags_data)<- c("data_mz", "data_dz", "n_mz", "n_dz") 
    jags <- jags.model(jags_file_sumscores_ade, jags_data, inits, n.chains = 1, quiet=FALSE)
    update(jags, n_burnin)
    out <- jags.samples(jags, c("tau_a", "tau_c", "beta0", "beta1"), n_iter)
    
    #==========================================================
    # III. Organize results  
    #==========================================================
    
    samples_var_a = 1/out$tau_a[,,1]; samples_var_c = 1/out$tau_c[,,1]
    samples_beta0 = out$beta0[,,1]; samples_beta1 = out$beta1[,,1]
    
    #Calculate mean values: 
    var_a = mean(samples_var_a); var_c = mean(samples_var_c)
    beta0 = mean(samples_beta0); beta1 = mean(samples_beta1)
    
    #Calculate SDs: 
    sd_var_a = sd(samples_var_a); sd_var_c = sd(samples_var_c)
    sd_beta0 = sd(samples_beta0); sd_beta1 = sd(samples_beta1)
    
    #Calculate HPD: 
    hpd_var_a = HPD(samples_var_a, 0.95); hpd_var_c = HPD(samples_var_c, 0.95)
    hpd_beta0 = HPD(samples_beta0, 0.95); hpd_beta1 = HPD(samples_beta1, 0.95)
    
    #Put results in a table 
    results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, beta0, sd_beta0, beta1, sd_beta1), 2, 4)
    colnames(results) <- c("varA","varC","beta0", "beta1")
    rownames(results) <- c("Posterior mean","Posterior standard deviation")
    results = as.table(results) 
    
    #And save output in a list: 
    output = list(results = results, samples_var_a = samples_var_a, samples_var_c = samples_var_c, 
                  samples_beta0 = samples_beta0, samples_beta1 = samples_beta1,
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
setwd("/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
source("simulate_twin_data.R")
source("HPD.R")
source("plot.samples.R")
data = simulate_twin_data(50, 20, n_items = 0)
data_mz = data$y_mz
data_dz = data$y_dz
n_burnin = 200; n_iter = 200
xx = sumscores_ge(data_mz = data_mz, data_dz = data_dz, n_burnin = 200, n_iter = 200)
xx
plot(xx$samples_var_a)'

# Raw data (example): 
list(mz = structure(.Data=c( 4, 3, 	2, 3,  4, 5,	 4,3, 	1, 2,	 6, 5,	 4, 5	, 3, 2  ), .Dim=c(8, 2)), dz = structure(.Data=c( 3, 5, 2, 4, 2, 2, 4,2, 2, 4, 6, 3, 3, 5, 5, 4 ), .Dim=c(8, 2)),   Singleton = c(4, 5, 2))
