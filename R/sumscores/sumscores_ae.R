#==========================================================
# sumscores_ae.R
# Analysis of twin data under the ACE model using sum scores.
# 
# This function writes a JAGS script for the model and 
# draws from the posterior distribution and then organizes 
# the output. Estimation of genotype by environment 
# interaction is optional (by defining ge = TRUE) when 
# calling the master function twinUniv)
# BayesTwin package
#==========================================================

sumscores_ae <- function(data_mz, data_dz, n_burnin, n_iter, ge,
                         var_prior){
  
    #Determine number of twin pairs
    n_mz <- nrow(data_mz); n_dz <- nrow(data_dz)
  
    #==========================================================
    # I. Write JAGS model file
    #==========================================================
    jags_model_sumscores_ae <- paste("model{
        ##MZ twins
        for (fam in 1:n_mz){
            a_mz[fam] ~ dnorm(mu, tau_a)
        ",ifelse(ge,"
            tau_e[fam] <- 1/(exp(beta0 + (beta1*a_mz[fam])))
                for (twin in 1:2){
                    data_mz[fam,twin] ~ dnorm(a_mz[fam],tau_e[fam])   
                }
        ",      "for (twin in 1:2){
                    data_mz[fam,twin] ~ dnorm(a_mz[fam],tau_e)   
                }"),"
            }
  
        ##DZ twins
        for (fam in 1:n_dz){
            a0_dz[fam] ~ dnorm(mu, doubletau_a)
        
            for (twin in 1:2){
                a2_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
            }

        ",ifelse(ge,"
                for (twin in 1:2){        								
                    tau_e_dz[fam,twin] <- 1/(exp(beta0 + (beta1*a2_dz[fam,twin])))
                    data_dz[fam,twin] ~ dnorm(a2_dz[fam,twin], tau_e_dz[fam,twin])
                }",
                
                "for (twin in 1:2){										
                    data_dz[fam,twin] ~ dnorm(a2_dz[fam,twin], tau_e)
                }"),"
            }
  
            doubletau_a <- 2*tau_a

            #Priors
            mu ~ dnorm(0,.1)
            ",ifelse(INV_GAMMA,"
            tau_a ~ dgamma(1,1) 
            tau_e ~ dgamma(1,1) #not used when ge = TRUE
            ","
            tau_a ~ dunif(0,100)
            tau_e ~ dunif(0,100) #not used when ge = TRUE
            "),"

            ",ifelse(ge,"
            beta0 ~ dnorm(-1,.5)
            beta1 ~ dnorm(0,.1)",
            ""),"
        }")
  
    jags_file_sumscores_ae <- tempfile(fileext=".txt")
    write(jags_model_sumscores_ae, jags_file_sumscores_ae)
  
    #==========================================================
    # II. Run JAGS analysis
    #==========================================================
    inits = list(tau_a = 2)
    jags_data <- list(data_mz, data_dz, n_mz, n_dz)
    names(jags_data)<- c("data_mz", "data_dz", "n_mz", "n_dz") 
    jags <- jags.model(jags_file_sumscores, jags_data, inits, n.chains = 1, quiet=FALSE)
    update(jags, n_burnin)
  
    if(ge == FALSE){
        out <- jags.samples(jags, c("tau_a", "tau_e", "mu"), n_iter)
    } else{
        out <- jags.samples(jags, c("tau_a", "beta0", "beta1", "mu"), n_iter)
    }

    #==========================================================
    # III. Organize results  
    #==========================================================
    if(ge == FALSE){
        #Save samples
        samples_var_a = 1/out$tau_a[,,1]; 
        samples_var_e = 1/out$tau_e[,,1]; samples_mu = out$mu[,,1]
        
        #Calculate mean values: 
        var_a = mean(samples_var_a); var_e = mean(samples_var_e)
        mu = mean(samples_mu)
        
        #Calculate SDs: 
        sd_var_a = sd(samples_var_a); sd_var_e = sd(samples_var_e)
        sd_mu = sd(samples_mu)
        
        #Calculate HPD: 
        hpd_var_a = HPD(samples_var_a, 0.95); 
        hpd_var_e = HPD(samples_var_e, 0.95); hpd_mu = HPD(samples_mu, 0.95)
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, var_e, sd_var_e, mu, sd_mu), 2, 3)
        colnames(results) <- c("varA","varE", "Mu")
        rownames(results) <- c("Posterior mean","Posterior standard deviation")
        results = as.table(results) 
        
        #And save output in a list: 
        output = list(results = results, samples_var_a = samples_var_a, 
                      samples_var_e = samples_var_e, samples_mu = samples_mu,
                      hpd_var_a = hpd_var_a, hpd_var_e = hpd_var_e, hpd_mu = hpd_mu)
        
        #Change class of objects in order to use right plot method 
        class(output) <- c("sumscores","list")    
        class(output$samples_var_a) <- "samples"
        class(output$samples_var_e) <- "samples"
        class(output$samples_mu) <- "samples"
        
    } else {
        #Save samples
        samples_var_a = 1/out$tau_a[,,1];
        samples_beta0 = out$beta0[,,1]; samples_beta1 = out$beta1[,,1]; samples_mu = out$mu[,,1]
        
        #Calculate mean values: 
        var_a = mean(samples_var_a);
        beta0 = mean(samples_beta0); beta1 = mean(samples_beta1)
        mu = mean(samples_mu)
        
        #Calculate SDs: 
        sd_var_a = sd(samples_var_a); 
        sd_beta0 = sd(samples_beta0); sd_beta1 = sd(samples_beta1)
        sd_mu = sd(samples_mu)
        
        #Calculate HPD: 
        hpd_var_a = HPD(samples_var_a, 0.95); 
        hpd_beta0 = HPD(samples_beta0, 0.95); hpd_beta1 = HPD(samples_beta1, 0.95)
        hpd_mu = HPD(samples_mu, 0.95)
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, beta0, sd_beta0, beta1, sd_beta1, mu, sd_mu), 2, 4)
        colnames(results) <- c("varA", "beta0", "beta1", "Mu")
        rownames(results) <- c("Posterior mean","Posterior standard deviation")
        results = as.table(results) 
        
        #And save output in a list: 
        output = list(results = results, samples_var_a = samples_var_a, 
                      samples_beta0 = samples_beta0, samples_beta1 = samples_beta1, samples_mu = samples_mu,
                      hpd_var_a = hpd_var_a, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1,
                      hpd_mu = hpd_mu)
        
        #Change class of objects in order to use right plot method 
        class(output) <- c("ge_irt","list")    
        class(output$samples_var_a) <- "samples"
        class(output$samples_beta0) <- "samples"
        class(output$samples_beta1) <- "samples"     
        class(output$samples_mu) <- "samples"
    }
  
    return(output)
}