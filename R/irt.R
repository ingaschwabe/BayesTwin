#==========================================================
# Analysis of twin data under the ACE model using item data
# Possible IRT models: 1PL, 2PL, GPCM, PCM
# 
# This function writes a JAGS script for the model and 
# draws from the posterior distribution and then organizes 
# the output. Estimation of genotype by environment 
# interaction is optional (by defining ge = TRUE) when 
# calling the master function twinUniv)
#==========================================================

irt <- function(data_mz, data_dz, n_burnin, n_iter, ge){
    
    #Determine number of twin pairs
    n_mz <- nrow(data_mz) ; n_dz<- nrow(data_dz)
    
    # determine number of phenotypic items
    n_var <- ncol(data_mz)/2
    
    #==========================================================
    # I. Write JAGS model file
    #==========================================================
    jags_model_irt <- paste("model{
        ##MZ twins
        for (fam in 1:n_mz){
            c_mz[fam] ~ dnorm(mu, tau_c)
            f_mz[fam] ~ dnorm(c_mz[fam], tau_a) 
          
        ",ifelse(ge,"
            a_mz[fam] <- f_mz[fam] - c_mz[fam]
            tau_e[fam] <- 1/(exp(beta0 + (beta1*a_mz[fam])))
            
            for (twin in 1:2){
                pheno_mz[fam,twin] ~ dnorm(f_mz[fam],tau_e[fam])   
            }
        ","
            for (twin in 1:2){
                pheno_mz[fam,twin] ~ dnorm(f_mz[fam],tau_e)   
            }"),"

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
    
        ",ifelse(ge,"
                for (twin in 1:2){            							
                    f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
                    a_dz[fam,twin] <- f_dz[fam,twin] - c_dz[fam]
                    tau_e_dz[fam,twin] <- 1/(exp(beta0 + (beta1*a_dz[fam,twin])))
                    pheno_dz[fam,twin] ~ dnorm(f_dz[fam,twin], tau_e_dz[fam,twin])
                }","
                for (twin in 1:2){										
                    f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
                    pheno_dz[fam,twin] ~ dnorm(f_dz[fam,twin], tau_e)
                }"),"

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
    doubletau_a <- 2*tau_a

    #Priors
    tau_c ~ dgamma(1,1)
    tau_a ~ dgamma(1,1)   

    for (i in 1:n_var){
      b[i] ~ dnorm(0,.1)
    }

    ",ifelse(ge,"
    beta0 ~ dnorm(-1,.5)
    beta1 ~ dnorm(0,.1)",
    "tau_e ~ dgamma(1,1)"),"
    }")

    jags_file_irt <- tempfile(fileext=".txt")
    write(jags_model_irt,jags_file_irt)
    
    #==========================================================
    # II. Run JAGS analysis
    #==========================================================
    inits = list(tau_a = 2, tau_c = 5)
    jags_data <- list(data_mz, data_dz, n_mz, n_dz, n_var)
    names(jags_data)<- c("data_mz", "data_dz", "n_mz", "n_dz", "n_var") 
    jags <- jags.model(jags_file_irt, jags_data, inits, n.chains = 1, quiet=FALSE)
    update(jags, n_burnin)
    
    if (ge == FALSE){
        out <- jags.samples(jags, c("tau_a", "tau_c", "tau_e", "b"), n_iter)
    } else {
        out <- jags.samples(jags, c("tau_a", "beta0", "beta1", "b"), n_iter)
    }
    
    #==========================================================
    # III. Organize results  
    #==========================================================
    if(ge == FALSE){
        #Save samples
        samples_var_a = 1/out$tau_a[,,1]; samples_var_c = 1/out$tau_c[,,1]
        samples_var_e = 1/out$tau_e[,,1]; samples_item_b = out$b[,,1]
        
        #Calculate mean values: 
        var_a = mean(samples_var_a); var_c = mean(samples_var_c); var_e = mean(samples_var_e)
        item_b = apply(out$b[,,1], 1, mean)
        
        #Calculate SDs: 
        sd_var_a = sd(samples_var_a); sd_var_c = sd(samples_var_c); sd_var_e = sd(samples_var_e)
        sd_item_b = apply(out$b[,,1], 1, sd)
        
        #Calculate HPD: 
        hpd_var_a = HPD(samples_var_a, 0.95); hpd_var_c = HPD(samples_var_c, 0.95)
        hpd_var_e = HPD(samples_var_e, 0.95); 
        hpd_item_b = apply(samples_item_b, 1, function (x) HPD(x, 0.95))
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, var_e, sd_var_e), 2, 3)
        colnames(results) <- c("varA","varC","varE")
        rownames(results) <- c("Posterior mean","Posterior standard deviation")
        results = as.table(results) 
        
        #And save output in a list: 
        output = list(results = results, samples_var_a = samples_var_a, samples_var_c = samples_var_c, 
                      samples_var_e = samples_var_e, samples_item_b = samples_item_b,
                      hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_var_e = hpd_var_e, 
                      hpd_item_b = hpd_item_b)
        
        #Change class of objects in order to use right plot method 
        class(output) <- c("sumscores","list")    
        class(output$samples_var_a) <- "samples"
        class(output$samples_var_c) <- "samples"
        class(output$samples_var_e) <- "samples"
        class(output$samples_item_b) <- "samples"
        
    } else {
        #Save samples
        samples_var_a = 1/out$tau_a[,,1]; samples_var_c = 1/out$tau_c[,,1]
        samples_beta0 = out$beta0[,,1]; samples_beta1 = out$beta1[,,1]; samples_item_b = out$b[,,1]
        
        #Calculate mean values: 
        var_a = mean(samples_var_a); var_c = mean(samples_var_c)
        beta0 = mean(samples_beta0); beta1 = mean(samples_beta1)
        item_b = apply(samples_item_b, 1, mean)
        
        #Calculate SDs: 
        sd_var_a = sd(samples_var_a); sd_var_c = sd(samples_var_c)
        sd_beta0 = sd(samples_beta0); sd_beta1 = sd(samples_beta1)
        sd_item_b = apply(samples_item_b, 1, sd)
        
        #Calculate HPD: 
        hpd_var_a = HPD(samples_var_a, 0.95); hpd_var_c = HPD(samples_var_c, 0.95)
        hpd_beta0 = HPD(samples_beta0, 0.95); hpd_beta1 = HPD(samples_beta1, 0.95)
        hpd_item_b = apply(samples_item_b, 1, function (x) HPD(x, 0.95))
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, beta0, sd_beta0, beta1, sd_beta1), 2, 4)
        colnames(results) <- c("varA","varC","beta0", "beta1")
        rownames(results) <- c("Posterior mean","Posterior standard deviation")
        results = as.table(results) 
        
        #And save output in a list: 
        output = list(results = results, samples_var_a = samples_var_a, samples_var_c = samples_var_c, 
                      samples_beta0 = samples_beta0, samples_beta1 = samples_beta1, samples_item_b = samples_item_b,
                      hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1,
                      hpd_item_b = hpd_item_b)
        
        #Change class of objects in order to use right plot method 
        class(output) <- c("ge_irt","list")    
        class(output$samples_var_a) <- "samples"
        class(output$samples_var_c) <- "samples"
        class(output$samples_beta0) <- "samples"
        class(output$samples_beta1) <- "samples"     
        class(output$samples_item_b) <- "samples"
    }
    
    return(output)
}
