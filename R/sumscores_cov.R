#==========================================================
# Analysis of twin data using sum scores under the ACE model.
# Inclusion of covariates is possible.  
# function sumscores_cov.R for BayesTwin package
# 
# This function writes a JAGS script for an ACE model 
# with covaristes. A prior probability distribution is 
# used for missing data in the covariates. 
#==========================================================

sumscores_cov <- function(data_mz, data_dz, 
                          X_mz_twin1, X_mz_twin2, 
                          X_dz_twin1, X_dz_twin2,
                          n_burnin, n_iter, ge, 
                          cont_cov, dich_cov, 
                          dich = dich, cont = cont){

    #Determine number of twin pairs
    n_mz <- nrow(data_mz); n_dz <- nrow(data_dz)
    
    #first dich/cont item and last dich/cont item: 
    if (dich == TRUE){
        b_d = dich_cov[1] 
        n_d = length(dich_cov)
    }
    
    if(cont == TRUE){
        b_c = cont_cov[1]
        n_c = length(cont_cov)
    }
  
    #Determine number of environmental covariates: 
    N = ncol(X_mz_twin1)

    #==========================================================
    # I. Write JAGS model file
    #==========================================================
    jags_model_sumscores_cov <- paste("model{
        for (fam in 1:n_mz){
            c_mz[fam] ~ dnorm(mu, tau_c)
            f_mz[fam] ~ dnorm(c_mz[fam], tau_a) 
        ",ifelse(ge,"
            a_mz[fam] <- f_mz[fam] - c_mz[fam]
            tau_e[fam] <- 1/(exp(beta0 + (beta1*a_mz[fam])))
            y_mz[fam,1] ~ dnorm(f_mz[fam] + inprod(X_mz_twin1[fam,], b),tau_e[fam])
            y_mz[fam,2] ~ dnorm(f_mz[fam] + inprod(X_mz_twin2[fam,], b),tau_e[fam]) 
        ","    
            y_mz[fam,1] ~ dnorm(f_mz[fam] + inprod(X_mz_twin1[fam,], b),tau_e) 
            y_mz[fam,2] ~ dnorm(f_mz[fam] + inprod(X_mz_twin2[fam,], b),tau_e) 
        "),"
        }
        
        ##DZ twins
        for (fam in 1:n_dz){ #for each DZ family: 
            c_dz[fam] ~ dnorm(mu, tau_c) 
            f0_dz[fam] ~ dnorm(c_dz[fam], doubletau_a)
                
        ", ifelse(ge,"
            for (twin in 1:2){        								
                f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
                a_dz[fam,twin] <- f_dz[fam,twin] - c_dz[fam]
                tau_e_dz[fam,twin] <- 1/(exp(beta0 + (beta1*a_dz[fam,twin])))
                }

            #Response variables: 
            y_dz[fam,1] ~ dnorm(f_dz[fam,1] + inprod(X_dz_twin1[fam,], b), tau_e_dz[fam,1]) 
            y_dz[fam,2] ~ dnorm(f_dz[fam,2] + inprod(X_dz_twin2[fam,], b), tau_e_dz[fam,2])
            
        ","
            for (twin in 1:2){            							
                f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
                
            }
        
            #Response variables: 
            y_dz[fam,1] ~ dnorm(f_dz[fam,1] + inprod(X_dz_twin1[fam,], b), tau_e) 
            y_dz[fam,2] ~ dnorm(f_dz[fam,2] + inprod(X_dz_twin2[fam,], b), tau_e)
                
        "),"
        }
        
            #For DZ twins:
            doubletau_a <- 2*tau_a    
        
            #Priors:
            tau_a ~ dgamma(1,1) 
            tau_c ~ dgamma(1,1)
        
            ",ifelse(ge,"        
            beta0 ~ dnorm(-1,0.5)
            beta1 ~ dnorm(0,.1)
            ","
            tau_e ~ dgamma(1,1)
            "),"
        
            #If response variable is standardized to have a mean value of 0, 
            #we can also fix its value to 0 (mu <- 0)
            mu ~ dnorm(0, .1)
        
            #Prior for regression coefficients: 
            #Multivariate normal 
            b[1:N] ~ dmnorm(mu_b[1:N], tau_b[,]) 
        
        
            ##Latent variable psi for MZ and DZ twins
            for(fam in 1:n_mz){
                psi_mz[fam, 1:N] ~ dmnorm(gamma_cov[1:N],tau_b_mz) 
            }
        
            for(fam in 1:n_dz){
                psi_dz[fam, 1:N] ~ dmnorm(gamma_cov[1:N],tau_b_dz) 
            }
        
            #For the rest (continous covariates), we use a normal distribution as prior
            #For the covariates that we were standardized, the expected value can also be fixed to 0. 
        
            for (i in 1:N){
                gamma_cov[i] ~  dnorm(0, .1)
            }
        
            #As a prior distribution for tau b, we use a wishart distribution. 
            tau_b_mz[1:N, 1:N] ~ dwish(omega_tau_b[,], N)  
            tau_b_dz[1:N, 1:N] ~ dwish(omega_tau_b[,], N)
            
        
            #Final prior distribution for the covariates: 
            for (fam in 1:n_mz){
                exp_mz_twin1[fam,1:N] ~ dmnorm(psi_mz[fam,], tau_w_mz)
                exp_mz_twin2[fam,1:N] ~ dmnorm(psi_mz[fam,], tau_w_mz)
            }
        
            for (fam in 1:n_dz){
                exp_dz_twin1[fam,1:N] ~ dmnorm(psi_dz[fam,], tau_w_dz)
                exp_dz_twin2[fam,1:N] ~ dmnorm(psi_dz[fam,], tau_w_dz)
            }
        
            #As a prior distribution for tau w, we use a wishart distribution. 
            tau_w_mz[1:N, 1:N] ~ dwish(omega_tau_w[,], N)  
            tau_w_dz[1:N, 1:N] ~ dwish(omega_tau_w[,], N)
            
        
            #For the dummy variables that are coded as 0  and 1, we use exp_mz_twin1, exp_mz_twin2, 
            #exp_dz_twin1 and exp_dz_twin2 as liabilites and use a bernoulli distribution as prior: 
        
            #To identify the model, the threshold t is fixed to 0. 
            ",ifelse(dich,"    
            t <- 0 
        
            for(fam in 1:n_mz){
                for (covariate in b_d:n_d){
                    V_mz_twin1[fam,covariate] <- step(exp_mz_twin1[fam,covariate] - t)
                    V_mz_twin2[fam,covariate] <- step(exp_mz_twin2[fam,covariate] - t)
                    X_mz_twin1[fam,covariate] ~ dbern(ifelse(V_mz_twin1[fam,covariate]==1,0.999,0.001))
                    X_mz_twin2[fam,covariate] ~ dbern(ifelse(V_mz_twin2[fam,covariate]==1,0.999,0.001))
                }
            }	
        
            for(fam in 1:n_dz){
                for (covariate in b_d:n_d){
                    V_dz_twin1[fam,covariate]<-step(exp_dz_twin1[fam,covariate] - t)
                    V_dz_twin2[fam,covariate]<-step(exp_dz_twin2[fam,covariate] - t)
                    X_dz_twin1[fam,covariate] ~ dbern(ifelse(V_dz_twin1[fam,covariate]==1,0.999,0.001))
                    X_dz_twin2[fam,covariate] ~ dbern(ifelse(V_dz_twin2[fam,covariate]==1,0.999,0.001))
                }
            }	
            ","     "), "
        
            #As JAGS does not allow to define a MVN distribution for partly missing data, we place
            #a distribution with exp_mz_twin1/2, exp_dz_twin1/2 as expected values on the continous 
            #covariates. We use a very small residual variance of 0.01. 
            ",ifelse(cont,"  
            for (fam in 1:n_mz){
                for (covariate in b_c:n_c){
                    X_mz_twin1[fam,covariate] ~ dnorm(exp_mz_twin1[fam,covariate],100)
                    X_mz_twin2[fam,covariate] ~ dnorm(exp_mz_twin2[fam,covariate],100)
                }
            }
        
            for (fam in 1:n_dz){
                for (covariate in b_c:n_c){
                    X_dz_twin1[fam,covariate] ~ dnorm(exp_dz_twin1[fam,covariate],100)
                    X_dz_twin2[fam,covariate] ~ dnorm(exp_dz_twin2[fam,covariate],100) 
                }
            }
            ","     "), "
        
        }")
    jags_file_sumscores_cov <- tempfile(fileext=".txt")
    write(jags_model_sumscores_cov,jags_file_sumscores_cov)
    #writeLines(jags_model_sumscores_cov, con = "file.txt", sep = "\n", useBytes = FALSE)
    
    #==========================================================
    # II. Run JAGS analysis
    #==========================================================
    inits = list(tau_a = 2, tau_c = 5)
    
    if(cont == TRUE && dich == TRUE){
        jags_data <- list("y_mz" = data_mz, "y_dz" = data_dz, 
                          "X_mz_twin1" = X_mz_twin1, "X_mz_twin2" = X_mz_twin2, 
                          "X_dz_twin1" = X_dz_twin1, "X_dz_twin2" = X_dz_twin2,
                          "n_mz" = n_mz, "n_dz" = n_dz, "N" = N, 
                          "b_d" = b_d, "n_d" = n_d, "b_c" = b_c, "n_c" = n_c,   ###Data and variables
        
                          "mu_b" = rep(0,N), "tau_b" = diag(1,N),
                          "omega_tau_b" = diag(1,N), "omega_tau_w" = diag(1,N)) ### for prior distributions
    
    } else if (cont == TRUE && dich == FALSE){
        jags_data <- list("y_mz" = data_mz, "y_dz" = data_dz, 
                          "X_mz_twin1" = X_mz_twin1, "X_mz_twin2" = X_mz_twin2, 
                          "X_dz_twin1" = X_dz_twin1, "X_dz_twin2" = X_dz_twin2,
                          "n_mz" = n_mz, "n_dz" = n_dz, "N" = N, 
                          "b_c" = b_c, "n_c" = n_c,                             ###Data and variables
                          
                          "mu_b" = rep(0,N), "tau_b" = diag(1,N),
                          "omega_tau_b" = diag(1,N), "omega_tau_w" = diag(1,N)) ### for prior distributions
    } else {
        jags_data <- list("y_mz" = data_mz, "y_dz" = data_dz, 
                          "X_mz_twin1" = X_mz_twin1, "X_mz_twin2" = X_mz_twin2, 
                          "X_dz_twin1" = X_dz_twin1, "X_dz_twin2" = X_dz_twin2,
                          "n_mz" = n_mz, "n_dz" = n_dz, "N" = N, 
                          "b_d" = b_d, "n_d" = n_d, ###Data and variables
                          
                          "mu_b" = rep(0,N), "tau_b" = diag(1,N),
                          "omega_tau_b" = diag(1,N), "omega_tau_w" = diag(1,N)) ### for prior distributions
        
    }
           
    jags <- jags.model(jags_file_sumscores_cov, jags_data, inits, n.chains = 1, quiet=FALSE)    
    update(jags, n_burnin)
    
    if(ge == FALSE){
        out <- jags.samples(jags, c("tau_a", "tau_c", "tau_e", "b", "mu"), n_iter)
    } else{
        out <- jags.samples(jags, c("tau_a", "tau_c", "beta0", "beta1", "b", "mu"), n_iter)
    }
        
    #==========================================================
    # III. Organize results  
    #==========================================================
    
    if(ge == FALSE){
        #Save samples: 
        samples_var_a = 1/out$tau_a[,,1]; samples_var_c = 1/out$tau_c[,,1]
        samples_var_e = 1/out$tau_e[,,1]; samples_b = out$b[,,1]; samples_mu = out$mu[,,1]
        
        #Calculate mean values: 
        var_a = mean(samples_var_a); var_c = mean(samples_var_c); var_e = mean(samples_var_e)
        mu = mean(samples_mu); b = apply(samples_b, 1, mean)
  
        #Calculate SDs: 
        sd_var_a = sd(samples_var_a); sd_var_c = sd(samples_var_c); sd_var_e = sd(samples_var_e)
        sd_mu = sd(samples_mu); sd_b = apply(samples_b, 1, sd)
  
        #Calculate HPD:   
        hpd_var_a = HPD(samples_var_a, 0.95); hpd_var_c = HPD(samples_var_c, 0.95)
        hpd_var_e = HPD(samples_var_e, 0.95); hpd_mu = HPD(samples_mu, 0.95)
        hpd_b = apply(samples_b, 1, function (x) HPD(x, 0.95))
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, var_e, sd_var_e, mu, sd_mu), 2, 4)
        colnames(results) <- c("varA","varC","varE", "Mu")
        rownames(results) <- c("Posterior mean","Posterior standard deviation")
        results = as.table(results) 
        
        colnames_reg <- rep(NA, length(b))
        for (i in 1:length(b)){
            colnames_reg[i] <- paste("b", i)
        }    
        results_reg = matrix(c(b, sd_b), 2, length(b))
        colnames(results_reg) <- colnames_reg
        rownames(results_reg) <- c("Posterior mean","Posterior standard deviation")    
        results_reg = as.table(results_reg)
        
        #And save output in a list: 
        output = list(results = results, results_reg = results_reg, samples_var_a = samples_var_a, samples_var_c = samples_var_c, 
                      samples_var_e = samples_var_e, samples_b = samples_b, hpd_b = hpd_b,
                      hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_var_e = hpd_var_e)
  
        #Change class of objects in order to use right plot method 
        class(output) <- c("sumscores","list")    
        class(output$samples_var_a) <- "samples"
        class(output$samples_var_c) <- "samples"
        class(output$samples_var_e) <- "samples"
        class(output$samples_b) <- "samples"
    } else {
        #Save samples
        samples_var_a = 1/out$tau_a[,,1]; samples_var_c = 1/out$tau_c[,,1]
        samples_beta0 = out$beta0[,,1]; samples_beta1 = out$beta1[,,1]; samples_b = out$b[,,1]
        
        #Calculate mean values: 
        var_a = mean(samples_var_a); var_c = mean(samples_var_c)
        beta0 = mean(samples_beta0); beta1 = mean(samples_beta1); b = apply(samples_b, 1, mean)
        
        #Calculate SDs: 
        sd_var_a = sd(samples_var_a); sd_var_c = sd(samples_var_c)
        sd_beta0 = sd(samples_beta0); sd_beta1 = sd(samples_beta1)
        sd_b = apply(samples_b, 1, sd)
        
        #Calculate HPD: 
        hpd_var_a = HPD(samples_var_a, 0.95); hpd_var_c = HPD(samples_var_c, 0.95)
        hpd_beta0 = HPD(samples_beta0, 0.95); hpd_beta1 = HPD(samples_beta1, 0.95)
        hpd_b = apply(samples_b, 1, function (x) HPD(x, 0.95))
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, beta0, sd_beta0, beta1, sd_beta1), 2, 4)
        colnames(results) <- c("varA","varC","beta0", "beta1")
        rownames(results) <- c("Posterior mean","Posterior standard deviation")
        results = as.table(results) 
        
        colnames_reg <- rep(NA, length(b))
        for (i in 1:length(b)){
            colnames_reg[i] <- paste("b", i)
        }    
        results_reg = matrix(c(b, sd_b), 2, length(b))
        colnames(results_reg) <- colnames_reg
        rownames(results_reg) <- c("Posterior mean","Posterior standard deviation")    
        results_reg = as.table(results_reg)
        
        #And save output in a list: 
        output = list(results = results, samples_var_a = samples_var_a, samples_var_c = samples_var_c, 
                      samples_beta0 = samples_beta0, samples_beta1 = samples_beta1, samples_b = samples_b,
                      hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1,
                      hpd_b = hpd_b)
        
        #Change class of objects in order to use right plot method 
        class(output) <- c("ge_irt","list")    
        class(output$samples_var_a) <- "samples"
        class(output$samples_var_c) <- "samples"
        class(output$samples_beta0) <- "samples"
        class(output$samples_beta1) <- "samples"    
        class(output$samples_b) <- "samples"    
    }
  
  return(output)
}