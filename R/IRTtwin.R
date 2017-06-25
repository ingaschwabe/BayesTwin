#==========================================================
# IRTtwin.R
#
# Master function that calls subroutines for univariate
# analyses with one phenotype. 
#
# BayesTwin package
#==========================================================
IRTtwin = function(data_mz, data_dz, 
                   twin1_datacols_p, twin2_datacols_p,
                   twin1_datacols_cov = NA, twin2_datacols_cov = NA, 
                   decomp_model = "ACE",
                   irt_model = "1PL",
                   ge = FALSE,
                   n_iter = 5000, n_burnin = 7000,
                   n_chains = 1, fit_stats = FALSE, 
                   var_prior = "INV_GAMMA",
                   N_cov = 0, inits = NA, Nk = 0){
    #==========================================================
    # Error messages
    #========================================================== 
    cat("Caution! To use this package, you need to install the MCMC program JAGS
        which can be freely obtained at http://mcmc-jags.sourceforge.net. \n")
    
    if(is.matrix(data_mz) == FALSE || is.matrix(data_dz) == FALSE){
        stop("The phenotypic data has to be stored in matrix form!")
    }
    
    if(var(c(data_dz[,c(twin1_datacols_p, twin2_datacols_p)],
             data_mz[,c(twin1_datacols_p, twin2_datacols_p)]), na.rm = TRUE) <= 0.05){
        warning("It seems that the total phenotypic variance is rather low (<= 0.05).
                If you are not already doing so, consider using univariate prior distributions!")
    }
    
    if (irt_model == "GPCM" && Nk == 0){
        stop("When you want to analyse the data under the GPCM, please specify the number of categories of your 
             phenotypic data (Nk)!")
    }
    
    if (irt_model == "PCM" && Nk == 0){
        stop("When you want to analyse the data under the PCM, please specify the number of categories of your 
             phenotypic data (Nk)!")
    }
    
    if(irt_model == "GPCM" && "0" %in% data_mz[,c(twin1_datacols_p, twin2_datacols_p)] || 
       irt_model == "PCM" && "0" %in% data_mz[,c(twin1_datacols_p, twin2_datacols_p)]){
        stop("It looks like you are trying to estimate either a GPCM or PCM and your phenotypic item data 
             contains zeros. When estimating a GPCM or PCM, item data needs to be coded as e.g. 1 2 3 
             for 3 answer categories instead of 0 1 2. Please recode your item data and run the analysis
             again")
    }
    
    if(irt_model == "GPCM" && "0" %in% data_dz[,c(twin1_datacols_p,twin2_datacols_p)] || 
       irt_model == "PCM" && "0" %in% data_dz[,c(twin1_datacols_p,twin2_datacols_p)]){
        stop("It looks like you are trying to estimate either a GPCM or PCM and your phenotypic item data 
             contains zeros. When estimating a GPCM or PCM, item data needs to be coded as e.g. 1 2 3 
             for 3 answer categories instead of 0 1 2. Please recode your item data and run the analysis
             again")
    }
    
    if(length(twin1_datacols_p) == 1 || length(twin2_datacols_p) == 1){
        stop("You need phenotypic data at item level to run the analysis! It is not possible to analyse sum scores")
    } 
    
    if( n_chains < 2 && fit_stats == TRUE){
        stop("If you want fit statisitcs, you need at least 2 Markov chains! Use n_chains to specify the number of chains")
    } 
    
    if(is.na(twin1_datacols_cov) == FALSE && is.na(twin2_datacols_cov) == TRUE || 
       is.na(twin1_datacols_cov) == TRUE && is.na(twin2_datacols_cov) == FALSE){
        stop("It looks like you have specified covariates for the first twin, but not for the second twin of a family. 
             This is not possible! When using covariates, these should be measured in both first and second twin of each family. 
             Please use twin1_datacols_cov to specifiy the columns of the first twin and twin2_datacols_cov 
             for the second twin of each family.")
    }
    
    #==========================================================
    # For JAGS analysis
    #==========================================================
    #Choose initial values if not specified otherwise by user: 
    if (is.na(inits) && ge == FALSE){
        if (decomp_model == "ACE"){
            inits = list(tau_a = 2, tau_c = 3, tau_e = 5)
        } else if (decomp_model == "ADE"){
            inits = list(tau_a = 2, tau_d = 3, tau_e = 5)
        } else {
            inits = list(tau_a = 2, tau_e = 2)
        }
        
    } 
    
    if (is.na(inits) && ge == TRUE){
        if (decomp_model == "ACE"){
            inits = list(tau_a = 2, tau_c = 3, beta0 = log(0.2), beta1 = 0)
        } else if (decomp_model == "ADE"){
            inits = list(tau_a = 2, tau_d = 3, beta0 = log(0.2), beta1 = 0)
        } else {
            inits = list(tau_a = 2, beta0 = log(0.2), beta1 = 0)
        }
    } 
    
    #==========================================================
    # Covariates
    #==========================================================  
    #Test if covariates are used:
    covariates = NA
    
    if (any(is.na(twin1_datacols_cov)) == FALSE){
        covariates = TRUE
        
        #Use only data with complete cases, meaning that cases with missing values
        #on covariate data will be ommited from the data anlysis, even if that means
        #that we do not use all known phenotypic values.
        if (any(is.na(data_mz[,c(twin1_datacols_cov, twin2_datacols_cov)])) == TRUE){
            data_mz = data_mz[complete.cases(data_mz[,c(twin1_datacols_cov, twin2_datacols_cov)]),] 
        }
        
        if (any(is.na(data_dz[,c(twin1_datacols_cov, twin2_datacols_cov)])) == TRUE){
            data_dz = data_dz[complete.cases(data_dz[,c(twin1_datacols_cov, twin2_datacols_cov)]),] 
        }

        #Select covariate data:
        X_mz_twin1 = as.matrix(data_mz[,twin1_datacols_cov])
        X_mz_twin2 = as.matrix(data_mz[,twin2_datacols_cov])
        X_dz_twin1 = as.matrix(data_dz[,twin1_datacols_cov])
        X_dz_twin2 = as.matrix(data_dz[,twin2_datacols_cov])
        
        #Print warning message when categorical data is used: 
        if(is.categorical(X_mz_twin1) == TRUE || is.categorical(X_mz_twin2) == TRUE || 
           is.categorical(X_dz_twin1) == TRUE || is.categorical(X_dz_twin2) == TRUE){
            warning("It looks like (part) of your covariate data is of categorical nature. 
                    Consider computing dummy variables! (this is not done automatically by the BayesTwin package")
        }
        
    } else {covariates = FALSE}
    
    
    #Select phenotypic data
    data_mz = data_mz[,c(twin1_datacols_p, twin2_datacols_p)]
    data_dz = data_dz[,c(twin1_datacols_p, twin2_datacols_p)]
    
    #==========================================================
    #==========================================================
    # Subroutines for AE model
    #==========================================================
    #I. With covariates: 
    if(covariates == TRUE && decomp_model == "AE"){
        output = irt_ae_cov(data_mz = data_mz, data_dz = data_dz, 
                            X_mz_twin1 = X_mz_twin1, 
                            X_mz_twin2 = X_mz_twin2, 
                            X_dz_twin1 = X_dz_twin1,
                            X_dz_twin2 = X_dz_twin2,
                            n_burnin = n_burnin, n_iter = n_iter, 
                            ge = ge, irt_model = irt_model, 
                            N_cov = N_cov, var_prior = var_prior,
                            n_chains = n_chains, fit_stats = fit_stats,
                            inits = inits, Nk = Nk)
    }
    
    #II. Without covariates: 
    if(covariates == FALSE && decomp_model == "AE"){
        output = irt_ae(data_mz = data_mz, data_dz = data_dz, 
                        n_burnin = n_burnin, n_iter = n_iter,
                        ge = ge, irt_model = irt_model,
                        var_prior = var_prior,
                        n_chains = n_chains, fit_stats = fit_stats,
                        inits = inits, Nk = Nk)              
    } 
    
    
    #==========================================================
    # Subroutines for ACE model
    #==========================================================
    #I. With covariates: 
    if(covariates == TRUE && decomp_model == "ACE"){
        output = irt_ace_cov(data_mz = data_mz, data_dz = data_dz, 
                             X_mz_twin1 = X_mz_twin1, 
                             X_mz_twin2 = X_mz_twin2, 
                             X_dz_twin1 = X_dz_twin1,
                             X_dz_twin2 = X_dz_twin2,
                             n_burnin = n_burnin, n_iter = n_iter, 
                             ge = ge, irt_model = irt_model, 
                             N_cov = N_cov,
                             var_prior = var_prior,
                             n_chains = n_chains, fit_stats = fit_stats,
                             inits = inits, Nk = Nk)
    }
    
    #II. Without covariates: 
    if(covariates == FALSE && decomp_model == "ACE"){
        output = irt_ace(data_mz = data_mz, data_dz = data_dz, 
                         n_burnin = n_burnin, n_iter = n_iter,
                         ge = ge, irt_model = irt_model,
                         var_prior = var_prior,
                         n_chains = n_chains, fit_stats = fit_stats,
                         inits = inits, Nk = Nk)              
    } 
    
    
    #==========================================================
    #==========================================================
    # Subroutines for ADE model
    #==========================================================
    #I. With covariates: 
    if(covariates == TRUE && decomp_model == "ADE"){
        output = irt_ade_cov(data_mz = data_mz, data_dz = data_dz, 
                             X_mz_twin1 = X_mz_twin1, 
                             X_mz_twin2 = X_mz_twin2, 
                             X_dz_twin1 = X_dz_twin1,
                             X_dz_twin2 = X_dz_twin2,
                             n_burnin = n_burnin, n_iter = n_iter, 
                             ge = ge, irt_model = irt_model, 
                             N_cov = N_cov, var_prior = var_prior,
                             n_chains = n_chains, fit_stats = fit_stats,
                             inits = inits, Nk = Nk)
    }
    
    #II. Without covariates: 
    if(covariates == FALSE && decomp_model == "ADE"){
        output = irt_ade(data_mz = data_mz, data_dz = data_dz, 
                         n_burnin = n_burnin, n_iter = n_iter,
                         var_prior = var_prior,
                         ge = ge, irt_model = irt_model,
                         n_chains = n_chains, fit_stats = fit_stats, 
                         inits = inits, Nk = Nk)              
    } 

    #==========================================================
    # Output
    #==========================================================
    if(ge == TRUE){
        cat("\n Analysis completed.") 
        cat("\n After a burn-in period of", n_burnin, "iterations,")
        cat("the characterisation of the posterior distribution for \n", 
            "the model parameters was based on an additional", n_iter, "iterations from", n_chains, "Markov chain(s).")
        cat("\n The analysis was run using an", decomp_model, "model", "with an integrated", irt_model, "IRT model.")
        cat("\n The model was analysed with genotype-environment interaction effect.\n")
        
    } else {
        cat("\n Analysis completed.") 
        cat("\n After a burn-in period of", n_burnin, "iterations,")
        cat("the characterisation of the posterior distribution for \n", 
            "the model parameters was based on an additional", n_iter, "iterations from", n_chains, "Markov chain(s).")
        cat("\n The analysis was run using an", decomp_model, "model", "with an integrated", irt_model, "IRT model.")
        cat("\n The model was analysed without genotype-environment interaction effect.\n")
    }
    cat("=================================================================")
    
    #Remind user of convergence-issue
    cat("\n WARNING! \n")
    cat("It is important to check that the MCMC algorithm has converged to the posterior distribution! \n")
    cat("To check convergence, you can use the posterior samples and the plot function with type = 'trace'. \n",
        "For example, call ACE_model$samples_var_a to obtain the samples for VAR(A) when the output of the analysis was saved in the object ACE_model. \n",
        "You can then use plot(ACE_model$samples_vara_a, type = 'trace') to check the convergence of the VAR(A) parameter.\n")
    cat("=================================================================")
    cat("\n For posterior means and standard deviations of all variance components, see below or use the summary function \n",
        "For example: summary(ACE_model$results) when the Ouptut of the analysis was saved in the object ACE_model).\n")
    cat("=================================================================")
    cat("\n")
    
    #Print results on the fly:
    if (covariates == TRUE){print(output$results_b)}
    print(output$results)
    
    #Return output, but do not print it
    class(output) = "bayestwin"
    output
}