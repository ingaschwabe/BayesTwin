#==========================================================
# irt_ae_cov.R
# Subroutine for master function twinUniv()
#
# Analysis of twin data under the ACE model using item 
# (phenotypic) data and environmental covariates
# Possible IRT models: 1PL, 2PL, GPCM, PCM
# 
# This function writes a JAGS script for the model and 
# draws from the posterior distribution and then organizes 
# the output. Estimation of genotype by environment 
# interaction is optional (by defining ge = TRUE) when 
# calling the master function twinUniv)
# BayesTwin package
#==========================================================

irt_ae_cov <- function(data_mz, data_dz, n_burnin, n_iter, ge, irt_model, N_cov,
                       var_prior, n_chains, fit_stats){
    
    #Make boolean variable to create model string with the 
    #right IRT model 
    PL_1 = FALSE; PL_2 = FALSE; GPCM = FALSE; PCM = FALSE
    
    if(irt_model == "1PL"){
        PL_1 = TRUE
    } else if (irt_model == "2PL"){
        PL_2 = TRUE
    } else if (irt_model == "GPCM"){
        GPCM = TRUE
    } else {
        PCM = TRUE
    }
    
    #Make boolean variable to create model string with the right prior
    INV_GAMMA = FALSE
    if(var_prior == "INV_GAMMA"){
        INV_GAMMA = TRUE
    }
    
    #Determine number of twin pairs
    n_mz <- nrow(data_mz) ; n_dz <- nrow(data_dz)
    
    # determine number of phenotypic items
    n_items <- ncol(data_mz)/2
    Nk <- length(unique(data_mz[,1]))
    
    #==========================================================
    # I. Write JAGS model file
    #==========================================================
    jags_model_irt_ae_cov <- paste("model{
        ##MZ twins
        for (fam in 1:n_mz){
            a_mz[fam] ~ dnorm(mu, tau_a)
                                
        ",ifelse(ge,"
            tau_e_mz[fam] <- 1/(exp(beta0 + (beta1*a_mz[fam])))
    	    pheno_mz[fam,1] ~ dnorm(a_mz[fam] + inprod(X_mz_twin1[fam,], b[1:N_cov]), tau_e_mz[fam]) 
		    pheno_mz[fam,2] ~ dnorm(a_mz[fam] + inprod(X_mz_twin2[fam,], b[1:N_cov]), tau_e_mz[fam])
         ","
            pheno_mz[fam,1] ~ dnorm(a_mz[fam] + inprod(X_mz_twin1[fam,], b[1:N_cov]), tau_e) 
		    pheno_mz[fam,2] ~ dnorm(a_mz[fam] + inprod(X_mz_twin2[fam,], b[1:N_cov]), tau_e)
        "),"
        ",ifelse(PL_1,"
            #1pl model twin1
            for (j in 1:n_items){
                logit(p[fam,j]) <- pheno_mz[fam,1] - item_b[j]
                data_mz[fam,j] ~ dbern(p[fam,j])
            }          
                                         
            #1pl model twin2
            for (j in (n_items+1):(2*n_items)){
                logit(p[fam,j]) <- pheno_mz[fam,2] - item_b[j-n_items]
                data_mz[fam,j] ~ dbern(p[fam,j])
            }
        ","
        "),"
        ",ifelse(PL_2,"
            #2pl model twin1
            for (j in 1:n_items){
                logit(p[fam,j]) <- alpha[j] * (pheno_mz[fam,1] - item_b[j])
                data_mz[fam,j] ~ dbern(p[fam,j])
            }          
                                         
            #1pl model twin2
            for (j in (n_items+1):(2*n_items)){
                logit(p[fam,j]) <- alpha[j - n_items] * (pheno_mz[fam,2] - item_b[j-n_items])
                data_mz[fam,j] ~ dbern(p[fam,j])
            }
        ","
        "),"
        ",ifelse(GPCM,"
            #GPCM for twin 1
            for (j in 1:n_items){
                for (k in 1:Nk){
                    eta[fam,j,k] <- alpha[j] * (pheno_mz[fam,1] - item_b[j,k])
                    psum [fam,j,k] <- sum(eta[fam,j,1:k])
                    exp.psum[fam,j,k] <- exp(psum[fam,j,k])
                    prob [fam,j,k] <- exp.psum [fam,j,k]/sum(exp.psum[fam,j,1:Nk])
                } 
            }
                                         
            #GPCM for twin 2
            for (j in (n_items+1):(2*n_items)){
                for (k in 1:Nk){
                    eta[fam,j,k] <- alpha [j-n_items] * (pheno_mz[fam,2] - item_b[j-n_items,k])
                    psum[fam,j,k] <- sum(eta[fam,j,1: k])
                    exp.psum[fam,j,k] <- exp(psum[fam,j,k])
                    prob[fam,j,k] <- exp.psum[fam,j,k]/sum(exp.psum[fam,j,1:Nk])
                } 
            }
            
            for (j in 1:(2*n_items)){
                data_mz[fam,j] ~ dcat(prob[fam,j,1:Nk]) # multinomial distr of data
            }
        ","
        "),"
        ",ifelse(PCM,"
            #PCM for twin 1
            for (j in 1:n_items){
                for (k in 1:Nk){
                    eta[fam,j,k] <- pheno_mz[fam,1] - item_b[j,k]
                    psum [fam,j,k] <- sum(eta[fam,j,1:k])
                    exp.psum[fam,j,k] <- exp(psum[fam,j,k])
                    prob [fam,j,k] <- exp.psum [fam,j,k]/sum(exp.psum[fam,j,1:Nk])
                } 
            }
                                         
            #PCM for twin 2
            for (j in (n_items+1):(2*n_items)){
                for (k in 1:Nk){
                    eta[fam,j,k] <- pheno_mz[fam,2] - item_b [j-n_items,k]
                    psum[fam,j,k] <- sum(eta[fam,j,1: k])
                    exp.psum [fam,j,k] <- exp(psum[fam,j,k])
                    prob[fam,j,k] <- exp.psum[fam,j,k]/sum(exp.psum[fam,j,1:Nk])
                } 
            }
            
            for (j in 1:(2*n_items)){
                data_mz[fam,j] ~ dcat(prob[fam,j,1:Nk]) # multinomial distr of data
            }
        ",""),"
        }
        ##DZ twins
        for(fam in 1:n_dz){
            a0_dz[fam] ~ dnorm(mu, doubletau_a)
        ",ifelse(ge,"
            for (twin in 1:2){            							
                a_dz[fam,twin] ~ dnorm(a0_dz[fam], doubletau_a)
                tau_e_dz[fam,twin] <- 1/(exp(beta0 + (beta1*a_dz[fam,twin])))
            }
            pheno_dz[fam,1] ~ dnorm(a_dz[fam,1] + inprod(X_dz_twin1[fam,], b[1:N_cov]), tau_e_dz[fam,1]) 
            pheno_dz[fam,2] ~ dnorm(a_dz[fam,2] + inprod(X_dz_twin2[fam,], b[1:N_cov]), tau_e_dz[fam,2])      
        ","
            for (twin in 1:2){										
                a_dz[fam,twin] ~ dnorm(a0_dz[fam], doubletau_a)
            }
            pheno_dz[fam,1] ~ dnorm(a_dz[fam,1] + inprod(X_dz_twin1[fam,], b[1:N_cov]), tau_e) 
            pheno_dz[fam,2] ~ dnorm(a_dz[fam,2] + inprod(X_dz_twin2[fam,], b[1:N_cov]), tau_e)
        "),"
                                
        ",ifelse(PL_1,"
        #1pl model twin1 (DZ)
        for (j in 1:n_items){
            logit(p2[fam,j]) <- pheno_dz[fam,1] - item_b[j]
            data_dz[fam,j] ~ dbern(p2[fam,j])
        }
        #1pl model twin2 (DZ)
        for (j in (n_items+1):(2*n_items)){
            logit(p2[fam,j]) <- pheno_dz[fam,2] - item_b[j-n_items]
            data_dz[fam,j] ~ dbern(p2[fam,j])
        }
        ","
        "),"
        ",ifelse(PL_2,"
        #2pl model twin1 (DZ)
        for (j in 1:n_items){
            logit(p2[fam,j]) <- alpha[j]*(pheno_dz[fam,1] - item_b[j])
                                         data_dz[fam,j] ~ dbern(p2[fam,j])
        }
                                         
        #2pl model twin2 (DZ)
        for (j in (n_items+1):(2*n_items)){
            logit(p2[fam,j]) <- alpha[j - n_items]*(pheno_dz[fam,2] - item_b[j-n_items])
            data_dz[fam,j] ~ dbern(p2[fam,j])
        }
        ","
        "),"
        ",ifelse(GPCM,"
        for (j in 1:n_items){
            for (k in 1:Nk){
                etadz[fam,j,k] <- alpha[j]*(pheno_dz[fam,1] - item_b[j,k])
                psumdz [fam,j,k] <- sum(etadz[fam,j,1:k])
                exp.psumdz[fam,j,k] <- exp(psumdz[fam,j,k])
                probdz[fam,j,k] <- exp.psumdz[fam,j,k]/sum(exp.psumdz[fam,j,1:Nk])
            } 
        }
        for (j in (n_items+1):(2*n_items)){
            for (k in 1:Nk){
                etadz[fam,j,k] <- alpha[j-n_items] * (pheno_dz[fam,2] - item_b[j-n_items,k])
                psumdz[fam,j,k] <- sum(etadz[fam,j,1:k])
                exp.psumdz[fam,j,k] <- exp(psumdz[fam,j,k])
                probdz[fam,j,k] <- exp.psumdz[fam,j,k]/sum(exp.psumdz[fam,j,1:Nk])
            } 
        }
        for (j in 1:(2*n_items)){
            data_dz[fam,j] ~ dcat(probdz[fam,j,1:Nk]) # multinomial distr of data
        }
        ","
        "),"
        ",ifelse(PCM,"
        for (j in 1:n_items){
            for (k in 1:Nk){
                etadz[fam,j,k] <- pheno_dz[fam,1] - item_b[j,k]
                psumdz [fam,j,k] <- sum(etadz[fam,j,1:k])
                exp.psumdz[fam,j,k] <- exp(psumdz[fam,j,k])
                probdz[fam,j,k] <- exp.psumdz[fam,j,k]/sum(exp.psumdz[fam,j,1:Nk])
            } 
        }
        for (j in (n_items+1):(2*n_items)){
            for (k in 1:Nk){
                etadz[fam,j,k] <- pheno_dz[fam,2] - item_b[j-n_items,k]
                psumdz[fam,j,k] <- sum(etadz[fam,j,1: k])
                exp.psumdz[fam,j,k] <- exp(psumdz[fam,j,k])
                probdz[fam,j,k] <- exp.psumdz[fam,j,k]/sum(exp.psumdz[fam,j,1:Nk])
            } 
        }
        for (j in 1:(2*n_items)){
            data_dz[fam,j] ~ dcat(probdz[fam,j,1:Nk]) # multinomial distr of data
        }
        ","
        "),"        
        }
        #Priors
        mu <- 0 #to identify scale 
        doubletau_a <- 2*tau_a
        b[1:N_cov] ~ dmnorm(mu_b[1:N_cov], tau_b[,]) #in R: mu_b = rep(0,N_cov), tau_b = diag(1, N_cov)
                                
        #Priors
        ",ifelse(INV_GAMMA,"
        tau_a ~ dgamma(1,1) 
        tau_e ~ dgamma(1,1) #not used when ge = TRUE
        ","
        tau_aa ~ dunif(0,100)
        tau_e ~ dunif(0,100) #not used when ge = TRUE
        "),"
        
        ",ifelse(PL_1,"
        for (j in 1:n_items){
            item_b[j] ~ dnorm(0,.1)
        }
        "," "),"
        ",ifelse(PL_2,"
            alpha[1] <- 1 #fix first item to identify scale
            for (j in 2:n_items){
                alpha[j] ~ dlnorm(0, .1)
            }
            for (j in 1:n_items){
                item_b[j] ~ dnorm(0,.1)
            }
        "," "),"
        ",ifelse(GPCM,"
            alpha[1] <- 1
            for (j in 2:n_items){
                alpha[j] ~ dlnorm(0, .1)
            }
            for (j in 1:n_items){
                item_b[j,1] <- 0.0 #first threshold = 0
                                         
                for (k in 2:Nk){
                    item_b[j,k] ~ dnorm (0, .1)
                }
            }
        "," "),"
        ",ifelse(PCM,"
        for (j in 1:n_items){
            item_b[j,1] <- 0.0 #first threshold = 0
            for (k in 2:Nk){
                item_b[j,k] ~ dnorm (0, .1)
            }
        }
        "," "),"
        ",ifelse(ge,"
        beta0 ~ dnorm(-1,.5)
        beta1 ~ dnorm(0,.1)",
        ""),"
}")

    jags_file_irt_ae_cov <- tempfile(fileext=".txt")
    write(jags_model_irt_ae_cov,jags_file_irt_ae_cov)
    #writeLines(jags_model_irt_ae_cov, con = "file.txt", sep = "\n", useBytes = FALSE)
    
    #==========================================================
    # II. Run JAGS analysis
    #==========================================================   
    inits = list(tau_a = 2)
    
    if (PCM == TRUE || GPCM == TRUE){
        jags_data <- list(data_mz, data_dz, n_mz, n_dz, n_items, Nk, rep(0, N_cov), diag(1,N_cov), N_cov,
                          X_mz_twin1, X_mz_twin2, X_dz_twin1, X_dz_twin2)
        names(jags_data)<- c("data_mz", "data_dz", "n_mz", "n_dz", "n_items", "Nk", "mu_b", "tau_b", "N_cov",
                             "X_mz_twin1", "X_mz_twin2", "X_dz_twin1", "X_dz_twin2") 
    } else {
        jags_data <- list(data_mz, data_dz, n_mz, n_dz, n_items, rep(0,N_cov), diag(1, N_cov), N_cov,
                          X_mz_twin1, X_mz_twin2, X_dz_twin1, X_dz_twin2)
        names(jags_data)<- c("data_mz", "data_dz", "n_mz", "n_dz", "n_items", "mu_b", "tau_b", "N_cov",
                             "X_mz_twin1", "X_mz_twin2", "X_dz_twin1", "X_dz_twin2")
    }
    
    jags <- jags.model(jags_file_irt_ae_cov, jags_data, inits, n_chains = 1, quiet=FALSE)
    update(jags, n_burnin)

    #Output, dependent on fit_stats, GE and IRT model
    if (fit_stats == FALSE){
        if (ge == FALSE && PL_1 == TRUE || ge == FALSE && PCM == TRUE){
            out <- jags.samples(jags, c("tau_a", "tau_e", "item_b", "b"), n_iter)
        } else if (ge == TRUE && PL_1 == TRUE || ge == TRUE && PCM == TRUE){
            out <- jags.samples(jags, c("tau_a", "beta0", "beta1", "item_b", "b"), n_iter)
        } else if (ge == FALSE && PL_2 == TRUE || ge == FALSE && GPCM == TRUE){
            out <- jags.samples(jags, c("tau_a", "tau_e", "item_b", "alpha", "b"), n_iter)
        } else {
            out <- jags.samples(jags, c("tau_a", "beta0", "beta1", "item_b", "alpha", "b"), n_iter)
        }
    } else {
        if (ge == FALSE && PL_1 == TRUE || ge == FALSE && PCM == TRUE){
            out <- jags.samples(jags, c("tau_a", "tau_e", "item_b", "b"), n_iter)
            out_dic <- dic.samples(jags, c("tau_a", "tau_e", "item_b", "b"), n_iter)
        } else if (ge == TRUE && PL_1 == TRUE || ge == TRUE && PCM == TRUE){
            out <- jags.samples(jags, c("tau_a", "beta0", "beta1", "item_b", "b"), n_iter)
            out_dic <- dic.samples(jags, c("tau_a", "beta0", "beta1", "item_b", "b"), n_iter)
        } else if (ge == FALSE && PL_2 == TRUE || ge == FALSE && GPCM == TRUE){
            out <- jags.samples(jags, c("tau_a", "tau_e", "item_b", "alpha", "b"), n_iter)
            out_dic <- dic.samples(jags, c("tau_a", "tau_e", "item_b", "alpha", "b"), n_iter)
        } else {
            out <- jags.samples(jags, c("tau_a", "beta0", "beta1", "item_b", "alpha", "b"), n_iter)
            out_dic <- dic.samples(jags, c("tau_a", "beta0", "beta1", "item_b", "alpha", "b"), n_iter)
        }
    }
    
    #==========================================================
    # III. Organize results  
    #==========================================================
    ##### Beta (regression coefficients) parameters: 
    list_b = list()
    
    for (i in 1:N_cov){
        #Fill list spaces with parameters
        list_b[[((i-1)*4+1)]] = out$b[i,,1]
        list_b[[((i-1)*4+2)]] = mean(out$b[i,,1])
        list_b[[((i-1)*4+3)]] = sd(out$b[i,,1])
        list_b[[((i-1)*4+4)]] = HPD(out$b[i,,1], 0.95)
        
        #Change parameter names: 
        names(list_b)[((i-1)*4+1)] = paste("samples_b_",i,sep="")
        names(list_b)[((i-1)*4+2)] = paste("b_",i,sep="")
        names(list_b)[((i-1)*4+3)] = paste("sd_b_",i,sep="")
        names(list_b)[((i-1)*4+4)] = paste("hpd_b_",i,sep="")
        
        #Change class of sample object:
        class(list_b)[((i-1)*4+1)] <- "bayestwin"
    }
    
    if(ge == FALSE){
        #Save samples
        samples_var_a = 1/out$tau_a[,,1]; 
        samples_var_e = 1/out$tau_e[,,1]; samples_item_b = out$item_b[,,1]
        if(PL_2 == TRUE || GPCM == TRUE){samples_alpha = out$alpha[,,1]}
        
        #Calculate mean values: 
        var_a = mean(samples_var_a); 
        item_b = apply(samples_item_b, 1, mean)
        if(PL_2 == TRUE || GPCM == TRUE){alpha = apply(samples_alpha, 1, mean)}
        
        #Calculate SDs: 
        sd_var_a = sd(samples_var_a); 
        sd_item_b = apply(samples_item_b, 1, sd)
        if(PL_2 == TRUE || GPCM == TRUE){alpha = apply(samples_alpha, 1, mean)}
        
        #Calculate HPD: 
        hpd_var_a = HPD(samples_var_a, 0.95); 
        hpd_var_e = HPD(samples_var_e, 0.95); 
        hpd_item_b = apply(samples_item_b, 1, function (x) HPD(x, 0.95))
        if(PL_2 == TRUE || GPCM == TRUE){hpd_alpha = apply(samples_alpha, 1, function (x) HPD(x, 0.95))}
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, var_e, sd_var_e), 2, 2)
        colnames(results) <- c("varA", "varE")
        rownames(results) <- c("Posterior mean","Posterior standard deviation")
        results = as.table(results) 
        
        #And save output in a list: 
        if (fit_stats == FALSE){
            if(PL_1 == TRUE || PCM == TRUE){
                output = list(results = results, samples_var_a = samples_var_a,  
                              samples_var_e = samples_var_e, samples_item_b = samples_item_b,
                              hpd_var_a = hpd_var_a, hpd_var_e = hpd_var_e, 
                              hpd_item_b = hpd_item_b)
                
                output = c(output, list_b)
            } else {
                output = list(results = results, samples_var_a = samples_var_a,  
                              samples_var_e = samples_var_e, samples_item_b = samples_item_b,
                              hpd_var_a = hpd_var_a, hpd_var_e = hpd_var_e, 
                              hpd_item_b = hpd_item_b, samples_alpha = samples_alpha, hpd_alpha = hpd_alpha)
                output = c(output, list_b)
                class(output$samples_alpha) <- "bayestwin"
            }
        } else {
            if(PL_1 == TRUE || PCM == TRUE){
                output = list(results = results, samples_var_a = samples_var_a,  
                              samples_var_e = samples_var_e, samples_item_b = samples_item_b,
                              hpd_var_a = hpd_var_a, hpd_var_e = hpd_var_e, 
                              hpd_item_b = hpd_item_b, dic = out_dic)
                
                output = c(output, list_b)
            } else {
                output = list(results = results, samples_var_a = samples_var_a,  
                              samples_var_e = samples_var_e, samples_item_b = samples_item_b,
                              hpd_var_a = hpd_var_a, hpd_var_e = hpd_var_e, 
                              hpd_item_b = hpd_item_b, samples_alpha = samples_alpha, hpd_alpha = hpd_alpha,
                              dic = out_dic)
                output = c(output, list_b)
                class(output$samples_alpha) <- "bayestwin"
            }
        }
        
        #Change class of objects in order to use right plot method 
        class(output) <- c("bayestwin","list")    
        class(output$samples_var_a) <- "bayestwin"
        class(output$samples_var_e) <- "bayestwin"
        class(output$samples_item_b) <- "bayestwin"
        
    } else {
        #Save samples
        samples_var_a = 1/out$tau_a[,,1]; 
        samples_beta0 = out$beta0[,,1]; samples_beta1 = out$beta1[,,1]; samples_item_b = out$item_b[,,1]
        if(PL_2 == TRUE || GPCM == TRUE){samples_alpha = out$alpha[,,1]}
        
        #Calculate mean values: 
        var_a = mean(samples_var_a); 
        beta0 = mean(samples_beta0); beta1 = mean(samples_beta1)
        item_b = apply(samples_item_b, 1, mean)
        if(PL_2 == TRUE || GPCM == TRUE){alpha = apply(samples_alpha, 1, mean)}
        
        #Calculate SDs: 
        sd_var_a = sd(samples_var_a); 
        sd_beta0 = sd(samples_beta0); sd_beta1 = sd(samples_beta1)
        sd_item_b = apply(samples_item_b, 1, sd)
        if(PL_2 == TRUE || GPCM == TRUE){alpha = apply(samples_alpha, 1, mean)}
        
        #Calculate HPD: 
        hpd_var_a = HPD(samples_var_a, 0.95); 
        hpd_beta0 = HPD(samples_beta0, 0.95); hpd_beta1 = HPD(samples_beta1, 0.95)
        hpd_item_b = apply(samples_item_b, 1, function (x) HPD(x, 0.95))
        if(PL_2 == TRUE || GPCM == TRUE){hpd_alpha = apply(samples_alpha, 1, function (x) HPD(x, 0.95))}
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, beta0, sd_beta0, beta1, sd_beta1), 2, 3)
        colnames(results) <- c("varA","beta0", "beta1")
        rownames(results) <- c("Posterior mean","Posterior standard deviation")
        results = as.table(results) 
        
        #And save output in a list: 
        if (fit_stats == FALSE){
            if(PL_1 == TRUE || PCM == TRUE){
                output = list(results = results, samples_var_a = samples_var_a,  
                              samples_beta0 = samples_beta0, samples_beta1 = samples_beta1, samples_item_b = samples_item_b,
                              hpd_var_a = hpd_var_a, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1,
                              hpd_item_b = hpd_item_b)
                output = c(output, list_b)
            } else {
                output = list(results = results, samples_var_a = samples_var_a,  
                              samples_beta0 = samples_beta0, samples_beta1 = samples_beta1, samples_item_b = samples_item_b,
                              hpd_var_a = hpd_var_a, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1,
                              hpd_item_b = hpd_item_b, hpd_alpha = hpd_alpha, samples_alpha = samples_alpha)
                
                output$samples_alpha <- "bayestwin"
            }    
        } else {
            if(PL_1 == TRUE || PCM == TRUE){
                output = list(results = results, samples_var_a = samples_var_a,  
                              samples_beta0 = samples_beta0, samples_beta1 = samples_beta1, samples_item_b = samples_item_b,
                              hpd_var_a = hpd_var_a, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1,
                              hpd_item_b = hpd_item_b, dic = out_dic)
                output = c(output, list_b)
            } else {
                output = list(results = results, samples_var_a = samples_var_a,  
                              samples_beta0 = samples_beta0, samples_beta1 = samples_beta1, samples_item_b = samples_item_b,
                              hpd_var_a = hpd_var_a, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1,
                              hpd_item_b = hpd_item_b, hpd_alpha = hpd_alpha, samples_alpha = samples_alpha, dic = out_dic)
                
                output$samples_alpha <- "bayestwin"
            }
        }
        
        #Change class of objects in order to use right plot method 
        class(output) <- c("bayestwin","list")    
        class(output$samples_var_a) <- "bayestwin"
        class(output$samples_beta0) <- "bayestwin"
        class(output$samples_beta1) <- "bayestwin"     
        class(output$samples_item_b) <- "bayestwin"        
    }
    
    return(output)
}