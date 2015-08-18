#==========================================================
# twinUniv.R
# Master function that calls subroutines for univariate
# analyses with one phenotype. 
# BayesTwin package
#==========================================================
#setwd("C:/Users/schwabei/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("/Users/stephanievandenberg/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("C:/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")

twinUniv = function(data_mz, data_dz, 
                    twin1_datacols_p, twin2_datacols_p,
                    twin1_datacols_cov = NA, twin2_datacols_cov = NA, #####KOMT LATER
                    decomp_model = "ACE",
                    irt_model = "1PL",
                    cov_imputation = FALSE,  ##### KOMT LATER
                    ge = FALSE,
                    n_iter = 10000, n_burnin = 5000,
                    var_prior = "INV_GAMMA"){
    #Load HPD function: 
    source("HPD.R")                                                         
    
    #Calculate number of MZ and DZ twins: 
    n_mz = nrow(data_mz)
    n_dz = nrow(data_dz)
    
    #Make sure that first columns are only twin 1, next the twin 2 columns
    data_mz = data_mz[,c(twin1_datacols_p, twin2_datacols_p)]
    data_dz = data_dz[,c(twin1_datacols_p, twin2_datacols_p)]
    
    #Test if its item data: 
    if(length(twin1_datacols_p) == 1){
        print("Error! You need phenotypic data on item level to run the analysis!")
    } 
    
    #When cov_imputation = FALSE, use only complete cases: 
    if(is.na(twin1_datacols_cov) && cov_imputation == FALSE){
        data_mz = na.omit(data_mz)
        data_dz = na.omit(data_dz)
    }
    
    #Test if covariates are used:
    covariates = NA
    
    #Option I: Use covaritaes without missing value imputation: 
    if (is.na(twin1_datacols_cov) == FALSE && cov_imputation == TRUE){
        X_mz_twin1 = data_mz[,twin1_datacols_cov]
        X_mz_twin2 = data_mz[,twin2_datacols_cov]
        X_dz_twin1 = data_dz[,twin1_datacols_cov]
        X_dz_twin2 = data_dz[,twin2_datacols_cov]
        covariates = TRUE
        
        #We have to know which covariates are dichotomous and which continous: 
        is.binary = function(x){
            length(unique(x)) <= 3 #allow for values 0,1 and NA. 
        }
        
        dich = NA; cont = NA
        
                
        #Apply on all covariate data: (assuming that the same covariates cannot 
        #be dichotomous for twin 1 and continous for twin 2 + assuming that 
        #the same covaraites are dichotomous for MZ and DZ twins
        if(length(twin1_datacols_cov) > 1){
            dich_cov = which(apply(X_mz_twin1, 2, is.binary) == TRUE)
            cont_cov = which(apply(X_mz_twin1, 2, is.binary) == FALSE)
            if(length(dich_cov) >= 1){
                dich = TRUE
            } else{
                dich = FALSE
            }
            
            if(length(cont_cov) >= 1){
                cont = TRUE
            } else {
                cont = FALSE
            }
                
        } else {
            dich_cov = is.binary(X_mz_twin1)
            if(dich_cov == FALSE){
                dich = FALSE; cont = TRUE 
                cont_cov = twin1_datacols_cov
            } else {
                dich = TRUE
                dich_cov = twin1_datacols_cov
            }
        }
                        
    #Option II: Use covariates with missing value imputation:
    } else if (is.na(twin1_datacols_cov) == FALSE && cov_imputation == FALSE){
        
        #Use only data with complete cases, meaning that cases with missing values
        #on covariate data will be ommited from the data anlysis, even if that means
        #that we do not use all known phenotypic values.
        
        #Which cases have missing values on covariate data: 
        data_mz[1:n_mz, ncol(data_mz)] = 1:n_mz            
        without_missings = na.omit(data_mz[,c(twin1_datacols_cov, twin2_datacols_cov)])
        data_mz = data_mz[without_missings[,ncol(data_mz)],] 
        
        data_dz[1:n_mz, ncol(data_dz)] = 1:n_dz         
        without_missings = na.omit(data_dz[,c(twin1_datacols_cov, twin2_datacols_cov)])
        data_dz = data_dz[without_missings[,ncol(data_dz)],] 
        
        #Select only covariate data on basis of new dataset without missings: 
        X_mz_twin1 = data_mz[,twin1_datacols_cov]
        X_mz_twin2 = data_mz[,twin2_datacols_cov]
        X_dz_twin1 = data_dz[,twin1_datacols_cov]
        X_dz_twin2 = data_dz[,twin2_datacols_cov]
        covariates = TRUE
    
    } else {covariates = FALSE}
    
    #==========================================================
    #==========================================================
    # Subroutines for AE model
    #==========================================================
    #Option for IRT model (1pl/2pl or (G)PCM) passed to function
    #with irt_model object and option fo genotype by environment
    #interaction by ge = FALSE/TRUE
    
    #I. With covariates: 
    if(covariates == TRUE && decomp_model == "AE"){
        source("irt_ae_cov.R")
        output = irt_ae_cov(data_mz = data_mz, data_dz = data_dz, 
                            n_burnin = n_burnin, n_iter = n_iter, 
                            ge = ge, irt_model = irt_model, 
                            N_cov = N_cov, var_prior = var_prior)
    }
    
    #II. Without covariates: 
    if(covariates == FALSE && decomp_model == "AE"){
        source("irt_ae.R")
        output = irt_ae(data_mz = data_mz, data_dz = data_dz, 
                        n_burnin = n_burnin, n_iter = n_iter,
                        ge = ge, irt_model = irt_model,
                        var_prior = var_prior)              
    } 
    
    #==========================================================
    # Subroutines for ACE model
    #==========================================================
    #Option for IRT model (1pl/2pl or (G)PCM) passed to function
    #with irt_model object and option fo genotype by environment
    #interaction by ge = FALSE/TRUE
    
    #I. With covariates: 
    if(covariates == TRUE && decomp_model == "ACE"){
        source("irt_ace_cov.R")
        output = irt_ace_cov(data_mz = data_mz, data_dz = data_dz, 
                            n_burnin = n_burnin, n_iter = n_iter, 
                            ge = ge, irt_model = irt_model, 
                            N_cov = N_cov,
                            var_prior = var_prior)
    }
    
    #II. Without covariates: 
    if(covariates == FALSE && decomp_model == "ACE"){
        source("irt_ace.R")
        output = irt_ace(data_mz = data_mz, data_dz = data_dz, 
                        n_burnin = n_burnin, n_iter = n_iter,
                        ge = ge, irt_model = irt_model,
                        var_prior = var_prior)              
    } 
    

    #==========================================================
    #==========================================================
    # Subroutines for ADE model
    #==========================================================
    #Option for IRT model (1pl/2pl or (G)PCM) passed to function
    #with irt_model object and option fo genotype by environment
    #interaction by ge = FALSE/TRUE
    
    #I. With covariates: 
    if(covariates == TRUE && decomp_model == "ADE"){
        source("irt_ade_cov.R")
        output = irt_ade_cov(data_mz = data_mz, data_dz = data_dz, 
                             n_burnin = n_burnin, n_iter = n_iter, 
                             ge = ge, irt_model = irt_model, 
                             N_cov = N_cov, var_prior = var_prior)
    }
    
    #II. Without covariates: 
    if(covariates == FALSE && decomp_model == "ADE"){
        source("irt_ade.R")
        output = irt_ade(data_mz = data_mz, data_dz = data_dz, 
                         n_burnin = n_burnin, n_iter = n_iter,
                         ge = ge, irt_model = irt_model)              
    } 

    
    #==========================================================
    # Output
    #==========================================================
    source("plot.samples.R") #load method for generic function (hoeft later niet meer)
    
    #Remind user of convergence-issue
    cat("\n Here give some advice over convergence")
    cat("\n Here explain how to get results etc")
    
    #Print results on the fly: 
    cat("\n") 
    print(output$results) #table that is made in subroutine (always give same name!)
    return(output) 
}

library(rjags)

#Run function: 
#Simulate data: 
source("simulate_twin_data.R")
data <- simulate_twin_data(n_mz = 10, n_dz = 20, var_a = 0.5, var_c = 0.3,  var_e = 0.2,
                           model = "AE", n_items = 8)

data$y_mz
outtwin = twinUniv(data_mz = data$y_mz, data_dz = data$y_dz, 
                   twin1_datacols_p = c(1:8), twin2_datacols_p = c(9:16),
                   decomp_model = "AE",
                   irt_model = "1PL",
                   cov_imputation = TRUE,  ##### KOMT LATER
                   ge = FALSE,
                   n_iter = 10000, n_burnin = 5000,
                   var_prior = "INV_GAMMA")

names(outtwin)
head(outtwin$samples_var_a)
head(outtwin$samples_var_e)
outtwin$results

#IRT + GE
#twin_analysis2 = twin_analysis(data_mz = data_mz, data_dz = data_dz, 
#                               twin1_datacols = 1:8, twin2_datacols = 9:16, 
#                               ordinal = T, model = "ACE", common = T,
#                               n_iter = 100, n_burnin = 100, ge = TRUE)
#twin_analysis2
#plot(twin_analysis2$samples_var_a, type = "Sampling plot")