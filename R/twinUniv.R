#==========================================================
# Master function that calls subroutines for univariate
# analyses with one phenotype. 
# BayesTwin package - Bayesian Analysis of Twin Data 
#==========================================================
#setwd("C:/Users/schwabei/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("/Users/stephanievandenberg/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("C:/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")

twinUniv = function(data_mz, data_dz, 
                    twin1_datacols_p, twin2_datacols_p,
                    twin1_datacols_cov = NA, twin2_datacols_cov = NA,
                    decomp_model = "ACE",
                    irt_model = "1PL",
                    #cov_imputation = TRUE, 
                    ge = FALSE,
                    n_iter = 10000, n_burnin = 5000){
    #Load HPD function: 
    source("HPD.R")                                                         
    
    #Calculate number of MZ and DZ twins: 
    n_mz = nrow(data_mz)
    n_dz = nrow(data_dz)
    
    #Make sure that first columns are only twin 1, next the twin 2 columns
    data_mz = data_mz[,c(twin1_datacols_p, twin2_datacols_p)]
    data_dz = data_dz[,c(twin1_datacols_p, twin2_datacols_p)]
    
    #Test if its item data: 
    itemdata = NA
    if(length(twin1_datacols) == 1){
        itemdata = FALSE
    } else {itemdata = TRUE}
    
    #When data_imputation = FALSE, use only complete cases: 
    if(is.na(twin1_datacols_cov) && data_imputation == FALSE){
        data_mz = na.omit(data_mz)
        data_dz = na.omit(data_dz)
    }
    
    #Test if covariates are used:
    covariates = NA
    
    #Option I: Use covaritaes without missing value imputation: 
    if (is.na(twin1_datacols_cov) == FALSE && data_imputation == TRUE){
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
    } else if (is.na(twin1_datacols_cov) == FALSE && data_imputation == FALSE){
        
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
    # Subroutines for ACE model
    #==========================================================

    #==========================================================
    # Call subroutines for analyses based on sumscores
    #==========================================================
    #### I. With covariates: #####
    #whether GE is estimated will be determined by ge = FALSE/TRUE
    if(itemdata == FALSE && covariates == TRUE && model == "ACE"){
        source("sumscores_cov.R")
        output = sumscores_cov(data_mz = data_mz, data_dz = data_dz,
                               X_mz_twin1 = X_mz_twin1, X_mz_twin2 = X_mz_twin2,
                               X_dz_twin1 = X_dz_twin1, X_dz_twin2 = X_dz_twin2,
                               dich_cov = dich_cov, cont_cov = cont_cov, 
                               n_burnin = n_burnin, n_iter = n_iter,
                               ge = ge, dich = dich) 
    }
    
    #### II. Without covariates: #####
    #whether GE is estimated will be determined by ge = FALSE/TRUE
    if(itemdata == FALSE && covariates == FALSE && model == "ACE"){
        source("sumscores.R")
        output = sumscores(data_mz = data_mz, data_dz = data_dz,
                           ge = ge, 
                           n_burnin = n_burnin, n_iter = n_iter)            
    } 
    
    #==========================================================
    # Call subroutines for item data analyses
    #==========================================================
    #I. With covariates: 
    #Option for IRT model (1pl/2pl or (G)PCM) passed to function
    #with irt_model object and option fo genotype by environment
    #interaction by ge = FALSE/TRUE
    if(itemdata == TRUE && covariates == TRUE && model == "ACE"){
        source("irt_cov.R")
        output = irt(data_mz = data_mz, data_dz = data_dz, 
                     n_burnin = n_burnin, n_iter = n_iter, 
                     ge = ge, irt_model = irt_model)
    }

    #II. Without covariates: 
    if(itemdata == TRUE && covariates == FALSE && model == "ACE"){
        source("irt.R")
        output = irt(data_mz = data_mz, data_dz = data_dz, 
                     n_burnin = n_burnin, n_iter = n_iter,
                     ge = ge, irt_model = irt_model)              
    } 
    

    #==========================================================
    #==========================================================
    # Subroutines for ADE model
    #==========================================================
    
    #==========================================================
    # Call subroutines for analyses based on sumscores
    #==========================================================
    #### I. With covariates: #####
    #whether GE is estimated will be determined by ge = FALSE/TRUE
    #if(itemdata == FALSE && covariates == TRUE && model == "ACE"){
    #source("sumscores_ade_cov.R")
    #output = sumscores_ade_cov(data_mz = data_mz, data_dz = data_dz,
    #                           X_mz_twin1 = X_mz_twin1, X_mz_twin2 = X_mz_twin2,
    #                           X_dz_twin1 = X_dz_twin1, X_dz_twin2 = X_dz_twin2,
    #                           dich_cov = dich_cov, cont_cov = cont_cov, 
    #                           n_burnin = n_burnin, n_iter = n_iter
    #                           ge = ge, dich = dich) 
    # }

    #### II. Without covariates: #####
    #whether GE is estimated will be determined by ge = FALSE/TRUE
    if(itemdata == FALSE && covariates == FALSE && model == "ADE"){
        source("sumscores_ade.R")
        output = sumscores_ade(data_mz = data_mz, data_dz = data_dz,
                               ge = ge, 
                               n_burnin = n_burnin, n_iter = n_iter)            
    } 

    
    #==========================================================
    # Call subroutines for item data analyses
    #==========================================================
    #I. With covariates: 
    #Option for IRT model (1pl/2pl or (G)PCM) passed to function
    #with irt_model object and option fo genotype by environment
    #interaction by ge = FALSE/TRUE

    #II. Without covariates: 
    if(itemdata == TRUE && covariates == FALSE && model == "ADE"){
        source("irt_ADE.R")
        output = irt(data_mz = data_mz, data_dz = data_dz, 
                    n_burnin = n_burnin, n_iter = n_iter,
                    ge = ge, irt_model = irt_model)              
    } 

    
    #==========================================================
    # Output
    #==========================================================
    #source("plot.samples.R") #load method for generic function (hoeft later niet meer)
    
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
#source("simulate_twin_data.R")
#data <- simulate_twin_data(nmz = 10, ndz = 20, var_a = 0.5, var_c = 0.3,  var_e = 0.2,
#                           model = "ACE", n_items = 8, n_var = 3)
#IRT + GE
#twin_analysis2 = twin_analysis(data_mz = data_mz, data_dz = data_dz, 
#                               twin1_datacols = 1:8, twin2_datacols = 9:16, 
#                               ordinal = T, model = "ACE", common = T,
#                               n_iter = 100, n_burnin = 100, ge = TRUE)
#twin_analysis2
#plot(twin_analysis2$samples_var_a, type = "Sampling plot")