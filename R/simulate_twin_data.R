#==========================================================
# Simulation of twin data under the AE, ADE or ACE model, 
# using sum  scores or item scores, optional with GxE 
# and/or environmental covariates
#==========================================================

#setwd("/Users/Inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
simulate_twin_data <- function(n_mz, n_dz, var_a = 0.5, var_c = 0.3,  var_e = 0.2, var_d = 0, 
                               model = "ACE", n_items = 20, n_var = 0, ge = FALSE, ge_beta0 = 
                               log(0.5), ge_beta1 = 1.5){
    
    #==========================================================
    # Error messages
    #==========================================================    
    if(n_items < 0){
        cat("Error: \n
             The number of test items is negative. Something went wrong here!")
    }
    
    
    if(n_items == 0){
        cat("Error: \n
             The number of test items used for the simulation is equal to 0. Without any items, the
             twin data cannot be simulated!")
    }
    
    if (ge == TRUE && ge_beta0 == 0){
        cat("Error: \n
                Please specify a value for beta0 when simulating data with genotype by environment interaction.
                The parameter beta0 is defined as average environmental variance (i.e., when A = 0) 
                For more information see Schwabe & van den Berg (2014), Behavior Genetics, 44 (4), 394-406.")
    }
    
    if (var_d > 0 && var_c > 0){
        cat("Warning: \n 
                You cannot specify variance due to shared-environmental effects (C) and variance due to 
                dominance effects (D) at the same time. This model is not identified! \n
                Depending on your model-choice (ACE/ADE),
                either var_c or var_d is used for the data simulation.")
    }
    
    if(n_dz == 0){
        cat("Warning: \n
                Unless you want to simulate data under the AE model, please chose a total number of DZ twins above 0.")
    }
    
    if(n_mz == 0){
        cat("Error: \n
                You specified a total number of MZ twins of 0. This won't work! Please choose a number above 0.")
    }
    #==========================================================
    
    #==========================================================
    # Data simulation ACE,ADE,AE:
    #==========================================================
    if(model == "ACE"){
        print("Simulating data under an ACE model...")
        source("simulate_ACE.R")
        output = simulate_ACE(n_mz = n_mz, n_dz = n_dz, var_a = var_a, var_c = var_c, 
                              var_e = var_e, n_items = n_items, n_var = n_var, ge = ge,
                              ge_beta0 = ge_beta0, ge_beta1 = ge_beta1)
    } else if (model == "ADE"){
        print("Simulating data under an ADE model...")
        source("simulate_ADE.R")
        output = simulate_ADE(n_mz = n_mz, n_dz = n_dz, var_a = var_a, var_d = var_d, 
                              var_e = var_e, n_items = n_items, n_var = n_var, ge = ge,
                              ge_beta0 = ge_beta0, ge_beta1 = ge_beta1)
    } else {
        print("Simulating data under an AE model...")
        source("simulate_AE.R")
        output = simulate_AE(n_mz = n_mz, n_dz = n_dz, var_a = var_a, var_e = var_e, 
                             n_items = n_items, n_var = n_var, ge = ge,
                             ge_beta0 = ge_beta0, ge_beta1 = ge_beta1)
    }
    #==========================================================
    return(output)    
}