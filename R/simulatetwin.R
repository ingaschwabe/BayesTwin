#==========================================================
# simulatetwin.R
# Simulation of twin data under the AE, ADE or ACE model,
# using sum  scores or item scores, optional with GxE 
# and/or environmental covariates
# BayesTwin package
#==========================================================

simulatetwin <- function(n_mz = 140, n_dz = 360, var_a = 0.5, var_c = 0.3,  var_e = 0.2, var_d = 0, 
                         model = "ACE", n_items = 20, n_cat = 0, 
                         ge = FALSE, ge_beta0 = 0, ge_beta1 = 0, irt_model = "1PL"){
    
    #==========================================================
    # Error messages
    #==========================================================  
    
    #Errors regarding IRT models/item data
    if (irt_model == "PCM" && n_cat < 3){
        stop("You have chosen to simulate data under the partial credit model but the number of categories is less
             than 3. To simulate data under the partial credit model, you need at least three categories (n_cat = 3)!")
    }
    
    if (irt_model == "GPCM" && n_cat< 3){
        stop("You have chosen to simulate data under the generalized partial credit model but the number of categories is less
             than 3. To simulate data under the genearalized partial credit model, you need at least three categories (n_cat = 3)!")
    }
    
    if (irt_model == "PCM" && n_cat == 0){
        stop("You have chosen to simulate data under the partial credit model but the number of categories is equal
             to 0. To simulate data under the partial credit model, you need at least three categories (n_cat = 3)!")
    }
    
    if (irt_model == "GPCM" && n_cat == 0){
        stop("You have chosen to simulate data under the generalized partial credit model but the number of categories is equal
             to 0. To simulate data under the genearalized partial credit model, you need at least three categories (n_cat = 3)!")
    }
    
    
    ##Errors regarding variance components and/or GxE
    if (var_e > 0 && ge_beta0 > 0){
        cat("Warning! \n
             Both VAR(E) and beta0 are specified. If a model without GxE is choosen, then the value given for var_e
             will be used for data simulation. For a model with GxE, the beta0 parameter will be used.")
    }
    
    if (var_d > 0 && var_c > 0){
        cat("Warning! \n
             You cannot specify variance due to shared-environmental effects (C) and variance due to 
             dominance effects (D) at the same time. This model is not identified!
             Depending on your model-choice (ACE/ADE), either var_c or var_d is used for the data simulation.")
    }
    
    
    if (ge == TRUE && ge_beta0 == 0){
        stop("Please specify a value for beta0 when simulating data with genotype by environment interaction.
              The parameter beta0 is defined as average environmental variance (i.e., when A = 0) 
              For more information see Schwabe & van den Berg (2014), Behavior Genetics, 44 (4), 394-406.")
    }
    
    if(model == "AE" && ge == FALSE && (round(var_a + var_e,1) != 1)){
        stop("Total phenotypic variance must be equal to 1! Please respecify variance components!")
    }
    
    if(model == "AE" && ge == TRUE && (round(var_a + exp(ge_beta0),1) != 1)){
        stop("Total phenotypic variance must be equal to 1! Please respecify variance components!")
    }

    if(model == "ACE" && ge == FALSE && (round(var_a + var_c + var_e,1) != 1)){
        stop("Total phenotypic variance must be equal to 1! Please respecify variance components!")
    }
    
    if(model == "ACE" && ge == TRUE && (round(var_a + var_c + exp(ge_beta0),1) != 1)){
        stop("Total phenotypic variance must be equal to 1! Please respecify variance components!")
    }
    
    if(model == "ADE" && ge == FALSE && (round(var_a + var_d + var_e,1) != 1)){
        stop("Total phenotypic variance must be equal to 1! Please respecify variance components!")
    }
    
    if(model == "ADE" && ge == TRUE && (round(var_a + var_d + exp(ge_beta0),1) != 1)){
        stop("Total phenotypic variance must be equal to 1! Please respecify variance components!")
    }
    
    if(ge == TRUE && var_e > 0){
        cat("Warning! \n
             A value has been chosen for argument value_e although the model under which the data will be 
             simulated includes GxE. Under this model, the parameter ge_beta0 will be used which represents 
             average environmental variance (see the help function for more details).")
    }
    
    
    #Errors regarding total twins and/or items: 
    if(n_items < 0){
        stop("The number of test items is negative. Something went wrong here!
             Use a number > 0 for the parameter 'n_items'.")
    }
    
    
    if(n_items == 0){
        stop("The number of test items used for the simulation is equal to 0. Without any items, the
             twin data cannot be simulated! Use a number > 0 for the parameter 'n_items'.")
    }
    
    if(n_dz == 0){
        cat("Warning: \n
                Unless you want to simulate data under the AE model, please chose a total number of DZ twins above 0.")
    }
    
    if(n_mz == 0){
        stop("You specified a total number of MZ twins of 0. This won't work! Please choose a number above 0.")
    }
    #==========================================================
    
    #==========================================================
    # Data simulation ACE,ADE,AE:
    #==========================================================
    if(model == "ACE"){
        print("Simulating data under an ACE model...")
        output = simulate_ACE(n_mz = n_mz, n_dz = n_dz, var_a = var_a, var_c = var_c, 
                              var_e = var_e, n_items = n_items, n_cat = n_cat,
                              ge = ge, ge_beta0 = ge_beta0, ge_beta1 = ge_beta1,
                              irt_model = irt_model)
    } else if (model == "ADE"){
        print("Simulating data under an ADE model...")
        output = simulate_ADE(n_mz = n_mz, n_dz = n_dz, var_a = var_a, var_d = var_d, 
                              var_e = var_e, n_items = n_items, n_cat = n_cat,
                              ge = ge,ge_beta0 = ge_beta0, ge_beta1 = ge_beta1,
                              irt_model = irt_model)
    } else {
        print("Simulating data under an AE model...")
        output = simulate_AE(n_mz = n_mz, n_dz = n_dz, var_a = var_a, var_e = var_e, 
                             n_items = n_items, n_cat = n_cat,
                             ge = ge, ge_beta0 = ge_beta0, ge_beta1 = ge_beta1,
                             irt_model = irt_model)
    }
    #==========================================================
    return(output)    
}