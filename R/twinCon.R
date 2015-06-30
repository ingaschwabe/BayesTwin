#==========================================================
# Master function that does a concordance rate analysis
# BayesTwin package - Bayesian Analysis of Twin Data 
#==========================================================
#setwd("C:/Users/schwabei/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("/Users/stephanievandenberg/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("C:/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")

twinCon = function(data_mz, data_dz, 
                   twin1_datacols_p, twin2_datacols_p,
                   n_iter = 10000, n_burnin = 5000){                                                 
    
    #Calculate number of MZ and DZ twins: 
    n_mz = nrow(data_mz)
    n_dz = nrow(data_dz)
    
    #Make sure that first columns are only twin 1, next the twin 2 columns
    data_mz = data_mz[,c(twin1_datacols_p, twin2_datacols_p)]
    data_dz = data_dz[,c(twin1_datacols_p, twin2_datacols_p)]
    
    
}