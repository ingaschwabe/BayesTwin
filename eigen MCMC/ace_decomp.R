mcmc <- function(y_mz, y_dz, n_iterations, burnin, skip){
    
    #Install packages when necessary: 
    #if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
    #if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
    
    #Prior values variance components
    alpha = 1
    beta = 1
    
    #Variables
    n_mz = nrow(y_mz) #Number of MZ twin pairs
    n_dz = nrow(y_dz) #Number of DZ twin pairs
    n = n_mz + n_dz
    pb = txtProgressBar(min = 0, max = n_iterations, style = 3, width = 50)
    n_save = floor((n_iterations-burnin)/skip) #Number of samples to be stored
    
    #STORE-variables:
    store_a_mz = matrix(NA, n_mz, n_save)
    store_c_mz = matrix(NA, n_mz, n_save)
    store_a1_dz = matrix(NA, n_dz, n_save)
    store_a2_dz = array(NA, c(n_dz, 2, n_save))
    store_c_dz = matrix(NA, n_dz, n_save)
    store_var_a = rep(NA, n_save)
    store_var_c = rep(NA, n_save)
    store_var_e = rep(NA, n_save)
    
    #Define start values:
    var_a = runif(1,0.1, 0.9)
    var_c = runif(1,0.1, 0.9) 
    var_e = runif(1,0.1, 0.9)
    a_mz = rnorm(n_mz, 0, sqrt(var_a))
    c_mz = rnorm(n_mz, 0, sqrt(var_c))
    a1_dz = rnorm(n_dz, 0, 0.5 * sqrt(var_a))
    a2_dz = cbind(rnorm(n_dz, a1_dz, sqrt(0.5 * var_a)),
                  rnorm(n_dz, a1_dz, sqrt(0.5 * var_a)))
    c_dz = rnorm(n_dz, 0, sqrt(var_c))
    
    #For each parameter, store start value(s) at first position of store vector:
    store_a_mz[,1] = a_mz 
    store_c_mz[,1] = c_mz
    store_a1_dz[,1] = a1_dz
    store_a2_dz[,,1] = a2_dz
    store_c_dz[,1] = c_dz
    store_var_a[1] = var_a
    store_var_c[1] = var_c
    store_var_e[1] = var_e
    
    #Begin with sampling procedure 
    tt = 0
    
    for (t in 1:n_iterations){
        
        #Save processing time and create progress bar: 
        ptm = proc.time()
        setTxtProgressBar(pb, t)
        
        #1. Sample C_i and A_i for MZ twins:
        e_c_mz = (1/var_a + 1/var_c)^-1 * (a_mz/var_a) 
        c_mz = rnorm(n_mz, e_c_mz, sqrt( (1/var_c + 1/var_a)^-1) )
        
        e_a_mz = (2/var_e + 1/var_a)^-1 * (    (y_mz[,1] + y_mz[,2])/var_e + c_mz/var_a) 
        a_mz = rnorm(n_mz, e_a_mz, sqrt((1/var_a + 2/var_e)^-1))
        
        #2. Sample C_i, A1_i and A2_ij for DZ twins: 
        e_c_dz = (2/var_a + 1/var_c)^-1 * ( (2*a1_dz)/var_a )
        c_dz = rnorm(n_dz, e_c_dz, sqrt((2/var_a + 1/var_c)^-1))        
        
        a1_dz = rnorm(n_dz, 1/3 * ( 2 * apply(a2_dz, 1, mean) + c_dz), sqrt(var_a/6))
        
        e_a2_dz_twin1 = (1/var_e + 2/var_a)^-1 * ( y_dz[,1]/var_e + (2*a1_dz)/var_a) 
        e_a2_dz_twin2 = (1/var_e + 2/var_a)^-1 * ( y_dz[,2]/var_e + (2*a1_dz)/var_a) 
        a2_dz = cbind(rnorm(n_dz, e_a2_dz_twin1, sqrt((2/var_a + 1/var_e)^-1)),
                      rnorm(n_dz, e_a2_dz_twin2, sqrt((2/var_a + 1/var_e)^-1)))
        
        #3. Sample variance components: 
        alpha_tilde_c = alpha + n/2
        beta_tilde_c = beta + 0.5 *  sum( (c(c_mz, c_dz)^2) )
        var_c = rinvgamma(1, alpha_tilde_c, beta_tilde_c)
        
        alpha_tilde_e = alpha + n
        beta_tilde_e =  beta + (0.5 *  sum( (y_dz - a2_dz)^2)) + (0.5 *  sum( (y_mz - a_mz)^2))
        var_e = rinvgamma(1, alpha_tilde_e, beta_tilde_e)
        
        alpha_tilde_a = alpha + 1.5*n_dz + 0.5*n_mz
        beta_tilde_a = beta + (sum((a1_dz - a2_dz[,1])^2)) + (sum((a1_dz - a2_dz[,2])^2))  + (0.5 * sum( (a_mz - c_mz)^2)) + 
                sum( (a1_dz - c_dz)^2) 
        var_a = rinvgamma(1, alpha_tilde_a, beta_tilde_a)
        
        #Set rest of parameters back to true parameter: 
        #         var_a = 0.5; 
        #         var_c = 0.3; 
        #         var_e = 0.2
        
        #       c_mz = true_c_mz
        #          c_dz = true_c_dz
        #     a_mz = true_a_mz
        #     a1_dz = true_a1_dz
        #     a2_dz = true_a2_dz
        
        
        #Store parameter estimates: 
        if( t>burnin && (t-burnin)%%skip == 0 ) {
            #print(paste("Iteration Number",t))     
            tt = tt+1
            #store_a_mz[,tt] = a_mz
            #store_c_mz[,tt] = c_mz
            #store_a1_dz[,tt] = a1_dz
            #store_a2_dz[,,tt] = a2_dz
            #store_c_dz[,tt] = c_dz
            store_var_a[tt] = var_a
            store_var_c[tt] = var_c
            store_var_e[tt] = var_e
        }
        
        
    }# n_iterations
    
    return(  list(
        #a_mz = store_a_mz, 
        #c_mz = store_c_mz,
        #a1_dz = store_a1_dz, 
        #a2_dz = store_a2_dz, 
        #c_dz = store_c_dz, 
        var_a = store_var_a,
        var_c = store_var_c,
        var_e = store_var_e
    ) )
}