
require(tmvmixnorm)
require(Matrix)
require(MASS)

source("post_inference_wo_equality.R")

remove_equality_condition = function(X, R, b, l) {

    # X   : [n x p] design matrix
    # R   : [m x p] constraint matrix, including both equality and inequality constraints
    # b   : [m]     constraint vector

    # n          : number of observations    
    # p          : number of parameters
    
    # m = l + k  : number of total constraints      
    # l          : number of equality constraints
    # k          : number of inequality constraints

    n = nrow(X)
    p = ncol(X)
    m = nrow(R)
    k = m-l

    if (l == 0) 
        return(list(
            new_R = R, 
            new_b = b, 
            new_X = X, 
            alpha = rep(0, n), 
            recover_mat = diag(1, p, p), 
            recover_bias = rep(0, p)
        ))
    
    X_1 = X[, 1:l]                                      # X_1    : [n x l]     matrix
    X_2 = X[, (l+1):p]                                  # X_2    : [n x (p-l)] matrix
     
    R_11 = R[1:l, 1:l]                                  # R_11   : [l x l]     matrix
    R_12 = R[1:l, (l+1):p]                              # R_12   : [l x (p-l)] matrix
    R_21 = R[(l+1):m, 1:l]                              # R_21   : [k x l]     matrix
    R_22 = R[(l+1):m, (l+1):p]                          # R_22   : [k x (p-l)] matrix

    b_1 = b[1:l]                                        # b_1    : [l]         vector  
    b_2 = b[(l+1):m]                                    # b_2    : [k]         vector
    
    inv_R_11 = MASS::ginv(R_11)                         
    
    # Constraints
    new_R = R_22 - R_21 %*% inv_R_11 %*% R_12           # new_R  : [k x (p-l)] matrix
    new_b = as.vector(b_2 - R_21 %*% inv_R_11 %*% b_1)  # new_b  : [k]         vector
    
    # Design matrix
    new_X = X_2 - X_1 %*% inv_R_11 %*% R_12             # new_X  : n x (p-l) matrix
    alpha = X_1 %*% inv_R_11 %*% b_1                    # alpha  : n         vector
    
    # Recovering matrix
    recover_mat = rbind(
        -inv_R_11 %*% R_12, 
        diag(1, p-l, p-l)
    )
    recover_bias = c(
        as.vector(inv_R_11 %*% b_1), 
        rep(0, p-l)
    )

    return(list(
        new_R = new_R, 
        new_b = new_b, 
        new_X = new_X, 
        alpha = alpha, 
        recover_mat = recover_mat, 
        recover_bias = recover_bias
    ))
}

lm_post_with_equality = function(
    X, y, R, b, l=0, 
    beta_init = NULL, sigma_sq_init = 1, 
    mu_0 = NULL, Sigma_0 = NULL, c_0 = 3, 
    nu_0 = 0.02, s_0 = 0.02, empirical_bayes = FALSE, 
    n_sample = 1000, thin = 1, burn = 1000, beta_burn = 10
) {
    n = nrow(X)
    p = ncol(X)
    m = nrow(R)

    if (l == 0) {
        return(lm_post_wo_equality(
            X = X, y = y, R = R, b = b, alpha = NULL, 
            beta_init = beta_init, sigma_sq_init = sigma_sq_init, 
            mu_0 = mu_0, Sigma_0 = Sigma_0, c_0 = c_0, 
            nu_0 = nu_0, s_0 = s_0, empirical_bayes = empirical_bayes,
            n_sample = n_sample, thin = thin, burn = burn, beta_burn = beta_burn
        ))
    }

    remove_list = remove_equality_condition(X, R, b, l)

    new_beta_init = NULL
    new_mu_0 = NULL
    new_Sigma_0 = NULL
    if (!is.null(beta_init)) new_beta_init = beta_init[(l+1):p]
    if (!is.null(mu_0)) new_mu_0 = mu_0[(l+1):p]
    if (!is.null(Sigma_0)) new_Sigma_0 = Sigma_0[(l+1):p, (l+1):p]

    transformed_samples = lm_post_wo_equality(
        X = remove_list$new_X, y = y, R = remove_list$new_R, b = remove_list$new_b, alpha = remove_list$alpha, 
        beta_init = new_beta_init, sigma_sq_init = sigma_sq_init, 
        mu_0 = new_mu_0, Sigma_0 = new_Sigma_0, c_0 = c_0, 
        nu_0 = nu_0, s_0 = s_0, empirical_bayes = empirical_bayes, 
        n_sample = n_sample, thin = thin, burn = burn, beta_burn = beta_burn
    )

    recovered_beta_samples = remove_list$recover_mat %*% transformed_samples$beta_samples + remove_list$recover_bias

    return(list(
        beta_samples = recovered_beta_samples, 
        sigma_sq_samples = transformed_samples$sigma_sq_samples
    ))
}

poisson_log_post_with_equality = function(
    X, y, R, b, l, 
    beta_init = NULL, 
    mu_0 = NULL, Sigma_0 = NULL, c_0 = 3, empirical_bayes = FALSE, 
    n_sample = 1000, thin = 1, burn = 1000, beta_burn = 10
) {
    n = nrow(X)
    p = ncol(X)
    m = nrow(R)

    if (l == 0) {
        return(poisson_log_post_wo_equality(
            X = X, y = y, R = R, b = b, alpha = NULL, 
            beta_init = beta_init, 
            mu_0 = mu_0, Sigma_0 = Sigma_0, c_0 = c_0, empirical_bayes = empirical_bayes, 
            n_sample = n_sample, thin = thin, burn = burn, beta_burn = 10
        ))
    }

    remove_list = remove_equality_condition(X, R, b, l)

    new_beta_init = NULL
    new_mu_0 = NULL
    new_Sigma_0 = NULL
    if (!is.null(beta_init)) new_beta_init = beta_init[(l+1):p]
    if (!is.null(mu_0)) new_mu_0 = mu_0[(l+1):p]
    if (!is.null(Sigma_0)) new_Sigma_0 = Sigma_0[(l+1):p, (l+1):p]

    transformed_samples = poisson_log_post_wo_equality(
        X = remove_list$new_X, y = y, R = remove_list$new_R, b = remove_list$new_b, alpha = remove_list$alpha, 
        beta_init = new_beta_init, 
        mu_0 = new_mu_0, Sigma_0 = new_Sigma_0, c_0 = c_0, empirical_bayes = empirical_bayes, 
        n_sample = n_sample, thin = thin, burn = burn, beta_burn = beta_burn
    )

    recovered_beta_samples = remove_list$recover_mat %*% transformed_samples + remove_list$recover_bias

    return(recovered_beta_samples)
}
