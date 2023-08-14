
require(tmvmixnorm)
require(Matrix)

lm_post_wo_equality = function(
    X, y, R, b, 
    beta_init = NULL, sigma_sq_init = 1, 
    mu_0 = NULL, Sigma_0 = NULL, c_0 = 3, 
    nu_0 = 0.02, s_0 = 0.02, empirical_bayes = FALSE, 
    n_sample = 1000, thin = 1, burn = 1000, beta_burn = 10
) {
    # X : [n x p] design matrix
    # y : [n]     response vector
    # R : [m x p] constraint matrix
    # b : [m]     constraint vector

    ##### Linear inequality constraints #####
    # R beta >= b

    ##### Bayesian model #####
    # model : y | X, beta, sigma_sq ~ N_n (X beta, sigma_sq I_n)
    # prior :                  beta ~ TN_p (mu_0, Sigma_0, R, b)
    #       :              sigma_sq ~ Inv-Gamma(nu_0 / 2, s_0 / 2)

    # Here TN_p is p-variate truncated normal distribution. 

    ##### Posterior inference #####
    # Gibbs sampler with conditional posterior
    #     beta | X, y, sigma_sq ~ TN_p (mu_tilde, Signa_tilde, R, b)
    # sigma_sq | X, y, beta     ~ Inv_Gamma(nu_n / 2, s_0 / 2), where

    # Sigma_tilde = (Sigma_0^{-1} + X^T X / sigma_sq)^{-1}, 
    # mu_tilde    = Sigma_tilde (Sigma_0^{-1} mu_0 + X^T y / sigma_sq), 
    # nu_n        = nu_0 + n, 
    # s_n         = s_0 + ||y - X beta||^2. 

    n = nrow(X)
    p = ncol(X)
    m = nrow(R)

    N = burn + n_sample * thin

    # initialization
    if (empirical_bayes) {
        freq_model = lm(y ~ X - 1)
        
        # beta_init = coef(temp_res)
        sigma_sq_init = var(freq_model$residuals)

        mu_0 = coef(freq_model)
        Sigma_0 = vcov(freq_model)
    }  else {
        if (is.null(mu_0)) mu_0 = rep(0, p)
        if (is.null(Sigma_0)) Sigma_0 = diag(c_0^2, p, p)
    }

    beta_samples = matrix(0, ncol = n_sample, nrow = p)
    sigma_sq_samples = rep(0, n_sample)

    cholesky_Sigma_0 = Matrix::chol(Sigma_0)
    inv_Sigma_0 = Matrix::chol2inv(cholesky_Sigma_0)
    xtx = t(X) %*% X
    xty = t(X) %*% y

    nu_n = nu_0 + n

    beta_post = beta_init
    sigma_sq_post = sigma_sq_init

    # Posterior inference ####
    for (i in 1:N) {
        # Gibbs sampler for beta | X, y, sigma_sq
        tilde_Sigma = Matrix::chol2inv(Matrix::chol(
            inv_Sigma_0 + xtx / sigma_sq_post
        ))
        tilde_mu = tilde_Sigma %*% (inv_Sigma_0 %*% mu_0 + xty / sigma_sq_post)

        beta_post = tmvmixnorm::rtmvn(
            n=1, Mean=tilde_mu, Sigma=tilde_Sigma, 
            D=R, lower = b, upper = rep(Inf, m), 
            int=beta_init, burn=beta_burn
        )

        # Gibbs sampler for sigma_sq | X, y, beta
        s_n = s_0 + sum((y - X %*% beta_post)^2)

        sigma_sq_post = 1 / rgamma(1, nu_n / 2, s_n / 2)

        # burning & thinning
        if(i > burn && (i - burn) %% thin ==0) {
            ind = (i - burn) / thin
            beta_samples[, ind] = beta_post
            sigma_sq_samples[ind] = sigma_sq_post
        }   
    }

    return(list(
        beta_samples = beta_samples, 
        sigma_sq_samples = sigma_sq_samples
    ))
}


poisson_log_post_wo_equality = function(
    X, y, R, b, 
    beta_init = NULL, 
    mu_0 = NULL, Sigma_0 = NULL, c_0 = 3, empirical_bayes = FALSE, 
    n_sample = 1000, thin = 1, burn = 1000, beta_burn = 10
) {
    # X   : [n x p] design matrix
    # y   : [n]     response vector
    # R   : [m x p] constraint matrix
    # b   : [m]     constraint vector

    # x_i : [p]     n-th covariate vector 
    # y_i :         n-th responsse scalar 

    # X = [x_1^T, ..., x_n^T]^T
    # y = [y_1,   ..., y_n  ]^T

    ##### Linear inequality constraints #####
    # R beta >= b

    ##### Bayesian model #####
    # model : y_i | x_i, beta ~ Poisson ( exp (x_i^T beta))
    # prior :            beta ~ TN_p (mu_0, Sigma_0, R, b)

    # Here TN_p is p-variate truncated normal distribution. 

    ##### Posterior inference #####
    # Gibbs sampler with conditional posterior by introducing auxiliary varibles u = (u_1, ... u_n)^T
    #         u_i | x_i, beta ~ Uniform (0, exp (-exp (x_i^T beta)))
    #        beta | X, y, u   ~ TN_p (tilde_mu, Sigma_0, tilde_R, tilde_b), where
    # tilde_mu = mu_0 + Sigma_0 X^T y    : [p]         vector 
    # tilde_R  = (R^T, -X^T)^T           : [(m+n) x p] matrix
    # tilde_b  = (b^T, log(-log(u))^T)^T : [(m+n)]     vector

    n = nrow(X)
    p = ncol(X)
    m = nrow(R)

    N = burn + n_sample * thin

    # initialization
    if (empirical_bayes) {
        freq_model = glm(y ~ X - 1, family = poisson(link = "log"))
        
        mu_0 = coef(freq_model)
        Sigma_0 = vcov(freq_model)
    }  else {
        if (is.null(mu_0)) mu_0 = rep(0, p)
        if (is.null(Sigma_0)) Sigma_0 = diag(c_0^2, p, p)
    }

    beta_samples = matrix(0, ncol = n_sample, nrow = p)

    beta_post = beta_init

    tilde_mu = mu_0 + as.vector(Sigma_0 %*% t(X) %*% y)

    # Posterior inference ####
    for (i in 1:N) {
        # Auxilary variable u_i | x_i, y_i, beta
        c = as.vector(sapply(1:n, function(i) {
            y_i = y[i]
            if (y_i == 0) return(0)
            
            x_i = X[i, ]
            lambda_i = sum(x_i * beta)
            u_i = runif(1) * lambda_i^y_i
            return(u_i^(1/y_i))
        }))

        # Gibbs sampler for beta | X, y, u
        tilde_R = rbind(R, X)
        tilde_b = c(b, c)

        beta_post = tmvmixnorm::rtmvn(
            n=1, Mean=tilde_mu, Sigma=Sigma_0, 
            D=tilde_R, lower = tilde_b, upper = rep(Inf, n+m), 
            int=beta_init, burn=beta_burn
        )


        # burning & thinning
        if(i > burn && (i - burn) %% thin ==0) {
            ind = (i - burn) / thin
            beta_samples[, ind] = beta_post
        }   
    }

    return(beta_samples)
}


poisson_iden_post_wo_equality = function(
    X, y, R, b, 
    beta_init = NULL, 
    mu_0 = NULL, Sigma_0 = NULL, c_0 = 3, empirical_bayes = FALSE, 
    n_sample = 1000, thin = 1, burn = 1000, beta_burn = 10
) {
    # X   : [n x p] design matrix
    # y   : [n]     response vector
    # R   : [m x p] constraint matrix
    # b   : [m]     constraint vector

    # x_i : [p]     n-th covariate vector 
    # y_i :         n-th responsse scalar 

    # X = [x_1^T, ..., x_n^T]^T
    # y = [y_1,   ..., y_n  ]^T

    ##### Linear inequality constraints #####
    # R beta >= b

    ##### Bayesian model #####
    # model : y_i | x_i, beta ~ Poisson (x_i^T beta)
    # prior :            beta ~ TN_p (mu_0, Sigma_0, R, b)

    # Here TN_p is p-variate truncated normal distribution. 

    ##### Posterior inference #####
    # Gibbs sampler with conditional posterior by introducing auxiliary varibles u = (u_1, ... u_n)^T
    #         u_i | x_i, y_i, beta ~ Uniform (0, (x_i^T beta)^y_i)
    #        beta | X, y, u        ~ TN_p    (tilde_mu, Sigma_0, tilde_R, tilde_b),   where
    
    # tilde_mu = mu_0 - Sigma_0 X^T 1_n  : [p]         vector 
    # tilde_R  = (R^T, X^T)^T            : [(m+n) x p] matrix
    # tilde_b  = (b^T, c^T)^T            : [(m+n)]     vector, where
    
    # c        = (c_1, ... c_n)^T        : [n]         vector is defined by
    # c_i      = u_i^{1/y_i}  if y_i != 0, 
    #            0            otherwise. 

    n = nrow(X)
    p = ncol(X)
    m = nrow(R)

    N = burn + n_sample * thin

    # initialization
    if (empirical_bayes) {
        freq_model = glm(y ~ X - 1, family = poisson(link = "identity"))
        
        mu_0 = coef(freq_model)
        Sigma_0 = vcov(freq_model)
    }  else {
        if (is.null(mu_0)) mu_0 = rep(0, p)
        if (is.null(Sigma_0)) Sigma_0 = diag(c_0^2, p, p)
    }

    beta_samples = matrix(0, ncol = n_sample, nrow = p)

    beta_post = beta_init
    u = rep(.5, n)

    tilde_mu = mu_0 - as.vector(Sigma_0 %*% t(X) %*% rep(1,n))

    # Posterior inference ####
    for (i in 1:N) {
        # Auxilary variable u_i | x_i, y, beta
        u = runif(n) * as.vector(exp(-exp(X %*% beta_post)))


        # Gibbs sampler for beta | X, y, u
        tilde_R = rbind(R, -X)
        tilde_b = c(b, log(-log(u)))

        beta_post = tmvmixnorm::rtmvn(
            n=1, Mean=tilde_mu, Sigma=Sigma_0, 
            D=tilde_R, lower = tilde_b, upper = rep(Inf, n+m), 
            int=beta_init, burn=beta_burn
        )


        # burning & thinning
        if(i > burn && (i - burn) %% thin ==0) {
            ind = (i - burn) / thin
            beta_samples[, ind] = beta_post
        }   
    }

    return(beta_samples)
}

remove_equality_constraints = function(
    X, R, b, q, l
) {
    # n   : number of observations
    
    # p   : number of parameters
    # q   : number of constrained parameters (q <= p)
    
    # m   : number of total constraints      (m = l + k)
    # l   : number of equality constraints
    # k   : number of inequality constraints

    n = nrow(X)
    p = ncol(X)
    m = nrow(R)
    k = m-l
    
    b_s = b_constraint[1:l]                             # b_s   : l         vector
    g_s = b_constraint[(l+1):m]                         # g_s   : k         vector
    
    R_0 = R[,1:q]                                       # R     : m x q     matrix
    
    R_11 = R[1:l, 1:l]                                  # R_11  : l x l     matrix
    R_12 = R[1:l, (l+1):q]                              # R_12  : l x (q-l) matrix
    R_21 = R[(l+1):m, 1:l]                              # R_21  : k x l     matrix
    R_22 = R[(l+1):m, (l+1):q]                          # R_22  : k x (q-l) matrix
    
    inv_R_11 = MASS::ginv(R_11)                         
    
    D_0 = R_22 - R_21 %*% inv_R_11 %*% R_12             # D_0   : k x (q-l) matrix
    w = as.vector(g_s - R_21 %*% inv_R_11 %*% b_s)      # w     : k         vector
    
    D = cbind(
        matrix(0, nrow = k, ncol = p-q), 
        D_0
    )                                                   # D     : k x (p-l) matrix 
    
    X_11 = X[, 1:l]                                     # X_11  : n x l     matrix
    X_12 = X[, (l+1):q]                                 # X_12  : n x (q-l) matrix
    X_2 = X[, (q+1):p]                                  # X_2   : n x (p-q) matrix
    
    U = cbind(
        X_2, 
        X_12 - X_11 %*% inv_R_11 %*% R_12
    )                                                   # U     : n x (p-l) matrix
    
    delta = X_11 %*% inv_R_11 %*% b_s                   # delta : (p-l)     vector
    
    A = 
    
    nu = c(
        as.vector(inv_R_11 %*% b_s), 
        rep(0, p-l)
    )                                                   # nu    : p         vector
    
    list(
        U = U, 
        delta = delta
    )
}