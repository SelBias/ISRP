
require(tmvmixnorm)
require(Matrix)
require(MASS)

lm_post_wo_equality = function(
    X, y, R, b, alpha = NULL,
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
    # model : y | X, beta, sigma_sq ~ N_n       (alpha + X beta, sigma_sq I_n)
    # prior :                  beta ~ TN_p      (mu_0, Sigma_0, R, b)
    #       :              sigma_sq ~ Inv-Gamma (nu_0 / 2, s_0 / 2)

    # Here TN_p is p-variate truncated normal distribution. 

    ##### Posterior inference #####
    # Gibbs sampler with conditional posterior
    #     beta | X, y, sigma_sq ~ TN_p (mu_tilde, Signa_tilde, R, b)
    # sigma_sq | X, y, beta     ~ Inv_Gamma(nu_n / 2, s_0 / 2), where

    # Sigma_tilde = (Sigma_0^{-1} + X^T X / sigma_sq)^{-1}, 
    # mu_tilde    = Sigma_tilde (Sigma_0^{-1} mu_0 + X^T (y - alpha) / sigma_sq), 
    # nu_n        = nu_0 + n, 
    # s_n         = s_0 + ||y - alpha - X beta||^2. 

    n = nrow(X)
    p = ncol(X)
    m = nrow(R)

    N = burn + n_sample * thin


    # initialization
    if (is.null(alpha)) alpha = rep(0, n)

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
    xty = t(X) %*% (y - alpha)

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
        s_n = s_0 + sum((y - alpha - X %*% beta_post)^2)

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
    X, y, R, b, alpha = NULL, 
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
    # model : y_i | x_i, beta ~ Poisson (exp(alpah + x_i^T beta))
    # prior :            beta ~ TN_p    (mu_0, Sigma_0, R, b)

    # Here TN_p is p-variate truncated normal distribution. 

    ##### Posterior inference #####
    # Gibbs sampler with conditional posterior by introducing auxiliary varibles u = (u_1, ... u_n)^T
    #         u_i | x_i, beta ~ Uniform (0, exp(-exp(alpha + x_i^T beta)))
    #        beta | X, y, u   ~ TN_p    (tilde_mu, Sigma_0, tilde_R, tilde_b), where
    # tilde_mu = mu_0 + Sigma_0 X^T y              : [p]         vector 
    # tilde_R  = (R^T, -X^T)^T                     : [(m+n) x p] matrix
    # tilde_b  = (b^T, alpha^T - log(-log(u))^T)^T : [(m+n)]     vector

    n = nrow(X)
    p = ncol(X)
    m = nrow(R)

    N = burn + n_sample * thin


    # initialization
    if (is.null(alpha)) alpha = rep(0, n)

    if (empirical_bayes) {
        freq_model = glm(y ~ X - 1, alpha = alpha, family = poisson(link = "log"))
        
        mu_0 = coef(freq_model)
        Sigma_0 = vcov(freq_model)
    }  else {
        if (is.null(mu_0)) mu_0 = rep(0, p)
        if (is.null(Sigma_0)) Sigma_0 = diag(c_0^2, p, p)
    }

    beta_samples = matrix(0, ncol = n_sample, nrow = p)

    beta_post = beta_init

    tilde_mu = mu_0 + as.vector(Sigma_0 %*% t(X) %*% y)
    tilde_R = rbind(R, -X)

    # Posterior inference ####
    for (i in 1:N) {
        # Auxilary variable u_i | x_i, y, beta
        u = runif(n) * as.vector(exp(-exp(alpha + X %*% beta_post)))

        # Gibbs sampler for beta | X, y, u
        tilde_b = c(b, alpha - log(-log(u)))

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
    X, y, R, b, alpha = NULL, 
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
    # model : y_i | x_i, beta ~ Poisson (alpha + x_i^T beta)
    # prior :            beta ~ TN_p    (mu_0, Sigma_0, R, b)

    # Here TN_p is p-variate truncated normal distribution. 

    ##### Posterior inference #####
    # Gibbs sampler with conditional posterior by introducing auxiliary varibles u = (u_1, ... u_n)^T
    #         u_i | x_i, y_i, beta ~ Uniform (0, (alpha + x_i^T beta)^y_i)
    #        beta | X, y, u        ~ TN_p    (tilde_mu, Sigma_0, tilde_R, tilde_b),   where
    
    # tilde_mu = mu_0 - Sigma_0 X^T 1_n  : [p]         vector 
    # tilde_R  = (R^T, X^T)^T            : [(m+n) x p] matrix
    # tilde_b  = (b^T, -a^T + c^T)^T            : [(m+n)]     vector, where
    
    # c        = (c_1, ... c_n)^T        : [n]         vector is defined by
    # c_i      = u_i^{1/y_i}  if y_i != 0, 
    #            0            otherwise. 

    n = nrow(X)
    p = ncol(X)
    m = nrow(R)

    N = burn + n_sample * thin

    # initialization
    if (is.null(alpha)) alpha = rep(0, n)

    if (empirical_bayes) {
        freq_model = glm(y ~ X - 1, alpha = alpha, family = poisson(link = "identity"))
        
        mu_0 = coef(freq_model)
        Sigma_0 = vcov(freq_model)
    }  else {
        if (is.null(mu_0)) mu_0 = rep(0, p)
        if (is.null(Sigma_0)) Sigma_0 = diag(c_0^2, p, p)
    }

    beta_samples = matrix(0, ncol = n_sample, nrow = p)
    beta_post = beta_init

    tilde_mu = mu_0 - as.vector(Sigma_0 %*% t(X) %*% rep(1,n))
    tilde_R = rbind(R, X)

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
        tilde_b = c(b, -alpha + c)

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