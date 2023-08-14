library(rstudioapi)
library(MASS)
library(Matrix)
library(tmvmixnorm)
library(restriktor)

library(ggplot2)
library(reshape2)
# library(ggpubr)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set dir
source("./posterior_inference.R")

#################################
#####     Simulation A1     #####
#################################

repitition = 100
sample_size = 10
set.seed(1)

# Hyperparameter for BICLS
beta_init = c(-0.5, -0.5, 0)
sigma_sq_init = 1
c_0 = 5
nu_0 = 0.02
s_0 = 0.02

# Set beta_true
beta_true = c(-1, -1, 1)
p = length(beta_true)

rho = 0.5
X_mean = rep(0, p)
X_precision = matrix(0, nrow=p, ncol=p)
for (i in 1:p) {
    for (j in 1:p) 
        X_precision[i,j] = rho^(abs(i-j))
}
X_cov = MASS::ginv(X_precision)
rm(list = c("i", "j", "rho"))

# 100 time repeated data generation

X_list = lapply(
    rep(sample_size, repitition), 
    MASS::mvrnorm, X_mean, X_cov
)
eps_list = lapply(
    rep(sample_size, repitition), 
    rnorm, 0, 3
)
y_list = lapply(
    1:repitition, 
    function(i, beta) {
        as.vector(X_list[[i]] %*% beta + eps_list[[i]])
    }, beta_true
)
rm(list = c("X_mean", "X_cov", "X_precision", "eps_list"))

# Linear inequality constraints
R = matrix(
    c( 1,-2,0, 
      -1, 0,0), 
    nrow=2, ncol=3, byrow=TRUE
)
b = rep(0,2)

################################
##### Parameter estimation #####
################################

# OLS estimates
fit_OLS = lapply(1:repitition, function(n) {
    X = X_list[[n]]
    y = y_list[[n]]
    lm(y ~ X - 1)
})

beta_OLS = sapply(fit_OLS, coef)
rownames(beta_OLS) = paste0("beta", rep(1:3))

beta_diff_OLS = beta_OLS - beta_true
(sd_OLS = apply(beta_diff_OLS, 1, sd))

# ICLS estimates
fit_ICLS = lapply(
    fit_OLS, 
    restriktor::restriktor, 
    constraints=R, rhs=b, neq=0
)
beta_ICLS = sapply(fit_ICLS, coef)
rownames(beta_ICLS) = paste0("beta", rep(1:3))

beta_diff_ICLS = beta_ICLS - beta_true
(sd_ICLS = apply(beta_diff_ICLS, 1, sd))

# BICLS estimates

fit_BICLS = lapply(
    1:repitition, function(
        n, 
        R, b, 
        beta_init = NULL, sigma_sq_init = 1, 
        mu_0 = NULL, Sigma_0 = NULL, c_0 = 5, 
        nu_0 = 0.02, s_0 = 0.02, empirical_bayes = FALSE, 
        n_sample = 1000, thin = 1, burn = 1000, beta_burn = 10
    ) {
        X = X_list[[n]]
        y = y_list[[n]]
        
        lm_post_wo_equality(
            X, y, R, b, 
            beta_init, sigma_sq_init, 
            mu_0, Sigma_0, c_0, 
            nu_0, s_0, empirical_bayes, 
            n_sample, thin, burn, beta_burn
        )$beta_samples
    }, 
    R, b, 
    beta_init=beta_init, sigma_sq_init=sigma_sq_init, c_0=c_0, 
    nu_0=nu_0, s_0=s_0
)

beta_BICLS = sapply(fit_BICLS, rowMeans)
rownames(beta_BICLS) = paste0("beta", rep(1:3))

beta_diff_BICLS = beta_BICLS - beta_true
(sd_BICLS = apply(beta_diff_BICLS, 1, sd))

###########################
##### Results summary #####
###########################

# Table 1 
table_1 = cbind(
    rowMeans(beta_diff_BICLS^2),
    rowMeans(beta_diff_ICLS^2),
    rowMeans(beta_diff_ICLS^2) / rowMeans(beta_diff_BICLS^2),
    rowMeans(beta_diff_OLS^2),
    rowMeans(beta_diff_OLS^2) / rowMeans(beta_diff_BICLS^2), 
    apply(beta_diff_BICLS, 1, var),
    apply(beta_diff_ICLS, 1, var),
    apply(beta_diff_ICLS, 1, var) / apply(beta_diff_BICLS, 1, var),
    apply(beta_diff_OLS, 1, var),
    apply(beta_diff_OLS, 1, var) / apply(beta_diff_BICLS, 1, var)
)
colnames(table_1) = c("BICLS_MSE", "ICLS_MSE", "ratio", "OLS_MSE", "ratio", "BICLS_var", "ICLS_var", "ratio", "OLS_var", "ratio")

# Table S1
table_s1 = cbind(sd_BICLS, sd_ICLS, sd_OLS)

#########################
#####  Save results #####
#########################

write.csv(table_1, file = paste0("./data/table_1_", sample_size, ".csv"), fileEncoding = "UTF-8")
write.csv(table_s1, file = paste0("./data/table_s1_", sample_size, ".csv"), fileEncoding= "UTF-8")
save(list = ls(all.names=TRUE), file = paste0("./data/A1_", sample_size, ".Rdata"))

