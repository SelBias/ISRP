# library(rstudioapi)
library(MASS)
library(Matrix)
library(tmvmixnorm)
library(restriktor)
library(glmnet)

library(ggplot2)
# library(reshape2)
# library(ggpubr)

source("./posterior_inference.R")

#################################
#####     Simulation A1     #####
#################################

repitition = 100
sample_size = 50
set.seed(7)

# Hyperparameter for BICLS
beta_init = c(-0.5, -0.5, 0.5, 0.5, 0.5, rep(0, 25))
sigma_sq_init = 1
nu_0 = 0.02
s_0 = 0.02
rho = 0.0

# Set beta_true
beta_true = c(rep(-1, 2), rep(1, 28))
p = length(beta_true)

X_mean = rep(0, p)
X_precision = matrix(0, nrow=p, ncol=p)
for (i in 1:p) {
    for (j in 1:p) 
        X_precision[i,j] = rho^(abs(i-j))
}
X_cov = MASS::ginv(X_precision)
rm(list = c("i", "j"))

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
R = matrix(0, nrow = 5, ncol = p)
R[1,1] = R[3,3] = R[4,4] = R[5,5] = 1
R[1,2] = -2
R[2,1] = -1

b = rep(0,5)

################################
##### Parameter estimation #####
################################

# ridge estimates
fit_ridge = lapply(1:repitition, function(n) {
    X = X_list[[n]]
    y = y_list[[n]]

    lamb = cv.glmnet(X, y, intercept=FALSE, alpha=0)$lambda.min

    glmnet(X, y, intercept=F, alpha=0, lambda=lamb)
})

beta_ridge = sapply(fit_ridge, function(fit) {
    as.vector(fit$beta)
})
rownames(beta_ridge) = paste0("beta", rep(1:30))

beta_diff_ridge = beta_ridge - beta_true

# BICLS estimates

fit_BICLS = lapply(
    1:repitition, function(
        n, 
        R, b, 
        beta_init = NULL, sigma_sq_init = 1, 
        mu_0 = NULL, Sigma_0 = NULL, 
        nu_0 = 0.02, s_0 = 0.02, empirical_bayes = FALSE, 
        n_sample = 1000, thin = 1, burn = 1000, beta_burn = 10
    ) {
        X = X_list[[n]]
        y = y_list[[n]]

        lamb = cv.glmnet(X, y, intercept=FALSE, alpha=0)$lambda.min
        
        lm_post_wo_equality(
            X, y, R, b, 
            beta_init, sigma_sq_init, 
            mu_0, Sigma_0, c_0 = lamb, 
            nu_0, s_0, empirical_bayes, 
            n_sample, thin, burn, beta_burn
        )$beta_samples
    }, 
    R, b, 
    beta_init=beta_init, sigma_sq_init=sigma_sq_init, 
    nu_0=nu_0, s_0=s_0
)

beta_BICLS = sapply(fit_BICLS, rowMeans)
rownames(beta_BICLS) = paste0("beta", rep(1:p))

beta_diff_BICLS = beta_BICLS - beta_true

###########################
##### Results summary #####
###########################

# Table 2
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

# fig 2
beta1 = data.frame(cbind(beta_diff_BICLS[1,], beta_diff_ICLS[1,], beta_diff_OLS[1,]))
colnames(beta1) = c("BICLS", "ICLS", "OLS")
beta2 = data.frame(cbind(beta_diff_BICLS[2,], beta_diff_ICLS[2,], beta_diff_OLS[2,]))
colnames(beta2) = c("BICLS", "ICLS", "OLS")
beta3 = data.frame(cbind(beta_diff_BICLS[3,], beta_diff_ICLS[3,], beta_diff_OLS[3,]))
colnames(beta3) = c("BICLS", "ICLS", "OLS")

plot1 = ggplot(data = melt(beta1), aes(x=variable, y=value)) + 
    theme_bw() + 
    theme(text = element_text(size = 20)) + 
    ylab("Beta1") + 
    xlab("") + 
    geom_boxplot(aes()) + 
    geom_hline(aes(yintercept = 0), col = 'red')

plot2 = ggplot(data = melt(beta2), aes(x=variable, y=value)) + 
    theme_bw() + 
    theme(text = element_text(size = 20)) + 
    ylab("Beta2") + 
    xlab("") + 
    geom_boxplot(aes()) + 
    geom_hline(aes(yintercept = 0), col = 'red')

plot3 = ggplot(data = melt(beta3), aes(x=variable, y=value)) + 
    theme_bw() + 
    theme(text = element_text(size = 20)) +  
    ylab("Beta3") + 
    xlab("") + 
    geom_boxplot(aes()) + 
    geom_hline(aes(yintercept = 0), col = 'red')


fig_1 = ggpubr::ggarrange(plot1, plot2, plot3, nrow = 1)

#########################
#####  Save results #####
#########################

# ggsave(plot1, filename="./fig/fig1_beta1.png", width = 3.5, height = 7)
# ggsave(plot2, filename="./fig/fig1_beta2.png", width = 3.5, height = 7)
# ggsave(plot3, filename="./fig/fig1_beta3.png", width = 3.5, height = 7)
# ggsave(fig_1, filename="./fig/fig1.png", width = 10.5, height = 7)

# write.csv(table_1, file = paste0("./data/table_1_", sample_size, ".csv"), fileEncoding = "UTF-8")
# write.csv(table_s1, file = paste0("./data/table_s1_", sample_size, ".csv"), fileEncoding= "UTF-8")
# save(list = ls(all.names=TRUE), file = paste0("./data/A1_", sample_size, ".Rdata"))