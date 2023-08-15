# library(rstudioapi)
library(MASS)
library(Matrix)
library(tmvmixnorm)
library(restriktor)
library(glmnet)

library(ggplot2)
library(reshape2)
# library(ggpubr)

source("./post_inference_wo_equality.R")

#################################
#####     Simulation A1     #####
#################################

repitition = 100
sample_size = 50
set.seed(7)

# Hyperparameter for BICLS
beta_init = c(-0.5, -0.5, 0.5, 0.5, 0.5, rep(0, 25))
sigma_sq_init = 1
# c_0 = 2 / sqrt(2 * lambda)
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

lamb_list = unlist(lapply(1:repitition, function(n) {
    X = X_list[[n]]
    y = y_list[[n]]

    cv.glmnet(X, y, intercept=FALSE, alpha=0)$lambda.min
}))

# ridge estimates
fit_ridge = lapply(1:repitition, function(n) {
    X = X_list[[n]]
    y = y_list[[n]]

    glmnet(X, y, intercept=F, alpha=0, lambda=lamb_list[n])
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


        
        lm_post_wo_equality(
            X, y, R, b, 
            beta_init, sigma_sq_init, 
            mu_0, Sigma_0, c_0 = 2 / sqrt(2 * lamb_list[n]), 
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
table_2 = cbind(
    rowMeans(beta_diff_BICLS^2),
    rowMeans(beta_diff_ridge^2),
    rowMeans(beta_diff_ridge^2) / rowMeans(beta_diff_BICLS^2)
)

# View(beta_diff_BICLS[1:5,])
# View(beta_diff_ridge[1:5,])

sse_ridge = apply(
    beta_diff_ridge, 2, function(beta_diff) {
        sum(beta_diff^2)
})

sse_BICLS = apply(
    beta_diff_BICLS, 2, function(beta_diff) {
        sum(beta_diff^2)
})

# fig 2
sse = data.frame(cbind(sse_BICLS, sse_ridge))
colnames(sse) = c("BICLS", "Ridge")

fig_2 = ggplot(data = melt(sse), aes(x=variable, y=value)) + 
    theme_bw() + 
    theme(text = element_text(size = 20)) + 
    ylab("Sum of squared errors") + 
    xlab("") + 
    ylim(5, 20) + 
    geom_boxplot(aes())

fig_2

#########################
#####  Save results #####
#########################

ggsave(fig_2, filename=paste0("./fig/fig2_rho_", rho, ".png"), width = 6, height = 7)
write.csv(table_2, file = paste0("./data/table_2_rho_", rho, ".csv"), fileEncoding = "UTF-8")
save(list = ls(all.names=TRUE), file = paste0("./data/A2_rho_", rho, ".Rdata"))
