library(MASS)
library(Matrix)
library(tmvmixnorm)
library(restriktor)

library(ggplot2)
library(reshape2)
library(ggpubr)

source("./post_inference_wo_equality.R")
source("./post_inference_with_equality.R")

#################################
#####     Simulation C      #####
#################################

repitition = 100
sample_size = 100
set.seed(10)

# Hyperparameter for BICGLM
beta_init = c(rep(1.5, 10), -3)
sigma_sq_init = 1
empirical_bayes = TRUE

# Set beta_true
beta_true = c(rep(1,10), 2)
p = length(beta_true)

# 100 time repeated data generation
X_list = lapply(
    rep(sample_size, repitition), function(n) {
        matrix(
            runif(n * p) - 0.5, 
            nrow = n, ncol = p
        )
    }

)

y_list = lapply(
    X_list, 
    function(X, beta) {
        n = nrow(X)
        rpois(n, exp(X %*% beta))
    }, beta_true
)


# Linear inequality constraints
R = matrix(0, nrow = 11, ncol = p)
R[1,] = rep(1, 11)
for(i in 1:10) R[i+1, i] = 1
b = c(12, rep(0.9, 10))
l = 1

################################
##### Parameter estimation #####
################################

# GLM estimates
fit_GLM = lapply(1:repitition, function(n) {
    X = X_list[[n]]
    y = y_list[[n]]
    glm(y ~ X - 1, family = poisson(link = "log"))
})


beta_GLM = sapply(fit_GLM, coef)
rownames(beta_GLM) = paste0("beta", 1:p)

beta_diff_GLM = beta_GLM - beta_true
(sd_GLM = apply(beta_diff_GLM, 1, sd))

# ICGLM estimates
fit_ICGLM = lapply(
    fit_GLM, 
    restriktor::restriktor, 
    constraints=R, rhs=b, neq=l
)


beta_ICGLM = sapply(fit_ICGLM, coef)
rownames(beta_ICGLM) = paste0("beta", 1:p)

beta_diff_ICGLM = beta_ICGLM - beta_true
(sd_ICGLM = apply(beta_diff_ICGLM, 1, sd))

# BICGLM estimates

fit_BICGLM = lapply(
    1:repitition, function(
        n, 
        R, b, l = 0, 
        beta_init = NULL, 
        mu_0 = NULL, Sigma_0 = NULL, c_0 = 3, empirical_bayes = FALSE, 
        n_sample = 1000, thin = 1, burn = 1000, beta_burn = 10
    ) {
        X = X_list[[n]]
        y = y_list[[n]]
        
        poisson_log_post_with_equality(
            X, y, R, b, l, 
            beta_init, 
            mu_0, Sigma_0, c_0, empirical_bayes, 
            n_sample, thin, burn, beta_burn
        )
    }, 
    R, b, l,
    beta_init=beta_init, empirical_bayes=TRUE
)

beta_BICGLM = sapply(fit_BICGLM, rowMeans)
rownames(beta_BICGLM) = paste0("beta", 1:p)

beta_diff_BICGLM = beta_BICGLM - beta_true
(sd_BICGLM = apply(beta_diff_BICGLM, 1, sd))

###########################
##### Results summary #####
###########################

# Table 4 
table_4 = cbind(
    rowMeans(beta_diff_BICGLM^2),
    rowMeans(beta_diff_ICGLM^2),
    rowMeans(beta_diff_ICGLM^2) / rowMeans(beta_diff_BICGLM^2),
    rowMeans(beta_diff_GLM^2),
    rowMeans(beta_diff_GLM^2) / rowMeans(beta_diff_BICGLM^2), 
    apply(beta_diff_BICGLM, 1, var),
    apply(beta_diff_ICGLM, 1, var),
    apply(beta_diff_ICGLM, 1, var) / apply(beta_diff_BICGLM, 1, var),
    apply(beta_diff_GLM, 1, var),
    apply(beta_diff_GLM, 1, var) / apply(beta_diff_BICGLM, 1, var)
)
colnames(table_4) = c("BICGLM_MSE", "ICGLM_MSE", "ratio", "GLM_MSE", "ratio", "BICGLM_var", "ICGLM_var", "ratio", "GLM_var", "ratio")

# Table S4
table_s4 = cbind(sd_BICGLM, sd_ICGLM, sd_GLM)

# fig s3
beta_list = lapply(
    1:p, function(i) {
        df = data.frame(cbind(beta_diff_BICGLM[i,], beta_diff_ICGLM[i,], beta_diff_GLM[i,]))
        colnames(df) = c("BICGLM", "ICGLM", "GLM")
        return (df)
    }
)

plot_list = lapply(
    1:p, 
    function(i) {
        ggplot(data = melt(beta_list[[i]]), aes(x=variable, y=value)) + 
            theme_bw() + 
            theme(text = element_text(size = 20)) +  
            ylab(paste0("Beta", i)) + 
            xlab("") + 
            geom_boxplot(aes()) + 
            geom_hline(aes(yintercept = 0), col = 'red')
    }
)

fig_s3 = ggpubr::ggarrange(plotlist = plot_list, nrow = 4, ncol = 3)
fig_s3
#########################
#####  Save results #####
#########################

ggsave(fig_s3, filename="./fig/fig_s3.png", width = 12, height = 16)

write.csv(table_4, file = "./data/table_4.csv", fileEncoding = "UTF-8")
write.csv(table_s4, file = paste0("./data/table_s4.csv"), fileEncoding= "UTF-8")
save(list = ls(all.names=TRUE), file = paste0("./data/B.Rdata"))
