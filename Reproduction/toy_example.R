library(MASS)
library(Matrix)
library(tmvmixnorm)
library(restriktor)

library(ggplot2)
library(reshape2)
library(ggpubr)

source("./post_inference_wo_equality.R")
source("./post_inference_with_equality.R")

repitition = 100
sample_size = 5
set.seed(10)

# Hyperparameter for BICGLM
beta_init = c(2,2,0)
empirical_bayes = TRUE

# Set beta_true
beta_true = c(1,1,-2)
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

X = X_list[[1]]
y = y_list[[1]]

# Linear inequality constraints
R = diag(1, 2, 3)
b = rep(0, 2)
l = 0

Bayes.icon.glm<-function(y,X,R,b,delta,n.samples=5000){
    n<-length(y)
    p<-ncol(X)
    m<-length(b)
    #Get Glm Betahat
    glmod<-glm(y~-1+X,offset = delta,family = poisson(link = "log"))   
    betaglm<-as.numeric(glmod$coefficients)
    #function for calculating information matrix
    Ibetafunc<-function(beta) {
        gamma2<-c()
        for(i in 1:n)
        {temp=crossprod(X[i,],beta)+delta[i]
        gamma2[i]=exp(temp)  #depends on the \psi function, exp for Poisson
        }
        GamaMat2<-diag(gamma2)
        Ibetahat<-t(X)%*%GamaMat2%*%X
        return(Ibetahat)
    }
    Ibetahat<-Ibetafunc(betaglm)
    mu1= betaglm
    Sigma1= chol2inv(chol(Ibetahat))      
    Sigmapost<-Sigma1
    mupost<-mu1+Sigma1%*%t(X)%*%y
    beta<- betaglm
    samples <- matrix(0,n.samples,p)
    Rstar<-rbind(R,-X)

    print(mu1)
    print(Sigmapost)
    #Start the MCMC sampler:
    for(i in 1:n.samples){
        #update U_i
        u<-c()
        for(l in 1:n) {
            up<- as.numeric(exp(-exp(t(X[l,])%*%beta+delta[l])))
            u[l]<-runif(1,0,up)
        }

        #update beta:     
        library(tmvmixnorm)
        bstar2<--log(-log(u))+delta
        bstar<-c(b,bstar2)
        mstar<-length(bstar)
        beta<-rtmvn(n=1, Mean=mupost, Sigmapost, D=Rstar, lower=bstar,
                    upper=rep(Inf,mstar),int=beta, burn=10)

        #store results:
        samples[i,]  <- c(beta)
        #print(i)
    }

    #return a list with the posterior samples:
    return(samples)
}

res = Bayes.icon.glm(y,X,R,b,delta = rep(0, sample_size))

source("./post_inference_wo_equality.R")
(res2 = poisson_log_post_wo_equality(X,y,R,b,beta_init = c(1,1,0), empirical_bayes=TRUE, thin = 1))

# source("./post_inference_wo_equality.R")
# poisson_log_post_wo_equality(
#     X, y, R, b, 
#     beta_init = beta_init, empirical_bayes = TRUE, 
#     n_sample = 100, burn = 1
# )
