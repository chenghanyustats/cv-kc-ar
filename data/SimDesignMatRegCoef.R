################################################################################
# Complex-valued fMRI Simulation Study: IID Design Matrix and Reg Coeff        #
# ~/Research/code/Bezener/Analysis/Data/SimDesignMatRegCoef.R                  #
# Cheng-Han Yu, UC Santa Cruz                                                  #
################################################################################
library(Matrix)

# Design Matrix
########################
# generate the covariate and reference function (No trend)
# =========================================================
x <- cbind(rep(1, n), canonical)
X <- as.matrix(bdiag(x, x))
II <- (diag(n) - rep(1, n) %*% solve(t(rep(1, n)) %*% rep(1, n)) %*% t(rep(1, n)))
Xstar <- as.matrix(bdiag(rep(1, n), rep(1, n)))
IIc <- (diag(2 * n) - Xstar %*% solve(t(Xstar) %*% Xstar) %*% t(Xstar))
x.mod <- II %*% canonical
X.cplx <- IIc %*% as.matrix(bdiag(canonical, canonical))

# Regression coefficients
#########################
phi <- pi / 4
alpha1 <- cos(phi)
alpha2 <- sin(phi)

sig2 <- 3
SNR <- c(1, 300) / sqrt(sig2)
CNR <- c(2, 5) / sqrt(sig2)

beta0 <- sqrt(sig2) * SNR
beta2 <- sqrt(sig2) * CNR

Beta.re <- Beta.im <- Beta.mag <- list() 

for (i in 1:length(SNR)) {
    for(j in 1:length(CNR)) {
        Beta.re[[(i - 1) * length(CNR) + j]] <- c(beta0[i], beta2[j]) * alpha1 
        Beta.im[[(i - 1) * length(CNR) + j]] <- c(beta0[i], beta2[j]) * alpha2
    }
}

for (i in 1:length(SNR)) {
    for(j in 1:length(CNR)) {
        Beta.mag[[(i - 1) * length(CNR) + j]] <- c(beta0[i], beta2[j])
        # Beta.im[[(i - 1) * length(CNR) + j]] <- c(beta0[i], beta2[j]) * alpha2
    }
}

# for (i in 1:length(SNR)) {
#     for(j in 1:length(CNR)) {
#         Beta.mag[[(i - 1) * length(CNR) + j]] <- c(beta0[i], beta2[j])
#         # Beta.im[[(i - 1) * length(CNR) + j]] <- c(beta0[i], beta2[j]) * alpha2
#     }
# }

