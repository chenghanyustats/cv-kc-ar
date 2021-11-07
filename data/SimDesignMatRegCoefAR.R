################################################################################
# Complex-valued fMRI Simulation Study: AR(1) Design Matrix and Reg Coeff      #
# ./data/SimDesignMatRegCoefAR.R #
# Cheng-Han Yu, Marquette University                                           #
################################################################################
# Need sig2_old, SNR, CNR vector
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
X_cplx <- X.cplx
x_mod <- x.mod
# Regression coefficients with different (SNR, CNR) combinations
########################
phi <- pi / 4
alpha1 <- cos(phi)
alpha2 <- sin(phi)

# Case 2: fixed beta, move sigma
# =================================
beta0 <- sqrt(sig2_old) * SNR
beta2 <- sqrt(sig2_old) * CNR

Beta.re <- Beta.im <- list()

for (i in 1:length(SNR)) {
    for(j in 1:length(CNR)) {
        Beta.re[[(i - 1) * length(CNR) + j]] <- c(beta0[i], beta2[j]) * alpha1 
        Beta.im[[(i - 1) * length(CNR) + j]] <- c(beta0[i], beta2[j]) * alpha2
    }
}











