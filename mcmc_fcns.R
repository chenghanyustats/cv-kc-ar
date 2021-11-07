################################################################################
# Self-created functions for MCMC algorithm in the CV-KC-AR paper              #
# ./mcmc_fcns.R                                                                #
# Cheng-Han Yu, Marquette University                                           #
################################################################################
# library(emulator)

# Bezier Kernel 
dist_mat <- function(coord, grid) {
    N.voxels <- dim(coord)[1]  # number of voxels
    G <- dim(grid)[1]
    dist.mat <- matrix(0, nrow = N.voxels, ncol = G)
    for (i in 1:N.voxels) {
        for (j in 1:G) {
            dist.mat[i, j] <- sqrt(sum((coord[i, ] - grid[j, ]) ^ 2))
        }
    }
    return(dist.mat)
} 

bezier2 <- function(dist.mat, nu, phi) {
    # 2-d spherical Bezier kernel
    # nu: smoothness
    # range parameter 
    # a bend matrix
    # coord: resolution of data (voxel observation site)
    # grid: resolution of process convolution (region site)
    K <- array(0, dim = dim(dist.mat))
    less.mat <- dist.mat < phi
    K[less.mat] <- (1 - (dist.mat[less.mat] / phi) ^ 2) ^ nu
    return(K)
}

normalized_bezier2 <- function(dist.mat, nu, phi) {
    K <- bezier2(dist.mat, nu, phi)
    rowsum <- rowSums(K)
    cond <- rowsum != 0
    K[cond, ] <- K[cond, ] / rowsum[cond]
    return(K)
}


################################################################################
# Functions to implement MCMC algorithm in Bezener (2015)                      #
# Cheng-Han Yu                                                                 #
# UC Santa Cruz                                                                #
# 10/1/2016                                                                    #
################################################################################
# This script contains functions used for MCMC algorithm created by Bezener (2015)
# There is no modifications of their algorithms. 
# For modified algorithm and old version of code, check MCMCfcns.R
library(emulator)

# Useful fucntions
################################################################################
# Compute hermitian of a matrix
herm <- function(cplxmat) {
    return (t(Conj(cplxmat)))
}

# Compute power exponential covariance function used in 
power_exp_cov_fcn <- function(centroid_mat, r) {
    # compute power exponential covariance function
    # C(dist) = sig2 * exp(-|dist / phi| ^ nu)
    # grid <- expand.grid(idx, idx)
    distance <- as.matrix(dist(centroid_mat))
    # powerexp <- exp(-(distance / r))
    # return(powerexp)
    return(exp(-(distance / r)))
}

# Compute the centroid of each region
compute_centroid <- function(region_mat) {
    G <- max(region_mat)
    centroid_mat <- matrix(0, nrow = G, ncol = 2)
    for (g in 1:G) {
        xycoord <- which(region_mat == g, arr.ind = TRUE)
        # centroid_mat[g, ] <- apply(xycoord, 2, mean)
        centroid_mat[g, ] <- colMeans(xycoord)
    }
    return(centroid_mat)
}

################################################################################
# MCMC GP methods
################################################################################
logSGaS <- function(Ga, S) {
    # quad <- quad.form.inv(Ga, S)
    return(-(length(S) / 2) * log(quad.form.inv(Ga, S)))
}


logS_r <- function(Ga, S) {
    return((-0.5) * determinant(Ga, logarithm = TRUE)$modulus + 
               logSGaS(Ga, S))
}

# MH ratio of conditional of S
MHratio_cond_S <- function(S, S_new, r, ga, Ga, g, regionlabel, ag.vec) {
    # conditional of S (Metropolis ratio)
    # ============================================
    # S is S_(j): G by 1 spatial random effect
    # r is r_j: range parameter
    # ga is ga_(j): N by 1 indicator
    # Ga is Ga_j: G by G covariance function
    # g: region g
    # region: region matrix showing number g
    # (Need to compute Ga)
    
    # G <- length(S)
    ga_g <- ga[regionlabel == g]  # all ga in region g
    # len_ga <- length(ga_g)
    p_old <- 1 / (1 + exp(-(ag.vec[g] + S[g])))
    # q_old <- 1 - p_old
    p_new <- 1 / (1 + exp(-(ag.vec[g] + S_new[g])))
    # q_new <- 1 - p_new
    # no <- sum(ga_g)
    
    #     if (p_new == 0) print(paste("p_new = 0, g=", g, "Snew =", S_new[g]))
    #     if (q_new == 0) print(paste("q_new = 0, g=", g, "Snew =", S_new[g]))
    #     if (p_old == 0) print(paste("p_old = 0, g=", g, "Snew =", S_new[g]))
    #     if (q_old == 0) print(paste("q_old = 0, g=", g, "Snew =", S_new[g]))
    
    # Avoid -Inf and Inf when take log
    if (p_new == 0) p_new <- 1e-8
    if (p_new == 1) p_new <- 1 - 1e-8
    if (p_old == 0) p_old <- 1e-8
    if (p_old == 1) p_old <- 1 - 1e-8
    # if (p_new == 1e-8 && p_old == 1e-8)
    # A <- no * log(p_new / p_old) + (len_ga - no) * log(q_new / q_old)
    A <- sum(ga_g) * log(p_new / p_old) + 
        (length(ga_g) - sum(ga_g)) * log((1 - p_new) / (1 - p_old))
    # is.nan(A)
    # Ga_inv <- chol2inv(chol(Ga))
    #     if (t(S_new) %*% Ga_inv %*% S_new < 0) {
    #         print(paste("t(S)*gainv*S_new", (t(S_new) %*% Ga_inv %*% S_new), "g=", g))
    #     }
    #     if (t(S) %*% Ga_inv %*% S < 0) {
    #         print(paste("t(S)*gainv*S", t(S) %*% Ga_inv %*% S, "g=", g))
    #     }
    # cholG <- t(chol(Ga))
    # quad_new <- quad.form.inv(Ga, S_new)
    # quad <- quad.form.inv(Ga, S) 
    # quad_new <- quad.form(Ga_inv, S_new)
    # quad <- quad.form(Ga_inv, S)    
    # B <- -(G / 2) * log(quad_new / quad)
    B <- -(length(S) / 2) * log(quad.form.inv(Ga, S_new) / quad.form.inv(Ga, S))
    #     B <- -(G / 2) * log((t(S_new) %*% Ga_inv %*% S_new) / 
    #                             t(S) %*% Ga_inv %*% S)
    #     if (is.nan(A)) {
    #         print(paste("p_new=", p_new))
    #         print(paste("q_new=", q_new))
    #         print(paste("p_old=", p_old))
    #         print(paste("q_old=", q_old))
    #     }
    #     if (is.nan(B)) {
    #         print(paste("quad_new=", quad_new))
    #         print(paste("quad=", quad))
    #     }
    return(A + B)
}

# marginal of y | ga integrating beta and sig out
marginal_yga <- function(Yvec.v, XtY.v, XtX) {
    M_ga0 <- crossprod(Yvec.v)
    M_ga1 <- M_ga0 - quad.form.inv(XtX, XtY.v)
    #     if (!is.null(XtX)) {
    #         M_ga1 <- M_ga0 - quad.form.inv(XtX, XtY[, v])
    #     } else {
    #         M_ga1 <- M_ga0 - quad.form(XtXinv, XtY[, v])
    #     }
    #     print(paste("M0", M_ga0))
    #     print(paste("M1", M_ga1))
    return(list(M0 = M_ga0, M1 = M_ga1))
}

# update binary variable gamma (one at a time) function
update_ga <- function(S, regionlabel, v, M0, M1, ag.vec, is.cplx) {
    g <- regionlabel[v]
    p_ber <- 1 / (1 + exp(-(ag.vec[g] + S[g])))
    # print(paste("p_ber", p_ber))
    # q_ber <- 1 - p_ber
    po <- n / 2
    if (is.cplx) {
        C <- 1 + n/2
    } else {
        C <- sqrt(1 + n)
    }
    # p_star <- 1 / (1 + (q_ber / p_ber) * C * (M1 / M0) ^ pow)
    p_star <- 1 / (1 + ((1 - p_ber) / p_ber) * C * ((M1 / M0) ^ po))
    # p_star = 0.5
    # print(paste("p-star", p_star))
    #     if (v == 400) {
    #         print(paste("p-ber", p_ber))
    #         print(paste("p-star", p_star))
    #     }
    #     return(list(p_ber = p_ber, p_star = p_star, gav = rbinom(1, 1, p_star),
    #                 pow = (M1 / M0) ^ po))
    return(rbinom(1, 1, p_star))
}

# update_ga <- function(S, regionlabel, v, M0, M1, ag.vec) {
#     g <- regionlabel[v]
#     p_ber <- 1 / (1 + exp(-(ag.vec[g] + S[g])))
#     # print(paste("p_ber", p_ber))
#     # q_ber <- 1 - p_ber
#     
#     C <- 1 + (n /2)
#     po <- n / 2
#     # p_star <- 1 / (1 + (q_ber / p_ber) * C * (M1 / M0) ^ pow)
#     p_star <- 1 / (1 + ((1 - p_ber) / p_ber) * C * ((M1 / M0) ^ po))
#     # p_star = 0.5
#     # print(paste("p-star", p_star))
#     #     if (v == 400) {
#     #         print(paste("p-ber", p_ber))
#     #         print(paste("p-star", p_star))
#     #     }
# #     return(list(p_ber = p_ber, p_star = p_star, gav = rbinom(1, 1, p_star),
# #                 pow = (M1 / M0) ^ po))
#     return(rbinom(1, 1, p_star))
# }

logcond_w <- function(Ga, S, w, df) {
    # need to compute Ga
    # r ~ chi-sq(8) 
    # G <- length(S)
    # r <- exp(w)
    return(logS_r(Ga, S) + (df / 2) * w - (exp(w) / 2))
}

# main M-H algorithm
# include is.cplx and adpative tuning
mcmc_gp <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                    name.par, adapt = TRUE, tune.lst, keep.lst, tune.r, keep.r,
                    tune.len = 50, regionlabel, region_mat, df = 8, 
                    ag.vec = rep(0, max(regionlabel)),
                    target.accept.rate, is.cplx = TRUE) {
    ############################################################################
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # n.mcmc: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    #################################################################
    
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  
        # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    
    #-------------------------------
    # Adaptive tuning
    #------------------------------- 
    
    # for S
    keep <- keep.lst
    keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
    
    # for r
    # keep.psi <- list(psi = 0)
    keep.tmp.r <- keep.r
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    p <- 1
    G <- max(regionlabel)
    XtX <- crossprod(Xr)
    # XtXinv <- chol2inv(chol(XtX))
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    n <- ncol(Yvec)
    
    #-------------------------------
    # Storage
    #-------------------------------
    #     Gamma <- array(NA, c(n.mcmc, N.voxels, p))
    #     SS <- array(NA, c(n.mcmc, G, p))
    #     R <- array(NA, c(n.mcmc, p))
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    P_star <- matrix(0, nrow = n.mcmc, ncol = N.voxels)
    
    #-------------------------------
    # Starting values
    #-------------------------------
    S <- start.lst$S  # S is a COLUMN vector
    ga <- start.lst$ga
    r <- start.lst$r
    w <- log(r)
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    
    for (i in 1:n.mcmc) { 
        if (i %% 1000 == 0) cat("iter:", i, "\r")
        flush.console()
        
        
        # Update tuning parameter of r
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.r$r <- keep.tmp.r$r / Tb
            tune.r$r <- get.tune(tune.r$r, keep.tmp.r$r, i)
            keep.tmp$r <- 0
        }
        
        # Update parameters
        #----------------------------
        
        # ========== sample S =============
        Ga_cov <- power_exp_cov_fcn(centroid_mat = centroid, r = r)
        
        for (g in 1:G) {
            # Update tuning parameter
            #------------------------
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp[[g]] <- keep.tmp[[g]] / Tb
                tune.lst[[g]] <- get.tune(tune.lst[[g]], keep.tmp[[g]], i)
                keep.tmp[[g]] <- 0
            }
            
            #----------------------------
            # Random walk Normal proposal 
            #----------------------------
            # use random walk proposal newS = S + N(0, c) to update S
            
            # Sg_new <- rnorm(1, S[g], tune.lst[[g]])
            S_new <- S
            # S_new[g] <- Sg_new 
            S_new[g] <- rnorm(1, S[g], tune.lst[[g]])
            if (log(runif(1)) < MHratio_cond_S(S, S_new, r, ga, Ga_cov, g, 
                                               regionlabel, ag.vec)) {
                S[g] <- S_new[g]
                keep[[g]] <- keep[[g]] + 1
                keep.tmp[[g]] <- keep.tmp[[g]] + 1
            }
        }
        
        # ============ sample ga ==============
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            #             ga_res <- update_ga(S, regionlabel, v, MM$M0, MM$M1, ag.vec = ag.vec,
            #                       is.cplx = is.cplx)
            ga[v] <- update_ga(S, regionlabel, v, MM$M0, MM$M1, ag.vec = ag.vec,
                               is.cplx = is.cplx)
            # ga[v] <- update_ga(S, regionlabel, v, MM$M0, MM$M1, ag.vec = ag.vec)
            #             ga[v] <- ga_res$gav
            #             P_star[i, v] <- ga_res$p_star
        }
        
        # ============ sample w and hence r ==============
        new_w <- rnorm(1, w, tune.r$r)
        new_Ga_cov <- power_exp_cov_fcn(centroid_mat = centroid, 
                                        r = exp(new_w))
        
        if (log(runif(1)) < (logcond_w(new_Ga_cov, S, new_w, df) - 
                             logcond_w(Ga_cov, S, w, df))) {
            w <- new_w
            r <- exp(w)
            # accept_r <- accept_r + 1
            keep.r$r <- keep.r$r + 1
            keep.tmp.r$r <- keep.tmp.r$r + 1
        } 
        
        # store results
        #         SS[i, , 1] <- S
        #         Gamma[i, , 1] <- ga
        #         R[i, 1] <- r
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                draws[(i - burn) %/% thin, ] <- c(ga, S, r)
            } 
        }
        
    }
    
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.r$r <- keep.r$r / n.mcmc
    
    # Write output
    #--------------
    return(list(draws = draws, accept = c(keep, keep.r), start = start.lst, 
                tune = c(tune.lst, tune.r), burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx, P_star = P_star))
}




################################################################################
# Process Convolution Method 1: 
# kernel combination on logit transformation of the probability of activaion
################################################################################
# update binary variable (one at a time)
ga_update_pc1 <- function(v, M0, M1, Kw, n, regionlabel, ag.vec, is.cplx) {
    # g <- regionlabel[v]
    # n <- dim(Yvec)[2]
    # S <- kernel %*% w
    # s_v <- Kw[v, ]
    #     if(length(unique(ag.vec)) > 1) {
    #         p_ber <- 1 / (1 + exp(-(ag.vec[regionlabel[v]] + Kw[v])))
    #     } else {
    #         p_ber <- 1 / (1 + exp(-(ag.vec[1] + Kw[v])))
    #     }
    if(sum(ag.vec) != 0) {
        p_ber <- 1 / (1 + exp(-(ag.vec[regionlabel[v]] + Kw[v])))
    } else {
        p_ber <- 1 / (1 + exp(-Kw[v]))
    }
    if (p_ber == 0) p_ber <- 1e-8
    if (p_ber == 1) p_ber <- 1 - 1e-8
    
    # q_ber <- 1 - p_ber
    #     if (is.cplx) {
    #         C <- 1 + n
    #         pow <- n
    #     } else {
    #         C <- sqrt(1 + n)
    #         pow <- n / 2
    #     }
    po <- n / 2
    if (is.cplx) {
        C <- 1 + n / 2
    } else {
        C <- sqrt(1 + n)
    }
    p_star <- 1 / (1 + ((1 - p_ber) / p_ber) * C * (M1 / M0) ^ po)
    return(rbinom(1, 1, p_star))
}

# update latent process
logcond_latent <- function(w, wg_new, ga, a.tau, b.tau, K, Kw, g, ag.vec, 
                           is.old = TRUE, regionlabel) {
    # N.voxels <- length(Kw)
    #     if(length(unique(ag.vec)) > 1) {
    #         a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
    #     } else {
    #         a.vec <- rep(ag.vec[1], length(Kw))
    #     }
    #     if(sum(ag.vec) != 0) {
    #         a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
    #     } else {
    #         a.vec <- rep(0, length(Kw))
    #     }
    if(is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
        return(list(value = sumlogber - B))
    } else {
        Kw <- Kw + K[, g] * (wg_new - w[g])
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        w[g] <- wg_new
        B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
        return(list(value = sumlogber - B, Kw = Kw))
    }
} 

logcond_latent_block <- function(w, w_new, ga, a.tau, b.tau, K, Kw, g, ag.vec, 
                                 is.old = TRUE, regionlabel) {
    if(is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
        return(list(value = sumlogber - B))
    } else {
        Kw <- K %*% w_new
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        # w[g] <- wg_new
        B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w_new ^ 2) + b.tau)
        return(list(value = sumlogber - B, Kw = Kw))
    }
} 


logcond_latent_block2 <- function(w, w_new, ga, a.tau, b.tau, K, Kw, w.idx, ag.vec, 
                                  is.old = TRUE, regionlabel) {
    if(is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
        return(list(value = sumlogber - B))
    } else {
        w[w.idx] <- w_new
        Kw <- K %*% w
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        # w[g] <- wg_new
        B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
        return(list(value = sumlogber - B, Kw = Kw))
    }
} 


# update kernel parameter (gamma prior Ga(phi | a.phi, b.phi))
logcond_psi <- function(psi, w, Kw, a.phi, b.phi, dist.mat, nu, ga, ag.vec, 
                        is.old = TRUE, regionlabel) {
    # N.voxels <- length(Kw)
    #     if(length(unique(ag.vec)) > 1) {
    #         a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
    #     } else {
    #         a.vec <- rep(ag.vec[1], length(Kw))
    #     }
    #     if(sum(ag.vec) != 0) {
    #         a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
    #     } else {
    #         a.vec <- rep(0, length(Kw))
    #     }
    phi <- exp(psi)
    if (is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
        loggamma <- dgamma(1 / phi, a.phi, b.phi, log = TRUE) - 2 * log(phi)
        return(list(value = logsumber + loggamma + psi))
    } else {
        K <- bezier2(dist.mat = dist.mat, nu = nu, phi = phi) 
        Kw <- K %*% w
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
        loggamma <- dgamma(1 / phi, a.phi, b.phi, log = TRUE) - 2 * log(phi)
        return(list(value = logsumber + loggamma + psi, K = K, Kw = Kw))
    }
}


marginal_yga_mat <- function(Yvec, XtY, XtX) {
    library(emulator)
    Nvoxels <- nrow(Yvec)
    M_ga0.mat <- rep(0, Nvoxels)
    M_ga1.mat <- rep(0, Nvoxels)
    for (v in 1:Nvoxels) {
        M_ga0 <- crossprod(Yvec[v, ])
        M_ga1 <- M_ga0 - quad.form.inv(XtX, XtY[, v])
        #     if (!is.null(XtX)) {
        #         M_ga1 <- M_ga0 - quad.form.inv(XtX, XtY[, v])
        #     } else {
        #         M_ga1 <- M_ga0 - quad.form(XtXinv, XtY[, v])
        #     }
        M_ga0.mat[v] <- M_ga0
        M_ga1.mat[v] <- M_ga1
    }
    return(list(M0 = M_ga0.mat, M1 = M_ga1.mat))
}

logcond_ga_pc1 <- function(MM_mat, draws, regionlabel, dist.mat, nu, 
                           ag.vec, n, is.cplx) {
    po <- n / 2
    Mv_logsum <- 0
    Nvoxel <- length(regionlabel)
    G <- max(regionlabel)
    post_ga <- draws[, 1:Nvoxel]
    post_w <- draws[, (Nvoxel+1):(Nvoxel+G)]
    post_phi <- draws[, ncol(draws)]
    post_prob <- apply(post_ga, 2, mean)
    activation <- post_prob > 0.8722
    active_idx <- which(post_prob > 0.8722)
    active_number <- sum(activation)
    for (v in 1:Nvoxel) {
        if (activation[v]) {
            logMv <- log(MM_mat$M1[v])
            Mv_logsum <- Mv_logsum + logMv
            
        } else {
            logMv <- log(MM_mat$M0[v])
            Mv_logsum <- Mv_logsum + logMv
        }
    }
    if(is.cplx) {
        qlogC <- -active_number * log(1 + n / 2)
    } else {
        qlogC <- -(active_number / 2) * log(1 + n)
    }
    active_ga <- post_ga[, activation]
    prob_vec <- rep(0, active_number)
    regionlabel_act <- regionlabel[active_idx]
    # ag_active <- ag.vec[active_idx]
    
    postmean.w <- apply(post_w, 2, mean)
    postmean.phi <- mean(post_phi)
    Ker <- bezier2(dist.mat, nu, postmean.phi)
    post_Kw <- Ker %*% postmean.w
    post_Kw_active <- post_Kw[active_idx]
    if(sum(ag.vec) != 0) {
        for (i in 1:active_number) {
            prob_vec[i] <- 1 / (1 + exp(-(ag.vec[regionlabel_act[i]]  + 
                                              post_Kw_active[i])))
        }
    } else {
        prob_vec[i] <- 1 / (1 + exp(-post_Kw_active[i]))
    }
    value <- qlogC - po * Mv_logsum + sum(dbinom(active_ga, 1, prob = prob_vec,
                                                 log = TRUE))
    print(value)
    return(value)
}

logcond_ga_pc2 <- function(MM_mat, draws, regionlabel, dist.mat, nu, n, is.cplx) {
    po <- n / 2
    Mv_logsum <- 0
    Nvoxel <- length(regionlabel)
    G <- max(regionlabel)
    post_ga <- draws[, 1:Nvoxel]
    post_w <- draws[, (Nvoxel+1):(Nvoxel+G)]
    post_phi <- draws[, ncol(draws)]
    post_prob <- apply(post_ga, 2, mean)
    activation <- post_prob > 0.8722
    active_idx <- which(post_prob > 0.8722)
    active_number <- sum(activation)
    for (v in 1:Nvoxel) {
        if (activation[v]) {
            logMv <- log(MM_mat$M1[v])
            Mv_logsum <- Mv_logsum + logMv
            
        } else {
            logMv <- log(MM_mat$M0[v])
            Mv_logsum <- Mv_logsum + logMv
        }
    }
    if(is.cplx) {
        qlogC <- -active_number * log(1 + n / 2)
    } else {
        qlogC <- -(active_number / 2) * log(1 + n)
    }
    active_ga <- post_ga[, activation]
    prob_vec <- rep(0, active_number)
    regionlabel_act <- regionlabel[active_idx]
    # ag_active <- ag.vec[active_idx]
    
    postmean.w <- apply(post_w, 2, mean)
    postmean.phi <- mean(post_phi)
    Ker <- normalized_bezier2(dist.mat, nu, postmean.phi)
    post_Kw <- Ker %*% postmean.w
    post_Kw_active <- post_Kw[active_idx]
    value <- qlogC - po * Mv_logsum + sum(dbinom(active_ga, 1, 
                                                 prob = post_Kw_active,
                                                 log = TRUE))
    print(value)
    return(value)
}

update_wg <- function(w, ga, tau2.vec, K, Kw, g, ag.vec, is.old = TRUE) {
    if(is.old) {
        p_bern <- 1 / (1 + exp(-(ag.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        sumlognorm <- sum(dnorm(w, 0, sqrt(tau2.vec), log = TRUE))
        # B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
        return(list(value = sumlogber + sumlognorm))
    } else {
        Kw <- K %*% w
        p_bern <- 1 / (1 + exp(-(ag.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        sumlognorm <- sum(dnorm(w, 0, sqrt(tau2.vec), log = TRUE))
        return(list(value = sumlogber + sumlognorm, Kw = Kw))
    }
}

# update_tau
update_tau2 <- function(a.tau, b.tau, w) {
    a.shape <- 0.5 + a.tau
    b.scale <- 0.5 * w ^ 2 + b.tau
    1 / rgamma(length(w), a.shape, b.scale)
}




mcmc_pc1 <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                     name.par, adapt = TRUE, tune.lst, keep.lst, 
                     tune.psi, keep.psi,
                     tune.len = 30, regionlabel, region_mat, nu = 2, 
                     a.tau = 1, b.tau = 1,
                     a.phi = 1, b.phi = 1, 
                     ag.vec = rep(0, max(regionlabel)),
                     target.accept.rate, is.cplx) {
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #-------------------------------
    # for w
    keep <- keep.lst
    keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
    
    # for psi
    # keep.psi <- list(psi = 0)
    keep.tmp.psi <- keep.psi
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    p <- 1
    G <- max(regionlabel)
    XtX <- crossprod(Xr)
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    n <- ncol(Yvec)
    MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- start.lst$w
    ga <- start.lst$ga
    phi <- start.lst$phi
    # z <- log(w / (1 - w))
    psi <- log(phi)
    
    # create kernel
    Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    Kw <- Ker %*% w
    # tau2.vec <- rep(1, G)
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) { 
        if (i %% 1000 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
            tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
            keep.tmp.psi$psi <- 0
        } 	
        
        
        # Update parameters
        #----------------------------
        
        # ========== sample w =============
        for (g in 1:G) {
            
            # Update tuning parameter
            #-------------------------
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp[[g]] <- keep.tmp[[g]] / Tb
                tune.lst[[g]] <- get.tune(tune.lst[[g]], keep.tmp[[g]], i)
                keep.tmp[[g]] <- 0
            }
            
            #----------------------------
            # Random walk Normal proposal 
            #----------------------------
            
            # use random walk proposal newS = S + N(0, c) to update S
            wg_new <- rnorm(1, w[g], tune.lst[[g]])
            # w_new <- w
            # w_new[g] <- wg_new 
            #             logcon.new <- logcond_latent(w_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
            #                                          ag.vec, is.old = FALSE, 
            #                                          regionlabel = regionlabel)
            #             logcon.old <- logcond_latent(w, ga, a.tau, b.tau, K = Ker, Kw, g, 
            #                                          ag.vec, is.old = TRUE, 
            #                                          regionlabel = regionlabel)
            
            logcon.new <- logcond_latent(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
                                         ag.vec, is.old = FALSE, 
                                         regionlabel = regionlabel)
            logcon.old <- logcond_latent(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
                                         ag.vec, is.old = TRUE, 
                                         regionlabel = regionlabel)
            
            
            if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
                w[g] <- wg_new
                Kw <- logcon.new$Kw
                keep[[g]] <- keep[[g]] + 1
                keep.tmp[[g]] <- keep.tmp[[g]] + 1
            }
            
            #             wg_new <- rnorm(1, w[g], tune.lst[[g]])
            #             w_new <- w
            #             w_new[g] <- wg_new 
            #             logcon.new <- update_wg(w_new, ga, tau2.vec, K = Ker, Kw, g, 
            #                                          ag.vec, is.old = FALSE)
            #             logcon.old <- update_wg(w, ga, tau2.vec, K = Ker, Kw, g, 
            #                                          ag.vec, is.old = TRUE)
            #             if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            #                 w[g] <- wg_new
            #                 Kw <- logcon.new$Kw
            #                 keep[[g]] <- keep[[g]] + 1
            #                 keep.tmp[[g]] <- keep.tmp[[g]] + 1
            #             }
        }
        
        # ============sample tau2 =============
        
        # tau2.vec <- update_tau2(a.tau, b.tau, w)
        
        # ============ sample ga ==============
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            ga[v] <- ga_update_pc1(v, MM$M0, MM$M1, Kw, n, regionlabel, ag.vec, 
                                   is.cplx)
        }
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        
        logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
                                  dist.mat = dist.mat, nu = nu, ga = ga, 
                                  ag.vec, is.old = FALSE, 
                                  regionlabel = regionlabel)
        logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
                                  dist.mat = dist.mat, nu = nu, ga = ga,
                                  ag.vec, is.old = TRUE, 
                                  regionlabel = regionlabel)
        
        #         logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
        #                                   dist.mat = dist.mat, nu = nu, ga = ga, 
        #                                   ag.vec, is.old = FALSE, regionlabel)
        #         logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
        #                                   dist.mat = dist.mat, nu = nu, ga = ga,
        #                                   ag.vec, is.old = TRUE, regionlabel)
        
        
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        } 
        # cat("r=", r, "\n")
        
        # store results
        # W[i, , 1] <- w
        # Gamma[i, , 1] <- ga
        # Phi[i, 1] <- phi
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                draws[(i - burn) %/% thin, ] <- c(ga, w, phi)
            } 
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi$psi <- keep.psi$psi / n.mcmc
    
    # Write output
    #--------------
    return(list(draws = draws, accept = c(keep, keep.psi), start = start.lst, 
                tune = c(tune.lst, tune.psi), burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx))
    
    # return(list(W = W, Gamma = Gamma, Phi = Phi, tune = tune, tune.psi = tune.psi,
    #             accept = keep, accept_psi = keep.psi))
    
}

# whole block
mcmc_pc1_block <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                           name.par, adapt = TRUE, tune.w, keep.w, 
                           tune.psi, keep.psi,
                           tune.len = 30, regionlabel, region_mat, nu = 2, 
                           a.tau = 1, b.tau = 1,
                           a.phi = 1, b.phi = 1, 
                           ag.vec = rep(0, max(regionlabel)),
                           target.accept.rate, is.cplx) {
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(0.5, 1000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #-------------------------------
    # for w
    keep.tmp.w <- keep.w
    # keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
    
    # for psi
    # keep.psi <- list(psi = 0)
    keep.tmp.psi <- keep.psi
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    p <- 1
    G <- max(regionlabel)
    XtX <- crossprod(Xr)
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    n <- ncol(Yvec)
    MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    W <- matrix(NA, nrow = n.mcmc, ncol = G)
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    # z <- log(w / (1 - w))
    psi <- log(phi)
    
    # create kernel
    Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    Kw <- Ker %*% w
    # tau2.vec <- rep(1, G)
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) { 
        if (i %% 10 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
            tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
            keep.tmp.psi$psi <- 0
        } 	
        
        
        # Update parameters
        #----------------------------
        # Update tuning parameter of w
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {
            # Adaptive tuning
            keep.tmp.w$w <- keep.tmp.w$w / Tb
            tune.w$w <- get.tune(tune.w$w, keep.tmp.w$w, i)
            keep.tmp.w$w <- 0
        }
        
        # ========== sample w =============
        # block move using empirical covariance structure
        if (i > 1000) {
            Sigma.w <- cov(W[(i-1000):(i-1), ])
            Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w)$mat)
            # print(round(Sigma.w, 2))
            w_new <- as.vector(rmvn(1, w, tune.w$w * Sigma.w))
        } else {
            w_new <- as.vector(rmvn(1, w, tune.w$w * diag(G)))
        }
        # w_new <- as.vector(rmvn(1, w, tune.w$w * diag(G)))
        logcon.new <- logcond_latent_block(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
                                           ag.vec, is.old = FALSE, 
                                           regionlabel = regionlabel)
        logcon.old <- logcond_latent_block(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
                                           ag.vec, is.old = TRUE, 
                                           regionlabel = regionlabel)
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            w <- w_new
            Kw <- logcon.new$Kw
            keep.w$w <- keep.w$w + 1
            keep.tmp.w$w <- keep.tmp.w$w + 1
        }
        W[i, ] <- w
        #         for (g in 1:G) {
        #             
        #             # Update tuning parameter
        #             #-------------------------
        #             if(adapt == TRUE & i %% Tb == 0) {
        #                 # Adaptive tuning
        #                 keep.tmp[[g]] <- keep.tmp[[g]] / Tb
        #                 tune.lst[[g]] <- get.tune(tune.lst[[g]], keep.tmp[[g]], i)
        #                 keep.tmp[[g]] <- 0
        #             }
        #             
        #             #----------------------------
        #             # Random walk Normal proposal 
        #             #----------------------------
        #             
        #             # use random walk proposal newS = S + N(0, c) to update S
        #             wg_new <- rnorm(1, w[g], tune.lst[[g]])
        #             # w_new <- w
        #             # w_new[g] <- wg_new 
        #             #             logcon.new <- logcond_latent(w_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #             #                                          ag.vec, is.old = FALSE, 
        #             #                                          regionlabel = regionlabel)
        #             #             logcon.old <- logcond_latent(w, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #             #                                          ag.vec, is.old = TRUE, 
        #             #                                          regionlabel = regionlabel)
        #             
        #             logcon.new <- logcond_latent(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #                                          ag.vec, is.old = FALSE, 
        #                                          regionlabel = regionlabel)
        #             logcon.old <- logcond_latent(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #                                          ag.vec, is.old = TRUE, 
        #                                          regionlabel = regionlabel)
        #             
        #             
        #             if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
        #                 w[g] <- wg_new
        #                 Kw <- logcon.new$Kw
        #                 keep[[g]] <- keep[[g]] + 1
        #                 keep.tmp[[g]] <- keep.tmp[[g]] + 1
        #             }
        #             
        #             #             wg_new <- rnorm(1, w[g], tune.lst[[g]])
        #             #             w_new <- w
        #             #             w_new[g] <- wg_new 
        #             #             logcon.new <- update_wg(w_new, ga, tau2.vec, K = Ker, Kw, g, 
        #             #                                          ag.vec, is.old = FALSE)
        #             #             logcon.old <- update_wg(w, ga, tau2.vec, K = Ker, Kw, g, 
        #             #                                          ag.vec, is.old = TRUE)
        #             #             if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
        #             #                 w[g] <- wg_new
        #             #                 Kw <- logcon.new$Kw
        #             #                 keep[[g]] <- keep[[g]] + 1
        #             #                 keep.tmp[[g]] <- keep.tmp[[g]] + 1
        #             #             }
        #         }
        
        # ============sample tau2 =============
        
        # tau2.vec <- update_tau2(a.tau, b.tau, w)
        
        # ============ sample ga ==============
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            ga[v] <- ga_update_pc1(v, MM$M0, MM$M1, Kw, n, regionlabel, ag.vec, 
                                   is.cplx)
        }
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        
        logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
                                  dist.mat = dist.mat, nu = nu, ga = ga, 
                                  ag.vec, is.old = FALSE, 
                                  regionlabel = regionlabel)
        logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
                                  dist.mat = dist.mat, nu = nu, ga = ga,
                                  ag.vec, is.old = TRUE, 
                                  regionlabel = regionlabel)
        
        #         logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
        #                                   dist.mat = dist.mat, nu = nu, ga = ga, 
        #                                   ag.vec, is.old = FALSE, regionlabel)
        #         logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
        #                                   dist.mat = dist.mat, nu = nu, ga = ga,
        #                                   ag.vec, is.old = TRUE, regionlabel)
        
        
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        } 
        # cat("r=", r, "\n")
        
        # store results
        # W[i, , 1] <- w
        # Gamma[i, , 1] <- ga
        # Phi[i, 1] <- phi
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                # print(w)
                draws[(i - burn) %/% thin, ] <- c(ga, w, phi)
            } 
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    # keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi$psi <- keep.psi$psi / n.mcmc
    keep.w$w <- keep.w$w / n.mcmc
    # Write output
    #--------------
    return(list(draws = draws, accept = c(keep.w, keep.psi), start = start.lst, 
                tune = c(tune.w, tune.psi), burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx))
    
    # return(list(W = W, Gamma = Gamma, Phi = Phi, tune = tune, tune.psi = tune.psi,
    #             accept = keep, accept_psi = keep.psi))
    
}

# several blocks
mcmc_pc1_block2 <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                            name.par, adapt = TRUE, tune.w, keep.w, 
                            tune.psi, keep.psi,
                            tune.len = 30, regionlabel, region_mat, nu = 2, 
                            a.tau = 1, b.tau = 1,
                            a.phi = 1, b.phi = 1, 
                            ag.vec = rep(0, max(regionlabel)),
                            target.accept.rate, is.cplx) {
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(0.5, 1000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #-------------------------------
    # for w
    keep <- keep.w
    keep.tmp.w <- keep  # track MH accpetance rate for adaptive tuning
    # keep.tmp.w <- keep.w
    # keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
    
    # for psi
    # keep.psi <- list(psi = 0)
    keep.tmp.psi <- keep.psi
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    p <- 1
    G <- max(regionlabel)
    XtX <- crossprod(Xr)
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    n <- ncol(Yvec)
    MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    W <- matrix(NA, nrow = n.mcmc, ncol = G)
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    # z <- log(w / (1 - w))
    psi <- log(phi)
    
    # create kernel
    Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    Kw <- Ker %*% w
    # tau2.vec <- rep(1, G)
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) { 
        if (i %% 1000 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
            tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
            keep.tmp.psi$psi <- 0
        } 	
        
        
        # Update parameters
        #----------------------------
        # Update tuning parameter of w
        #------------------------
        #         if(adapt == TRUE & i %% Tb == 0) {
        #             # Adaptive tuning
        #             keep.tmp.w$w <- keep.tmp.w$w / Tb
        #             tune.w$w <- get.tune(tune.w$w, keep.tmp.w$w, i)
        #             keep.tmp.w$w <- 0
        #         }
        
        # ========== sample w =============
        for (j in 1:5) {
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp.w[[j]] <- keep.tmp.w[[j]] / Tb
                tune.w[[j]] <- get.tune(tune.w[[j]], keep.tmp.w[[j]], i)
                keep.tmp.w[[j]] <- 0
            }
            w.idx <- (5*(j-1) + 1:5)
            # block move using empirical covariance structure
            if (i > 1000) {
                Sigma.w <- cov(W[(i-1000):(i-1), w.idx])
                Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w)$mat)
                # print(round(Sigma.w, 2))
                w_new <- as.vector(mvnfast::rmvn(1, w[w.idx], tune.w[[j]] * Sigma.w))
            } else {
                w_new <- as.vector(mvnfast::rmvn(1, w[w.idx], tune.w[[j]] * diag(5)))
            }
            logcon.new <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
                                                ag.vec, is.old = FALSE, 
                                                regionlabel = regionlabel)
            logcon.old <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
                                                ag.vec, is.old = TRUE, 
                                                regionlabel = regionlabel)
            if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
                w[w.idx] <- w_new
                Kw <- logcon.new$Kw
                keep[[j]] <- keep[[j]] + 1
                keep.tmp.w[[j]] <- keep.tmp.w[[j]] + 1
            }
            W[i, w.idx] <- w_new
        }
        # 
        #         w_new <- as.vector(rmvn(1, w, tune.w$w * diag(G)))
        #         logcon.new <- logcond_latent_block(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #                                            ag.vec, is.old = FALSE, 
        #                                            regionlabel = regionlabel)
        #         logcon.old <- logcond_latent_block(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #                                            ag.vec, is.old = TRUE, 
        #                                            regionlabel = regionlabel)
        #         if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
        #             w <- w_new
        #             Kw <- logcon.new$Kw
        #             keep.w$w <- keep.w$w + 1
        #             keep.tmp.w$w <- keep.tmp.w$w + 1
        #         }
        # W[i, ] <- w
        #         for (g in 1:G) {
        #             
        #             # Update tuning parameter
        #             #-------------------------
        #             if(adapt == TRUE & i %% Tb == 0) {
        #                 # Adaptive tuning
        #                 keep.tmp[[g]] <- keep.tmp[[g]] / Tb
        #                 tune.lst[[g]] <- get.tune(tune.lst[[g]], keep.tmp[[g]], i)
        #                 keep.tmp[[g]] <- 0
        #             }
        #             
        #             #----------------------------
        #             # Random walk Normal proposal 
        #             #----------------------------
        #             
        #             # use random walk proposal newS = S + N(0, c) to update S
        #             wg_new <- rnorm(1, w[g], tune.lst[[g]])
        #             # w_new <- w
        #             # w_new[g] <- wg_new 
        #             #             logcon.new <- logcond_latent(w_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #             #                                          ag.vec, is.old = FALSE, 
        #             #                                          regionlabel = regionlabel)
        #             #             logcon.old <- logcond_latent(w, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #             #                                          ag.vec, is.old = TRUE, 
        #             #                                          regionlabel = regionlabel)
        #             
        #             logcon.new <- logcond_latent(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #                                          ag.vec, is.old = FALSE, 
        #                                          regionlabel = regionlabel)
        #             logcon.old <- logcond_latent(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
        #                                          ag.vec, is.old = TRUE, 
        #                                          regionlabel = regionlabel)
        #             
        #             
        #             if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
        #                 w[g] <- wg_new
        #                 Kw <- logcon.new$Kw
        #                 keep[[g]] <- keep[[g]] + 1
        #                 keep.tmp[[g]] <- keep.tmp[[g]] + 1
        #             }
        #             
        #             #             wg_new <- rnorm(1, w[g], tune.lst[[g]])
        #             #             w_new <- w
        #             #             w_new[g] <- wg_new 
        #             #             logcon.new <- update_wg(w_new, ga, tau2.vec, K = Ker, Kw, g, 
        #             #                                          ag.vec, is.old = FALSE)
        #             #             logcon.old <- update_wg(w, ga, tau2.vec, K = Ker, Kw, g, 
        #             #                                          ag.vec, is.old = TRUE)
        #             #             if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
        #             #                 w[g] <- wg_new
        #             #                 Kw <- logcon.new$Kw
        #             #                 keep[[g]] <- keep[[g]] + 1
        #             #                 keep.tmp[[g]] <- keep.tmp[[g]] + 1
        #             #             }
        #         }
        
        # ============sample tau2 =============
        
        # tau2.vec <- update_tau2(a.tau, b.tau, w)
        
        # ============ sample ga ==============
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            ga[v] <- ga_update_pc1(v, MM$M0, MM$M1, Kw, n, regionlabel, ag.vec, 
                                   is.cplx)
        }
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        
        logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
                                  dist.mat = dist.mat, nu = nu, ga = ga, 
                                  ag.vec, is.old = FALSE, 
                                  regionlabel = regionlabel)
        logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
                                  dist.mat = dist.mat, nu = nu, ga = ga,
                                  ag.vec, is.old = TRUE, 
                                  regionlabel = regionlabel)
        
        #         logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
        #                                   dist.mat = dist.mat, nu = nu, ga = ga, 
        #                                   ag.vec, is.old = FALSE, regionlabel)
        #         logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
        #                                   dist.mat = dist.mat, nu = nu, ga = ga,
        #                                   ag.vec, is.old = TRUE, regionlabel)
        
        
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        } 
        # cat("r=", r, "\n")
        
        # store results
        # W[i, , 1] <- w
        # Gamma[i, , 1] <- ga
        # Phi[i, 1] <- phi
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                # print(w)
                draws[(i - burn) %/% thin, ] <- c(ga, w, phi)
            } 
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi$psi <- keep.psi$psi / n.mcmc
    # keep.w$w <- keep.w$w / n.mcmc
    # keep <- lapply(keep, function(x) x / n.mcmc)
    # Write output
    #--------------
    return(list(draws = draws, accept = c(keep, keep.psi), start = start.lst, 
                tune = c(tune.w, tune.psi), burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx))
    
    # return(list(W = W, Gamma = Gamma, Phi = Phi, tune = tune, tune.psi = tune.psi,
    #             accept = keep, accept_psi = keep.psi))
    
}

# update region-specific a_g

logcond_a <- function(a.vec, ga, Kw, mu, Sig) {
    p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
    sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
    logNormal <- mvnfast::dmvn(a.vec, mu, sigma = Sig, log = TRUE)
    return(sumlogber + logNormal)
}



mcmc_pc1_connect <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                             name.par, adapt = TRUE, tune.w, keep.w, 
                             tune.psi, keep.psi, tune.a, keep.a,
                             tune.len = 30, regionlabel, region_mat, nu = 2, 
                             a.tau = 1, b.tau = 1,
                             a.phi = 1, b.phi = 1, 
                             ag.vec = rep(0, max(regionlabel)),
                             target.accept.rate, is.cplx) {
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(0.5, 1000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #-------------------------------
    # for w
    keep <- keep.w
    keep.tmp.w <- keep  # track MH accpetance rate for adaptive tuning
    
    # for psi
    keep.tmp.psi <- keep.psi
    
    # for a
    keep.tmp.a <- keep.a
    
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    p <- 1
    G <- max(regionlabel)
    XtX <- crossprod(Xr)
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    n <- ncol(Yvec)
    MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    W <- matrix(NA, nrow = n.mcmc, ncol = G)
    A <- matrix(NA, nrow = n.mcmc, ncol = G)
    Sig.array <- array(0, dim = c(n.mcmc, G, G))
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    avec <- start.lst$a
    Sig <- start.lst$Sig
    
    # create kernel
    Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    Kw <- Ker %*% w
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) { 
        if (i %% 10 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
            tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
            keep.tmp.psi$psi <- 0
        } 	
        
        # Update tuning parameter of a
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.a$a <- keep.tmp.a$a / Tb
            tune.a$a <- get.tune(tune.a$a, keep.tmp.a$a, i)
            keep.tmp.a$a <- 0
        } 	
        
        # ========== sample w =============
        for (j in 1:5) {
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp.w[[j]] <- keep.tmp.w[[j]] / Tb
                tune.w[[j]] <- get.tune(tune.w[[j]], keep.tmp.w[[j]], i)
                keep.tmp.w[[j]] <- 0
            }
            w.idx <- (5*(j-1) + 1:5)
            # block move using empirical covariance structure
            if (i > 500) {
                Sigma.w <- cov(W[(i-Tb):(i-1), w.idx])
                Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w)$mat)
                # print(round(Sigma.w, 2))
                w_new <- as.vector(rmvn(1, w[w.idx], tune.w[[j]] * Sigma.w))
            } else {
                w_new <- as.vector(rmvn(1, w[w.idx], tune.w[[j]] * diag(5)))
            }
            logcon.new <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
                                                ag.vec = avec, is.old = FALSE, 
                                                regionlabel = regionlabel)
            logcon.old <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
                                                ag.vec = avec, is.old = TRUE, 
                                                regionlabel = regionlabel)
            if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
                w[w.idx] <- w_new
                Kw <- logcon.new$Kw
                keep[[j]] <- keep[[j]] + 1
                keep.tmp.w[[j]] <- keep.tmp.w[[j]] + 1
            }
            W[i, w.idx] <- w[w.idx]
        }
        
        # ============ sample ga ==============
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            ga[v] <- ga_update_pc1(v, MM$M0, MM$M1, Kw, n, regionlabel, 
                                   ag.vec = avec, 
                                   is.cplx)
        }
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        
        logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
                                  dist.mat = dist.mat, nu = nu, ga = ga, 
                                  ag.vec = avec, is.old = FALSE, 
                                  regionlabel = regionlabel)
        logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
                                  dist.mat = dist.mat, nu = nu, ga = ga,
                                  ag.vec = avec, is.old = TRUE, 
                                  regionlabel = regionlabel)
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        } 
        
        # ============ sample a ==============
        if (i > 500) {
            Sigma.a <- cov(A[(i-Tb):(i-1), ])
            Sigma.a <- as.matrix(Matrix::nearPD(Sigma.a)$mat)
            # print(round(Sigma.w, 2))
            new_avec <- as.vector(rmvn(1, avec, tune.a$a * Sigma.a))
            # new_avec <- as.vector(rmvn(1, avec, tune.a$a * diag(G)))
        } else {
            new_avec <- as.vector(rmvn(1, avec, tune.a$a * diag(G)))
        }
        
        logcond_a_old <- logcond_a(avec, ga, Kw, rep(0, G), Sig)
        logcond_a_new <- logcond_a(new_avec, ga, Kw, rep(0, G), Sig)
        
        if (log(runif(1)) < (logcond_a_new - logcond_a_old)) {
            avec <- new_avec
            keep.a$a <- keep.a$a + 1
            keep.tmp.a$a <- keep.tmp.a$a + 1
        } 
        
        A[i, ] <- avec
        
        
        # ============ sample Sigma ==============
        # prior IW(G, diag(G))
        
        Sig <- MCMCpack::riwish(v = G, S = avec %*% t(avec) + diag(G))
        
        Sig.array[i, , ] <- Sig
        
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                # print(w)
                draws[(i - burn) %/% thin, ] <- c(ga, w, phi)
            } 
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi$psi <- keep.psi$psi / n.mcmc
    keep.a$a <- keep.a$a / n.mcmc
    # keep.w$w <- keep.w$w / n.mcmc
    # keep <- lapply(keep, function(x) x / n.mcmc)
    # Write output
    #--------------
    return(list(draws = draws, 
                A = A,
                Sig.array = Sig.array,
                accept = c(keep, keep.psi, keep.a), 
                start = start.lst, 
                tune = c(tune.w, tune.psi, tune.a), 
                burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx))
}

# ====== (This is not working well: bad mixing, nonidentifiable?)
logcond_w <- function(w, a.vec, ga, K, Kw, Sig, regionlabel, is.old) {
    if(!is.old) Kw <- K %*% w
    if(sum(a.vec) != 0) {
        a.vec <- sapply(1:length(Kw), function(x) a.vec[regionlabel[x]])
        p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
    } else {
        p_bern <- 1 / (1 + exp(-Kw))
    }
    
    sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
    logNormal <- mvnfast::dmvn(w, rep(0, length(w)), sigma = Sig, 
                               log = TRUE)
    if(is.old) {
        return(list(value = sumlogber + logNormal))
    } else {
        return(list(value = sumlogber + logNormal, Kw = Kw))
    }
}

mcmc_pc1_connect_w <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                               name.par, adapt = TRUE, tune.w, keep.w, 
                               tune.psi, keep.psi,
                               tune.len = 30, regionlabel, region_mat, nu = 2, 
                               a.phi = 1, b.phi = 1, 
                               ag.vec = rep(0, max(regionlabel)),
                               target.accept.rate, is.cplx) {
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(0.5, 1000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #-------------------------------
    # for w
    keep <- keep.w
    keep.tmp.w <- keep  # track MH accpetance rate for adaptive tuning
    
    # for psi
    # keep.tmp.psi <- keep.psi
    
    # for a
    keep.tmp.a <- keep.a
    
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    p <- 1
    G <- max(regionlabel)
    XtX <- crossprod(Xr)
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    n <- ncol(Yvec)
    MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    W <- matrix(NA, nrow = n.mcmc, ncol = G)
    A <- matrix(NA, nrow = n.mcmc, ncol = G)
    Sig.lst <- list()
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    Sig <- start.lst$Sig
    
    # create kernel
    Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    Kw <- Ker %*% w
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) { 
        if (i %% 10 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        #         if(adapt == TRUE & i %% Tb == 0) {  
        #             # Adaptive tuning
        #             keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
        #             tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
        #             keep.tmp.psi$psi <- 0
        #         } 	
        
        # Update tuning parameter of w
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {
            # Adaptive tuning
            keep.tmp.w$w <- keep.tmp.w$w / Tb
            tune.w$w <- get.tune(tune.w$w, keep.tmp.w$w, i)
            keep.tmp.w$w <- 0
        }	
        
        # ========== sample w =============
        # block move using empirical covariance structure
        if (i > 500) {
            Sigma.w <- cov(W[(i-Tb):(i-1), ])
            Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w)$mat)
            # print(round(Sigma.w, 2))
            w_new <- as.vector(rmvn(1, w, tune.w$w * Sigma.w))
        } else {
            w_new <- as.vector(rmvn(1, w, tune.w$w * diag(G)))
        }
        logcon.new <- logcond_w(w_new, ag.vec, ga, K = Ker, Kw, 
                                Sig, regionlabel, is.old = FALSE)
        logcon.old <- logcond_w(w, ag.vec, ga, K = Ker, Kw, 
                                Sig, regionlabel, is.old = TRUE)
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            w <- w_new
            Kw <- logcon.new$Kw
            keep.w$w <- keep.w$w + 1
            keep.tmp.w$w <- keep.tmp.w$w + 1
        }
        W[i, ] <- w
        
        # ============ sample ga ==============
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            ga[v] <- ga_update_pc1(v, MM$M0, MM$M1, Kw, n, regionlabel, 
                                   ag.vec = avec, 
                                   is.cplx)
        }
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        
        logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
                                  dist.mat = dist.mat, nu = nu, ga = ga, 
                                  ag.vec = avec, is.old = FALSE, 
                                  regionlabel = regionlabel)
        logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
                                  dist.mat = dist.mat, nu = nu, ga = ga,
                                  ag.vec = avec, is.old = TRUE, 
                                  regionlabel = regionlabel)
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        } 
        
        # ============ sample Sigma ==============
        # prior IW(G, diag(G))
        
        Sig <- MCMCpack::riwish(v = G, S = w %*% t(w) + diag(G))
        
        Sig.lst[[i]] <- Sig
        
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                # print(w)
                draws[(i - burn) %/% thin, ] <- c(ga, w, phi)
            } 
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi$psi <- keep.psi$psi / n.mcmc
    # keep.a$a <- keep.a$a / n.mcmc
    keep.w$w <- keep.w$w / n.mcmc
    # keep <- lapply(keep, function(x) x / n.mcmc)
    # Write output
    #--------------
    return(list(draws = draws, 
                # A = A,
                Sig.lst = Sig.lst,
                accept = c(keep.w, keep.psi), 
                start = start.lst, 
                tune = c(tune.w, tune.psi), 
                burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx))
}
# =====
mcmc_pc1_connect_mu <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                                name.par, adapt = TRUE, tune.w, keep.w, 
                                tune.psi, keep.psi, tune.a, keep.a,
                                tune.len = 30, regionlabel, region_mat, nu = 2, 
                                a.phi = 1, b.phi = 1, 
                                a.tau = 1, b.tau = 1,
                                ag.vec = rep(0, max(regionlabel)),
                                target.accept.rate, is.cplx) {
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(0.5, 1000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #-------------------------------
    # for w
    keep <- keep.w
    keep.tmp.w <- keep  # track MH accpetance rate for adaptive tuning
    
    # for psi
    keep.tmp.psi <- keep.psi
    
    # for a
    keep.tmp.a <- keep.a
    
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    p <- 1
    G <- max(regionlabel)
    XtX <- crossprod(Xr)
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    n <- ncol(Yvec)
    MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    W <- matrix(NA, nrow = n.mcmc, ncol = G)
    A <- matrix(NA, nrow = n.mcmc, ncol = G)
    Mu <- matrix(NA, nrow = n.mcmc, ncol = G)
    Sig.array <- array(0, dim = c(n.mcmc, G, G))
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    avec <- start.lst$a
    mu <- start.lst$mu
    Sig <- start.lst$Sig
    
    # create kernel
    Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    Kw <- Ker %*% w
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) { 
        if (i %% 10 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
            tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
            keep.tmp.psi$psi <- 0
        } 	
        
        # Update tuning parameter of a
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.a$a <- keep.tmp.a$a / Tb
            tune.a$a <- get.tune(tune.a$a, keep.tmp.a$a, i)
            keep.tmp.a$a <- 0
        } 	
        
        # ========== sample w =============
        for (j in 1:5) {
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp.w[[j]] <- keep.tmp.w[[j]] / Tb
                tune.w[[j]] <- get.tune(tune.w[[j]], keep.tmp.w[[j]], i)
                keep.tmp.w[[j]] <- 0
            }
            w.idx <- (5*(j-1) + 1:5)
            # block move using empirical covariance structure
            if (i > 500) {
                Sigma.w <- cov(W[(i-Tb):(i-1), w.idx])
                Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w)$mat)
                # print(round(Sigma.w, 2))
                w_new <- as.vector(rmvn(1, w[w.idx], tune.w[[j]] * Sigma.w))
            } else {
                w_new <- as.vector(rmvn(1, w[w.idx], tune.w[[j]] * diag(5)))
            }
            logcon.new <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
                                                ag.vec = avec, is.old = FALSE, 
                                                regionlabel = regionlabel)
            logcon.old <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
                                                ag.vec = avec, is.old = TRUE, 
                                                regionlabel = regionlabel)
            if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
                w[w.idx] <- w_new
                Kw <- logcon.new$Kw
                keep[[j]] <- keep[[j]] + 1
                keep.tmp.w[[j]] <- keep.tmp.w[[j]] + 1
            }
            W[i, w.idx] <- w_new
        }
        
        # ============ sample ga ==============
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            ga[v] <- ga_update_pc1(v, MM$M0, MM$M1, Kw, n, regionlabel, 
                                   ag.vec = avec, 
                                   is.cplx)
        }
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        
        logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
                                  dist.mat = dist.mat, nu = nu, ga = ga, 
                                  ag.vec = avec, is.old = FALSE, 
                                  regionlabel = regionlabel)
        logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
                                  dist.mat = dist.mat, nu = nu, ga = ga,
                                  ag.vec = avec, is.old = TRUE, 
                                  regionlabel = regionlabel)
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        } 
        
        # ============ sample a ==============
        if (i > 500) {
            Sigma.a <- cov(A[(i-Tb):(i-1), ])
            Sigma.a <- as.matrix(Matrix::nearPD(Sigma.a)$mat)
            # print(round(Sigma.w, 2))
            new_avec <- as.vector(rmvn(1, avec, tune.a$a * Sigma.a))
            # new_avec <- as.vector(rmvn(1, avec, tune.a$a * diag(G)))
        } else {
            new_avec <- as.vector(rmvn(1, avec, tune.a$a * diag(G)))
        }
        
        logcond_a_old <- logcond_a(avec, ga, Kw, mu, Sig)
        logcond_a_new <- logcond_a(new_avec, ga, Kw, mu, Sig)
        
        if (log(runif(1)) < (logcond_a_new - logcond_a_old)) {
            avec <- new_avec
            keep.a$a <- keep.a$a + 1
            keep.tmp.a$a <- keep.tmp.a$a + 1
        } 
        
        A[i, ] <- avec
        
        
        # ============ sample mu ==============
        Sig.inv <- chol2inv(chol(Sig))
        mu.var <- chol2inv(chol(Sig.inv + diag(G)))
        mu <- rmvn(1, mu.var %*% Sig.inv %*% avec, mu.var)
        Mu[i, ] <- mu
        
        # ============ sample Sigma ==============
        # prior IW(G, diag(G))
        S1 <- (avec - as.vector(mu)) %*% t(avec - as.vector(mu)) + diag(G)
        Sig <- MCMCpack::riwish(v = G, S = S1)
        
        Sig.array[i, , ] <- Sig
        
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                # print(w)
                draws[(i - burn) %/% thin, ] <- c(ga, w, phi)
            } 
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi$psi <- keep.psi$psi / n.mcmc
    keep.a$a <- keep.a$a / n.mcmc
    # keep.w$w <- keep.w$w / n.mcmc
    # keep <- lapply(keep, function(x) x / n.mcmc)
    # Write output
    #--------------
    return(list(draws = draws, 
                A = A,
                Mu = Mu,
                Sig.array = Sig.array,
                accept = c(keep, keep.psi, keep.a), 
                start = start.lst, 
                tune = c(tune.w, tune.psi, tune.a), 
                burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx))
}




mcmc_pc1_multires <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                              name.par, adapt = TRUE, tune.lst, keep.lst, 
                              tune.lstf, keep.lstf, 
                              tune.psi, keep.psi,
                              tune.len = 30, 
                              regionlabel, region_mat, 
                              regionlabel.finer, region_mat_finer, 
                              nu = 2, 
                              a.tau = 1, b.tau = 1,
                              a.phi = 1, b.phi = 1, 
                              a.tauf = 1, b.tauf = 1,
                              a.phif = 1, b.phif = 1, 
                              ag.vec = rep(0, max(regionlabel)),
                              target.accept.rate, is.cplx) {
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate,
                         a = min(1, 10000 / sqrt(k))) {  # adaptive tuning
        # a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #-------------------------------
    # for w
    keep <- keep.lst
    keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
    
    # for wf
    keepf <- keep.lstf
    keep.tmpf <- keepf  # track MH accpetance rate for adaptive tuning
    
    # for psi
    # keep.psi <- list(psi = 0)
    keep.tmp.psi <- keep.psi
    
    
    
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    N <- sqrt(N.voxels)
    p <- 1
    G <- max(regionlabel)
    # regionlabelf <- relabel(regionlabel.finer)  # define relabel function!
    # Gf <- max(regionlabelf)
    region.idx <- regionlabel.finer != 0
    Gf <- length(unique(regionlabel.finer[region.idx]))
    XtX <- crossprod(Xr)
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    centroidf <- compute_centroid(region_mat_finer)
    centroidf2 <- centroidf[unique(regionlabel.finer[region.idx]), ]
    Coord2 <- Coord[region.idx, ]
    dist.matf <- dist_mat(coord = Coord2, grid = centroidf2)
    
    n <- ncol(Yvec)
    MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- start.lst$w
    wf <- start.lst$wf
    ga <- start.lst$ga
    phi <- start.lst$phi
    phif <- start.lst$phif
    # z <- log(w / (1 - w))
    psi <- log(phi)
    psif <- log(phif)
    
    # create kernel
    Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    Kerf <- bezier2(dist.matf, nu = nu, phi = phif)
    Kw <- Ker %*% w
    Kwf <- matrix(0, N.voxels, 1)
    Kwf[region.idx, ] <- Kerf %*% wf
    
    # tau2.vec <- rep(1, G)
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) { 
        if (i %% 1000 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
            tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
            keep.tmp.psi$psi <- 0
        } 	
        
        if(adapt == TRUE & i %% Tb == 0) {  
            # Adaptive tuning
            keep.tmp.psi$psif <- keep.tmp.psi$psif / Tb
            tune.psi$psif <- get.tune(tune.psi$psif, keep.tmp.psi$psif, i)
            keep.tmp.psi$psif <- 0
        } 
        
        # Update parameters
        #----------------------------
        
        # ========== sample w =============
        for (g in 1:G) {
            
            # Update tuning parameter
            #-------------------------
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp[[g]] <- keep.tmp[[g]] / Tb
                tune.lst[[g]] <- get.tune(tune.lst[[g]], keep.tmp[[g]], i)
                keep.tmp[[g]] <- 0
            }
            
            #----------------------------
            # Random walk Normal proposal 
            #----------------------------
            
            # use random walk proposal newS = S + N(0, c) to update S
            wg_new <- rnorm(1, w[g], tune.lst[[g]])
            # print(paste("wg_new = ", wg_new))
            #             logcond_wc <- function(w, wg_new, ga, a.tau, b.tau, K, Kw, Kwf, g, ag.vec, 
            #                                    is.old = TRUE, regionlabel)
            logcon.new <- logcond_wc(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, Kwf, g, 
                                     ag.vec, is.old = FALSE, 
                                     regionlabel = regionlabel)
            # print(paste("logcon.new = ", logcon.new))
            logcon.old <- logcond_wc(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, Kwf, g, 
                                     ag.vec, is.old = TRUE, 
                                     regionlabel = regionlabel)
            # print(paste("logcon.old = ", logcon.old))
            
            if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
                w[g] <- wg_new
                Kw <- logcon.new$Kw
                keep[[g]] <- keep[[g]] + 1
                keep.tmp[[g]] <- keep.tmp[[g]] + 1
            }
        }
        
        
        
        # ========== sample wf =============
        for (g in 1:Gf) {
            
            # Update tuning parameter
            #-------------------------
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmpf[[g]] <- keep.tmpf[[g]] / Tb
                tune.lstf[[g]] <- get.tune(tune.lstf[[g]], keep.tmpf[[g]], i)
                keep.tmpf[[g]] <- 0
            }
            
            #----------------------------
            # Random walk Normal proposal 
            #----------------------------
            
            # use random walk proposal newS = S + N(0, c) to update S
            wfg_new <- rnorm(1, wf[g], tune.lstf[[g]])
            # print(paste("wfg_new = ", wfg_new))
            logcon.newf <- logcond_wf(wf, wfg_new, ga, a.tauf, b.tauf, Kf = Kerf, Kw, Kwf, g, 
                                      ag.vec, is.old = FALSE, 
                                      regionlabel = regionlabel, region.idx)
            # print(paste("logcon.newf = ", logcon.newf))
            logcon.oldf <- logcond_wf(wf, wfg_new, ga, a.tauf, b.tauf, Kf = Kerf, Kw, Kwf, g, 
                                      ag.vec, is.old = TRUE, 
                                      regionlabel = regionlabel, region.idx)
            # print(paste("logcon.oldf = ", logcon.oldf))
            if (log(runif(1)) < (logcon.newf$value - logcon.oldf$value)) {
                wf[g] <- wfg_new
                Kwf <- logcon.newf$Kwf
                keepf[[g]] <- keepf[[g]] + 1
                keep.tmpf[[g]] <- keep.tmpf[[g]] + 1
            }
        }
        
        
        # ============sample tau2 =============
        
        # tau2.vec <- update_tau2(a.tau, b.tau, w)
        
        # ============ sample ga ==============
        Kw_combine <- Kw + Kwf
        
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            ga[v] <- ga_update_pc1(v, MM$M0, MM$M1, Kw_combine, n, regionlabel, ag.vec, 
                                   is.cplx)
        }
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        
        logcon.new <- logcond_psic(new_psi, w, Kw, Kwf, a.phi, b.phi, 
                                   dist.mat = dist.mat, nu = nu, ga = ga, 
                                   ag.vec, is.old = FALSE, 
                                   regionlabel = regionlabel)
        logcon.old <- logcond_psic(psi, w, Kw, Kwf, a.phi, b.phi,
                                   dist.mat = dist.mat, nu = nu, ga = ga,
                                   ag.vec, is.old = TRUE, 
                                   regionlabel = regionlabel)
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        } 
        
        # ============ sample psif and hence phif ==============
        new_psif <- rnorm(1, psif, tune.psi$psif)
        # new_psif <- rnorm(1, psif, 1)
        # print(paste("new_psif = ", new_psif))
        logcon.newf <- logcond_psif(new_psif, wf, Kw, Kwf, a.phif, b.phif, 
                                    dist.matf = dist.matf, nu = nu, ga = ga, 
                                    ag.vec, is.old = FALSE, 
                                    regionlabel = regionlabel, region.idx)
        # print(paste("logcon.newf = ", logcon.newf$value))
        logcon.oldf <- logcond_psif(psif, wf, Kw, Kwf, a.phif, b.phif,
                                    dist.matf = dist.matf, nu = nu, ga = ga,
                                    ag.vec, is.old = TRUE, 
                                    regionlabel = regionlabel, region.idx)
        # print(paste("logcon.oldf = ", logcon.oldf$value))
        
        if (log(runif(1)) < (logcon.newf$value - logcon.oldf$value)) {
            psif <- new_psif
            phif <- exp(psif)
            Kerf <- logcon.newf$Kf
            Kwf <- logcon.newf$Kwf
            keep.psi$psif <- keep.psi$psif + 1
            keep.tmp.psi$psif <- keep.tmp.psi$psif + 1
        } 
        
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                draws[(i - burn) %/% thin, ] <- c(ga, w, wf, phi, phif)
            } 
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi <- lapply(keep.psi, function(x) x / n.mcmc)
    keepf <- lapply(keepf, function(x) x / n.mcmc)
    
    # Write output
    #--------------
    return(list(draws = draws, accept = c(keep, keepf, keep.psi), start = start.lst, 
                tune = c(tune.lst, tune.lstf, tune.psi), burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx))
    
    # return(list(W = W, Gamma = Gamma, Phi = Phi, tune = tune, tune.psi = tune.psi,
    #             accept = keep, accept_psi = keep.psi))
    
} 


# update latent process

logcond_wc <- function(w, wg_new, ga, a.tau, b.tau, K, Kw, Kwf, g, ag.vec, 
                       is.old = TRUE, regionlabel) {
    Kw_combine <- Kw + Kwf
    if(is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            
            p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw_combine))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
        return(list(value = sumlogber - B))
    } else {
        # print(paste("g=", g))
        # print(paste("dim K=", dim(K)))
        # print(paste("wg_new=", wg_new))
        # print(paste("w[g]=", w[g]))
        # print(paste("length Kw=", length(Kw)))
        Kw <- Kw + K[, g] * (wg_new - w[g])
        Kw_combine <- Kw + Kwf
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw_combine))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        w[g] <- wg_new
        B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
        return(list(value = sumlogber - B, Kw = Kw))
    }
} 
logcond_wf <- function(wf, wfg_new, ga, a.tauf, b.tauf, Kf, Kw, Kwf, g, ag.vec, 
                       is.old = TRUE, regionlabel, region.idx) {
    Kw_combine <- Kw + Kwf
    if(is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            
            p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw_combine))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        B <- (length(wf) / 2 + a.tauf) * log((0.5) * sum(wf ^ 2) + b.tauf)
        return(list(value = sumlogber - B))
    } else {
        #         b = Kf[, g] * (wfg_new - wf[g])
        #         bb = Kwf[region.idx, ]
        #         print(paste("b len = ", length(b)))
        #         print(paste("bb len = ", length(bb)))
        Kwf[region.idx, ] <- Kwf[region.idx, ] + Kf[, g] * (wfg_new - wf[g])
        # print(paste("Kwf[region.idx, ] len = ", length(Kwf[region.idx, ])))
        # Kwf[region.idx, ]
        # Kwf[region.idx, ] <- Kwf[region.idx, ] + b
        Kw_combine <- Kw + Kwf
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw_combine))
        }
        # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        # G <- dim(w)[1]
        wf[g] <- wfg_new
        B <- (length(wf) / 2 + a.tauf) * log((0.5) * sum(wf ^ 2) + b.tauf)
        return(list(value = sumlogber - B, Kwf = Kwf))
    }
} 

logcond_psic <- function(psi, w, Kw, Kwf, a.phi, b.phi, dist.mat, nu, ga, ag.vec, 
                         is.old = TRUE, regionlabel) {
    Kw_combine <- Kw + Kwf
    phi <- exp(psi)
    if (is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw_combine))
        }
        logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
        return(list(value = logsumber + loggamma + psi))
    } else {
        K <- bezier2(dist.mat = dist.mat, nu = nu, phi = phi) 
        Kw <- K %*% w
        Kw_combine <- Kw + Kwf
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw_combine))
        }
        logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
        return(list(value = logsumber + loggamma + psi, K = K, Kw = Kw))
    }
}

logcond_psif <- function(psif, wf, Kw, Kwf, a.phif, b.phif, dist.matf, nu, ga, ag.vec, 
                         is.old = TRUE, regionlabel, region.idx) {
    Kw_combine <- Kw + Kwf
    phif <- exp(psif)
    if (is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw_combine))
        }
        logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        loggamma <- dgamma(phif, a.phif, b.phif, log = TRUE)
        return(list(value = logsumber + loggamma + psif))
    } else {
        Kf <- bezier2(dist.mat = dist.matf, nu = nu, phi = phif) 
        Kwf <- matrix(0, length(Kw), 1)
        Kwf[region.idx, ] <- Kf %*% wf
        Kw_combine <- Kw + Kwf
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw_combine))
        }
        logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        loggamma <- dgamma(phif, a.phif, b.phif, log = TRUE)
        return(list(value = logsumber + loggamma + psif, Kf = Kf, Kwf = Kwf))
    }
}


# Functions including Gaussian kernels
################################################################################
# library(kernlab)
GaussKer <- function(dist.mat, phi = 1) {
    # 2-d Gaussian Kernel (radial basis function, RBF)
    # defined as exp(-||x - x'|| ^ 2 / (2 * sigma ^ 2) )
    # phi = 1 / sigma ^ 2, i.e., precision
    # coord: resolution of data (voxel observation site)
    # grid: resolution of process convolution (region site)
    ker <- exp(- dist.mat ^ 2 / (2 * phi))
    return(ker)
}


logcond_psi_type <- function(psi, w, Kw, a.phi, b.phi, dist.mat, nu, ga, ag.vec, 
                             is.old = TRUE, regionlabel, kerneltype) {
    phi <- exp(psi)
    if (is.old) {
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
        return(list(value = logsumber + loggamma + psi))
    } else {
        if (kerneltype == "bezier") {
            K <- bezier2(dist.mat = dist.mat, nu = nu, phi = phi) 
        } else {
            K <- GaussKer(dist.mat = dist.mat, phi = phi) 
        }
        Kw <- K %*% w
        if(sum(ag.vec) != 0) {
            a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
            p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
        } else {
            p_bern <- 1 / (1 + exp(-Kw))
        }
        logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
        loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
        return(list(value = logsumber + loggamma + psi, K = K, Kw = Kw))
    }
}



mcmc_pc1_type <- function(Yvec, Xr, start.lst, n.mcmc, burn, thin, 
                          name.par, adapt = TRUE, tune.lst, keep.lst, 
                          tune.psi, keep.psi,
                          tune.len = 30, regionlabel, region_mat, nu = 2, 
                          a.tau = 1, b.tau = 1,
                          a.phi = 1, b.phi = 1, 
                          ag.vec = rep(0, max(regionlabel)),
                          target.accept.rate, is.cplx, kerneltype = "bezier") {
    # Yvec: N*N by n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    # kerneltype: bezier, gaussian
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #-------------------------------
    # for w
    keep <- keep.lst
    keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
    
    # for psi
    # keep.psi <- list(psi = 0)
    keep.tmp.psi <- keep.psi
    Tb <- tune.len  # frequency of adaptive tuning
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N.voxels <- nrow(Yvec)
    p <- 1
    G <- max(regionlabel)
    XtX <- crossprod(Xr)
    XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    n <- ncol(Yvec)
    MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- start.lst$w
    ga <- start.lst$ga
    phi <- start.lst$phi
    # z <- log(w / (1 - w))
    psi <- log(phi)
    
    # create kernel
    if (kerneltype == "bezier") {
        Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    } else if (kerneltype == "gaussian") {
        Ker <- GaussKer(dist.mat, phi = phi)
    }
    # Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
    Kw <- Ker %*% w
    # tau2.vec <- rep(1, G)
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) { 
        if (i %% 1000 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        if(kerneltype == "bezier") {
            if(adapt == TRUE & i %% Tb == 0) {  
                # Adaptive tuning
                keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
                tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
                keep.tmp.psi$psi <- 0
            }
        }
        
        # Update parameters
        #----------------------------
        
        # ========== sample w =============
        for (g in 1:G) {
            
            # Update tuning parameter
            #-------------------------
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp[[g]] <- keep.tmp[[g]] / Tb
                tune.lst[[g]] <- get.tune(tune.lst[[g]], keep.tmp[[g]], i)
                keep.tmp[[g]] <- 0
            }
            
            #----------------------------
            # Random walk Normal proposal 
            #----------------------------
            
            # use random walk proposal newS = S + N(0, c) to update S
            wg_new <- rnorm(1, w[g], tune.lst[[g]])
            # w_new <- w
            # w_new[g] <- wg_new 
            #             logcon.new <- logcond_latent(w_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
            #                                          ag.vec, is.old = FALSE, 
            #                                          regionlabel = regionlabel)
            #             logcon.old <- logcond_latent(w, ga, a.tau, b.tau, K = Ker, Kw, g, 
            #                                          ag.vec, is.old = TRUE, 
            #                                          regionlabel = regionlabel)
            
            logcon.new <- logcond_latent(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
                                         ag.vec, is.old = FALSE, 
                                         regionlabel = regionlabel)
            logcon.old <- logcond_latent(w, wg_new, ga, a.tau, b.tau, K = Ker, Kw, g, 
                                         ag.vec, is.old = TRUE, 
                                         regionlabel = regionlabel)
            
            
            if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
                w[g] <- wg_new
                Kw <- logcon.new$Kw
                keep[[g]] <- keep[[g]] + 1
                keep.tmp[[g]] <- keep.tmp[[g]] + 1
            }
            
            #             wg_new <- rnorm(1, w[g], tune.lst[[g]])
            #             w_new <- w
            #             w_new[g] <- wg_new 
            #             logcon.new <- update_wg(w_new, ga, tau2.vec, K = Ker, Kw, g, 
            #                                          ag.vec, is.old = FALSE)
            #             logcon.old <- update_wg(w, ga, tau2.vec, K = Ker, Kw, g, 
            #                                          ag.vec, is.old = TRUE)
            #             if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            #                 w[g] <- wg_new
            #                 Kw <- logcon.new$Kw
            #                 keep[[g]] <- keep[[g]] + 1
            #                 keep.tmp[[g]] <- keep.tmp[[g]] + 1
            #             }
        }
        
        # ============sample tau2 =============
        
        # tau2.vec <- update_tau2(a.tau, b.tau, w)
        
        # ============ sample ga ==============
        for (v in 1:N.voxels) {
            MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
            ga[v] <- ga_update_pc1(v, MM$M0, MM$M1, Kw, n, regionlabel, ag.vec, 
                                   is.cplx)
        }
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        logcon.new <- logcond_psi_type(new_psi, w, Kw, a.phi, b.phi, 
                                       dist.mat = dist.mat, nu = nu, ga = ga, 
                                       ag.vec, is.old = FALSE, 
                                       regionlabel = regionlabel, kerneltype)
        logcon.old <- logcond_psi_type(psi, w, Kw, a.phi, b.phi,
                                       dist.mat = dist.mat, nu = nu, ga = ga,
                                       ag.vec, is.old = TRUE, 
                                       regionlabel = regionlabel, kerneltype)
        #         if(kerneltype == "bezier") {
        #             new_psi <- rnorm(1, psi, tune.psi$psi)
        #             logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
        #                                       dist.mat = dist.mat, nu = nu, ga = ga, 
        #                                       ag.vec, is.old = FALSE, 
        #                                       regionlabel = regionlabel)
        #             logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
        #                                       dist.mat = dist.mat, nu = nu, ga = ga,
        #                                       ag.vec, is.old = TRUE, 
        #                                       regionlabel = regionlabel)
        #         } else if (kerneltype == "gaussian") {
        #             new_psi <- rnorm(1, psi, tune.psi$psi)
        #             logcon.new.gaussian <- logcond_psi_gaussian(new_psi, w, Kw, a.phi, b.phi, 
        #                                       dist.mat = dist.mat, nu = nu, ga = ga, 
        #                                       ag.vec, is.old = FALSE, 
        #                                       regionlabel = regionlabel)
        #             logcon.old.gaussian <- logcond_psi(psi, w, Kw, a.phi, b.phi,
        #                                       dist.mat = dist.mat, nu = nu, ga = ga,
        #                                       ag.vec, is.old = TRUE, 
        #                                       regionlabel = regionlabel)
        #         }
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        } 
        
        
        #         logcon.new <- logcond_psi(new_psi, w, Kw, a.phi, b.phi, 
        #                                   dist.mat = dist.mat, nu = nu, ga = ga, 
        #                                   ag.vec, is.old = FALSE, regionlabel)
        #         logcon.old <- logcond_psi(psi, w, Kw, a.phi, b.phi,
        #                                   dist.mat = dist.mat, nu = nu, ga = ga,
        #                                   ag.vec, is.old = TRUE, regionlabel)
        
        
        
        
        # cat("r=", r, "\n")
        
        # store results
        # W[i, , 1] <- w
        # Gamma[i, , 1] <- ga
        # Phi[i, 1] <- phi
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                draws[(i - burn) %/% thin, ] <- c(ga, w, phi)
            } 
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi$psi <- keep.psi$psi / n.mcmc
    
    # Write output
    #--------------
    return(list(draws = draws, accept = c(keep, keep.psi), start = start.lst, 
                tune = c(tune.lst, tune.psi), burn = burn, thin = thin, 
                n.mcmc = n.mcmc, sampleidx = sampleidx))
    
    # return(list(W = W, Gamma = Gamma, Phi = Phi, tune = tune, tune.psi = tune.psi,
    #             accept = keep, accept_psi = keep.psi))
    
}



