################################################################################
# Self-created cpp functions for MCMC algorithms in the CV-KC-AR paper         #
# ./mcmc_cpp_fcns.R                                                            #
# Cheng-Han Yu, Marquette University                                           #
################################################################################

library(compiler)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(microbenchmark)

# for KC methods
sourceCpp("./rcpp/dist_map_cpp.cpp")
sourceCpp("./rcpp/bezier2_cpp.cpp")

# CV-KC
sourceCpp("./rcpp/ga_update_pc1_cpp.cpp")
sourceCpp("./rcpp/logcond_psi_cpp.cpp")
sourceCpp("./rcpp/logcond_latent_cpp.cpp")
sourceCpp("./rcpp/mcmc_pc1_main_cpp.cpp")
sourceCpp("./rcpp/mcmc_kcs_main_cpp.cpp")

# CV-GP
sourceCpp("./rcpp/MHratio_cond_S_cpp.cpp")
sourceCpp("./rcpp/ga_update_gp_cpp.cpp")
sourceCpp("./rcpp/logcond_w_cpp.cpp")
sourceCpp("./rcpp/mcmc_gp_main_cpp.cpp")


# # CV-KC-phi
# sourceCpp("./rcpp/mcmc_pc1_phi_main_cpp.cpp")


# CV-KC-multiresolution
sourceCpp("./rcpp/mcmc_pc1_multires_cpp.cpp")

# CV-KC-multisubjects
sourceCpp("./rcpp/mcmc_multisubject_main_cpp.cpp")

###########
# CV-GP-AR
###########


mcmc_gp_cpp_full <- function(Yvec, Xr, start.lst, nmcmc, burn, thin, 
                             name.par, adapt = TRUE, tunelst, keeplst, 
                             tuner, keepr,
                             tunelen = 50, regionlabel, region_mat, df = 8, 
                             agvec = rep(0, 25),
                             target.accept.rate, iscplx) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    centroid <- compute_centroid(region_mat)
    
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = nmcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    
    #-------------------------------
    # Starting values
    #-------------------------------
    S <- start.lst$S  # S is a COLUMN vector
    ga <- start.lst$ga
    r <- start.lst$r
    w <- log(r)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    mcmc_gp_main_cpp(Yvec, Xr, S, ga, w, r, nmcmc, burn, thin, 
                     adapt, tunelst, keeplst, tuner, keepr, tunelen, 
                     regionlabel, df, agvec, centroid, iscplx, 
                     get_tune, power_exp_cov_fcn, update_ga, draws)
}


mcmc_gp_ar_cpp_full <- function(Yvec, Xr, arcoefvec, start.lst, nmcmc, burn, thin, 
                                name.par, adapt = TRUE, tunelst, keeplst, 
                                tuner, keepr,
                                tunelen = 50, regionlabel, region_mat, df = 8, 
                                agvec = rep(0, 25),
                                target.accept.rate, iscplx) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    centroid <- compute_centroid(region_mat)
    
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = nmcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    
    #-------------------------------
    # Starting values
    #-------------------------------
    S <- start.lst$S  # S is a COLUMN vector
    ga <- start.lst$ga
    r <- start.lst$r
    w <- log(r)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    mcmc_gp_ar_main_cpp(Yvec, Xr, arcoefvec, S, ga, w, r, nmcmc, burn, thin, 
                        adapt, tunelst, keeplst, tuner, keepr, tunelen, 
                        regionlabel, df, agvec, centroid, iscplx, 
                        get_tune, power_exp_cov_fcn, update_ga, draws)
}

###########
# CV-KC-AR
###########
mcmc_kcs_cpp_full <- function(Yvec, Xr, start.lst, nmcmc, burn, thin, 
                              name.par, adapt = TRUE, tunelst, keeplst, 
                              tunepsi, keeppsi,
                              tunelen = 50, regionlabel, region_mat, nu = 2, 
                              atau = 1/2, btau = 1/2,
                              aphi = 1/2, bphi = 1/2, 
                              agvec = rep(0, 25),
                              target.accept.rate, iscplx, kerneltype) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N <- sqrt(nrow(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = nmcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    mcmc_kcs_main_cpp(Yvec, Xr, w, ga, phi, psi, nmcmc, burn, thin, 
                      adapt, tunelst, keeplst, tunepsi, keeppsi, tunelen, 
                      regionlabel, nu, atau, btau, aphi, bphi, 
                      agvec, distmat, iscplx, get_tune, draws, kerneltype)
}


mcmc_kcs_ar_cpp_full <- function(Yvec, Xr, arcoefvec, start.lst, nmcmc, burn, thin, 
                                 name.par, adapt = TRUE, tunelst, keeplst, 
                                 tunepsi, keeppsi,
                                 tunelen = 50, regionlabel, region_mat, nu = 2, 
                                 atau = 1/2, btau = 1/2,
                                 aphi = 1/2, bphi = 1/2, 
                                 agvec = rep(0, 25),
                                 target.accept.rate, iscplx, kerneltype) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N <- sqrt(nrow(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = nmcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    mcmc_kcs_ar_main_cpp(Yvec, Xr, arcoefvec, w, ga, phi, psi, nmcmc, burn, thin, 
                         adapt, tunelst, keeplst, tunepsi, keeppsi, tunelen, 
                         regionlabel, nu, atau, btau, aphi, bphi, 
                         agvec, distmat, iscplx, get_tune, draws, kerneltype)
}



mcmc_kcs_block_cpp_full <- function(Yvec, Xr, start.lst, nmcmc, burn, thin, 
                                    name.par, adapt = TRUE, tunelst, keeplst, 
                                    tunepsi, keeppsi,
                                    tunelen = 50, regionlabel, region_mat, nu = 2, 
                                    atau = 1/2, btau = 1/2,
                                    aphi = 1/2, bphi = 1/2, 
                                    agvec = rep(0, 25),
                                    target.accept.rate, iscplx, kerneltype) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N <- sqrt(nrow(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = nmcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    mcmc_kcs_block_main_cpp(Yvec, Xr, w, ga, phi, psi, nmcmc, burn, thin, 
                            adapt, tunelst, keeplst, tunepsi, keeppsi, tunelen, 
                            regionlabel, nu, atau, btau, aphi, bphi, 
                            agvec, distmat, iscplx, get_tune, draws, kerneltype)
}

mcmc_kcs_ar_block_cpp_full <- function(Yvec, Xr, arcoefvec, start.lst, nmcmc, burn, thin, 
                                       name.par, adapt = TRUE, tunelst, keeplst, 
                                       tunepsi, keeppsi,
                                       tunelen = 50, regionlabel, region_mat, nu = 2, 
                                       atau = 1/2, btau = 1/2,
                                       aphi = 1/2, bphi = 1/2, 
                                       agvec = rep(0, 25),
                                       target.accept.rate, iscplx, kerneltype) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    N <- sqrt(nrow(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = nmcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    mcmc_kcs_ar_block_main_cpp(Yvec, Xr, arcoefvec, w, ga, phi, psi, nmcmc, burn, thin, 
                               adapt, tunelst, keeplst, tunepsi, keeppsi, tunelen, 
                               regionlabel, nu, atau, btau, aphi, bphi, 
                               agvec, distmat, iscplx, get_tune, draws, kerneltype)
}




################################################################################

###########
# CV-KC-AR-phi
###########

# full cpp function mcmc_pc1_phi_main_cpp()
mcmc_pc1_phi_cpp_full <- function(Yvec, Xr, start.lst, nmcmc, burn, thin, 
                                  name.par, adapt = TRUE, tunelst, keeplst, 
                                  tunepsilst, keeppsilst,
                                  tunelen = 50, regionlabel, region_mat, nu = 2, 
                                  atau = 1, btau = 1,
                                  aphi = 1, bphi = 1, 
                                  agvec = rep(0, 25),
                                  target.accept.rate, iscplx) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    mcmc_pc1_phi_main_cpp(Yvec, Xr, w, ga, phi, psi, nmcmc, burn, thin, 
                          adapt, tunelst, keeplst, tunepsilst, keeppsilst, tunelen, 
                          regionlabel, nu, atau, btau, aphi, bphi, 
                          agvec, distmat, iscplx, get_tune, draws)
}


# full cpp function mcmc_pc1_phi_block_main_cpp()
mcmc_pc1_phi_block_cpp_full <- function(Yvec, Xr, start.lst, nmcmc, burn, thin, 
                                        name.par, adapt = TRUE, tunelst, keeplst, 
                                        tunepsilst, keeppsilst,
                                        tunepsilst_block, keeppsilst_block,
                                        tunelen = 50, regionlabel, region_mat, nu = 2, 
                                        atau = 1/2, btau = 1/2,
                                        aphi = 1/2, bphi = 1/2, 
                                        agvec = rep(0, 25),
                                        target.accept.rate, iscplx) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    # List mcmc_pc1_phi_block_main_cpp(arma::mat Yvec, arma::mat Xr, 
    #                                  arma::vec w, arma::vec ga, arma::vec phi, arma::vec psi,
    #                                  int nmcmc, int burn, int thin, 
    #                                  bool adapt, List tunelst, List keeplst, 
    #                                  List tunepsilst, List keeppsilst, int tunelen, 
    #                                  arma::vec regionlabel, int nu, 
    #                                  double atau, double btau,
    #                                  double aphi, double bphi, 
    #                                  arma::vec agvec, arma::mat distmat, bool iscplx,
    #                                  Function get_tune, arma::mat draws)
    mcmc_pc1_phi_block_main_cpp(Yvec, Xr, w, ga, phi, psi, nmcmc, burn, thin, 
                                adapt, tunelst, keeplst, tunepsilst, keeppsilst, 
                                tunepsilst_block, keeppsilst_block, tunelen, 
                                regionlabel, nu, atau, btau, aphi, bphi, 
                                agvec, distmat, iscplx, get_tune, draws)
}

mcmc_pc1_phi_block_cpp_full_vec <- function(Yvec, Xr, start.lst, nmcmc, burn, thin, 
                                            name.par, adapt = TRUE, tunevec, keepvec, 
                                            tunepsivec, keeppsivec,
                                            tunepsilst_block, keeppsilst_block,
                                            tunelen = 50, regionlabel, region_mat, nu = 2, 
                                            atau = 1/2, btau = 1/2,
                                            aphi = 1/2, bphi = 1/2, 
                                            agvec = rep(0, 25),
                                            target.accept.rate, iscplx) {
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    # List mcmc_pc1_phi_block_main_cpp(arma::mat Yvec, arma::mat Xr, 
    #                                  arma::vec w, arma::vec ga, arma::vec phi, arma::vec psi,
    #                                  int nmcmc, int burn, int thin, 
    #                                  bool adapt, List tunelst, List keeplst, 
    #                                  List tunepsilst, List keeppsilst, int tunelen, 
    #                                  arma::vec regionlabel, int nu, 
    #                                  double atau, double btau,
    #                                  double aphi, double bphi, 
    #                                  arma::vec agvec, arma::mat distmat, bool iscplx,
    #                                  Function get_tune, arma::mat draws)
    mcmc_pc1_phi_block_main_cpp_vec(Yvec, Xr, w, ga, phi, psi, nmcmc, burn, thin, 
                                    adapt, tunevec, keepvec, tunepsivec, keeppsivec, 
                                    tunepsilst_block, keeppsilst_block, tunelen, 
                                    regionlabel, nu, atau, btau, aphi, bphi, 
                                    agvec, distmat, iscplx, get_tune, draws)
}

################################################################################
# CV-KC Multi-resolution 
################################################################################
# full cpp function mcmc_pc1_main_cpp()

mcmc_pc1_multires_cpp_full <- function(Yvec, Xr, start.lst, nmcmc, burn, thin, 
                                       name.par, adapt = TRUE, tunelst, keeplst, 
                                       tunelstf, keeplstf, 
                                       tunepsi, keeppsi,
                                       tunelen = 50, regionlabel, region_mat, 
                                       regionlabelfiner, region_mat_finer,
                                       nu = 2, 
                                       atau = 1/2, btau = 1/2,
                                       aphi = 1/2, bphi = 1/2, 
                                       atauf = 1/2, btauf = 1/2,
                                       aphif = 1/2, bphif = 1/2, 
                                       agvec = rep(0, max(regionlabel)),
                                       target.accept.rate, iscplx) {
    #################################################################
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
    get_tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        a <- min(1, 10000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Here create some objects that will be used later.
    #-------------------------------
    # N <- sqrt(nrow(Yvec))
    # centroid <- compute_centroid(region_mat)
    # Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    # distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    
    # N.voxels <- nrow(Yvec)
    # p <- 1
    # G <- max(regionlabel)
    # XtX <- crossprod(Xr)
    # XtY <- crossprod(Xr, t(Yvec))
    # centroid <- compute_centroid(region_mat)
    # Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    # dist.mat <- dist_mat(coord = Coord, grid = centroid)
    # n <- ncol(Yvec)
    # MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    
    
    N.voxels <- nrow(Yvec)
    N <- sqrt(N.voxels)
    # p <- 1
    # G <- max(regionlabel)
    # regionlabelf <- relabel(regionlabel.finer)  # define relabel function!
    # Gf <- max(regionlabelf)
    regionidx <- which(regionlabelfiner != 0)
    # regionidx <- regionidx - 1
    Gf <<- length(unique(regionlabel.finer[region.idx]))
    # XtX <- crossprod(Xr)
    # XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    distmat <- dist_mat_cpp(coord = Coord, grid = centroid)
    centroidf <- compute_centroid(region_mat_finer)
    centroidf2 <- centroidf[unique(regionlabelfiner[regionidx]), ]
    Coord2 <- Coord[regionidx, ]
    distmatf <- dist_mat_cpp(coord = Coord2, grid = centroidf2)
    
    # n <- ncol(Yvec)
    # MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    lensam <- length(sampleidx)
    draws <- matrix(NA, nrow = lensam, ncol = length(name.par))
    colnames(draws) <- name.par
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga <- start.lst$ga
    phi <- start.lst$phi
    psi <- log(phi)
    wf <- as.vector(start.lst$wf)
    phif <- start.lst$phif
    psif <- log(phif)
    
    #-------------------------------
    # Main Rcpp mcmc function
    #-------------------------------
    mcmc_pc1_multires_cpp_main(Yvec, Xr, w, ga, phi, psi, wf, phif, psif,
                               nmcmc, burn, thin,
                               adapt, tunelst, keeplst, 
                               tunelstf, keeplstf,
                               tunepsi, keeppsi, 
                               tunelen, 
                               regionlabel, regionlabelfiner, regionidx,
                               nu, atau, btau, aphi, bphi, 
                               atauf, btauf,
                               aphif, bphif,
                               agvec, distmat, distmatf,
                               iscplx, get_tune, draws)
}

################################################################################
# CV-KC Multi-subject
################################################################################
# full cpp function mcmc_multisubject_main_cpp()
mcmc_pc1_multi_prec_cpp_full <- function(Yvec_multi, Xr, start.lst, nmcmc, burn, thin,
                                         name.par, adapt = TRUE, tunew, keepw,
                                         tune.psi, keep.psi, tune.a, keep.a,
                                         tune.len = 30, regionlabel, region_mat,
                                         region_conn, nu = 2,
                                         a.tau = 1, b.tau = 1,
                                         a.phi = 1, b.phi = 1,
                                         target.accept.rate, is.cplx) {
    #################################################################
    # Yvec_multi: ns x N*N x n COMPLEX-valued data
    # Xr: n by p design matrix (p = 1)
    # init_ga: initial indicator variable
    # init_S: initial spatial effect (ROW factor)
    # init_r: initial range parameter
    # m: number of iterations
    # regionlabel: N*N by 1 region vector
    # region_mat: N by N region matrix
    # region_conn: N*N by 1 region vector for connectivity
    #################################################################
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    
    if(missing(region_conn)) {
        region_conn <- regionlabel
    }
    
    get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
        aa <- min(0.5, 1000 / sqrt(k))
        # a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - aa, log(tune) + aa))
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
    ns <- dim(Yvec_multi)[1]
    N.voxels <- dim(Yvec_multi)[2]
    # print(paste("N.voxels", N.voxels))
    n <- dim(Yvec_multi)[3]
    N <- sqrt(N.voxels)
    p <- 1
    G <- max(regionlabel)
    D <- max(region_conn)
    XtX <- crossprod(Xr)
    # XtY <- crossprod(Xr, t(Yvec))
    centroid <- compute_centroid(region_mat)
    Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
    dist.mat <- dist_mat(coord = Coord, grid = centroid)
    
    # MM_mat <- marginal_yga_multi(Yvec_multi, XtX, Xr)
    # MM_mat <- marginal_yga_mat(Yvec, XtY, XtX)
    Yvec_multi_cpp <- aperm(Yvec_multi, c(2, 3, 1))
    MM_mat <- marginal_yga_multi_cpp(Yvec_multi_cpp, XtX, Xr)
    M0_mat <- MM_mat$M0_mat
    M1_mat <- MM_mat$M1_mat
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    W <- matrix(NA, nrow = n.mcmc, ncol = G)
    A <- array(NA, dim = c(n.mcmc, ns, D))
    Ga_mat <- array(NA, dim = c(n.mcmc, ns, N.voxels))
    Prec_array <- array(0, dim = c(n.mcmc, D, D))
    
    #-------------------------------
    # Starting values
    #-------------------------------
    w <- as.vector(start.lst$w)
    ga_mat <- start.lst$ga_mat
    phi <- start.lst$phi
    psi <- log(phi)
    a_mat <- start.lst$a_mat
    Prec <- start.lst$Prec
    
    # create kernel
    Ker <- bezier2_cpp(dist.mat, nu, phi)
    # Ker <- bezier2(dist.mat, nu, phi)
    Kw <- Ker %*% w
    
    #-------------------------------
    # M-H algorithm
    #-------------------------------
    for (i in 1:n.mcmc) {
        if (i %% 100 == 0) cat("iter:", i, "\r")
        
        # Update tuning parameter of psi
        #------------------------
        if(adapt == TRUE & i %% Tb == 0) {
            # Adaptive tuning
            keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
            tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
            keep.tmp.psi$psi <- 0
        }
        
        #         # Update tuning parameter of a
        #         #------------------------
        #         if(adapt == TRUE & i %% Tb == 0) {
        #             # Adaptive tuning
        #             keep.tmp.a$a <- keep.tmp.a$a / Tb
        #             tune.a$a <- get.tune(tune.a$a, keep.tmp.a$a, i)
        #             keep.tmp.a$a <- 0
        #         }
        
        # ========== sample w =============
        for (j in 1:5) {
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp.w[[j]] <- keep.tmp.w[[j]] / Tb
                tune.w[[j]] <- get.tune(tune.w[[j]], keep.tmp.w[[j]], i)
                keep.tmp.w[[j]] <- 0
            }
            w.idx <- (5 * (j - 1) + 1:5)
            # block move using empirical covariance structure
            if (i > 1000) {
                Sigma.w <- cov(W[(i - (1000)):(i - 1), w.idx])
                Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w, doSym = TRUE)$mat)
                # print(round(Sigma.w, 2))
                w_new <- as.vector(rmvn(1, w[w.idx], tune.w[[j]] * Sigma.w))
            } else {
                w_new <- as.vector(rmvn(1, w[w.idx], tune.w[[j]] * diag(5)))
            }
            # logcond_wg_block_multisub <- function(w, w_new, ga_mat, a.tau, b.tau, K, Kw,
            # w.idx, ag_vec, is.old = TRUE, regionlabel)
            # logcond_wg_block_multi_cpp(arma::vec w, arma::vec w_new, arma::mat ga_mat, 
            #                            double atau, double btau, arma::mat K, arma::mat Kw, 
            #                            arma::uvec widx, arma::mat a_mat, 
            #                            bool isold, arma::vec region_conn)
            # print("Before Check")
            logcon.new <- logcond_wg_block_multi_cpp(w, w_new, ga_mat,
                                                     atau = a.tau, btau = b.tau,
                                                     K = Ker, Kw, widx = w.idx - 1,
                                                     a_mat = a_mat, isold = FALSE,
                                                     region_conn = region_conn)
            logcon.old <- logcond_wg_block_multi_cpp(w, w_new, ga_mat,
                                                     atau = a.tau, btau = b.tau,
                                                     K = Ker, Kw, widx = w.idx - 1,
                                                     a_mat = a_mat, isold = TRUE,
                                                     region_conn = region_conn)
            # logcon.new <- logcond_wg_block_multi(w, w_new, ga_mat, a.tau, b.tau,
            #                                      K = Ker, Kw, w.idx,
            #                                      a_mat = a_mat, is.old = FALSE,
            #                                      region_conn = region_conn)
            # logcon.old <- logcond_wg_block_multi(w, w_new, ga_mat, a.tau, b.tau,
            #                                      K = Ker, Kw, w.idx,
            #                                      a_mat = a_mat, is.old = TRUE,
            #                                      region_conn = region_conn)
            if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
                w[w.idx] <- w_new
                Kw <- logcon.new$Kw
                keep[[j]] <- keep[[j]] + 1
                keep.tmp.w[[j]] <- keep.tmp.w[[j]] + 1
                # W[i, w.idx] <- w_new
            }
            W[i, w.idx] <- w[w.idx]
        }
        # ============ sample ga ==============
        for (s in 1:ns) {
            for (v in 1:N.voxels) {
                # MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
                # ga_update_multisub <- function(s, v, ga_mat, M0_mat, M1_mat, Kw, n, regionlabel,
                # ag_vec, is.cplx)
                
                # ga_update_multi <- function(s, v, ga_mat, M0_mat, M1_mat, Kw, n,
                #                             region_conn, a_mat, is.cplx)
                # ga_mat[s, v] <- ga_update_multi(s, v, ga_mat, MM_mat$M0,
                #                                 MM_mat$M1, Kw, n, region_conn,
                #                                 a_mat = a_mat, is.cplx)
                # ga_mat[s, v] <- ga_update_multi(s, v, ga_mat, M0_mat,
                #                                 M1_mat, Kw, n, region_conn,
                #                                 a_mat = a_mat, is.cplx)
                ga_mat[s, v] <- ga_update_multi_cpp(s - 1, v - 1, ga_mat,
                                                    M0_mat, M1_mat, Kw, n,
                                                    region_conn, a_mat,
                                                    iscplx = is.cplx)
            }
        }
        
        Ga_mat[i, , ] <- ga_mat
        
        # ============ sample psi and hence phi ==============
        new_psi <- rnorm(1, psi, tune.psi$psi)
        #         logcond_psi_multisub <- function(psi, w, Kw, a_phi, b_phi, dist.mat, nu,
        #                                          ga_mat, ag_vec, is.old = TRUE, regionlabel)
        # logcond_psi_multi <- function(psi, w, Kw, a_phi, b_phi, dist.mat, nu,
        #                               ga_mat, a_mat, is.old = TRUE, region_conn)
        
        # logcond_psi_multi_cpp(double psi, arma::vec w, arma::mat Kw, 
        #                       double aphi, double bphi, arma::mat distmat, 
        #                       double nu, arma::mat ga_mat, arma::mat a_mat, 
        #                       bool isold, arma::vec region_conn)
        
        logcon.new <- logcond_psi_multi_cpp(new_psi, w, Kw, aphi = a.phi,
                                            bphi = b.phi,
                                            distmat = dist.mat, nu = nu,
                                            ga_mat = ga_mat,
                                            a_mat = a_mat, isold = FALSE,
                                            region_conn = region_conn)
        logcon.old <- logcond_psi_multi_cpp(psi, w, Kw, aphi = a.phi,
                                            bphi = b.phi,
                                            distmat = dist.mat, nu = nu,
                                            ga_mat = ga_mat,
                                            a_mat = a_mat, isold = TRUE,
                                            region_conn = region_conn)
        
        
        
        # logcon.new <- logcond_psi_multi(new_psi, w, Kw, a.phi, b.phi,
        #                                 dist.mat = dist.mat, nu = nu,
        #                                 ga_mat = ga_mat,
        #                                 a_mat = a_mat, is.old = FALSE,
        #                                 region_conn = region_conn)
        # logcon.old <- logcond_psi_multi(psi, w, Kw, a.phi, b.phi,
        #                                 dist.mat = dist.mat, nu = nu,
        #                                 ga_mat = ga_mat,
        #                                 a_mat = a_mat, is.old = TRUE,
        #                                 region_conn = region_conn)
        
        if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            psi <- new_psi
            phi <- exp(psi)
            Ker <- logcon.new$K
            Kw <- logcon.new$Kw
            keep.psi$psi <- keep.psi$psi + 1
            keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
        }
        
        # ============ sample a ==============
        for (s in 1:ns) {
            # Update tuning parameter of a
            #------------------------
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp.a[[s]] <- keep.tmp.a[[s]] / Tb
                tune.a[[s]] <- get.tune(tune.a[[s]], keep.tmp.a[[s]], i)
                keep.tmp.a[[s]] <- 0
            }
            # A <- array(NA, dim = c(n.mcmc, ns, D))
            if (i > 1000) {
                Sigma.a <- cov(A[(i - (1000)):(i - 1), s, ])
                Sigma.a <- as.matrix(Matrix::nearPD(Sigma.a)$mat)
                # print(round(Sigma.w, 2))
                # a_new <- as.vector(rmvn(1, a_mat[s, ], tune.a[[s]] * Sigma.a))
                a_new <- rmvn(1, a_mat[s, ], tune.a[[s]] * Sigma.a)
                # new_avec <- as.vector(rmvn(1, avec, tune.a$a * diag(G)))
            } else {
                a_new <- rmvn(1, a_mat[s, ], tune.a[[s]] * diag(D))
            }
            # logcond_a_multi(s, a_mat, a_new, ga_mat, Kw, priorSig, region_conn,
            #                 is.old = TRUE)
            # logcond_a_multi <- function(ag_vec, ag_new, ga_mat, Kw, priorSig, g,
            #                              is.old = TRUE)
            # List logcond_a_multi_prec_cpp(int s, arma::mat a_mat, arma::mat a_new,
            #                               arma::mat ga_mat, arma::mat Kw, arma::mat Prec,
            #                               arma::vec region_conn, bool isold)
            logcond_a_old <- logcond_a_multi_prec_cpp(s-1, a_mat, a_new, ga_mat, Kw, Prec,
                                                      region_conn, isold = TRUE)
            logcond_a_new <- logcond_a_multi_prec_cpp(s-1, a_mat, a_new, ga_mat, Kw, Prec,
                                                      region_conn, isold = FALSE)
            # logcond_a_old <- logcond_a_multi_prec(s, a_mat, a_new, ga_mat, Kw, Prec,
            #                                       region_conn, is.old = TRUE)
            # logcond_a_new <- logcond_a_multi_prec(s, a_mat, a_new, ga_mat, Kw, Prec,
            #                                       region_conn, is.old = FALSE)
            if (log(runif(1)) < (logcond_a_new$value - logcond_a_old$value)) {
                a_mat[s, ] <- a_new
                keep.a[[s]] <- keep.a[[s]] + 1
                keep.tmp.a[[s]] <- keep.tmp.a[[s]] + 1
            }
            
            A[i, s, ] <- a_mat[s, ]
            # for (j in 1:5) {
            #     # if(adapt == TRUE & i %% Tb == 0) {
            #     #     # Adaptive tuning
            #     #     keep.tmp.a[[j]] <- keep.tmp.a[[j]] / Tb
            #     #     tune.a[[j]] <- get.tune(tune.a[[j]], keep.tmp.a[[j]], i)
            #     #     keep.tmp.a[[j]] <- 0
            #     # }
            #     a.idx <- (5 * (j - 1) + 1:5)
            #     if (i > 500) {
            #         Sigma.a <- cov(A[(i-Tb):(i-1), a.idx])
            #         Sigma.a <- as.matrix(Matrix::nearPD(Sigma.a)$mat)
            #         # print(round(Sigma.w, 2))
            #         new_avec <- as.vector(rmvn(1, ag_vec[a.idx], tune.a[[j]] * Sigma.a))
            #         # new_avec <- as.vector(rmvn(1, avec, tune.a$a * diag(G)))
            #     } else {
            #         new_avec <- as.vector(rmvn(1, ag_vec[a.idx], tune.a[[j]] * diag(G)))
            #     }
            #     #             logcond_ag_multisub <- function(ag_vec, ag_new, ga_mat, Kw, priorSig, g,
            #     #                                             is.old = TRUE)
            #     logcond_a_old <- logcond_ag_multisub(ag_vec, new_avec, ga_mat, Kw,
            #                                          rep(0, G), a.idx, TRUE)
            #     logcond_a_new <- logcond_ag_multisub(ag_vec, new_avec, ga_mat, Kw,
            #                                          rep(0, G), a.idx, FALSE)
            #
            #     if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
            #         ag_vec[a.idx] <- new_avec
            #         keep.a[[j]] <- keep.a[[j]] + 1
            #         keep.tmp.a[[j]] <- keep.tmp.a[[j]] + 1
            #         # A[i, a.idx] <- new_avec
            #     }
            #     A[i, a.idx] <- ag_vec[a.idx]
            # }
            
        }
        
        
        #         # ============ sample Sigma ==============
        #         # prior IW(G, diag(G))
        #
        
        # Sig <- MCMCpack::riwish(v = ns + r0, S = (t(a_mat) %*% a_mat + r0 * diag(D))^-1)
        r0 <- D
        Covar <- t(a_mat) %*% a_mat + r0 * diag(D)
        Sigma <- chol2inv(chol(Covar))
        Prec <- stats::rWishart(1, df = ns + r0, Sigma = Sigma)
        Prec <- Prec[, , 1]
        Prec_array[i, , ] <- Prec
        
        #  Save samples
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                # print(w)
                draws[(i - burn) %/% thin, ] <- c(as.vector(ga_mat), w,
                                                  as.vector(a_mat), phi)
                # (v = 1, s = 1), ..., (v = 1, s = ns), ...,
                # (v = N.voxels, s = 1), ..., (v = N.voxels, s = ns)
            }
        }
        
        
    }
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    keep.psi$psi <- keep.psi$psi / n.mcmc
    keep.a <- lapply(keep.a, function(x) x / n.mcmc)
    # keep.a$a <- keep.a$a / n.mcmc
    # keep.w$w <- keep.w$w / n.mcmc
    # keep <- lapply(keep, function(x) x / n.mcmc)
    # Write output
    #--------------
    return(list(draws = draws,
                A = A,
                Prec_array = Prec_array,
                Ga_mat = Ga_mat,
                # W = W,
                accept = c(keep, keep.psi, keep.a),
                start = start.lst,
                tune = c(tune.w, tune.psi, tune.a),
                burn = burn, thin = thin,
                n.mcmc = n.mcmc, sampleidx = sampleidx))
}











