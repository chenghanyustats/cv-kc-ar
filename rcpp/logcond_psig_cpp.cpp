//logcond_psig_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat bezier2_phi_cpp(arma::mat distmat, double nu, arma::vec phivec) {
    int N = distmat.n_rows;
    int M = distmat.n_cols;
    arma::mat K = zeros(N, M);
    for (int i=0; i < N; i++) {
        for (int j=0; j < M; j++) {
            if (distmat(i, j) < phivec[j]) {
                K(i, j) = pow(1 - pow(distmat(i, j) / phivec[j], 2), nu);
            }
        }
    }
    return (K);
}


// # log conditional posterior for psi_g
// logcond_psig <- function(psi, w, Kw, a.phi, b.phi, dist.mat, nu, ga, ag.vec, 
//                          is.old = TRUE, regionlabel, g) {
//     phi <- exp(psi)
//     loggamma <- dgamma(phi[g], a.phi, b.phi, log = TRUE)
//     if (is.old) {
//         if(sum(ag.vec) != 0) {
//             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
//             p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         } else {
//             p_bern <- 1 / (1 + exp(-Kw))
//         }
//         logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
//             return(list(value = logsumber + loggamma + psi[g]))
//     } else {
//         K <- bezier2_phi(dist.mat = dist.mat, nu = nu, phi.vec = phi) 
//         Kw <- K %*% w
//         if(sum(ag.vec) != 0) {
//             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
//             p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         } else {
//             p_bern <- 1 / (1 + exp(-Kw))
//         }
//         logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
//             return(list(value = logsumber + loggamma + psi[g], K = K, Kw = Kw))
//     }
// }


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_psig_cpp(arma::vec psi, arma::vec w, arma::mat Kw, 
                     double aphi, double bphi, arma::mat distmat, 
                     double nu, arma::vec ga, arma::vec agvec, 
                     bool isold, arma::vec regionlabel, int gidx) {
    
    arma::vec Kw_vec = vectorise(Kw);
    int N = Kw.n_rows;
    arma::vec phi = exp(psi);
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double logsumber;
    double loggamma;
    arma::vec logber(N);
    arma::mat K;
    
    loggamma = R::dgamma(phi[gidx], aphi, 1/bphi, 1);

    if (isold) { 
        if(sum(agvec) != 0) {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi[gidx]
        );
    } else {
        K = bezier2_phi_cpp(distmat, nu, phi);
        Kw = K * w;
        Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi[gidx],
            Rcpp::Named("K")     = K,
            Rcpp::Named("Kw")    = Kw
        );
    }
}







