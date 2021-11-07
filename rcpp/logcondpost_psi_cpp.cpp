//logcondpost_psi_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat bezier2_cpp(arma::mat distmat, double nu, double phi) {
    int N = distmat.n_rows;
    int M = distmat.n_cols;
    arma::mat K = zeros(N, M);
    for (int i=0; i < N; i++) {
        for (int j=0; j < M; j++) {
            if (distmat(i, j) < phi) {
                K(i, j) = pow(1 - pow(distmat(i, j) / phi, 2), nu);
            }
        }
    }
    return (K);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma:: vec rowSumsC(arma::mat x) {
    int nrow = x.n_rows, ncol = x.n_cols;
    arma::vec out(nrow);
    for (int i = 0; i < nrow; i++) {
        double total = 0;
        for (int j = 0; j < ncol; j++) {
            total += x(i, j);
        }
        out[i] = total;
    }
    return out;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat normalized_bezier2_cpp(arma::mat distmat, double nu, double phi) {
    arma::mat K = bezier2_cpp(distmat, nu, phi);
    arma::vec rowsum = rowSumsC(K);
    int nrow = K.n_rows;
    int ncol = K.n_cols;
    for (int i = 0; i < nrow; i++) {
        if (rowsum[i] != 0) {
            K(i, span(0, ncol - 1)) = K(i, span(0, ncol - 1)) / rowsum[i];
        }
    }
    return (K);
}


// logcondpost_psi <- function(psi, w, Kw, a.phi, b.phi, dist.mat, nu, ga, is.old = TRUE) {
//     phi <- exp(psi)
//     loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
//     if (is.old) {
//         logsumber <- sum(dbinom(ga, size = 1, prob = Kw, log = TRUE))
//         return(list(value = logsumber + loggamma + psi))
//     } else {
//         K <- normalized_bezier2(dist.mat = dist.mat, nu = nu, phi = phi) 
//         Kw <- K %*% w
//         logsumber <- sum(dbinom(ga, size = 1, prob = Kw, log = TRUE))
//         return(list(value = logsumber + loggamma + psi, K = K, Kw = Kw))
//     }
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcondpost_psi_cpp(double psi, arma::vec w, arma::mat Kw, 
                     double aphi, double bphi, arma::mat distmat, 
                     double nu, arma::vec ga, bool isold) {
    
    arma::vec Kw_vec = vectorise(Kw);
    int N = Kw.n_rows;
    double phi = exp(psi);
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double logsumber;
    double loggamma = R::dgamma(phi, aphi, 1/bphi, 1);;
    arma::vec logber(N);
    arma::mat K;

    // loggamma = R::dgamma(phi, aphi, 1/bphi, 1);
    
    if (isold) { 
        for (int i = 0; i < N; ++i) {
            logber[i] = R::dbinom(ga[i], 1, Kw_vec[i], 1);
        }
        logsumber = sum(logber);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi
        );
    } else {
        K = normalized_bezier2_cpp(distmat = distmat, nu = nu, phi = phi);
        Kw = K * w;
        Kw_vec = vec(Kw);
        for (int i = 0; i < N; ++i) {
            if (Kw_vec[i] >= 1) Kw_vec[i] = 1 - 1e-8;
            if (Kw_vec[i] <= 0) Kw_vec[i] = 1e-8;
            logber[i] = R::dbinom(ga[i], 1, Kw_vec[i], 1);
        }
        logsumber = sum(logber);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi,
            Rcpp::Named("K")     = K,
            Rcpp::Named("Kw")    = Kw
        );
    }
}
