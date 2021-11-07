//logcond_psi_cpp.cpp

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
    arma::mat K = zeros<mat>(N, M);
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
List logcond_psi_cpp(double psi, arma::vec w, arma::mat Kw, 
                       double aphi, double bphi, arma::mat distmat, 
                       double nu, arma::vec ga, arma::vec agvec, 
                       bool isold, arma::vec regionlabel) {
    
    arma::vec Kw_vec = vectorise(Kw);
    int N = Kw.n_rows;
    double phi = exp(psi);
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double logsumber;
    double loggamma;
    arma::vec logber(N);
    arma::mat K;
    
    
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
        loggamma = R::dgamma(phi, aphi, 1/bphi, 1);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi
        );
    } else {
        K = bezier2_cpp(distmat, nu, phi);
        Kw = K * w;
        Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; ++i) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
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
        loggamma = R::dgamma(phi, aphi, 1/bphi, 1);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi,
            Rcpp::Named("K")     = K,
            Rcpp::Named("Kw")    = Kw
        );
    }
}
