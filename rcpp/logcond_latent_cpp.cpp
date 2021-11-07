//logcond_latent_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_latent_cpp(arma::vec w, double wg_new, arma::vec ga, 
                        double atau, double btau, arma::mat K, arma::mat Kw, 
                        int gidx, arma::vec agvec, 
                        bool isold, arma::vec regionlabel) {
    arma::vec Kw_vec = vec(Kw);
    int G = w.n_elem;
    int N = Kw.n_rows;
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double sumlogber;
    arma::vec logber(N);

    if(isold) {
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        sumlogber = sum(logber);
        double B = (G / 2 + atau) * log((0.5) * sum(square(w)) + btau);
        return List::create(
            Rcpp::Named("value") = sumlogber - B
        );
    } else {
        Kw_vec = Kw_vec + K.col(gidx) * (wg_new - w(gidx));
        Kw = mat(Kw_vec);
//         w[gidx] = wg_new;
//         Kw = K * w;
        // Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        sumlogber = sum(logber);
        w[gidx] = wg_new;
        double B = (G / 2 + atau) * log((0.5) * sum(square(w)) + btau);
        return List::create(
            Rcpp::Named("value") = sumlogber - B,
            Rcpp::Named("Kw") = Kw
        );
    }
}



