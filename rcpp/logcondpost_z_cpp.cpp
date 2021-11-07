//logcondpost_z_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcondpost_z_cpp(double zg, double zg_new, arma::vec ga, 
                        double awvec_g, double bwvec_g, 
                        arma::vec Kg, arma::mat Kw,
                        bool isold) {
    arma::vec Kw_vec = vec(Kw);
    int N = Kw.n_rows;
    double p_ber;
    // arma::vec a_vec(N);
    double sumlogber;
    arma::vec logber(N);
    const double eps = 1e-8;
    
    double wg = exp(zg) / (1 + exp(zg));
    double wg_new = exp(zg_new) / (1 + exp(zg_new));
    if(isold) {
        double logbeta = R::dbeta(wg, awvec_g, bwvec_g, 1);

        for (int i=0; i < N; i++) {
            p_ber = Kw_vec[i];
            if (Kw_vec[i] < eps) p_ber = eps;
            if (Kw_vec[i] > (1 - eps)) p_ber = 1 - eps;
            logber[i] = R::dbinom(ga[i], 1, p_ber, 1);
        }
        sumlogber = sum(logber);
        double logja = zg - 2 * log(1 + exp(zg));
        return List::create(
            Rcpp::Named("value") = logbeta + sumlogber + logja
        );
    } else {
        double logbeta = R::dbeta(wg_new, awvec_g, bwvec_g, 1);
        Kw_vec = Kw_vec + Kg * (wg_new - wg);
        Kw = mat(Kw_vec);
        for (int i=0; i < N; i++) {
            p_ber = Kw_vec[i];
            if (Kw_vec[i] < eps) p_ber = eps;
            if (Kw_vec[i] > (1 - eps)) p_ber = 1 - eps;
            // Kw_vec[i] = Kw_vec[i] + Kg[i] * (wg_new - wg);
            // logber[i] = R::dbinom(ga[i], 1, Kw_vec[i], 1);
            logber[i] = R::dbinom(ga[i], 1, p_ber, 1);
            // Rcout << "\r" << " Kw_vec[i] " << Kw_vec[i] << endl;
            // Rcout << "\r" << " logber[i] " << logber[i] << endl;
        }
        sumlogber = sum(logber);
        double logja = zg_new - 2 * log(1 + exp(zg_new));
        // Rcout << "\r" << "logbeta " << logbeta << endl;
        // Rcout << "\r" << "sumlogber " << sumlogber << endl;
        // Rcout << "\r" << "logja " << logja << endl;
        return List::create(
            Rcpp::Named("value") = logbeta + sumlogber + logja,
            Rcpp::Named("Kw") = Kw
        );
    }
}








