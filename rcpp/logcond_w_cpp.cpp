//logcond_w_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double logSGaS_cpp(arma::mat Ga, arma::vec S) {
    
    int lenS = S.n_elem;
    arma::mat Gainv = inv_sympd(Ga);
    double quadS = (S.t() * Gainv * S).eval()(0, 0);

    return(-(lenS / 2) * log(quadS));
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double logS_r_cpp(arma::mat Ga, arma::vec S) {
    
    return((-0.5) * log(det(Ga)) + logSGaS_cpp(Ga, S));
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double logcond_w_cpp(arma::mat Ga, arma::vec S, double w, int df) {
    return(logS_r_cpp(Ga, S) + 0.5 * df * w - (exp(w) / 2));
}
