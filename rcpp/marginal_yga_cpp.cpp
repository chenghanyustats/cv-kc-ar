//marginal_yga_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List marginal_yga_cpp(arma::vec Yvecv, arma::vec XtYv, arma::mat XtX) {
    arma::vec M_ga0 = Yvecv.t() * Yvecv;
    arma::vec M_ga1 = M_ga0 - XtYv.t() * inv_sympd(XtX) * XtYv;
//     Rcout << "\r" << "M_ga0: " << M_ga0 << endl;
//     Rcout << "\r" << "M_ga1: " << M_ga1 << endl;
    return List::create(
        Rcpp::Named("M0") = M_ga0,
        Rcpp::Named("M1") = M_ga1
    );
}