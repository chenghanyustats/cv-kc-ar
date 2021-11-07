//ga_update_pc1_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int ga_update_pc1_cpp(int v, double M0, double M1, arma::mat Kw, int n, 
                      arma::vec regionlabel, arma::vec agvec, bool iscplx) {
    arma::vec Kw_vec = vectorise(Kw);
    int g = regionlabel[v];
    double p_ber;
    int C;
    int po;
    double p_star;
    
    if(sum(agvec) != 0) {
        p_ber = 1 / (1 + exp(-(agvec[g - 1] + Kw_vec[v])));
    } else {
        p_ber = 1 / (1 + exp(-Kw_vec[v]));
    }
    if (p_ber == 0) p_ber = 1e-8;
    if (p_ber == 1) p_ber = 1 - 1e-8;
    
//     if (iscplx) {
//         C = 1 + n;
//         po = n;
//     } else {
//         C = sqrt(1 + n);
//         po = n / 2;
//     }
    po = n / 2;
    if (iscplx) {
        C = 1 + n / 2;
    } else {
        C = sqrt(1 + n);
    }
    p_star = 1 / (1 + ((1 - p_ber) / p_ber) * C * pow(M1 / M0, po));
    return rbinom(1, 1, p_star)[0];
}


