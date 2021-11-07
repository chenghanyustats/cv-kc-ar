//ga_update_gp_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <random>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int ga_update_gp_cpp(int v, double M0, double M1, arma::vec S, double n,
                      arma::vec regionlabel, arma::vec agvec, bool iscplx) {
    regionlabel -= 1;
    int g = regionlabel[v];
    double p_ber = 1 / (1 + exp(-(agvec[g] + S[g])));
    // double p_ber = 1 / (1 + exp(-(S[g])));
    double C;
    double po;
    double p_star;
    po = n / 2;
    // po = n;
//     if(sum(agvec) != 0) {
//         p_ber = 1 / (1 + exp(-(agvec[regionlabel[v]] + Kw_vec[v])));
//     } else {
//         p_ber = 1 / (1 + exp(-Kw_vec[v]));
//     }
    if (iscplx) {
        C = 1 + n / 2;
        // C = 1 + n;
    } else {
        C = sqrt(1 + n);
    }
//     C = sqrt(1 + n);
//     po = n / 2;
    // p_star = 0.5;
    p_star = 1 / (1 + ((1 - p_ber) / p_ber) * C * std::pow(M1 / M0, po));
    // Rcout << "\r" << "p_star: " << p_star << endl;
    return (rbinom(1, 1, p_star)[0]);
//     if(v == 399) {
//         Rcout << "\r" << "p_ber: " << p_ber << endl;
//         Rcout << "\r" << "p_star: " << p_star << endl;
//         Rcout << "\r" << "p_ber: " << p_ber << endl;
//         Rcout << "\r" << "p_star: " << p_star << endl;
//     }
//     std::default_random_engine generator;
//     std::binomial_distribution<int> distribution(1, p_star);
//     return(distribution(generator));
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List ga_update_gp_cpp_lst(int v, double M0, double M1, arma::vec S, double n,
                     arma::vec regionlabel, arma::vec agvec, bool iscplx) {
    regionlabel -= 1;
    int g = regionlabel[v];
    double p_ber = 1 / (1 + exp(-(agvec[g] + S[g])));
    // Rcout << "\r" << "p_ber: " << p_ber << endl;
    // double p_ber = 1 / (1 + exp(-(S[g])));
    double C;
    double po;
    double p_star;
    if (iscplx) {
        C = 1 + n;
        // Rcout << "\r" << "C: " << C << endl;
        po = n;
        // Rcout << "\r" << "po: " << po << endl;
    } else {
        C = sqrt(1 + n);
        po = n / 2;
    }
    double denom = (1 + ((1 - p_ber) / p_ber) * C * std::pow(M1 / M0, po));
    p_star = pow(denom, -1);
    // p_star = 0.5;
    return List::create(
//         _["M0"] = M0,
//         _["M1"] = M1,
        _["pow"] = std::pow(M1 / M0, po),
        _["p_ber"] = p_ber, 
        _["p_star"] = p_star, 
        _["gav"] = rbinom(1, 1, p_star)[0]
    );
    // return (rbinom(1, 1, p_star)[0]);
}










