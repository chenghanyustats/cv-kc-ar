//MHratio_cond_S_cpp.cpp


// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double MHratio_cond_S_cpp(arma::vec S, arma::vec S_new, double r, 
                          arma::vec ga, arma::mat Ga, int g, arma::vec regionlabel,
                          arma::vec agvec) {
    
    const double eps = 1e-8;
    arma::vec gag = ga.elem(find(regionlabel == g + 1));
    double p_old = 1 / (1 + exp(-(agvec[g] + S[g])));
    double p_new = 1 / (1 + exp(-(agvec[g] + S_new[g])));
//     Rcout << "\r" << "p_old: " << p_old << endl;
//     Rcout << "\r" << "p_new: " << p_new << endl;
    // # Avoid -Inf and Inf when take log
    if (p_new < eps) p_new = eps;
    if (p_new > (1 - eps)) p_new = 1 - eps;
    if (p_old < eps) p_old = eps;
    if (p_old > (1 - eps)) p_old = 1 - eps;
    
//     Rcout << "\r" <<  "ga_g: " << gag << endl;
//     Rcout << "\r" << "sum(ga_g): " << sum(gag) << endl;
//     Rcout << "\r" << "ga_g.n_elem: " << gag.n_elem << endl;
//     Rcout << "\r" << "log(p_new / p_old): " << log(p_new / p_old) << endl;
    double A = sum(gag) * log(p_new / p_old) + 
        (gag.n_elem - sum(gag)) * log((1 - p_new) / (1 - p_old));
    // Rcout << "\r" << "A: " << A << endl;
    arma::mat Gainv = inv_sympd(Ga);
    double quadSnew = (S_new.t() * Gainv * S_new).eval()(0, 0);
    double quadS = (S.t() * Gainv * S).eval()(0, 0);
    double lenS = S.n_elem; 
    double B = -(lenS / 2) * log(quadSnew / quadS);
//     Rcout << "\r" << "quadSnew : " << quadSnew << endl;
//     Rcout << "\r" << "quadS: " << quadS << endl;
//     Rcout << "\r" << "log(quadSnew / quadS): " << log(quadSnew / quadS) << endl;
//     Rcout << "\r" << "S.n_elem: " << lenS << endl;
//     Rcout << "\r" << "-(S.n_elem / 2): " << -(lenS / 2) << endl;
//     Rcout << "\r" << "B: " << B << endl;
//     return List::create(
//         Rcpp::Named("ga_g") = gag,
//         Rcpp::Named("MHratio") = A + B
//     );
    return (A + B);
}

