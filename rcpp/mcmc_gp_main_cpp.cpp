#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
double callgettune(double tune, int keep, int k, Function f) {
    double res;
    res = as<double>(f(tune, keep, k));
    return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat callpowerexpcov(arma::mat centroid_mat, double r, Function f) {
    arma::mat res;
    res = as<arma::mat>(f(centroid_mat, r));
    return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int callupdatega(arma::vec S, arma::vec regionlabel, int v, double M0, 
                       double M1, arma::vec agvec, bool iscplx, Function f) {
    int res;
    res = as<int>(f(S, regionlabel, v, M0, M1, agvec, iscplx));
    return res;
}

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double MHratio_cond_S_cpp(arma::vec S, arma::vec S_new, double r,
//                           arma::vec ga, arma::mat Ga, int g, arma::vec regionlabel,
//                           arma::vec agvec) {
// 
//     const double eps = 1e-8;
//     arma::vec gag = ga.elem(find(regionlabel == g + 1));
//     double p_old = 1 / (1 + exp(-(agvec[g] + S[g])));
//     double p_new = 1 / (1 + exp(-(agvec[g] + S_new[g])));
//     //     Rcout << "\r" << "p_old: " << p_old << endl;
//     //     Rcout << "\r" << "p_new: " << p_new << endl;
//     // # Avoid -Inf and Inf when take log
//     // if (p_new < eps) p_new = eps;
//     // if (p_new > (1 - eps)) p_new = 1 - eps;
//     // if (p_old < eps) p_old = eps;
//     // if (p_old > (1 - eps)) p_old = 1 - eps;
// 
//     //     Rcout << "\r" <<  "ga_g: " << gag << endl;
//     //     Rcout << "\r" << "sum(ga_g): " << sum(gag) << endl;
//     //     Rcout << "\r" << "ga_g.n_elem: " << gag.n_elem << endl;
//     //     Rcout << "\r" << "log(p_new / p_old): " << log(p_new / p_old) << endl;
//     double A = sum(gag) * log(p_new / p_old) +
//         (gag.n_elem - sum(gag)) * log((1 - p_new) / (1 - p_old));
//     // Rcout << "\r" << "A: " << A << endl;
//     arma::mat Gainv = inv_sympd(Ga);
//     double quadSnew = (S_new.t() * Gainv * S_new).eval()(0, 0);
//     double quadS = (S.t() * Gainv * S).eval()(0, 0);
//     double lenS = S.n_elem;
//     double B = -(lenS / 2) * log(quadSnew / quadS);
//     //     Rcout << "\r" << "quadSnew : " << quadSnew << endl;
//     //     Rcout << "\r" << "quadS: " << quadS << endl;
//     //     Rcout << "\r" << "log(quadSnew / quadS): " << log(quadSnew / quadS) << endl;
//     //     Rcout << "\r" << "S.n_elem: " << lenS << endl;
//     //     Rcout << "\r" << "-(S.n_elem / 2): " << -(lenS / 2) << endl;
//     //     Rcout << "\r" << "B: " << B << endl;
//     //     return List::create(
//     //         Rcpp::Named("ga_g") = gag,
//     //         Rcpp::Named("MHratio") = A + B
//     //     );
//     return (A + B);
// }



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double MHratio_cond_S_cpp(arma::vec S, arma::vec S_new, double r,
                          arma::vec ga, arma::mat Gainv, int g, arma::vec regionlabel,
                          arma::vec agvec) {

    const double eps = 1e-8;
    arma::vec gag = ga.elem(find(regionlabel == g + 1));
    double p_old = 1 / (1 + exp(-(agvec[g] + S[g])));
    double p_new = 1 / (1 + exp(-(agvec[g] + S_new[g])));
    //     Rcout << "\r" << "p_old: " << p_old << endl;
    //     Rcout << "\r" << "p_new: " << p_new << endl;
    // # Avoid -Inf and Inf when take log
    // if (p_new < eps) p_new = eps;
    // if (p_new > (1 - eps)) p_new = 1 - eps;
    // if (p_old < eps) p_old = eps;
    // if (p_old > (1 - eps)) p_old = 1 - eps;

    //     Rcout << "\r" <<  "ga_g: " << gag << endl;
    //     Rcout << "\r" << "sum(ga_g): " << sum(gag) << endl;
    //     Rcout << "\r" << "ga_g.n_elem: " << gag.n_elem << endl;
    //     Rcout << "\r" << "log(p_new / p_old): " << log(p_new / p_old) << endl;
    double A = sum(gag) * log(p_new / p_old) +
        (gag.n_elem - sum(gag)) * log((1 - p_new) / (1 - p_old));
    // Rcout << "\r" << "A: " << A << endl;
    // arma::mat Gainv = inv_sympd(Ga);
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



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List marginal_yga_cpp(arma::mat Yvec, arma::mat XtY, arma::mat XtX) {
    int Nvoxels = Yvec.n_rows;
    arma::vec M_ga0(Nvoxels);
    arma::vec M_ga1(Nvoxels);
    arma::vec Yvecv;
    arma::vec XtYv;
    arma::mat inv_XtX = inv_sympd(XtX);
    // arma::vec XtYvt = XtYv.t();
    for (int v=0; v < Nvoxels; ++v) {
        Yvecv = vectorise(Yvec.row(v));
        XtYv = vectorise(XtY.col(v));
        M_ga0[v] = (Yvecv.t() * Yvecv).eval()(0,0);
        M_ga1[v] = M_ga0[v] -(XtYv.t() * inv_XtX * XtYv).eval()(0,0);
    }
    return List::create(
        Rcpp::Named("M0") = M_ga0,
        Rcpp::Named("M1") = M_ga1
    );
}


// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List get_ar_mat(double arcoef, int n) {
//     arma::mat Lambda(n, n);
//     arma::mat Lambda_inv;
//     for (int i=0; i < n; ++i) {
//         for (int j=0; j < n; ++j) {
//             Lambda(i, j) = pow(arcoef, abs(i-j));
//         }
//     }
//     Lambda_inv = inv_sympd(Lambda);
//     
//     // return Lambda;
//     return List::create(
//         Rcpp::Named("Lambda") = Lambda,
//         Rcpp::Named("Lambda_inv") = Lambda_inv)
//     ;
// }


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat get_ar1_mat(double arcoef, int n) {
    arma::mat Lambda(n, n);
    // arma::mat Lambda_inv;
    for (int i=0; i < n; ++i) {
        for (int j=0; j < n; ++j) {
            Lambda(i, j) = pow(arcoef, abs(i-j));
        }
    }
    // Lambda_inv = inv_sympd(Lambda);
    
    return Lambda;
    // return List::create(
    //     // Rcpp::Named("Lambda") = Lambda
        // Rcpp::Named("Lambda_inv") = Lambda_inv
    // )
    // ;
}

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::mat get_ar1_inv_mat(double arcoef, int n) {
//     // arma::mat Lambda(n, n);
//     double d = 1 / (1 - pow(arcoef, 2));
//     arma::mat L(n, n, fill::zeros);
//     for (int i=0; i < (n - 1); ++i) {
//         L(i + 1, i) = d * (-arcoef);
//     }
//     arma::mat D = L.t() + L;
//     D(0, 0) = d;
//     D(n - 1, n - 1) = d;
//     for (int i=1; i < (n - 1); ++i) {
//         D(i, i) = d * (1 + pow(arcoef, 2));
//     }
// 
//     return D;
// }

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::mat get_ar1_inv_mat_new(double arcoef, int n) {
//     // arma::mat Lambda(n, n);
//     double sq = pow(arcoef, 2);
//     double d = 1 / (1 - sq);
//     arma::mat L(n, n, fill::zeros);
//     for (int i=0; i < (n - 1); ++i) {
//         L(i + 1, i) = -d * arcoef;
//     }
//     arma::mat D = L.t() + L;
//     // D(0, 0) = d;
//     // D(n - 1, n - 1) = d;
//     for (int i=0; i < n; ++i) {
//         if (i == 0 || i == n - 1) {
//             D(i, i) = d;
//         } else {
//             D(i, i) = d * (1 + sq);
//         }
//     }
//     
//     return D;
// }


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat get_ar1_inv_mat_new2(double arcoef, int n, bool iscplx) {
    // arma::mat Lambda(n, n);
    int m;
    if(iscplx) {
        m = n / 2;
    } else {
        m = n;
    }
    double sq = pow(arcoef, 2);
    double d = 1 / (1 - sq);
    double q = -d * arcoef;
    // arma::mat D(n, n, fill::zeros);
    arma::mat D = zeros<mat>(m, m);
    for (int i=0; i < m; ++i) {
        if (i == 0) {
            D(i + 1, i) = q;
            D(i, i) = d;
        } else if (i == m - 1) {
            D(i - 1, i) = q;
            D(i, i) = d;
        } else {
            D(i - 1, i) = q;
            D(i + 1, i) = q;
            D(i, i) = d * (1 + sq);
        }
    }
    // for (int i=1; i < m - 1; ++i) {
    //     D(i - 1, i) = q;
    //     D(i + 1, i) = q;
    //     D(i, i) = d * (1 + sq);
    // }
    // D(1, 0) = q;
    // D(0, 0) = d;
    // D(m - 2, m - 1) = q;
    // D(m - 1, m - 1) = d;
    
    if (iscplx) {
        arma::mat B = zeros<mat>(2 * m, 2 * m);  
        // B.zeros(2 * m, 2 * m);
        B.submat(0, 0, m - 1, m - 1) = D;
        B.submat(m, m, 2 * m - 1, 2 * m - 1) = D;
        return B;
    } else {
        return D;
    }
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat get_ar1_inv_mat_new2_cv(double arcoef_re, double arcoef_im, int n) {
    // arma::mat Lambda(n, n);
    int m = n / 2;
    double sq_re = pow(arcoef_re, 2);
    double d_re = 1 / (1 - sq_re);
    double q_re = -d_re * arcoef_re;
    double sq_im = pow(arcoef_im, 2);
    double d_im = 1 / (1 - sq_im);
    double q_im = -d_im * arcoef_im;
    // arma::mat D(n, n, fill::zeros);
    arma::mat D_re = zeros<mat>(m, m);
    arma::mat D_im = zeros<mat>(m, m);
    for (int i=0; i < m; ++i) {
        if (i == 0) {
            D_re(i + 1, i) = q_re;
            D_re(i, i) = d_re;
            D_im(i + 1, i) = q_im;
            D_im(i, i) = d_im;
        } else if (i == m - 1) {
            D_re(i - 1, i) = q_re;
            D_re(i, i) = d_re;
            D_im(i - 1, i) = q_im;
            D_im(i, i) = d_im;
        } else {
            D_re(i - 1, i) = q_re;
            D_re(i + 1, i) = q_re;
            D_re(i, i) = d_re * (1 + sq_re);
            D_im(i - 1, i) = q_im;
            D_im(i + 1, i) = q_im;
            D_im(i, i) = d_im * (1 + sq_im);
        }
    }
    // for (int i=1; i < m - 1; ++i) {
    //     D(i - 1, i) = q;
    //     D(i + 1, i) = q;
    //     D(i, i) = d * (1 + sq);
    // }
    // D(1, 0) = q;
    // D(0, 0) = d;
    // D(m - 2, m - 1) = q;
    // D(m - 1, m - 1) = d;
    
    arma::mat B = zeros<mat>(n, n);  
    // B.zeros(2 * m, 2 * m);
    B.submat(0, 0, m - 1, m - 1) = D_re;
    B.submat(m, m, n - 1, n - 1) = D_im;
    return B;
}
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List marginal_yga_ar_cpp(arma::mat Yvec, arma::mat Xr, arma::vec arcoefvec) {
//     int Nvoxels = Yvec.n_rows;
//     int n = Yvec.n_cols;
//     arma::vec M_ga0(Nvoxels);
//     arma::vec M_ga1(Nvoxels);
//     arma::mat XtX_ar;
//     arma::mat XtY_ar;
//     arma::vec Yvecv;
//     arma::vec XtYv_ar;
//     arma::mat XtXinv_ar;
//     // arma::vec XtYvt = XtYv.t();
//     for (int v=0; v < Nvoxels; ++v) {
//         List Lam_list = get_ar_mat(arcoefvec[v], n);
//         arma::mat Lambda_inv = Lam_list["Lambda_inv"];
//         Yvecv = vectorise(Yvec.row(v));
//         XtX_ar = Xr.t() * Lambda_inv * Xr;
//         XtY_ar = Xr.t() * Lambda_inv * Yvec.t();
//         XtXinv_ar = inv_sympd(XtX_ar);
//         XtYv_ar = vectorise(XtY_ar.col(v));
//         M_ga0[v] = (Yvecv.t() * Lambda_inv * Yvecv).eval()(0,0);
//         M_ga1[v] = M_ga0[v] -(XtYv_ar.t() * XtXinv_ar * XtYv_ar).eval()(0,0);
//     }
//     return List::create(
//         Rcpp::Named("M0") = M_ga0,
//         Rcpp::Named("M1") = M_ga1
//     );
// }



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int ga_update_gp_cpp(int v, double M0, double M1, arma::vec S, int n,
                     arma::vec regionlabel, arma::vec agvec, bool iscplx) {
    int g = regionlabel[v];
    const double eps = 1e-8;
    double p_ber = 1 / (1 + exp(-(agvec[g - 1] + S[g - 1])));
    if (p_ber < eps) p_ber = eps;
    if (p_ber > (1 - eps)) p_ber = 1 - eps;
    int C;
    int po;
    double p_star;
    // if (iscplx) {
    //     C = 1 + n / 2;
    //     po = n;
    // } else {
    //     C = sqrt(1 + n);
    //     po = n / 2;
    // }
    po = n / 2;
    if (iscplx) {
        C = 1 + n / 2;
    } else {
        C = sqrt(1 + n);
    }
    p_star = 1 / (1 + ((1 - p_ber) / p_ber) * C * pow(M1 / M0, po));
    return rbinom(1, 1, p_star)[0];
}


// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double logSGaS_cpp(arma::mat Ga, arma::vec S) {
// 
//     int lenS = S.n_elem;
//     arma::mat Gainv = inv_sympd(Ga);
//     double quadS = (S.t() * Gainv * S).eval()(0, 0);
// 
//     return(-(lenS / 2) * log(quadS));
// }


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double logSGaS_cpp(arma::mat Gainv, arma::vec S) {

    int lenS = S.n_elem;
    // arma::mat Gainv = inv_sympd(Ga);
    double quadS = (S.t() * Gainv * S).eval()(0, 0);

    return(-(lenS / 2) * log(quadS));
}

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double logS_r_cpp(arma::mat Ga, arma::vec S) {
// 
//     return((-0.5) * log(det(Ga)) + logSGaS_cpp(Ga, S));
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double logS_r_cpp(arma::mat Ga, arma::mat Gainv, arma::vec S) {
    // double val;
    // double sign;
    return((-0.5) * log(det(Ga)) + logSGaS_cpp(Gainv, S));
}


// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double logcond_w_cpp(arma::mat Ga, arma::vec S, double w, int df) {
//     return(logS_r_cpp(Ga, S) + 0.5 * df * w - (exp(w) / 2));
// }


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double logcond_w_cpp(arma::mat Ga, arma::mat Gainv, arma::vec S, double w, int df) {
    return(logS_r_cpp(Ga, Gainv, S) + 0.5 * df * w - (exp(w) / 2));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_gp_main_cpp(arma::mat Yvec, arma::mat Xr, 
                       arma::vec S, arma::vec ga, double w, double r,
                       int nmcmc, int burn, int thin, 
                       bool adapt, List tunelst, List keeplst, 
                       List tuner, List keepr, int tunelen, 
                       arma::vec regionlabel, int df, 
                       arma::vec agvec, arma::mat centroid, bool iscplx,
                       Function get_tune, Function power_exp_cov_fcn,
                       Function update_ga,
                       arma::mat draws) {
    
    
    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec.n_rows;
    const int n = Yvec.n_cols;
    const int p = 1;
    int G = regionlabel.max();
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();

    List MM = marginal_yga_cpp(Yvec, XtY, XtX);
    arma::vec M0 = MM["M0"];
    arma::vec M1 = MM["M1"];
    
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    Rcpp::List keep = keeplst;
    Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    
    // # for psi
    // # keep.psi <- list(psi = 0)
    Rcpp::List keeptmpr = keepr;
    int Tb = tunelen;  //# frequency of adaptive tuning
    
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    if (iscplx) {
        Rcout << "CV-GP Begins" << endl;
    } else {
        Rcout << "MO-GP Begins" << endl;
    }

    for (int i = 0; i < nmcmc; ++i) {
        
        // print iteration
        
        
        // # Update tuning parameter of psi
        // #------------------------
        if(adapt == TRUE & fmod(i , Tb) == 0) {  
            // # Adaptive tuning
            // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
            keeptmpr["r"] = as<double>(keeptmpr["r"]) / Tb;
            tuner["r"] = callgettune(tuner["r"], keeptmpr["r"], i, get_tune);
            keeptmpr["r"] = 0;
        }
        
        // # Update parameters
        // #----------------------------
        //         
        // # ========== sample S =============
        arma::mat Ga_cov = callpowerexpcov(centroid, r, power_exp_cov_fcn);
        arma::mat inv_Ga_cov = inv_sympd(Ga_cov);
            
        for (int g = 0; g < G; ++g) {
            // Rcout << "\r" << "g: " << g + 1 << endl;
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmp[g] = as<double>(keeptmp[g]) / Tb;
                tunelst[g] = callgettune(tunelst[g], keeptmp[g], i, get_tune);
                keeptmp[g] = 0;
            }
            
            // #----------------------------
            // # Random walk Normal proposal 
            // #----------------------------
            arma::vec S_new = S;
            S_new[g] = rnorm(1, S[g], tunelst[g])[0];
            if (log(runif(1, 0, 1)[0]) < MHratio_cond_S_cpp(S, S_new, r, ga, 
                    inv_Ga_cov, g, regionlabel, agvec)) {
                    S[g] = S_new[g];
                    keep[g] = as<int>(keep[g]) + 1;
                    keeptmp[g] = as<int>(keeptmp[g]) + 1;
            }
        }
        
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            // Rcout << "\r" << "v: " << v + 1 << endl;
//             List MM = marginal_yga_cpp(vectorise(Yvec.row(v)), 
//                                        vectorise(XtY.col(v)), XtX);
//             ga[v] = ga_update_gp_cpp(v, as<double>(MM["M0"]), as<double>(MM["M1"]), 
//                                      S, n, regionlabel, agvec, iscplx);
//             ga[v] = callupdatega(S, regionlabel, v + 1, as<double>(MM["M0"]), 
//                                  as<double>(MM["M1"]), agvec, iscplx, update_ga);
            ga[v] = ga_update_gp_cpp(v, M0[v], M1[v], S, n, regionlabel, agvec, iscplx);
        }
        
        // # ============ sample w and hence r ==============
        double new_w = rnorm(1, w, tuner["r"])[0];
        arma::mat new_Ga_cov = callpowerexpcov(centroid, exp(new_w), 
                                               power_exp_cov_fcn);
        arma::mat new_Ga_inv = inv_sympd(new_Ga_cov);
        if (log(runif(1, 0, 1)[0]) < (logcond_w_cpp(new_Ga_cov, new_Ga_inv, S, new_w, df) -
            logcond_w_cpp(Ga_cov, inv_Ga_cov, S, w, df))) {
            w = new_w;
            r = exp(w);
            keepr["r"] = as<int>(keepr["r"]) + 1;
            keeptmpr["r"] = as<int>(keeptmpr["r"]) + 1;
        }
        // if (log(runif(1, 0, 1)[0]) < (logcond_w_cpp(new_Ga_cov, S, new_w, df) - 
        //     logcond_w_cpp(Ga_cov, S, w, df))) {
        //     w = new_w;
        //     r = exp(w);
        //     keepr["r"] = as<int>(keepr["r"]) + 1;
        //     keeptmpr["r"] = as<int>(keeptmpr["r"]) + 1;
        // } 
        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = S.t();
                draws(d - 1, draws.n_cols - 1) = r;
                //                 Gamma.row(d) = ga.t();
                //                 W.row(d) = w.t();
                //                 Phi.row(d) <- phi;
            } 
        }
        if (fmod(i , 1000) == 0) {
            Rcout << "\r" << "MCMC Iter: " << i + 1 << std::flush;
        }
    }
    // draws = join_rows(join_rows(Gamma, W), Phi);
    
    Rcout << "\n" << "Done! " << endl;
    
    // # Write output
    // #--------------
    return List::create(
        _["draws"] = draws, 
        _["accept_keep"] = keep, 
        _["accept_keep_r"] = keepr, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_gp_ar_main_cpp(arma::mat Yvec, arma::mat Xr, arma::vec arcoefvec, 
                      arma::vec S, arma::vec ga, double w, double r,
                      int nmcmc, int burn, int thin, 
                      bool adapt, List tunelst, List keeplst, 
                      List tuner, List keepr, int tunelen, 
                      arma::vec regionlabel, int df, 
                      arma::vec agvec, arma::mat centroid, bool iscplx,
                      Function get_tune, Function power_exp_cov_fcn,
                      Function update_ga,
                      arma::mat draws) {
    
    
    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec.n_rows;
    const int n = Yvec.n_cols;
    const int p = 1;
    int G = regionlabel.max();
    // arma::mat XtX = Xr.t() * Xr;
    // arma::mat XtY = Xr.t() * Yvec.t();
    
    // List MM = marginal_yga_cpp(Yvec, XtY, XtX);
    // arma::vec M0 = MM["M0"];
    // arma::vec M1 = MM["M1"];
    // int Nvoxels = Yvec.n_rows;
    // int n = Yvec.n_cols;
    
    // arma::vec M_ga0(Nvoxels);
    // arma::vec M_ga1(Nvoxels);
    double M_ga0;
    double M_ga1;
    arma::mat XtLaminv_ar;
    // arma::mat XtX_ar;
    arma::mat XtY_ar;
    arma::vec Yvecv;
    arma::vec XtYv_ar;
    arma::mat Lambda_inv;
    // arma::mat XtXinv_ar;
    // arma::vec XtYvt = XtYv.t();
    // for (int v=0; v < Nvoxels; ++v) {
    // 
    // }
    
    
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    Rcpp::List keep = keeplst;
    Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    
    // # for psi
    // # keep.psi <- list(psi = 0)
    Rcpp::List keeptmpr = keepr;
    int Tb = tunelen;  //# frequency of adaptive tuning
    

    
    arma::vec M_ga0_vec(Nvoxels);
    arma::vec M_ga1_vec(Nvoxels);
    for (int v = 0; v < Nvoxels; ++v) {
        Lambda_inv = get_ar1_inv_mat_new2(arcoefvec[v], n, iscplx);
        Yvecv = vectorise(Yvec.row(v));
        XtLaminv_ar = Xr.t() * Lambda_inv;
        XtYv_ar = XtLaminv_ar * Yvecv;
        M_ga0_vec[v] = (Yvecv.t() * Lambda_inv * Yvecv).eval()(0,0);
        M_ga1_vec[v] = M_ga0_vec[v] - (XtYv_ar.t() * inv_sympd(XtLaminv_ar * Xr) * XtYv_ar).eval()(0,0);
    }   
    
    // for (int v = 0; v < Nvoxels; ++v) {
    //     Lambda_inv = get_ar1_inv_mat_new2(arcoefvec[v], n);
    //     Yvecv = vectorise(Yvec.row(v));
    //     XtLaminv_ar = Xr.t() * Lambda_inv;
    //     XtYv_ar = XtLaminv_ar * Yvecv;
    //     M_ga0_vec[v] = (Yvecv.t() * Lambda_inv * Yvecv).eval()(0,0);
    //     M_ga1_vec[v] = M_ga0_vec[v] - (XtYv_ar.t() * inv_sympd(XtLaminv_ar * Xr) * XtYv_ar).eval()(0,0);
    // }
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    if (iscplx) {
        Rcout << "CV-GP-AR Begins" << endl;
    } else {
        Rcout << "MO-GP-AR Begins" << endl;
    }
    
    for (int i = 0; i < nmcmc; ++i) {
        
        // print iteration
        
        
        // # Update tuning parameter of psi
        // #------------------------
        if(adapt == TRUE & fmod(i , Tb) == 0) {  
            // # Adaptive tuning
            // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
            keeptmpr["r"] = as<double>(keeptmpr["r"]) / Tb;
            tuner["r"] = callgettune(tuner["r"], keeptmpr["r"], i, get_tune);
            keeptmpr["r"] = 0;
        }
        
        // # Update parameters
        // #----------------------------
        //         
        // # ========== sample S =============
        arma::mat Ga_cov = callpowerexpcov(centroid, r, power_exp_cov_fcn);
        arma::mat inv_Ga_cov = inv_sympd(Ga_cov);
        
        for (int g = 0; g < G; ++g) {
            // Rcout << "\r" << "g: " << g + 1 << endl;
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmp[g] = as<double>(keeptmp[g]) / Tb;
                tunelst[g] = callgettune(tunelst[g], keeptmp[g], i, get_tune);
                keeptmp[g] = 0;
            }
            
            // #----------------------------
            // # Random walk Normal proposal 
            // #----------------------------
            arma::vec S_new = S;
            S_new[g] = rnorm(1, S[g], tunelst[g])[0];
            if (log(runif(1, 0, 1)[0]) < MHratio_cond_S_cpp(S, S_new, r, ga, 
                    inv_Ga_cov, g, regionlabel, agvec)) {
                S[g] = S_new[g];
                keep[g] = as<int>(keep[g]) + 1;
                keeptmp[g] = as<int>(keeptmp[g]) + 1;
            }
        }
        
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            // Rcout << "\r" << "v: " << v + 1 << endl;
            //             List MM = marginal_yga_cpp(vectorise(Yvec.row(v)), 
            //                                        vectorise(XtY.col(v)), XtX);
            //             ga[v] = ga_update_gp_cpp(v, as<double>(MM["M0"]), as<double>(MM["M1"]), 
            //                                      S, n, regionlabel, agvec, iscplx);
            //             ga[v] = callupdatega(S, regionlabel, v + 1, as<double>(MM["M0"]), 
            //                                  as<double>(MM["M1"]), agvec, iscplx, update_ga);
            // List Lam_list = get_ar_mat(arcoefvec[v], n);
            // arma::mat Lambda_inv = Lam_list["Lambda_inv"];
            // arma::mat Lambda_inv = get_ar1_inv_mat(arcoefvec[v], n);
            // arma::mat Lambda_inv = get_ar1_inv_mat_new(arcoefvec[v], n);
            // Lambda_inv = get_ar1_inv_mat_new2(arcoefvec[v], n, iscplx);
            // Yvecv = vectorise(Yvec.row(v));
            // XtLaminv_ar = Xr.t() * Lambda_inv;
            // if(v == 0) {
            //     Rcout << "\r" << "Size of Lam: " << Lambda_inv.submat(0, 0, 49, 49) << endl;
            // }
            // Rcout << "\r" << "Size of Xr: " << Xr.size() << endl;
            // XtX_ar = Xr.t() * Lambda_inv * Xr;
            // XtX_ar = XtLaminv_ar * Xr;
            // Rcout << "XtX_ar:" << size(XtX_ar) << endl;
            // XtXinv_ar = inv_sympd(XtX_ar);
            // XtY_ar = Xr.t() * Lambda_inv * Yvec.t();
            // Rcout << "XtY_ar:" << size(XtY_ar) << endl;
            // XtYv_ar = vectorise(XtY_ar.col(v));
            // XtYv_ar = Xr.t() * Lambda_inv * Yvecv;
            // XtYv_ar = XtLaminv_ar * Yvecv;
            // Rcout << "XtYv_ar:" << size(XtYv_ar) << endl;
            // M_ga0[v] = (Yvecv.t() * Lambda_inv * Yvecv).eval()(0,0);
            // M_ga1[v] = M_ga0[v] - (XtYv_ar.t() * XtXinv_ar * XtYv_ar).eval()(0,0);
            // M_ga0 = (Yvecv.t() * Lambda_inv * Yvecv).eval()(0,0);
            // M_ga1 = M_ga0[v] - (XtYv_ar.t() * XtXinv_ar * XtYv_ar).eval()(0,0);
            // M_ga1 = M_ga0 - (XtYv_ar.t() * XtXinv_ar * XtYv_ar).eval()(0,0);
            // M_ga1 = M_ga0 - (XtYv_ar.t() * inv_sympd(XtX_ar) * XtYv_ar).eval()(0,0);
            // M_ga1 = M_ga0 - (XtYv_ar.t() * inv_sympd(XtLaminv_ar * Xr) * XtYv_ar).eval()(0,0);
            // M_ga1 = M_ga0 - (XtYv_ar.t() * solve(XtLaminv_ar * Xr, XtYv_ar)).eval()(0,0);
            // List MM = marginal_yga_ar_cpp(Yvec, Xr, arcoefvec);
            // arma::vec M0 = MM["M0"];
            // arma::vec M1 = MM["M1"];
            // ga[v] = ga_update_gp_cpp(v, M_ga0[v], M_ga1[v], S, n, regionlabel, agvec, iscplx);
            ga[v] = ga_update_gp_cpp(v, M_ga0_vec[v], M_ga1_vec[v], S, n, regionlabel, agvec, iscplx);
        }
        
        // # ============ sample w and hence r ==============
        double new_w = rnorm(1, w, tuner["r"])[0];
        arma::mat new_Ga_cov = callpowerexpcov(centroid, exp(new_w), 
                                               power_exp_cov_fcn);
        arma::mat new_Ga_inv = inv_sympd(new_Ga_cov);
        if (log(runif(1, 0, 1)[0]) < (logcond_w_cpp(new_Ga_cov, new_Ga_inv, S, new_w, df) -
            logcond_w_cpp(Ga_cov, inv_Ga_cov, S, w, df))) {
            w = new_w;
            r = exp(w);
            keepr["r"] = as<int>(keepr["r"]) + 1;
            keeptmpr["r"] = as<int>(keeptmpr["r"]) + 1;
        }
        // if (log(runif(1, 0, 1)[0]) < (logcond_w_cpp(new_Ga_cov, S, new_w, df) - 
        //     logcond_w_cpp(Ga_cov, S, w, df))) {
        //     w = new_w;
        //     r = exp(w);
        //     keepr["r"] = as<int>(keepr["r"]) + 1;
        //     keeptmpr["r"] = as<int>(keeptmpr["r"]) + 1;
        // } 
        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = S.t();
                draws(d - 1, draws.n_cols - 1) = r;
                //                 Gamma.row(d) = ga.t();
                //                 W.row(d) = w.t();
                //                 Phi.row(d) <- phi;
            } 
        }
        if (fmod(i , 1000) == 0) {
            Rcout << "\r" << "MCMC Iter: " << i + 1 << std::flush;
        }
    }
    // draws = join_rows(join_rows(Gamma, W), Phi);
    
    Rcout << "\n" << "Done! " << endl;
    
    // # Write output
    // #--------------
    return List::create(
        _["draws"] = draws, 
        _["accept_keep"] = keep, 
        _["accept_keep_r"] = keepr, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}












