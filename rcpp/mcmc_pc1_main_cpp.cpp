#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;


const double log2pi = std::log(2.0 * M_PI);

// struct add_multiple {
//     int incr;
//     int count;
//     add_multiple(int incr)
//         : incr(incr), count(0)
//     {}
//     inline int operator()(int d) {
//         return d + incr * count++;
//     }
// };
// 
// 
// 
// // [[Rcpp::export]]
// Rcpp::NumericVector rcpp_seq(double from_, double to_, double by_ = 1.0) {
//     int adjust = std::pow(10, std::ceil(std::log10(10 / by_)) - 1);
//     int from = adjust * from_;
//     int to = adjust * to_;
//     int by = adjust * by_;
//     
//     std::size_t n = ((to - from) / by) + 1;
//     Rcpp::IntegerVector res = Rcpp::rep(from, n);
//     add_multiple ftor(by);
//     
//     std::transform(res.begin(), res.end(), res.begin(), ftor);
//     return Rcpp::NumericVector(res) / adjust;
// }



// [[Rcpp::export]]
double callgettune(double tune, double keep, int k, Function f) {
    double res;
    res = as<double>(f(tune, keep, k));
    return res;
}


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
    return R::rbinom(1, p_star);
}


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


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_latent_block_cpp(arma::vec w, arma::vec w_new, arma::vec ga, 
                        double atau, double btau, arma::mat K, arma::mat Kw, 
                        arma::uvec widx, arma::vec agvec, 
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
        w.elem(widx) = w_new;
        // Rcout << "w_new_vec" << w_new_vec.n_elem << endl;
        
        Kw = K * w;
        Kw_vec = vec(Kw);
        // Kw = mat(Kw_vec);
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
        // w[gidx] = wg_new;
        double B = (G / 2 + atau) * log((0.5) * sum(square(w)) + btau);
        return List::create(
            Rcpp::Named("value") = sumlogber - B,
            Rcpp::Named("Kw") = Kw
        );
    }
}

// logcond_latent_block2 <- function(w, w_new, ga, a.tau, b.tau, K, Kw, w.idx, ag.vec, 
//                                   is.old = TRUE, regionlabel) {
//     if(is.old) {
//         if(sum(ag.vec) != 0) {
//             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
//             p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         } else {
//             p_bern <- 1 / (1 + exp(-Kw))
//         }
// # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
// # G <- dim(w)[1]
//             B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
//             return(list(value = sumlogber - B))
//     } else {
//         w[w.idx] <- w_new
//         Kw <- K %*% w
//         if(sum(ag.vec) != 0) {
//             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
//             p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         } else {
//             p_bern <- 1 / (1 + exp(-Kw))
//         }
// # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
// # G <- dim(w)[1]
// # w[g] <- wg_new
//             B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
//             return(list(value = sumlogber - B, Kw = Kw))
//     }
// } 
















// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List marginal_yga_cpp(arma::mat Yvec, arma::mat XtY, arma::mat XtX) {
    int Nvoxels = Yvec.n_rows;
    arma::vec M_ga0(Nvoxels);
    arma::vec M_ga1(Nvoxels);
    arma::vec Yvecv;
    arma::vec XtYv;
    arma::mat XtXinv = inv_sympd(XtX);
    // arma::vec XtYvt = XtYv.t();
    for (int v=0; v < Nvoxels; ++v) {
        Yvecv = vectorise(Yvec.row(v));
        XtYv = vectorise(XtY.col(v));
        M_ga0[v] = (Yvecv.t() * Yvecv).eval()(0,0);
        M_ga1[v] = M_ga0[v] -(XtYv.t() * XtXinv * XtYv).eval()(0,0);
    }
    return List::create(
        Rcpp::Named("M0") = M_ga0,
        Rcpp::Named("M1") = M_ga1
    );
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List get_ar_mat(double arcoef, int n) {
    arma::mat Lambda(n, n);
    arma::mat Lambda_inv;
    for (int i=0; i < n; ++i) {
        for (int j=0; j < n; ++j) {
            Lambda(i, j) = pow(arcoef, abs(i-j));
        }
    }
    Lambda_inv = inv_sympd(Lambda);
    
    return List::create(
        Rcpp::Named("Lambda") = Lambda,
        Rcpp::Named("Lambda_inv") = Lambda_inv
    );
}


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


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat get_ar1_inv_mat_new2(double arcoef, int n) {
    // arma::mat Lambda(n, n);
    double sq = pow(arcoef, 2);
    double d = 1 / (1 - sq);
    double q = -d * arcoef;
    // arma::mat D(n, n, fill::zeros);
    arma::mat D = zeros<mat>(n, n);
    for (int i=0; i < n; ++i) {
        if (i == 0) {
            D(i + 1, i) = q;
            D(i, i) = d;
        } else if (i == n - 1) {
            D(i - 1, i) = q;
            D(i, i) = d;
        } else {
            D(i - 1, i) = q;
            D(i + 1, i) = q;
            D(i, i) = d * (1 + sq);
        }
    }
    return D;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List marginal_yga_ar_cpp(arma::mat Yvec, arma::mat Xr, arma::vec arcoefvec) {
    int Nvoxels = Yvec.n_rows;
    int n = Yvec.n_cols;
    arma::vec M_ga0(Nvoxels);
    arma::vec M_ga1(Nvoxels);
    arma::mat XtX_ar;
    arma::mat XtY_ar;
    arma::vec Yvecv;
    arma::vec XtYv_ar;
    arma::mat XtXinv_ar = inv_sympd(XtX_ar);
    // arma::vec XtYvt = XtYv.t();
    for (int v=0; v < Nvoxels; ++v) {
        List Lam_list = get_ar_mat(arcoefvec[v], n);
        arma::mat Lambda_inv = Lam_list["Lambda_inv"];
        Yvecv = vectorise(Yvec.row(v));
        XtX_ar = Xr.t() * Lambda_inv * Xr;
        XtY_ar = Xr.t() * Lambda_inv * Yvec.t();
        XtXinv_ar = inv_sympd(XtX_ar);
        XtYv_ar = vectorise(XtY_ar.col(v));
        M_ga0[v] = (Yvecv.t() * Lambda_inv * Yvecv).eval()(0,0);
        M_ga1[v] = M_ga0[v] -(XtYv_ar.t() * XtXinv_ar * XtYv_ar).eval()(0,0);
    }
    return List::create(
        Rcpp::Named("M0") = M_ga0,
        Rcpp::Named("M1") = M_ga1
    );
}






// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_psi_cpp(double psi, arma::vec w, arma::mat Kw, 
                     double aphi, double bphi, arma::mat distmat, 
                     double nu, arma::vec ga, arma::vec agvec, 
                     bool isold, arma::vec regionlabel) {
    
    arma::vec Kw_vec = vec(Kw);
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
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        loggamma = R::dgamma(phi, aphi, 1/bphi, 1);
        // loggamma = R::dgamma(1 / phi, aphi, 1 / bphi, 1) - 2 * log(phi);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi
        );
    } else {
        K = bezier2_cpp(distmat, nu, phi);
        Kw = K * w;
        Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        loggamma = R::dgamma(phi, aphi, 1/bphi, 1);
        // loggamma = R::dgamma(1 / phi, aphi, 1 / bphi, 1) - 2 * log(phi);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi,
            Rcpp::Named("K")     = K,
            Rcpp::Named("Kw")    = Kw
        );
    }
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_main_cpp(arma::mat Yvec, arma::mat Xr, 
                       arma::vec w, arma::vec ga, double phi, double psi,
                       int nmcmc, int burn, int thin, 
                       bool adapt, List tunelst, List keeplst, 
                       List tunepsi, List keeppsi, int tunelen, 
                       arma::vec regionlabel, double nu, 
                       double atau, double btau,
                       double aphi, double bphi, 
                       arma::vec agvec, arma::mat distmat, bool iscplx,
                       Function get_tune, arma::mat draws) {
    

    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec.n_rows;
    const int n = Yvec.n_cols;
    // const int p = 1;
    int G = regionlabel.max();
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_cpp(distmat, nu, phi);
    arma::mat Kw = Ker * w;

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
    Rcpp::List keeptmppsi = keeppsi;
    int Tb = tunelen;  //# frequency of adaptive tuning
    arma::vec count = zeros<vec>(25);
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    if (iscplx) {
        Rcout << "CV-KC Begins" << endl;
    } else {
        Rcout << "MO-KC Begins" << endl;
    }
    
    for (int i = 0; i < nmcmc; ++i) {
        
        // print iteration

        
        // # Update tuning parameter of psi
        // #------------------------
        if(adapt == TRUE & fmod(i , Tb) == 0) {  
            // # Adaptive tuning
            // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
            keeptmppsi["psi"] = as<double>(keeptmppsi["psi"]) / Tb;
            tunepsi["psi"] = callgettune(tunepsi["psi"], keeptmppsi["psi"], i, get_tune);
            keeptmppsi["psi"] = 0;
        }
        
        // # Update parameters
        // #----------------------------
        //         
        // # ========== sample w =============
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
            
            // # use random walk proposal newS = S + N(0, c) to update S
            double wg_new = rnorm(1, w[g], tunelst[g])[0];
                
            List logcon_new_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau, 
                                                     Ker, Kw, g, agvec, 
                                                     FALSE, regionlabel);
            List logcon_old_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau,
                                                     Ker, Kw, g, agvec,
                                                     TRUE, regionlabel);
                
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
                as<double>(logcon_old_cpp["value"]))) {
                    w[g] = wg_new;
                    Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                    keep[g] = as<int>(keep[g]) + 1;
                    // Rcout << "\r" << "keep[g] " << as<int>(keep[g]) << endl;
                    count[g] = count[g] + 1;
                    // Rcout << "\r" << "count " << count << endl;
                    keeptmp[g] = as<int>(keeptmp[g]) + 1;
            }
        }
        
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            // Rcout << "\r" << "v: " << v + 1 << endl;
//             List MM = marginal_yga_cpp(vectorise(Yvec.row(v)), 
//                                        vectorise(XtY.col(v)), XtX);
//             ga[v] = ga_update_pc1_cpp(v, as<double>(MM["M0"]), 
//                           as<double>(MM["M1"]), 
//                           Kw, n, regionlabel, agvec, iscplx);

            ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw, n, regionlabel, agvec, iscplx);
        }
        
        // # ============ sample psi and hence phi ==============
        double new_psi = rnorm(1, psi, tunepsi["psi"])[0]; 
            
        List logcon_new_psi_cpp = logcond_psi_cpp(new_psi, w, Kw, aphi, bphi, 
                                              distmat, nu, ga, agvec, FALSE, 
                                              regionlabel);
        List logcon_old_psi_cpp = logcond_psi_cpp(psi, w, Kw, aphi, bphi, 
                                              distmat, nu, ga, agvec, TRUE, 
                                              regionlabel);
            
            
        if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) - 
            as<double>(logcon_old_psi_cpp["value"]))) {
            psi = new_psi;
            phi = exp(psi);
            Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
            Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
            keeppsi["psi"] = as<int>(keeppsi["psi"]) + 1;
            keeptmppsi["psi"] = as<int>(keeptmppsi["psi"]) + 1;
        } 

        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                draws(d - 1, draws.n_cols - 1) = phi;
//                 Gamma.row(d) = ga.t();
//                 W.row(d) = w.t();
//                 Phi.row(d) <- phi;
            } 
        }
        if (fmod(i , 1000) == 0) {
            Rcout << "\r" << "MCMC Iter: " << i + 1 << endl;
        }
    }
    // draws = join_rows(join_rows(Gamma, W), Phi);
        
    Rcout << "Done! " << endl;
    
    // # Write output
    // #--------------
    return List::create(
        _["draws"] = draws, 
        _["accept_keep"] = count, 
        _["accept_keep_psi"] = keeppsi, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_ar_main_cpp(arma::mat Yvec, arma::mat Xr, arma::vec arcoefvec, 
                       arma::vec w, arma::vec ga, double phi, double psi,
                       int nmcmc, int burn, int thin, 
                       bool adapt, List tunelst, List keeplst, 
                       List tunepsi, List keeppsi, int tunelen, 
                       arma::vec regionlabel, double nu, 
                       double atau, double btau,
                       double aphi, double bphi, 
                       arma::vec agvec, arma::mat distmat, bool iscplx,
                       Function get_tune, arma::mat draws) {
    
    
    // #-------------------------------
    // # Here create/declare some objects/variables that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec.n_rows;
    const int n = Yvec.n_cols;
    // const int p = 1;
    int G = regionlabel.max();
    double M_ga0;
    double M_ga1;
    arma::mat XtLaminv_ar;
    // arma::mat XtX_ar;
    arma::mat XtY_ar;
    arma::vec Yvecv;
    arma::vec XtYv_ar;
    arma::mat Lambda_inv;
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_cpp(distmat, nu, phi);
    arma::mat Kw = Ker * w;

    
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    Rcpp::List keep = keeplst;
    Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    
    // # for psi
    // # keep.psi <- list(psi = 0)
    Rcpp::List keeptmppsi = keeppsi;
    int Tb = tunelen;  //# frequency of adaptive tuning
    arma::vec count = zeros<vec>(25);
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    if (iscplx) {
        Rcout << "CV-KC-AR Begins" << endl;
    } else {
        Rcout << "MO-KC-AR Begins" << endl;
    }
    
    for (int i = 0; i < nmcmc; ++i) {
        
        // print iteration
        
        
        // # Update tuning parameter of psi
        // #------------------------
        if(adapt == TRUE & fmod(i , Tb) == 0) {  
            // # Adaptive tuning
            // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
            keeptmppsi["psi"] = as<double>(keeptmppsi["psi"]) / Tb;
            tunepsi["psi"] = callgettune(tunepsi["psi"], keeptmppsi["psi"], i, get_tune);
            keeptmppsi["psi"] = 0;
        }
        
        // # Update parameters
        // #----------------------------
        //         
        // # ========== sample w =============
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
            
            // # use random walk proposal newS = S + N(0, c) to update S
            double wg_new = rnorm(1, w[g], tunelst[g])[0];
            
            List logcon_new_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau, 
                                                     Ker, Kw, g, agvec, 
                                                     FALSE, regionlabel);
            List logcon_old_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau,
                                                     Ker, Kw, g, agvec,
                                                     TRUE, regionlabel);
            
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
                as<double>(logcon_old_cpp["value"]))) {
                w[g] = wg_new;
                Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                keep[g] = as<int>(keep[g]) + 1;
                // Rcout << "\r" << "keep[g] " << as<int>(keep[g]) << endl;
                count[g] = count[g] + 1;
                // Rcout << "\r" << "count " << count << endl;
                keeptmp[g] = as<int>(keeptmp[g]) + 1;
            }
        }
        
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            // Rcout << "\r" << "v: " << v + 1 << endl;
            //             List MM = marginal_yga_cpp(vectorise(Yvec.row(v)), 
            //                                        vectorise(XtY.col(v)), XtX);
            //             ga[v] = ga_update_pc1_cpp(v, as<double>(MM["M0"]), 
            //                           as<double>(MM["M1"]), 
            //                           Kw, n, regionlabel, agvec, iscplx);

            // List MM = marginal_yga_ar_cpp(Yvec, Xr, arcoefvec);
            Lambda_inv = get_ar1_inv_mat_new2(arcoefvec[v], n);
            Yvecv = vectorise(Yvec.row(v));
            XtLaminv_ar = Xr.t() * Lambda_inv;
            XtYv_ar = XtLaminv_ar * Yvecv;
            M_ga0 = (Yvecv.t() * Lambda_inv * Yvecv).eval()(0,0);
            M_ga1 = M_ga0 - (XtYv_ar.t() * inv_sympd(XtLaminv_ar * Xr) * XtYv_ar).eval()(0,0);
            // arma::vec M0 = MM["M0"];
            // arma::vec M1 = MM["M1"];
            // ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw, n, regionlabel, agvec, iscplx);
            ga[v] = ga_update_pc1_cpp(v, M_ga0, M_ga1, Kw, n, regionlabel, agvec, iscplx);
        }
        
        // # ============ sample psi and hence phi ==============
        double new_psi = rnorm(1, psi, tunepsi["psi"])[0]; 
        
        List logcon_new_psi_cpp = logcond_psi_cpp(new_psi, w, Kw, aphi, bphi, 
                                                  distmat, nu, ga, agvec, FALSE, 
                                                  regionlabel);
        List logcon_old_psi_cpp = logcond_psi_cpp(psi, w, Kw, aphi, bphi, 
                                                  distmat, nu, ga, agvec, TRUE, 
                                                  regionlabel);
        
        
        if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) - 
            as<double>(logcon_old_psi_cpp["value"]))) {
            psi = new_psi;
            phi = exp(psi);
            Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
            Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
            keeppsi["psi"] = as<int>(keeppsi["psi"]) + 1;
            keeptmppsi["psi"] = as<int>(keeptmppsi["psi"]) + 1;
        } 
        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                draws(d - 1, draws.n_cols - 1) = phi;
                //                 Gamma.row(d) = ga.t();
                //                 W.row(d) = w.t();
                //                 Phi.row(d) <- phi;
            } 
        }
        if (fmod(i , 1000) == 0) {
            Rcout << "\r" << "MCMC Iter: " << i + 1 << endl;
        }
    }
    // draws = join_rows(join_rows(Gamma, W), Phi);
    
    Rcout << "\n" << "Done! " << endl;
    
    // # Write output
    // #--------------
    return List::create(
        _["draws"] = draws, 
        _["accept_keep"] = count, 
        _["accept_keep_psi"] = keeppsi, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_main_cpp_vec(arma::mat Yvec, arma::mat Xr, 
                       arma::vec w, arma::vec ga, double phi, double psi,
                       int nmcmc, int burn, int thin, 
                       bool adapt, arma::vec tunevec, arma::vec keepvec, 
                       double tunepsi, double keeppsi, int tunelen, 
                       arma::vec regionlabel, int nu, 
                       double atau, double btau,
                       double aphi, double bphi, 
                       arma::vec agvec, arma::mat distmat, bool iscplx,
                       Function get_tune, arma::mat draws) {
    
    
    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec.n_rows;
    const int n = Yvec.n_cols;
    // const int p = 1;
    int G = regionlabel.max();
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_cpp(distmat, nu, phi);
    arma::mat Kw = Ker * w;
    
    List MM = marginal_yga_cpp(Yvec, XtY, XtX);
    arma::vec M0 = MM["M0"];
    arma::vec M1 = MM["M1"];
    
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    // Rcpp::List keep = keeplst;
    // Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    // Rcpp::List keep = keepvec;
    arma::vec keep = keepvec;
    arma::vec keeptmp = keepvec;  // track MH accpetance rate for adaptive tuning
    // # for psi
    // # keep.psi <- list(psi = 0)
    double keeptmppsi = keeppsi;
    int Tb = tunelen;  //# frequency of adaptive tuning
    // arma::vec count = zeros<vec>(25);
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    Rcout << "KC-S Begins" << endl;
    
    for (int i = 0; i < nmcmc; ++i) {
        
        // print iteration
        
        
        // # Update tuning parameter of psi
        // #------------------------
        if(adapt == TRUE & fmod(i , Tb) == 0) {  
            // # Adaptive tuning
            // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
            keeptmppsi = keeptmppsi / Tb;
            tunepsi = callgettune(tunepsi, keeptmppsi, i, get_tune);
            keeptmppsi = 0;
        }
        
        // # Update parameters
        // #----------------------------
        //         
        // # ========== sample w =============
        for (int g = 0; g < G; ++g) {
            // Rcout << "\r" << "g: " << g + 1 << endl;
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmp[g] = keeptmp[g] / Tb;
                tunevec[g] = callgettune(tunevec[g], keeptmp[g], i, get_tune);
                keeptmp[g] = 0;
            }
            
            // #----------------------------
            // # Random walk Normal proposal 
            // #----------------------------
            
            // # use random walk proposal newS = S + N(0, c) to update S
            double wg_new = rnorm(1, w[g], tunevec[g])[0];
            
            List logcon_new_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau, 
                                                     Ker, Kw, g, agvec, 
                                                     FALSE, regionlabel);
            List logcon_old_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau,
                                                     Ker, Kw, g, agvec,
                                                     TRUE, regionlabel);
            
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
                as<double>(logcon_old_cpp["value"]))) {
                w[g] = wg_new;
                Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                keep[g] = keep[g] + 1;
                // Rcout << "\r" << "keep[g] " << as<int>(keep[g]) << endl;
                // count[g] = count[g] + 1;
                // Rcout << "\r" << "count " << count << endl;
                keeptmp[g] = keeptmp[g] + 1;
            }
        }
        
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            // Rcout << "\r" << "v: " << v + 1 << endl;
            //             List MM = marginal_yga_cpp(vectorise(Yvec.row(v)), 
            //                                        vectorise(XtY.col(v)), XtX);
            //             ga[v] = ga_update_pc1_cpp(v, as<double>(MM["M0"]), 
            //                           as<double>(MM["M1"]), 
            //                           Kw, n, regionlabel, agvec, iscplx);
            
            ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw, n, regionlabel, agvec, iscplx);
        }
        
        // # ============ sample psi and hence phi ==============
        double new_psi = rnorm(1, psi, tunepsi)[0]; 
        
        List logcon_new_psi_cpp = logcond_psi_cpp(new_psi, w, Kw, aphi, bphi, 
                                                  distmat, nu, ga, agvec, FALSE, 
                                                  regionlabel);
        List logcon_old_psi_cpp = logcond_psi_cpp(psi, w, Kw, aphi, bphi, 
                                                  distmat, nu, ga, agvec, TRUE, 
                                                  regionlabel);
        
        
        if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) - 
            as<double>(logcon_old_psi_cpp["value"]))) {
            psi = new_psi;
            phi = exp(psi);
            Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
            Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
            keeppsi = keeppsi + 1;
            keeptmppsi = keeptmppsi + 1;
        } 
        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                draws(d - 1, draws.n_cols - 1) = phi;
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
    
    Rcout << "Done! " << endl;
    
    // # Write output
    // #--------------
    return List::create(
        _["draws"] = draws, 
        _["accept_keep"] = keep, 
        _["accept_keep_psi"] = keeppsi, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    for (int i=0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
        out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
    }  
    
    if (logd == false) {
        out = exp(out);
    }
    return(out);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_block_main_cpp(arma::mat Yvec, arma::mat Xr, 
                       arma::vec w, arma::vec ga, double phi, double psi,
                       int nmcmc, int burn, int thin, 
                       bool adapt, List tunelst, List keeplst, 
                       List tunepsi, List keeppsi, int tunelen, 
                       arma::vec regionlabel, int nu, 
                       double atau, double btau,
                       double aphi, double bphi, 
                       arma::vec agvec, arma::mat distmat, bool iscplx,
                       Function get_tune, arma::mat draws) {
    
    
    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec.n_rows;
    const int n = Yvec.n_cols;
    // const int p = 1;
    int G = regionlabel.max();
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_cpp(distmat, nu, phi);
    arma::mat Kw = Ker * w;
    
    List MM = marginal_yga_cpp(Yvec, XtY, XtX);
    arma::vec M0 = MM["M0"];
    arma::vec M1 = MM["M1"];
    arma::mat Wmat(nmcmc, G);
    arma::mat w_new;
    arma::uvec widx;
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    Rcpp::List keep = keeplst;
    Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    
    // # for psi
    // # keep.psi <- list(psi = 0)
    Rcpp::List keeptmppsi = keeppsi;
    int Tb = tunelen;  //# frequency of adaptive tuning
    
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    // Begin algorithm
    if (iscplx) {
        Rcout << "CV-KC (block) Begins" << endl;
    } else {
        Rcout << "MO-KC (block) Begins" << endl;
    }
    
    for (uword i = 0; i < nmcmc; ++i) {
        
        // print iteration
        
        
        // # Update tuning parameter of psi
        // #------------------------
        if(adapt == TRUE & fmod(i , Tb) == 0) {  
            // # Adaptive tuning
            // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
            keeptmppsi["psi"] = as<double>(keeptmppsi["psi"]) / Tb;
            tunepsi["psi"] = callgettune(tunepsi["psi"], keeptmppsi["psi"], i, get_tune);
            keeptmppsi["psi"] = 0;
        }
        
        // # Update parameters
        // #----------------------------
        //         
        // # ========== sample w =============
        int sqrtG = sqrt(G);
        // int sqrtG = ceil(sqrt(G)); 
        for(int j = 0; j < sqrtG; ++j) {
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmp[j] = as<double>(keeptmp[j]) / Tb;
                tunelst[j] = callgettune(tunelst[j], keeptmp[j], i, get_tune);
                keeptmp[j] = 0;
            }
            
            // widx = sqrtG * j + linspace<uvec>(0, sqrtG - 1, sqrtG);
            widx = sqrtG * (j - 1) + linspace<uvec>(0, sqrtG - 1, sqrtG);
            // Rcout << " widx " << widx.n_elem << endl;
            if (i > 1000) {
                // cov(W[(i-1000):(i-1), w.idx])
                arma::mat Sigma_w = cov(Wmat(span(i - 1000, i - 1), span(widx[0], widx[sqrtG - 1])));
                // Sigma_w = (Sigma_w + Sigma_w.t()) / 2;
                arma::mat sigma = as<double>(tunelst[j]) * Sigma_w;
                w_new = mvrnormArma(1, w.elem(widx), sigma);
                // Rcout << "w_new " << w_new.n_elem << endl;
            } else {
                w_new = mvrnormArma(1, w.elem(widx), eye<mat>(sqrtG, sqrtG));
                // Rcout << "w_new " << w_new.n_elem << endl;
            }
            
            // arma::vec w, arma::vec w_new, arma::vec ga, 
            // double atau, double btau, arma::mat K, arma::mat Kw, 
            // arma::uvec widx, arma::vec agvec, 
            // bool isold, arma::vec regionlabel
            arma::vec w_new_vec = vectorise(w_new);
            // Rcout << "w_new_vec" << w_new_vec.n_elem << endl;
                
            List logcon_new_cpp = logcond_latent_block_cpp(w, w_new_vec, ga, atau, btau, 
                                                     Ker, Kw, widx, agvec, 
                                                     FALSE, regionlabel);
            List logcon_old_cpp = logcond_latent_block_cpp(w, w_new_vec, ga, atau, btau,
                                                     Ker, Kw, widx, agvec,
                                                     TRUE, regionlabel);
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
                as<double>(logcon_old_cpp["value"]))) {
                w.elem(widx) = w_new;
                Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                keep[j] = as<int>(keep[j]) + 1;
                keeptmp[j] = as<int>(keeptmp[j]) + 1;
            }
            Wmat(i, span(widx[0], widx[sqrtG - 1])) = mat(w.elem(widx)).t();
        }
        
        
        // // The rest of regions
        // if(adapt == TRUE & fmod(i , Tb) == 0) {
        //     // # Adaptive tuning
        //     keeptmp[sqrtG - 1] = keeptmp[sqrtG - 1] / Tb;
        //     tunevec[sqrtG - 1] = callgettune(tunevec[sqrtG - 1], 
        //                                      keeptmp[sqrtG - 1], i, get_tune);
        //     keeptmp[sqrtG - 1] = 0;
        // }
        // 
        // // #----------------------------
        // // # Random walk Normal proposal 
        // // #----------------------------
        // 
        // // # use random walk proposal newS = S + N(0, c) to update S
        // double wg_new = rnorm(1, w[sqrtG - 1], tunevec[sqrtG - 1])[0];
        // 
        // List logcon_new_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau, 
        //                                          Ker, Kw, g, agvec, 
        //                                          FALSE, regionlabel);
        // List logcon_old_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau,
        //                                          Ker, Kw, g, agvec,
        //                                          TRUE, regionlabel);
        // 
        // if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
        //     as<double>(logcon_old_cpp["value"]))) {
        //     w[g] = wg_new;
        //     Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
        //     keep[g] = keep[g] + 1;
        //     // Rcout << "\r" << "keep[g] " << as<int>(keep[g]) << endl;
        //     // count[g] = count[g] + 1;
        //     // Rcout << "\r" << "count " << count << endl;
        //     keeptmp[g] = keeptmp[g] + 1;
        // }

        
        
        // for (int g = 0; g < G; ++g) {
        //     // Rcout << "\r" << "g: " << g + 1 << endl;
        //     // # Update tuning parameter
        //     // #-------------------------
        //     if(adapt == TRUE & fmod(i , Tb) == 0) {
        //         // # Adaptive tuning
        //         keeptmp[g] = as<double>(keeptmp[g]) / Tb;
        //         tunelst[g] = callgettune(tunelst[g], keeptmp[g], i, get_tune);
        //         keeptmp[g] = 0;
        //     }
        //     
        //     // #----------------------------
        //     // # Random walk Normal proposal 
        //     // #----------------------------
        //     
        //     // # use random walk proposal newS = S + N(0, c) to update S
        //     double wg_new = rnorm(1, w[g], tunelst[g])[0];
        //     
        //     List logcon_new_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau, 
        //                                              Ker, Kw, g, agvec, 
        //                                              FALSE, regionlabel);
        //     List logcon_old_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau,
        //                                              Ker, Kw, g, agvec,
        //                                              TRUE, regionlabel);
        //     
        //     if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
        //         as<double>(logcon_old_cpp["value"]))) {
        //         w[g] = wg_new;
        //         Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
        //         keep[g] = as<int>(keep[g]) + 1;
        //         keeptmp[g] = as<int>(keeptmp[g]) + 1;
        //     }
        // }
//         for (j in 1:5) {
//             if(adapt == TRUE & i %% Tb == 0) {
// # Adaptive tuning
//                 keep.tmp.w[[j]] <- keep.tmp.w[[j]] / Tb
//                 tune.w[[j]] <- get.tune(tune.w[[j]], keep.tmp.w[[j]], i)
//                 keep.tmp.w[[j]] <- 0
//             }
//             w.idx <- (5*(j-1) + 1:5)
// # block move using empirical covariance structure
//                 if (i > 1000) {
//                     Sigma.w <- cov(W[(i-1000):(i-1), w.idx])
//                     Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w)$mat)
// # print(round(Sigma.w, 2))
//                     w_new <- as.vector(mvnfast::rmvn(1, w[w.idx], tune.w[[j]] * Sigma.w))
//                 } else {
//                     w_new <- as.vector(mvnfast::rmvn(1, w[w.idx], tune.w[[j]] * diag(5)))
//                 }
//                 logcon.new <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
//                                                     ag.vec, is.old = FALSE, 
//                                                     regionlabel = regionlabel)
//                     logcon.old <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
//                                                         ag.vec, is.old = TRUE, 
//                                                         regionlabel = regionlabel)
//                     if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
//                         w[w.idx] <- w_new
//                         Kw <- logcon.new$Kw
//                         keep[[j]] <- keep[[j]] + 1
//                         keep.tmp.w[[j]] <- keep.tmp.w[[j]] + 1
//                     }
//                     W[i, w.idx] <- w_new
//         }
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            // Rcout << "\r" << "v: " << v + 1 << endl;
            //             List MM = marginal_yga_cpp(vectorise(Yvec.row(v)), 
            //                                        vectorise(XtY.col(v)), XtX);
            //             ga[v] = ga_update_pc1_cpp(v, as<double>(MM["M0"]), 
            //                           as<double>(MM["M1"]), 
            //                           Kw, n, regionlabel, agvec, iscplx);
            
            ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw, n, regionlabel, agvec, iscplx);
        }
        
        // # ============ sample psi and hence phi ==============
        double new_psi = rnorm(1, psi, tunepsi["psi"])[0]; 
        
        List logcon_new_psi_cpp = logcond_psi_cpp(new_psi, w, Kw, aphi, bphi, 
                                                  distmat, nu, ga, agvec, FALSE, 
                                                  regionlabel);
        List logcon_old_psi_cpp = logcond_psi_cpp(psi, w, Kw, aphi, bphi, 
                                                  distmat, nu, ga, agvec, TRUE, 
                                                  regionlabel);
        
        
        if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) - 
            as<double>(logcon_old_psi_cpp["value"]))) {
            psi = new_psi;
            phi = exp(psi);
            Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
            Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
            keeppsi["psi"] = as<int>(keeppsi["psi"]) + 1;
            keeptmppsi["psi"] = as<int>(keeptmppsi["psi"]) + 1;
        } 
        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                draws(d - 1, draws.n_cols - 1) = phi;
                //                 Gamma.row(d) = ga.t();
                //                 W.row(d) = w.t();
                //                 Phi.row(d) <- phi;
            } 
        }
        if (fmod(i , 1000) == 0) {
            Rcout << "MCMC Iter: " << i + 1 << "\r" << std::flush;
        }
    }
    // draws = join_rows(join_rows(Gamma, W), Phi);
    
    Rcout << "\n" << "Done! " << endl;
    
    // # Write output
    // #--------------
    return List::create(
        _["draws"] = draws, 
        _["accept_keep"] = keep, 
        _["accept_keep_psi"] = keeppsi, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_ar_block_main_cpp(arma::mat Yvec, arma::mat Xr, arma::vec arcoefvec, 
                             arma::vec w, arma::vec ga, double phi, double psi,
                             int nmcmc, int burn, int thin, 
                             bool adapt, List tunelst, List keeplst, 
                             List tunepsi, List keeppsi, int tunelen, 
                             arma::vec regionlabel, int nu, 
                             double atau, double btau,
                             double aphi, double bphi, 
                             arma::vec agvec, arma::mat distmat, bool iscplx,
                             Function get_tune, arma::mat draws) {
    
    
    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec.n_rows;
    const int n = Yvec.n_cols;
    // const int p = 1;
    int G = regionlabel.max();
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();
    double M_ga0;
    double M_ga1;
    arma::mat XtLaminv_ar;
    // arma::mat XtX_ar;
    arma::mat XtY_ar;
    arma::vec Yvecv;
    arma::vec XtYv_ar;
    arma::mat Lambda_inv;
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_cpp(distmat, nu, phi);
    arma::mat Kw = Ker * w;
    
    // List MM = marginal_yga_cpp(Yvec, XtY, XtX);
    // arma::vec M0 = MM["M0"];
    // arma::vec M1 = MM["M1"];
    arma::mat Wmat(nmcmc, G);
    arma::mat w_new;
    arma::uvec widx;
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    Rcpp::List keep = keeplst;
    Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    
    // # for psi
    // # keep.psi <- list(psi = 0)
    Rcpp::List keeptmppsi = keeppsi;
    int Tb = tunelen;  //# frequency of adaptive tuning
    
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    // Begin algorithm
    if (iscplx) {
        Rcout << "CV-KC-AR (block) Begins" << endl;
    } else {
        Rcout << "MO-KC-AR (block) Begins" << endl;
    }
    
    for (uword i = 0; i < nmcmc; ++i) {
        
        // print iteration
        
        
        // # Update tuning parameter of psi
        // #------------------------
        if(adapt == TRUE & fmod(i , Tb) == 0) {  
            // # Adaptive tuning
            // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
            keeptmppsi["psi"] = as<double>(keeptmppsi["psi"]) / Tb;
            tunepsi["psi"] = callgettune(tunepsi["psi"], keeptmppsi["psi"], i, get_tune);
            keeptmppsi["psi"] = 0;
        }
        
        // # Update parameters
        // #----------------------------
        //         
        // # ========== sample w =============
        int sqrtG = sqrt(G); 
        for(int j = 0; j < sqrtG; ++j) {
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmp[j] = as<double>(keeptmp[j]) / Tb;
                tunelst[j] = callgettune(tunelst[j], keeptmp[j], i, get_tune);
                keeptmp[j] = 0;
            }
            
            widx = sqrtG * j + linspace<uvec>(0, sqrtG - 1, sqrtG);
            // Rcout << " widx " << widx.n_elem << endl;
            if (i > 1000) {
                // cov(W[(i-1000):(i-1), w.idx])
                arma::mat Sigma_w = cov(Wmat(span(i - 1000, i - 1), span(widx[0], widx[sqrtG - 1])));
                // Sigma_w = (Sigma_w + Sigma_w.t()) / 2;
                arma::mat sigma = as<double>(tunelst[j]) * Sigma_w;
                w_new = mvrnormArma(1, w.elem(widx), sigma);
                // Rcout << "w_new " << w_new.n_elem << endl;
            } else {
                w_new = mvrnormArma(1, w.elem(widx), eye<mat>(sqrtG, sqrtG));
                // Rcout << "w_new " << w_new.n_elem << endl;
            }
            
            // arma::vec w, arma::vec w_new, arma::vec ga, 
            // double atau, double btau, arma::mat K, arma::mat Kw, 
            // arma::uvec widx, arma::vec agvec, 
            // bool isold, arma::vec regionlabel
            arma::vec w_new_vec = vectorise(w_new);
            // Rcout << "w_new_vec" << w_new_vec.n_elem << endl;
            
            List logcon_new_cpp = logcond_latent_block_cpp(w, w_new_vec, ga, atau, btau, 
                                                           Ker, Kw, widx, agvec, 
                                                           FALSE, regionlabel);
            List logcon_old_cpp = logcond_latent_block_cpp(w, w_new_vec, ga, atau, btau,
                                                           Ker, Kw, widx, agvec,
                                                           TRUE, regionlabel);
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
                as<double>(logcon_old_cpp["value"]))) {
                w.elem(widx) = w_new;
                Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                keep[j] = as<int>(keep[j]) + 1;
                keeptmp[j] = as<int>(keeptmp[j]) + 1;
            }
            Wmat(i, span(widx[0], widx[sqrtG - 1])) = mat(w.elem(widx)).t();
        }
        
        
        
        
        
        // for (int g = 0; g < G; ++g) {
        //     // Rcout << "\r" << "g: " << g + 1 << endl;
        //     // # Update tuning parameter
        //     // #-------------------------
        //     if(adapt == TRUE & fmod(i , Tb) == 0) {
        //         // # Adaptive tuning
        //         keeptmp[g] = as<double>(keeptmp[g]) / Tb;
        //         tunelst[g] = callgettune(tunelst[g], keeptmp[g], i, get_tune);
        //         keeptmp[g] = 0;
        //     }
        //     
        //     // #----------------------------
        //     // # Random walk Normal proposal 
        //     // #----------------------------
        //     
        //     // # use random walk proposal newS = S + N(0, c) to update S
        //     double wg_new = rnorm(1, w[g], tunelst[g])[0];
        //     
        //     List logcon_new_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau, 
        //                                              Ker, Kw, g, agvec, 
        //                                              FALSE, regionlabel);
        //     List logcon_old_cpp = logcond_latent_cpp(w, wg_new, ga, atau, btau,
        //                                              Ker, Kw, g, agvec,
        //                                              TRUE, regionlabel);
        //     
        //     if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
        //         as<double>(logcon_old_cpp["value"]))) {
        //         w[g] = wg_new;
        //         Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
        //         keep[g] = as<int>(keep[g]) + 1;
        //         keeptmp[g] = as<int>(keeptmp[g]) + 1;
        //     }
        // }
        //         for (j in 1:5) {
        //             if(adapt == TRUE & i %% Tb == 0) {
        // # Adaptive tuning
        //                 keep.tmp.w[[j]] <- keep.tmp.w[[j]] / Tb
        //                 tune.w[[j]] <- get.tune(tune.w[[j]], keep.tmp.w[[j]], i)
        //                 keep.tmp.w[[j]] <- 0
        //             }
        //             w.idx <- (5*(j-1) + 1:5)
        // # block move using empirical covariance structure
        //                 if (i > 1000) {
        //                     Sigma.w <- cov(W[(i-1000):(i-1), w.idx])
        //                     Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w)$mat)
        // # print(round(Sigma.w, 2))
        //                     w_new <- as.vector(mvnfast::rmvn(1, w[w.idx], tune.w[[j]] * Sigma.w))
        //                 } else {
        //                     w_new <- as.vector(mvnfast::rmvn(1, w[w.idx], tune.w[[j]] * diag(5)))
        //                 }
        //                 logcon.new <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
        //                                                     ag.vec, is.old = FALSE, 
        //                                                     regionlabel = regionlabel)
        //                     logcon.old <- logcond_latent_block2(w, w_new, ga, a.tau, b.tau, K = Ker, Kw, w.idx, 
        //                                                         ag.vec, is.old = TRUE, 
        //                                                         regionlabel = regionlabel)
        //                     if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
        //                         w[w.idx] <- w_new
        //                         Kw <- logcon.new$Kw
        //                         keep[[j]] <- keep[[j]] + 1
        //                         keep.tmp.w[[j]] <- keep.tmp.w[[j]] + 1
        //                     }
        //                     W[i, w.idx] <- w_new
        //         }
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            // Rcout << "\r" << "v: " << v + 1 << endl;
            //             List MM = marginal_yga_cpp(vectorise(Yvec.row(v)), 
            //                                        vectorise(XtY.col(v)), XtX);
            //             ga[v] = ga_update_pc1_cpp(v, as<double>(MM["M0"]), 
            //                           as<double>(MM["M1"]), 
            //                           Kw, n, regionlabel, agvec, iscplx);
            Lambda_inv = get_ar1_inv_mat_new2(arcoefvec[v], n);
            Yvecv = vectorise(Yvec.row(v));
            XtLaminv_ar = Xr.t() * Lambda_inv;
            XtYv_ar = XtLaminv_ar * Yvecv;
            M_ga0 = (Yvecv.t() * Lambda_inv * Yvecv).eval()(0,0);
            M_ga1 = M_ga0 - (XtYv_ar.t() * inv_sympd(XtLaminv_ar * Xr) * XtYv_ar).eval()(0,0);
            // ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw, n, regionlabel, agvec, iscplx);
            ga[v] = ga_update_pc1_cpp(v, M_ga0, M_ga1, Kw, n, regionlabel, agvec, iscplx);
        }
        
        // # ============ sample psi and hence phi ==============
        double new_psi = rnorm(1, psi, tunepsi["psi"])[0]; 
        
        List logcon_new_psi_cpp = logcond_psi_cpp(new_psi, w, Kw, aphi, bphi, 
                                                  distmat, nu, ga, agvec, FALSE, 
                                                  regionlabel);
        List logcon_old_psi_cpp = logcond_psi_cpp(psi, w, Kw, aphi, bphi, 
                                                  distmat, nu, ga, agvec, TRUE, 
                                                  regionlabel);
        
        
        if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) - 
            as<double>(logcon_old_psi_cpp["value"]))) {
            psi = new_psi;
            phi = exp(psi);
            Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
            Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
            keeppsi["psi"] = as<int>(keeppsi["psi"]) + 1;
            keeptmppsi["psi"] = as<int>(keeptmppsi["psi"]) + 1;
        } 
        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                draws(d - 1, draws.n_cols - 1) = phi;
                //                 Gamma.row(d) = ga.t();
                //                 W.row(d) = w.t();
                //                 Phi.row(d) <- phi;
            } 
        }
        if (fmod(i , 1000) == 0) {
            Rcout << "MCMC Iter: " << i + 1 << "\r" << std::flush;
        }
    }
    // draws = join_rows(join_rows(Gamma, W), Phi);
    
    Rcout << "\n" << "Done! " << endl;
    
    // # Write output
    // #--------------
    return List::create(
        _["draws"] = draws, 
        _["accept_keep"] = keep, 
        _["accept_keep_psi"] = keeppsi, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}








