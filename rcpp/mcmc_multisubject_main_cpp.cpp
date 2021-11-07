#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>
#include <string.h>

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
double callgettune(double tune, int keep, int k, Function f) {
    double res;
    res = as<double>(f(tune, keep, k));
    return res;
}



// [[Rcpp::export]]
arma::mat callrWishart(int df, arma::mat Sigma, Function rWishart) {
    // stats::rWishart(1, df = ns + D, Sigma = Sigma)
    arma::mat res;
    res = as<mat>(rWishart(1, df, Sigma));
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
arma::mat GaussianKer_cpp(arma::mat distmat, double phi) {
    int N = distmat.n_rows;
    int M = distmat.n_cols;
    // arma::mat K = zeros(N, M);
    arma::mat K = exp(- pow(distmat, 2) / (2 * phi));
    return (K);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int ga_update_multi_cpp(int s, int v, arma::mat ga_mat, 
                        arma::mat M0_mat, arma::mat M1_mat, arma::mat Kw, int n, 
                      arma::vec region_conn, arma::mat a_mat, bool iscplx) {
    // # a_mat: ns x D
    // # ga_mat: ns x N.voxels
    // Rcout << "Kw.n_rows" << Kw.n_rows << endl;
    // Rcout << "a_mat.row(s)" << a_mat.row(s) << endl;
    arma::vec a_vec = vectorise(a_mat.row(s));
    // Rcout << "a_mat.row(s)" << a_mat.row(s) << endl;
    arma::vec Kw_vec = vectorise(Kw);
    int N = Kw.n_rows;
    // Rcout << "Kw.n_rows" << Kw.n_rows << endl;
    int ns = ga_mat.n_rows;
    // Rcout << "ga_mat.n_rows" << ga_mat.n_rows << endl;
    // int g = regionlabel[v];
    double p_ber;
    int C;
    int po;
    double p_star;
    // Rcout << "region_conn[v] - 1" << region_conn[v] - 1 << endl;
    p_ber = 1 / (1 + exp(-(a_vec[region_conn[v] - 1] + Kw_vec[v])));
    if (p_ber == 0) p_ber = 1e-8;
    if (p_ber == 1) p_ber = 1 - 1e-8;
    
    arma::mat Msv = M1_mat;
    arma::mat Msv1;
    arma::mat Msv0;
    // Rcout << "Kw.n_rows" << Kw.n_rows << endl;
    for (int i = 0; i < ns; i++) {
        for (int j = 0; j < N; j++) {
            if (ga_mat(i, j) == 0) {
                Msv(i, j) = M0_mat(i, j);
            }
        }
    }
    // Rcout << " " << endl;
    Msv1 = Msv;
    Msv0 = Msv;
    Msv1(s, v) = M1_mat(s, v);
    Msv0(s, v) = M0_mat(s, v);
    
    double sum1 = sum(Msv1.col(v));
    double sum0 = sum(Msv0.col(v));
    
    //     if (iscplx) {
    //         C = 1 + n;
    //         po = n;
    //     } else {
    //         C = sqrt(1 + n);
    //         po = n / 2;
    //     }
    po = ns * n / 2;
    if (iscplx) {
        C = 1 + n / 2;
    } else {
        C = sqrt(1 + n);
    }
    p_star = 1 / (1 + ((1 - p_ber) / p_ber) * C * pow(sum1 / sum0, po));
    return R::rbinom(1, p_star);
}




// ga_update_multi <- function(s, v, ga_mat, M0_mat, M1_mat, Kw, n,
//                             region_conn, a_mat, is.cplx) {
// # a_mat: ns x D
// # ga_mat: ns x N.voxels
// 
//     a_vec <- a_mat[s, ]
//     p_ber <- 1 / (1 + exp(-(a_vec[region_conn[v]] + Kw[v])))
//     if (p_ber == 0) p_ber <- 1e-8
//     if (p_ber == 1) p_ber <- 1 - 1e-8
// 
//     Msv <- ifelse(ga_mat, M1_mat, M0_mat)
//         Msv1 <- Msv
//         Msv0 <- Msv
//         Msv1[s, v] <- M1_mat[s, v]
//     Msv0[s, v] <- M0_mat[s, v]
// 
//     sum1 <- sum(Msv1[, v])
//         sum0 <- sum(Msv0[, v])
// 
//         po <- ns * n / 2
//     if (is.cplx) {
//         C <- 1 + n / 2
//     } else {
//         C <- sqrt(1 + n)
//     }
// 
//     p_star <- 1 / (1 + ((1 - p_ber) / p_ber) * C * (sum1 / sum0) ^ po)
//         return(rbinom(1, 1, p_star))
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_wg_block_multi_cpp(arma::vec w, arma::vec w_new, arma::mat ga_mat, 
                        double atau, double btau, arma::mat K, arma::mat Kw, 
                        arma::uvec widx, arma::mat a_mat, 
                        bool isold, arma::vec region_conn) {
    int G = w.n_elem;
    // double sumlogber;
    int ns = ga_mat.n_rows;
    // Rcout << "\r" << "check 0 " << endl;
    arma::vec Kw_vec = vec(Kw);
    // Rcout << "\r" << "check 00 " << endl;
    int N = Kw.n_rows;
    arma::mat p_bern_mat(ns, N);
    arma::mat logber_mat(ns, N);
    // arma::vec a_vec(N);
    double logsumber;
    // Rcout << "\r" << "check 1 " << endl;
    if(isold) {
        for (int s = 0; s < ns; s++) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern_mat(s, i) = 1 / (1 + exp(-(a_mat(s, region_conn[i] - 1) + Kw_vec[i])));
                logber_mat(s, i) = R::dbinom(ga_mat(s, i), 1, p_bern_mat(s, i), 1);
            }
        }
        logsumber = accu(logber_mat);
        double B = (G / 2 + atau) * log((0.5) * sum(square(w)) + btau);
        // Rcout << "\r" << "check 2 " << endl;
        return List::create(
            Rcpp::Named("value") = logsumber - B
        );
    } else {
        // be careful about index
        w.elem(widx) = w_new;
        // Rcout << "\r" << "check 3 " << endl;
        Kw = K * w;
        // Rcout << "\r" << "check 4 " << endl;
        Kw_vec = vec(Kw);
        for (int s = 0; s < ns; s++) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern_mat(s, i) = 1 / (1 + exp(-(a_mat(s, region_conn[i] - 1) + Kw_vec[i])));
                logber_mat(s, i) = R::dbinom(ga_mat(s, i), 1, p_bern_mat(s, i), 1);
            }
        }
        logsumber = accu(logber_mat);
        double B = (G / 2 + atau) * log((0.5) * sum(square(w)) + btau);
        // Rcout << "\r" << "check 5 " << endl;
        return List::create(
            Rcpp::Named("value") = logsumber - B,
            Rcpp::Named("Kw") = Kw
        );
    }
}

// logcond_wg_block_multi <- function(w, w_new, ga_mat, a.tau, b.tau, K, Kw, 
//                                    w.idx, a_mat, is.old = TRUE, region_conn) {
//     ns <- dim(ga_mat)[1]
//     N.voxels <- dim(ga_mat)[2]
//     if(is.old) {
//         a_mat_v <- sapply(1:length(Kw), function(x) a_mat[, region_conn[x]])
//         Kw_mat <- matrix(rep(Kw, each = ns), ns, N.voxels)
//         p_bern_mat <- 1 / (1 + exp(-(a_mat_v + Kw_mat)))
//         sumlogber <- sum(dbinom(ga_mat, size = 1, prob = p_bern_mat, log = TRUE))
//         B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
//         return(list(value = sumlogber - B))
//     } else {
//         w[w.idx] <- w_new
//         Kw <- K %*% w
//         a_mat_v <- sapply(1:length(Kw), function(x) a_mat[, region_conn[x]])
//         Kw_mat <- matrix(rep(Kw, each = ns), ns, N.voxels)
//         p_bern_mat <- 1 / (1 + exp(-(a_mat_v + Kw_mat)))
//         sumlogber <- sum(dbinom(ga_mat, size = 1, prob = p_bern_mat, log = TRUE))
//         B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
//         return(list(value = sumlogber - B, Kw = Kw))
//     }
// } 





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_wg_multi_cpp(arma::vec w, double wg_new, arma::mat ga_mat, 
                                double atau, double btau, arma::mat K, arma::mat Kw, 
                                int gidx, arma::mat a_mat, 
                                bool isold, arma::vec region_conn){
    int ns = ga_mat.n_rows;
    arma::vec Kw_vec = vec(Kw);
    int G = w.n_elem;
    int N = Kw.n_rows;
    arma::mat p_bern_mat(ns, N);
    arma::mat logber_mat(ns, N);
    // arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double sumlogber;
    // arma::vec logber(N);
    
    if(isold) {
        for (int s = 0; s < ns; s++) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern_mat(s, i) = 1 / (1 + exp(-(a_mat(s, region_conn[i] - 1) + Kw_vec[i])));
                logber_mat(s, i) = R::dbinom(ga_mat(s, i), 1, p_bern_mat(s, i), 1);
            }
        }
        sumlogber = accu(logber_mat);
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
        for (int s = 0; s < ns; s++) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern_mat(s, i) = 1 / (1 + exp(-(a_mat(s, region_conn[i] - 1) + Kw_vec[i])));
                logber_mat(s, i) = R::dbinom(ga_mat(s, i), 1, p_bern_mat(s, i), 1);
            }
        }
        sumlogber = accu(logber_mat);
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
List marginal_yga_multi_cpp(arma::cube Yvec_multi, arma::mat XtX, arma::mat Xr) {
    // dimension of Yvec_multi is Nvoxels x n x ns
    int Nvoxels = Yvec_multi.n_rows;
    int ns = Yvec_multi.n_slices;
    
    arma::mat M_ga0_mat(ns, Nvoxels);
    arma::mat M_ga1_mat(ns, Nvoxels);
    arma::vec Yvecv;
    arma::vec XtYv;
    arma::mat XtXinv = inv_sympd(XtX);

    for (int s = 0; s < ns; ++s) {
        arma::mat Yvec = Yvec_multi.slice(s);
        arma::mat XtYs = Xr.t() * Yvec.t();
        for (int v=0; v < Nvoxels; ++v) {
            Yvecv = vectorise(Yvec.row(v));
            XtYv = vectorise(XtYs.col(v));
            double M_ga0 = (Yvecv.t() * Yvecv).eval()(0,0);
            double M_ga1 = M_ga0 -(XtYv.t() * XtXinv * XtYv).eval()(0,0);
            M_ga0_mat(s, v) = M_ga0;
            M_ga1_mat(s, v) = M_ga1;
        }
    }
    return List::create(
        Rcpp::Named("M0_mat") = M_ga0_mat,
        Rcpp::Named("M1_mat") = M_ga1_mat
    );
}


// marginal_yga_multi <- function(Yvec_multi, XtX, Xr) {
//     library(emulator)
// # Yvec_multi: multi-subject data with dimension ns x Nvoxels (N * N) x n
// # ns: number of subjects
// # Nvoxels: number of voxels
// # n: number of time points
// # Return: M values for each subject and voxel
// ##########################
//     ns <- dim(Yvec_multi)[1]
//     Nvoxels <- dim(Yvec_multi)[2]
//     
//     M_ga0.mat <- matrix(0, ns, Nvoxels)
//     M_ga1.mat <- matrix(0, ns, Nvoxels)
//     
//     for (s in 1:ns) {
//         XtYs <- crossprod(Xr, t(Yvec_multi[s, , ]))
//         for(v in 1:Nvoxels) {
//             M_ga0 <- crossprod(Yvec_multi[s, v, ])
//             M_ga1 <- M_ga0 - quad.form.inv(XtX, XtYs[, v])
//             
//             M_ga0.mat[s, v] <- M_ga0
//             M_ga1.mat[s, v] <- M_ga1
//         }
//     }
//     return(list(M0 = M_ga0.mat, M1 = M_ga1.mat))
// }



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_psi_multi_cpp(double psi, arma::vec w, arma::mat Kw, 
                          double aphi, double bphi, arma::mat distmat, 
                          double nu, arma::mat ga_mat, arma::mat a_mat, 
                          bool isold, arma::vec region_conn) {
    int ns = ga_mat.n_rows;
    arma::vec Kw_vec = vec(Kw);
    int N = Kw.n_rows;
    double phi = exp(psi);
    arma::mat p_bern_mat(ns, N);
    arma::mat logber_mat(ns, N);
    // arma::vec a_vec(N);
    double logsumber;
    double loggamma;
    
    // arma::vec logber(N);
    arma::mat K;
    
    if (isold) { 
        for (int s = 0; s < ns; s++) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern_mat(s, i) = 1 / (1 + exp(-(a_mat(s, region_conn[i] - 1) + Kw_vec[i])));
                logber_mat(s, i) = R::dbinom(ga_mat(s, i), 1, p_bern_mat(s, i), 1);
            }
        }
        logsumber = accu(logber_mat);
        // loggamma = R::dgamma(phi, aphi, 1/bphi, 1);
        loggamma = R::dgamma(1 / phi, aphi, 1 / bphi, 1) - 2 * log(phi);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi
        );
    } else {
        K = bezier2_cpp(distmat = distmat, nu = nu, phi = phi);
        // K = bezier2_cpp(distmat = distmat, nu = nu, phi = phi);
        Kw = K * w;
        Kw_vec = vec(Kw);
        for (int s = 0; s < ns; s++) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern_mat(s, i) = 1 / (1 + exp(-(a_mat(s, region_conn[i] - 1) + Kw_vec[i])));
                logber_mat(s, i) = R::dbinom(ga_mat(s, i), 1, p_bern_mat(s, i), 1);
            }
        }
        logsumber = accu(logber_mat);
        // loggamma = R::dgamma(phi, aphi, 1/bphi, 1);
        loggamma = R::dgamma(1 / phi, aphi, 1 / bphi, 1) - 2 * log(phi);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi,
            Rcpp::Named("K")     = K,
            Rcpp::Named("Kw")    = Kw
        );
    }
}


// 
// logcond_psi_multi <- function(psi, w, Kw, a_phi, b_phi, dist.mat, nu, 
//                               ga_mat, a_mat, is.old = TRUE, region_conn) {
//     ns <- dim(ga_mat)[1]
//     N.voxels <- dim(ga_mat)[2]
//     phi <- exp(psi)
//     if (is.old) {
//         a_mat_v <- sapply(1:length(Kw), function(x) a_mat[, region_conn[x]])
//         Kw_mat <- matrix(rep(Kw, each = ns), ns, N.voxels)
//         p_bern_mat <- 1 / (1 + exp(-(a_mat_v + Kw_mat)))
//         logsumber <- sum(dbinom(ga_mat, size = 1, prob = p_bern_mat, log = TRUE))
//         loggamma <- dgamma(phi, a_phi, b_phi, log = TRUE)
//         return(list(value = logsumber + loggamma + psi))
//     } else {
//         K <- bezier2(dist.mat = dist.mat, nu = nu, phi = phi) 
//         Kw <- K %*% w
//         a_mat_v <- sapply(1:length(Kw), function(x) a_mat[, region_conn[x]])
//         Kw_mat <- matrix(rep(Kw, each = ns), ns, N.voxels)
//         p_bern_mat <- 1 / (1 + exp(-(a_mat_v + Kw_mat)))
//         logsumber <- sum(dbinom(ga_mat, size = 1, prob = p_bern_mat, log = TRUE))
//         loggamma <- dgamma(phi, a_phi, b_phi, log = TRUE)
//         return(list(value = logsumber + loggamma + psi, K = K, Kw = Kw))
//     }
// }


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



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_a_multi_prec_cpp(int s, arma::mat a_mat, arma::mat a_new,
                              arma::mat ga_mat, arma::mat Kw, arma::mat Prec,
                              arma::vec region_conn, bool isold) {
    // Rcout << "\r" << "check 1 " << endl;
    int ns = ga_mat.n_rows;
    int Nvoxels = ga_mat.n_cols;
    arma::vec Kw_vec = vec(Kw);
    // Rcout << "\r" << "check 1.1 " << endl;
    arma::vec p_bern(Nvoxels);
    arma::vec logber(Nvoxels);
    double logsumber;
    double loggamma;
    // Rcout << "\r" << "check 1.2 " << endl;
    arma::rowvec a_vec = a_mat.row(s);
    // Rcout << "\r" << "check 1.2.1 " << endl;
    arma::mat Sig = inv_sympd(Prec);
    // Rcout << "\r" << "check 2 " << endl;
    arma::mat K;
    
    if (isold) { 
        for (int i=0; i < Nvoxels; i++) {
            p_bern[i] = 1 / (1 + exp(-(a_vec[region_conn[i] - 1] + Kw_vec[i])));
            logber[i] = R::dbinom(ga_mat(s, i), 1, p_bern[i], 1);
            // p_bern_mat(s, i) = 1 / (1 + exp(-(a_mat(s, region_conn[i]) + Kw_vec[i])));
            // logber_mat(s, i) = R::dbinom(ga_mat(s, i), 1, p_bern_mat(s, i), 1);
        }
        // Rcout << "\r" << "check 3 " << endl;
        logsumber = accu(logber);
        // Rcout << "size(a_vec)" << size(a_vec) << endl;
        
        arma::rowvec a_mean = zeros<rowvec>(size(a_vec));
        
        arma::vec logNormal = dmvnrm_arma(mat(a_vec), a_mean, Sig, true);
        // Rcout << "\r" << "check 4 " << endl;
        return List::create(
            Rcpp::Named("value") = logsumber + logNormal
        );
    } else {
        a_vec = a_new;
        // Rcout << "\r" << "check 5 " << endl;
        for (int i=0; i < Nvoxels; i++) {
            p_bern[i] = 1 / (1 + exp(-(a_vec[region_conn[i] - 1] + Kw_vec[i])));
            logber[i] = R::dbinom(ga_mat(s, i), 1, p_bern[i], 1);
            // p_bern_mat(s, i) = 1 / (1 + exp(-(a_mat(s, region_conn[i]) + Kw_vec[i])));
            // logber_mat(s, i) = R::dbinom(ga_mat(s, i), 1, p_bern_mat(s, i), 1);
        }
        logsumber = accu(logber);
        
        arma::rowvec a_mean = zeros<rowvec>(size(a_vec));
        // Rcout << "\r" << "check 6 " << endl;
        arma::vec logNormal = dmvnrm_arma(mat(a_vec), a_mean, Sig, true);
        // Rcout << "\r" << "check 7 " << endl;
        a_mat.row(s) = a_new;
        
        return List::create(
            Rcpp::Named("value") = logsumber + logNormal,
            Rcpp::Named("a_mat") = a_mat
        );
    }
}

//     logcond_a_multi_prec <- function(s, a_mat, a_new, ga_mat, Kw, Prec, region_conn,
//                                      is.old = TRUE) {
//         ns <- dim(ga_mat)[1]
//         N.voxels <- dim(ga_mat)[2]
//         if (is.old) {
//             a_vec <- a_mat[s, ]
//             a_vec_v <- sapply(1:length(Kw), function(x) a_vec[region_conn[x]])
//             p_bern <- 1 / (1 + exp(-(a_vec_v + Kw)))
//             sumlogber <- sum(dbinom(ga_mat[s, ], size = 1, prob = p_bern, log = TRUE))
//             Sig <- chol2inv(chol(Prec))
//             logNormal <- mvnfast::dmvn(a_vec, mu = rep(0, length(a_vec)),
//                                        sigma = Sig, log = TRUE)
//             return(list(value = sumlogber + logNormal))
//         } else {
//             a_vec <- a_new
//             a_vec_v <- sapply(1:length(Kw), function(x) a_vec[region_conn[x]])
//             p_bern <- 1 / (1 + exp(-(a_vec_v + Kw)))
//             sumlogber <- sum(dbinom(ga_mat[s, ], size = 1, prob = p_bern, log = TRUE))
// # logNormal <- mvnfast::dmvn(a.vec, mu, sigma = Sig, log = TRUE)
//             Sig <- chol2inv(chol(Prec))
//             logNormal <- mvnfast::dmvn(a_vec, mu = rep(0, length(a_vec)),
//                                        sigma = Sig, log = TRUE)
// # logNormal <- dnorm(ag_new)
//             a_mat[s, ] <- a_new
//             return(list(value = sumlogber + logNormal, a_mat = a_mat))
//         }
//     }



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_multi_prec_cpp(arma::cube Yvec_multi, arma::mat Xr,
                       arma::vec w, arma::mat ga_mat, double phi, double psi,
                       arma::mat a_mat, arma::mat Prec,
                       int nmcmc, int burn, int thin,
                       bool adapt, List tunelst, List keeplst,
                       List tunepsi, List keeppsi,
                       List tunea, List keepa,
                       int tunelen,
                       arma::vec regionlabel, arma::mat region_mat,
                       arma::vec region_conn, int nu,
                       double atau, double btau,
                       double aphi, double bphi,
                       arma::vec agvec, arma::mat distmat, bool iscplx,
                       Function get_tune, Function rWishart, arma::mat draws) {

    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec_multi.n_rows;
    const int n = Yvec_multi.n_cols;
    const int ns = ga_mat.n_rows;
    const int p = 1;
    int G = regionlabel.max();
    arma::uvec uniquelabel = find_unique(mat(region_conn));
    int D = uniquelabel.size() - 1;  // exclude 0
    arma::mat XtX = Xr.t() * Xr;
    // arma::mat XtY = Xr.t() * Yvec.t();

    arma::mat Ker = bezier2_cpp(distmat, nu = nu, phi = phi);
    List MM = marginal_yga_multi_cpp(Yvec_multi, XtX, Xr);
    arma::mat M0_mat = MM["M0_mat"];
    arma::mat M1_mat = MM["M1_mat"];
    // arma::mat Ker = bezier2_cpp(distmat, nu = nu, phi = phi);
    arma::mat Kw = Ker * w;
    arma::mat w_new;
    arma::uvec widx;
    arma::mat Wmat(nmcmc, G);
    // arma::mat Amat(nmcmc, D);
    arma::cube Acube(nmcmc, D, ns);
    arma::cube Gacube(ns, Nvoxels, nmcmc);
    arma::cube Prec_cube(D, D, nmcmc);
    // Prec_array <- array(0, dim = c(n.mcmc, D, D))
    // Ga_mat <- array(NA, dim = c(n.mcmc, ns, N.voxels))
    Rcout << "\r" << "check 1 " << endl;



    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    Rcpp::List keep = keeplst;
    Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning

    // # for psi
    // # keep.psi <- list(psi = 0)
    Rcpp::List keeptmppsi = keeppsi;

    Rcpp::List keeptmpa = keepa;
    int Tb = tunelen;  //# frequency of adaptive tuning



    // W <- matrix(NA, nrow = n.mcmc, ncol = G)

    // Ga_mat <- array(NA, dim = c(n.mcmc, ns, N.voxels))
    // Sig_array <- array(0, dim = c(n.mcmc, D, D))

    Rcout << "\r" << "check 2 " << endl;
    // #-------------------------------
    // # M-H algorithm
    // #-------------------------------
    // Begin algorithm
    Rcout << "KC-S Multi-subject-connectvity Begins" << endl;

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
        
        for(int j = 0; j < 5; ++j) {
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmp[j] = as<double>(keeptmp[j]) / Tb;
                tunelst[j] = callgettune(tunelst[j], keeptmp[j], i, get_tune);
                keeptmp[j] = 0;
            }
            
            widx = 5 * j + linspace<uvec>(0, 4, 5);
            // Rcout << " widx " << widx.n_elem << endl;
            if (i > 1000) {
                // cov(W[(i-1000):(i-1), w.idx])
                arma::mat Sigma_w = cov(Wmat(span(i - 1000, i - 1), span(widx[0], widx[4])));
                Sigma_w = (Sigma_w + Sigma_w.t()) / 2;
                arma::mat sigma = as<double>(tunelst[j]) * Sigma_w;
                w_new = mvrnormArma(1, w.elem(widx), sigma);
                // Rcout << "w_new " << w_new.n_elem << endl;
            } else {
                w_new = mvrnormArma(1, w.elem(widx), eye<mat>(5,5));
                // Rcout << "w_new " << w_new.n_elem << endl;
            }
            
            // arma::vec w, arma::vec w_new, arma::vec ga, 
            // double atau, double btau, arma::mat K, arma::mat Kw, 
            // arma::uvec widx, arma::vec agvec, 
            // bool isold, arma::vec regionlabel
            arma::vec w_new_vec = vectorise(w_new);
            // Rcout << "w_new_vec" << w_new_vec.n_elem << endl;
            
            List logcon_new_cpp = logcond_wg_block_multi_cpp(w, w_new_vec, ga_mat, 
                                                             atau, btau,
                                                             Ker, Kw, widx,
                                                             a_mat, FALSE,
                                                             region_conn);
            List logcon_old_cpp = logcond_wg_block_multi_cpp(w, w_new_vec, ga_mat, 
                                                             atau, btau,
                                                             Ker, Kw, widx,
                                                             a_mat, TRUE,
                                                             region_conn);
            
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
                as<double>(logcon_old_cpp["value"]))) {
                w.elem(widx) = w_new;
                Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                keep[j] = as<int>(keep[j]) + 1;
                keeptmp[j] = as<int>(keeptmp[j]) + 1;
            }
            Wmat(i, span(widx[0], widx[4])) = mat(w.elem(widx)).t();
        }
        
        
        
        
        
        

        // # ============ sample ga ==============
        for (int s = 0; s < ns; ++s) {
            for (int v = 0; v < Nvoxels; ++v) {
                // ga_update_multi_cpp(int s, int v, arma::mat ga_mat,
                //                     arma::mat M0_mat, arma::mat M1_mat, arma::mat Kw, int n,
                //                     arma::vec region_conn, arma::mat a_mat, bool iscplx)

                ga_mat(s, v) = ga_update_multi_cpp(s, v, ga_mat, M0_mat, M1_mat,
                                                   Kw, n, region_conn, a_mat,
                                                   iscplx);
            }
        }
        // Gacube(ns, Nvoxels, nmcmc);
        Gacube.slice(i) = ga_mat;
        // Ga_mat(i, , )  ga_mat

        // # ============ sample psi and hence phi ==============
        double new_psi = rnorm(1, psi, tunepsi["psi"])[0];

        // List logcond_psi_multi_cpp(double psi, arma::vec w, arma::mat Kw,
        //                            double aphi, double bphi, arma::mat distmat,
        //                            double nu, arma::mat ga_mat, arma::mat a_mat,
        //                            bool isold, arma::vec region_conn)
        List logcon_new_psi_cpp = logcond_psi_multi_cpp(new_psi, w, Kw, aphi, bphi,
                                                        distmat, nu, ga_mat,
                                                        a_mat, FALSE, region_conn);
        List logcon_old_psi_cpp = logcond_psi_multi_cpp(psi, w, Kw, aphi, bphi,
                                                        distmat, nu, ga_mat,
                                                        a_mat, TRUE, region_conn);


        if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) -
            as<double>(logcon_old_psi_cpp["value"]))) {
            psi = new_psi;
            phi = exp(psi);
            Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
            Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
            keeppsi["psi"] = as<int>(keeppsi["psi"]) + 1;
            keeptmppsi["psi"] = as<int>(keeptmppsi["psi"]) + 1;
        }

        // # ============ sample a ==============
        for (int s = 0; s < ns; ++s) {
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmpa[s] = as<double>(keeptmpa[s]) / Tb;
                tunea[s] = callgettune(tunea[s], keeptmpa[s], i, get_tune);
                keeptmpa[s] = 0;
            }
            arma::vec mu = vec(a_mat.row(s));
            arma::mat a_new;
            if (i > 1000) {
                arma::mat Sigma_a = cov(Acube.slice(s).rows(i - 1000, i - 1));
                Sigma_a = (Sigma_a + Sigma_a.t()) / 2;
                arma::mat sigma = as<double>(tunea[s]) * Sigma_a;
                // Sigma.a <- as.matrix(Matrix::nearPD(Sigma.a)$mat)
                // mvrnormArma(int n, arma::vec mu, arma::mat sigma)
                a_new = mvrnormArma(1, mu, sigma);
                // arma::mat a_new = mvrnormArma(1, vec(a_mat.row(s)), tunea[s] * Sigma_a);
            } else {
                a_new = mvrnormArma(1, vec(a_mat.row(s)), as<double>(tunea[s]) * eye<mat>(D,D));
            }
            // logcond_a_old <- logcond_a_multi(s, a_mat, a_new, ga_mat, Kw, Sig,
            //                                  region_conn, is.old = TRUE)
            //     logcond_a_new <- logcond_a_multi(s, a_mat, a_new, ga_mat, Kw, Sig,
            //                                      region_conn, is.old = FALSE)
            List logcond_a_old = logcond_a_multi_prec_cpp(s, a_mat, a_new,
                                                          ga_mat, Kw, Prec,
                                                          region_conn, TRUE);
            List logcond_a_new = logcond_a_multi_prec_cpp(s, a_mat, a_new,
                                                          ga_mat, Kw, Prec,
                                                          region_conn, FALSE);
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcond_a_new["value"]) -
                as<double>(logcond_a_old["value"]))) {
                a_mat.row(s) = a_new;
                keepa[s] = as<int>(keepa[s]) + 1;
                keeptmpa[s] = as<int>(keeptmpa[s]) + 1;
            }
            
            Acube.slice(s).row(i) = a_mat.row(s);
            // arma::cube B(nmcmc, D, ns);
            // B(1, 2, 3) = 5;
            // B(i, span(0, D - 1), s) = vectorise(a_mat.row(s));
        }



        // # ============ sample Prec ==============
        int r0 = D;
        arma::mat Covar = a_mat.t() * a_mat + r0 * eye<mat>(D, D);
        arma::mat Sigma = inv_sympd(Covar);
        Prec = callrWishart(ns + r0, Sigma, rWishart);
        
        Prec_cube.slice(i) = Prec;
            // Prec <- Prec[, , 1]
        // Prec_array[i, , ] <- Prec



        // #  Save samples
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, ns * Nvoxels - 1)) = vectorise(ga_mat, 0).t();
                draws(d - 1, span(ns * Nvoxels, ns * Nvoxels + G - 1)) = w.t();
                draws(d - 1, span(ns * Nvoxels + G, ns * Nvoxels + G + ns * D - 1)) = vectorise(a_mat, 0).t();
                draws(d - 1, draws.n_cols - 1) = phi;
                //                 Gamma.row(d) = ga.t();
                //                 W.row(d) = w.t();
                //                 Phi.row(d) <- phi;
                // ns * Nvoxels
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
        _["A"] = Acube,
        _["Prec_array"] = Prec_cube,
        _["Ga_mat"] = Gacube,
        _["W"] = Wmat,
        _["accept_keep"] = keep,
        _["accept_keep_psi"] = keeppsi,
        _["burn"] = burn,
        _["thin"] = thin,
        _["nmcmc"] = nmcmc
    );
}

// 
// mcmc_pc1_multi_prec <- function(Yvec_multi, Xr, start.lst, n.mcmc, burn, thin, 
//                                 name.par, adapt = TRUE, tune.w, keep.w, 
//                                 tune.psi, keep.psi, tune.a, keep.a,
//                                 tune.len = 30, regionlabel, region_mat, 
//                                 region_conn, nu = 2, 
//                                 a.tau = 1, b.tau = 1,
//                                 a.phi = 1, b.phi = 1, 
//                                 ag.vec = rep(0, max(regionlabel)),
//                                 target.accept.rate, is.cplx) {
// #################################################################
// # Yvec_multi: ns x N*N x n COMPLEX-valued data
// # Xr: n by p design matrix (p = 1)
// # init_ga: initial indicator variable
// # init_S: initial spatial effect (ROW factor)
// # init_r: initial range parameter
// # m: number of iterations
// # regionlabel: N*N by 1 region vector
// # region_mat: N by N region matrix
// # region_conn: N*N by 1 region vector for connectivity
// #################################################################
// #-------------------------------
// # Functions used in the algorithm
// #-------------------------------
//     
//     if(missing(region_conn)) {
//         region_conn <- regionlabel
//     }
//     
//     get.tune <- function(tune, keep, k, target = target.accept.rate){  # adaptive tuning
//         aa <- min(0.5, 1000 / sqrt(k))
// # a <- min(0.025, 1 / sqrt(k))
//         exp(ifelse(keep < target, log(tune) - aa, log(tune) + aa))
//     }
//     
// #-------------------------------
// # Adaptive tuning
// #-------------------------------
// # for w
//     keep <- keep.w
//         keep.tmp.w <- keep  # track MH accpetance rate for adaptive tuning
//             
// # for psi
//             keep.tmp.psi <- keep.psi
//                 
// # for a
//                 keep.tmp.a <- keep.a
//                     
//                     Tb <- tune.len  # frequency of adaptive tuning
//                     
// #-------------------------------
// # Here create some objects that will be used later.
// #-------------------------------
//                     ns <- dim(Yvec_multi)[1]
//                     N.voxels <- dim(Yvec_multi)[2]
//                     n <- dim(Yvec_multi)[3]
//                     N <- sqrt(N.voxels)
//                         p <- 1
//                     G <- max(regionlabel)
//                         D <- max(region_conn)
//                         XtX <- crossprod(Xr)
// 
//                         centroid <- compute_centroid(region_mat)
//                         Coord <- cbind(rep(1:N, N), rep(1:N, each = N))
//                         dist.mat <- dist_mat(coord = Coord, grid = centroid)
//                         
//                         MM_mat <- marginal_yga_multi(Yvec_multi, XtX, Xr)
//                         
// #-------------------------------
// # Storage
// #-------------------------------
//                         sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
//                             draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
//                             colnames(draws) <- name.par
//                             W <- matrix(NA, nrow = n.mcmc, ncol = G)
//                             A <- array(NA, dim = c(n.mcmc, ns, D))
//                             Ga_mat <- array(NA, dim = c(n.mcmc, ns, N.voxels))
//                             Prec_array <- array(0, dim = c(n.mcmc, D, D))
//                             
// #-------------------------------
// # Starting values
// #-------------------------------
//                             w <- as.vector(start.lst$w)
//                                 ga_mat <- start.lst$ga_mat
//                                 phi <- start.lst$phi
//                                 psi <- log(phi)
//                                 a_mat <- start.lst$a_mat
//                                 Prec <- start.lst$Prec
//                                 
// # create kernel
//                                 Ker <- bezier2(dist.mat, nu = nu, phi = phi) 
//                                     Kw <- Ker %*% w
//                                     
// #-------------------------------
// # M-H algorithm
// #------------------------------- 
//                                     for (i in 1:n.mcmc) { 
//                                         if (i %% 1 == 0) cat("iter:", i, "\r")
//                                             
// # Update tuning parameter of psi
// #------------------------
//                                             if(adapt == TRUE & i %% Tb == 0) {  
// # Adaptive tuning
//                                                 keep.tmp.psi$psi <- keep.tmp.psi$psi / Tb
//                                                 tune.psi$psi <- get.tune(tune.psi$psi, keep.tmp.psi$psi, i)
//                                                 keep.tmp.psi$psi <- 0
//                                             } 	
//                                             
// #         # Update tuning parameter of a
// #         #------------------------
// #         if(adapt == TRUE & i %% Tb == 0) {  
// #             # Adaptive tuning
// #             keep.tmp.a$a <- keep.tmp.a$a / Tb
// #             tune.a$a <- get.tune(tune.a$a, keep.tmp.a$a, i)
// #             keep.tmp.a$a <- 0
// #         } 	
//                                             
// # ========== sample w =============
//                                             for (j in 1:5) {
//                                                 if(adapt == TRUE & i %% Tb == 0) {
// # Adaptive tuning
//                                                     keep.tmp.w[[j]] <- keep.tmp.w[[j]] / Tb
//                                                     tune.w[[j]] <- get.tune(tune.w[[j]], keep.tmp.w[[j]], i)
//                                                     keep.tmp.w[[j]] <- 0
//                                                 }
//                                                 w.idx <- (5 * (j - 1) + 1:5)
// # block move using empirical covariance structure
//                                                     if (i > 500) {
//                                                         Sigma.w <- cov(W[(i - (Tb+200)):(i - 1), w.idx])
//                                                         Sigma.w <- as.matrix(Matrix::nearPD(Sigma.w, doSym = TRUE)$mat)
// # print(round(Sigma.w, 2))
//                                                         w_new <- as.vector(rmvn(1, w[w.idx], tune.w[[j]] * Sigma.w))
//                                                     } else {
//                                                         w_new <- as.vector(rmvn(1, w[w.idx], tune.w[[j]] * diag(5)))
//                                                     }
// # logcond_wg_block_multisub <- function(w, w_new, ga_mat, a.tau, b.tau, K, Kw, 
// # w.idx, ag_vec, is.old = TRUE, regionlabel)
//                                                     logcon.new <- logcond_wg_block_multi(w, w_new, ga_mat, a.tau, b.tau, 
//                                                                                          K = Ker, Kw, w.idx, 
//                                                                                          a_mat = a_mat, is.old = FALSE, 
//                                                                                          region_conn = region_conn)
//                                                         logcon.old <- logcond_wg_block_multi(w, w_new, ga_mat, a.tau, b.tau, 
//                                                                                              K = Ker, Kw, w.idx, 
//                                                                                              a_mat = a_mat, is.old = TRUE, 
//                                                                                              region_conn = region_conn)
//                                                         if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
//                                                             w[w.idx] <- w_new
//                                                             Kw <- logcon.new$Kw
//                                                             keep[[j]] <- keep[[j]] + 1
//                                                             keep.tmp.w[[j]] <- keep.tmp.w[[j]] + 1
// # W[i, w.idx] <- w_new
//                                                         } 
//                                                         W[i, w.idx] <- w[w.idx]
//                                             }
// # ============ sample ga ==============
//                                             for (s in 1:ns) {
//                                                 for (v in 1:N.voxels) {
// # MM <- marginal_yga(Yvec[v, ], XtY[, v], XtX)
// # ga_update_multisub <- function(s, v, ga_mat, M0_mat, M1_mat, Kw, n, regionlabel, 
// # ag_vec, is.cplx)
//                                                     
// # ga_update_multi <- function(s, v, ga_mat, M0_mat, M1_mat, Kw, n,
// #                             region_conn, a_mat, is.cplx) 
//                                                     ga_mat[s, v] <- ga_update_multi(s, v, ga_mat, MM_mat$M0, 
//                                                                                     MM_mat$M1, Kw, n, region_conn, 
//                                                                                     a_mat = a_mat, is.cplx)
//                                                 }
//                                             }
//                                             
//                                             Ga_mat[i, , ] <- ga_mat
//                                                 
// # ============ sample psi and hence phi ==============
//                                                 new_psi <- rnorm(1, psi, tune.psi$psi)
// #         logcond_psi_multisub <- function(psi, w, Kw, a_phi, b_phi, dist.mat, nu, 
// #                                          ga_mat, ag_vec, is.old = TRUE, regionlabel) 
// # logcond_psi_multi <- function(psi, w, Kw, a_phi, b_phi, dist.mat, nu, 
// #                               ga_mat, a_mat, is.old = TRUE, region_conn) 
//                                                     logcon.new <- logcond_psi_multi(new_psi, w, Kw, a.phi, b.phi, 
//                                                                                     dist.mat = dist.mat, nu = nu, 
//                                                                                     ga_mat = ga_mat, 
//                                                                                     a_mat = a_mat, is.old = FALSE, 
//                                                                                     region_conn = region_conn)
//                                                     logcon.old <- logcond_psi_multi(psi, w, Kw, a.phi, b.phi,
//                                                                                     dist.mat = dist.mat, nu = nu, 
//                                                                                     ga_mat = ga_mat,
//                                                                                     a_mat = a_mat, is.old = TRUE, 
//                                                                                     region_conn = region_conn)
//                                                     
//                                                     if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
//                                                         psi <- new_psi
//                                                         phi <- exp(psi)
//                                                         Ker <- logcon.new$K
//                                                         Kw <- logcon.new$Kw
//                                                         keep.psi$psi <- keep.psi$psi + 1
//                                                         keep.tmp.psi$psi <- keep.tmp.psi$psi + 1
//                                                     } 
//                                                     
// # ============ sample a ==============
//                                                     for (s in 1:ns) {
// # Update tuning parameter of a
// #------------------------
//                                                         if(adapt == TRUE & i %% Tb == 0) {  
// # Adaptive tuning
//                                                             keep.tmp.a[[s]] <- keep.tmp.a[[s]] / Tb
//                                                             tune.a[[s]] <- get.tune(tune.a[[s]], keep.tmp.a[[s]], i)
//                                                             keep.tmp.a[[s]] <- 0
//                                                         } 
// # A <- array(NA, dim = c(n.mcmc, ns, D))
//                                                         if (i > 500) {
//                                                             Sigma.a <- cov(A[(i - (Tb + 200)):(i - 1), s, ])
//                                                             Sigma.a <- as.matrix(Matrix::nearPD(Sigma.a)$mat)
// # print(round(Sigma.w, 2))
//                                                             a_new <- as.vector(rmvn(1, a_mat[s, ], tune.a[[s]] * Sigma.a))
// # new_avec <- as.vector(rmvn(1, avec, tune.a$a * diag(G)))
//                                                         } else {
//                                                             a_new <- as.vector(rmvn(1, a_mat[s, ], tune.a[[s]] * diag(D)))
//                                                         }
// # logcond_a_multi(s, a_mat, a_new, ga_mat, Kw, priorSig, region_conn,
// #                 is.old = TRUE)
// # logcond_a_multi <- function(ag_vec, ag_new, ga_mat, Kw, priorSig, g,
// #                              is.old = TRUE) 
//                                                         logcond_a_old <- logcond_a_multi_prec(s, a_mat, a_new, ga_mat, Kw, Prec, 
//                                                                                               region_conn, is.old = TRUE)
//                                                             logcond_a_new <- logcond_a_multi_prec(s, a_mat, a_new, ga_mat, Kw, Prec, 
//                                                                                                   region_conn, is.old = FALSE)
//                                                             
//                                                             if (log(runif(1)) < (logcond_a_new$value - logcond_a_old$value)) {
//                                                                 a_mat[s, ] <- a_new
//                                                                 keep.a[[s]] <- keep.a[[s]] + 1
//                                                                 keep.tmp.a[[s]] <- keep.tmp.a[[s]] + 1
//                                                             } 
//                                                             
//                                                             A[i, s, ] <- a_mat[s, ]
// # for (j in 1:5) {
// #     # if(adapt == TRUE & i %% Tb == 0) {
// #     #     # Adaptive tuning
// #     #     keep.tmp.a[[j]] <- keep.tmp.a[[j]] / Tb
// #     #     tune.a[[j]] <- get.tune(tune.a[[j]], keep.tmp.a[[j]], i)
// #     #     keep.tmp.a[[j]] <- 0
// #     # }
// #     a.idx <- (5 * (j - 1) + 1:5)
// #     if (i > 500) {
// #         Sigma.a <- cov(A[(i-Tb):(i-1), a.idx])
// #         Sigma.a <- as.matrix(Matrix::nearPD(Sigma.a)$mat)
// #         # print(round(Sigma.w, 2))
// #         new_avec <- as.vector(rmvn(1, ag_vec[a.idx], tune.a[[j]] * Sigma.a))
// #         # new_avec <- as.vector(rmvn(1, avec, tune.a$a * diag(G)))
// #     } else {
// #         new_avec <- as.vector(rmvn(1, ag_vec[a.idx], tune.a[[j]] * diag(G)))
// #     }
// #     #             logcond_ag_multisub <- function(ag_vec, ag_new, ga_mat, Kw, priorSig, g,
// #     #                                             is.old = TRUE) 
// #     logcond_a_old <- logcond_ag_multisub(ag_vec, new_avec, ga_mat, Kw, 
// #                                          rep(0, G), a.idx, TRUE)
// #     logcond_a_new <- logcond_ag_multisub(ag_vec, new_avec, ga_mat, Kw, 
// #                                          rep(0, G), a.idx, FALSE)
// #     
// #     if (log(runif(1)) < (logcon.new$value - logcon.old$value)) {
// #         ag_vec[a.idx] <- new_avec
// #         keep.a[[j]] <- keep.a[[j]] + 1
// #         keep.tmp.a[[j]] <- keep.tmp.a[[j]] + 1
// #         # A[i, a.idx] <- new_avec
// #     }
// #     A[i, a.idx] <- ag_vec[a.idx]
// # }
//                                                             
//                                                     }
//                                                     
//                                                     
// #         # ============ sample Sigma ==============
// #         # prior IW(G, diag(G))
// #        
//                                                     
// # Sig <- MCMCpack::riwish(v = ns + D, S = t(a_mat) %*% a_mat + diag(D))
//                                                     Covar <- t(a_mat) %*% a_mat + diag(D)
//                                                         Sigma <- chol2inv(chol(Covar))
//                                                         Prec <- stats::rWishart(1, df = ns + D, Sigma = Sigma)
//                                                         Prec <- Prec[, , 1]
//                                                     Prec_array[i, , ] <- Prec
//                                                         
//                                                         
// #  Save samples 
// # -----------------------------------------------------
//                                                         if (i > burn) {
//                                                             if (i %% thin == 0) {
// # print(w)
//                                                                 draws[(i - burn) %/% thin, ] <- c(as.vector(ga_mat), w, 
//                                                                        as.vector(a_mat), phi)
// # (v = 1, s = 1), ..., (v = 1, s = ns), ..., 
// # (v = N.voxels, s = 1), ..., (v = N.voxels, s = ns)
//                                                             } 
//                                                         }
//                                                         
//                                                         
//                                     }
// # Acceptance Probability
// #----------------------------
//                                     keep <- lapply(keep, function(x) x / n.mcmc)
//                                         keep.psi$psi <- keep.psi$psi / n.mcmc
//                                         keep.a <- lapply(keep.a, function(x) x / n.mcmc)
// # keep.a$a <- keep.a$a / n.mcmc
// # keep.w$w <- keep.w$w / n.mcmc
// # keep <- lapply(keep, function(x) x / n.mcmc)
// # Write output
// #--------------
//                                         return(list(draws = draws, 
//                                                     A = A,
//                                                     Prec_array = Prec_array,
//                                                     Ga_mat = Ga_mat,
//                                                     W = W,
//                                                     accept = c(keep, keep.psi, keep.a), 
//                                                     start = start.lst, 
//                                                     tune = c(tune.w, tune.psi, tune.a), 
//                                                     burn = burn, thin = thin, 
//                                                     n.mcmc = n.mcmc, sampleidx = sampleidx))
// }
// 
// 
// 
// 








