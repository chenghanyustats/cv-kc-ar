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
List logcond_wc_cpp(arma::vec w, double wg_new, arma::vec ga, 
                        double atau, double btau, arma::mat K, arma::mat Kw, 
                        arma::mat Kwf, int gidx, arma::vec agvec, 
                        bool isold, arma::vec regionlabel) {
    arma::vec Kw_vec = vec(Kw);
    int G = w.n_elem;
    int N = Kw.n_rows;
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double sumlogber;
    arma::vec logber(N);
    
    arma::mat Kw_combine = Kw + Kwf;
    arma::vec Kw_combine_vec = vec(Kw_combine);
    // Rcout << "Kw_combine_vec " << Kw_combine_vec.n_elem << endl;
    if(isold) {
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                // Rcout << "i " << i << endl;
                p_bern[i] = 1 / (1 + exp(-(Kw_combine_vec[i])));
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
        Kw_combine = Kw + Kwf;
        arma::vec Kw_combine_vec = vec(Kw_combine);
        // Rcout << "Kw_combine_vec " << Kw_combine_vec.n_elem << endl;
        //         w[gidx] = wg_new;
        //         Kw = K * w;
        // Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                // Rcout << "i " << i << endl;
                p_bern[i] = 1 / (1 + exp(-(Kw_combine_vec[i])));
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
    
//     Kw_combine <- Kw + Kwf
//     if(is.old) {
//         if(sum(ag.vec) != 0) {
//             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
//             
//             p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
//         } else {
//             p_bern <- 1 / (1 + exp(-Kw_combine))
//         }
// # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
// # G <- dim(w)[1]
//             B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
//             return(list(value = sumlogber - B))
//     } else {
//         Kw <- Kw + K[, g] * (wg_new - w[g])
//         Kw_combine <- Kw + Kwf
//         if(sum(ag.vec) != 0) {
//             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
//             p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
//         } else {
//             p_bern <- 1 / (1 + exp(-Kw_combine))
//         }
// # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
// # G <- dim(w)[1]
//             w[g] <- wg_new
//             B <- (length(w) / 2 + a.tau) * log((0.5) * sum(w ^ 2) + b.tau)
//             return(list(value = sumlogber - B, Kw = Kw))
//     }
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_wf_cpp(arma::vec wf, double wfg_new, arma::vec ga, 
                    double atauf, double btauf, arma::mat Kf, arma::mat Kw, 
                    arma::mat Kwf, int gidx, arma::vec agvec, 
                    bool isold, arma::vec regionlabel, arma::uvec regionidx) {
    // arma::vec Kwf_vec = vec(Kwf);
    // NumericVector Kwf_vec = vec(Kwf);
    int Gf = wf.n_elem;
    int N = Kw.n_rows;
    // Kwf.rows(regionidx - 1) = Kf * wf;
    // Rcout << "Kw.n_rows" << Kw.n_rows << endl;
    // Rcout << "Kwf.n_rows" << Kwf.n_rows << endl;
    // Rcout << "Kwf.n_rows" << Kwf.rows(regionidx - 1) << endl;
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double sumlogber;
    arma::vec logber(N);
    arma::mat Kw_combine = Kw + Kwf;
    arma::vec Kw_combine_vec = vec(Kw_combine);
    // Rcout << "Kw_combine_vec " << Kw_combine_vec.n_elem << endl;
    if(isold) {
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        sumlogber = sum(logber);
        double B = (Gf / 2 + atauf) * log((0.5) * sum(square(wf)) + btauf);
        return List::create(
            Rcpp::Named("value") = sumlogber - B
        );
    } else {
        // Kwf_vec.elem(regionidx).fill(Kwf_vec.elem(regionidx) + Kf.col(gidx) * (wfg_new - wf(gidx)));
        // Kwf_vec(regionidx) = Kwf_vec(regionidx) + Kf.col(gidx) * (wfg_new - wf(gidx));
        // Kwf.rows(regionidx) = Kwf.rows(regionidx) + Kf.col(gidx) * (wfg_new - wf(gidx));
        Kwf.rows(regionidx - 1) = Kwf.rows(regionidx - 1) + Kf.col(gidx) * (wfg_new - wf(gidx));
        // Rcout << "Kwf.n_rows" << Kwf.n_rows << endl;
        // Kwf_vec.elem(regionidx) = Kwf_vec.elem(regionidx) + Kf.col(gidx) * (wfg_new - wf(gidx));
        // Kwf[region.idx, ] <- Kwf[region.idx, ] + Kf[, g] * (wfg_new - wf[g])
        // Kwf = mat(Kwf_vec);
        Kw_combine = Kw + Kwf;
        arma::vec Kw_combine_vec = vec(Kw_combine);
        // Rcout << "Kw_combine_vec " << Kw_combine_vec.n_elem << endl;
        //         w[gidx] = wg_new;
        //         Kw = K * w;
        // Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        sumlogber = sum(logber);
        wf[gidx] = wfg_new;
        double B = (Gf / 2 + atauf) * log((0.5) * sum(square(wf)) + btauf);
        return List::create(
            Rcpp::Named("value") = sumlogber - B,
            Rcpp::Named("Kwf") = Kwf
        );
    }
    
    
//     Kw_combine <- Kw + Kwf
//     if(is.old) {
//         if(sum(ag.vec) != 0) {
//             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
//             
//             p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
//         } else {
//             p_bern <- 1 / (1 + exp(-Kw_combine))
//         }
// # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
// # G <- dim(w)[1]
//             B <- (length(wf) / 2 + a.tauf) * log((0.5) * sum(wf ^ 2) + b.tauf)
//             return(list(value = sumlogber - B))
//     } else {
// #         b = Kf[, g] * (wfg_new - wf[g])
// #         bb = Kwf[region.idx, ]
// #         print(paste("b len = ", length(b)))
// #         print(paste("bb len = ", length(bb)))
//         Kwf[region.idx, ] <- Kwf[region.idx, ] + Kf[, g] * (wfg_new - wf[g])
// # print(paste("Kwf[region.idx, ] len = ", length(Kwf[region.idx, ])))
// # Kwf[region.idx, ]
// # Kwf[region.idx, ] <- Kwf[region.idx, ] + b
//         Kw_combine <- Kw + Kwf
//         if(sum(ag.vec) != 0) {
//             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
//             p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
//         } else {
//             p_bern <- 1 / (1 + exp(-Kw_combine))
//         }
// # p_bern <- 1 / (1 + exp(-(a.vec + Kw)))
//         sumlogber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
// # G <- dim(w)[1]
//             wf[g] <- wfg_new
//             B <- (length(wf) / 2 + a.tauf) * log((0.5) * sum(wf ^ 2) + b.tauf)
//             return(list(value = sumlogber - B, Kwf = Kwf))
//     }
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List marginal_yga_cpp(arma::mat Yvec, arma::mat XtY, arma::mat XtX) {
    int Nvoxels = Yvec.n_rows;
    arma::vec M_ga0(Nvoxels);
    arma::vec M_ga1(Nvoxels);
    arma::vec Yvecv;
    arma::vec XtYv;
    arma::mat XtXinv = inv_sympd(XtX);
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
List logcond_psic_cpp(double psi, arma::vec w, arma::mat Kw, arma::mat Kwf,
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

    
    arma::mat Kw_combine = Kw + Kwf;
    arma::vec Kw_combine_vec = vec(Kw_combine);
    if (isold) { 
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_combine_vec[i])));
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
        arma::mat K = bezier2_cpp(distmat, nu, phi);
        Kw = K * w;
        Kw_combine = Kw + Kwf;
        arma::vec Kw_combine_vec = vec(Kw_combine);
        // Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_combine_vec[i])));
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
    
    // logcond_psic <- function(psi, w, Kw, Kwf, a.phi, b.phi, dist.mat, nu, ga, ag.vec, 
    //                          is.old = TRUE, regionlabel) {
    //     Kw_combine <- Kw + Kwf
    //     phi <- exp(psi)
    //     if (is.old) {
    //         if(sum(ag.vec) != 0) {
    //             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
    //             p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
    //         } else {
    //             p_bern <- 1 / (1 + exp(-Kw_combine))
    //         }
    //         logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
    //             loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
    //             return(list(value = logsumber + loggamma + psi))
    //     } else {
    //         K <- bezier2(dist.mat = dist.mat, nu = nu, phi = phi) 
    //         Kw <- K %*% w
    //         Kw_combine <- Kw + Kwf
    //         if(sum(ag.vec) != 0) {
    //             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
    //             p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
    //         } else {
    //             p_bern <- 1 / (1 + exp(-Kw_combine))
    //         }
    //         logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
    //             loggamma <- dgamma(phi, a.phi, b.phi, log = TRUE)
    //             return(list(value = logsumber + loggamma + psi, K = K, Kw = Kw))
    //     }
    // }
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_psif_cpp(double psif, arma::vec wf, arma::mat Kw, arma::mat Kwf,
                      double aphif, double bphif, arma::mat distmatf, 
                      double nu, arma::vec ga, arma::vec agvec, 
                      bool isold, arma::vec regionlabel, arma::uvec regionidx) {
    // logcond_psif <- function(psif, wf, Kw, Kwf, a.phif, b.phif, dist.matf, nu, ga, ag.vec, 
    //                          is.old = TRUE, regionlabel, region.idx) 
    arma::vec Kwf_vec = vec(Kwf);
    int N = Kw.n_rows;
    double phif = exp(psif);
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double logsumber;
    double loggamma;
    arma::vec logber(N);
    // arma::mat K;
    
    arma::mat Kw_combine = Kw + Kwf;
    arma::vec Kw_combine_vec = vec(Kw_combine);
    if (isold) { 
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        loggamma = R::dgamma(phif, aphif, 1/bphif, 1);
        // loggamma = R::dgamma(1 / phi, aphi, 1 / bphi, 1) - 2 * log(phi);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psif
        );
    } else {
        arma::mat Kf = bezier2_cpp(distmatf, nu, phif);
        // Rcout << "Kwf.n_rows " << Kwf.n_rows << endl;
        // Rcout << "Kw.n_rows" << Kw.n_rows << endl;
        // Kwf = Kf * wf;
        // Kwf.rows(regionidx - 1) = Kf * wf;
        Kwf.rows(regionidx - 1) = Kf * wf;
        
        Kw_combine = Kw + Kwf;
        arma::vec Kw_combine_vec = vec(Kw_combine);
        // Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; i++) {
                // a_vec[i] = agvec(regionlabel(i));
                // p_bern[i] = 1 / (1 + exp(-(a_vec[i] + Kw_vec[i])));
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; i++) {
                p_bern[i] = 1 / (1 + exp(-(Kw_combine_vec[i])));
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        loggamma = R::dgamma(phif, aphif, 1/bphif, 1);
        // loggamma = R::dgamma(1 / phi, aphi, 1 / bphi, 1) - 2 * log(phi);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psif,
            Rcpp::Named("Kf")     = Kf,
            Rcpp::Named("Kwf")    = Kwf
        );
    }
    
    // logcond_psif <- function(psif, wf, Kw, Kwf, a.phif, b.phif, dist.matf, nu, ga, ag.vec, 
    //                          is.old = TRUE, regionlabel, region.idx) {
    //     Kw_combine <- Kw + Kwf
    //     phif <- exp(psif)
    //     if (is.old) {
    //         if(sum(ag.vec) != 0) {
    //             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
    //             p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
    //         } else {
    //             p_bern <- 1 / (1 + exp(-Kw_combine))
    //         }
    //         logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
    //             loggamma <- dgamma(phif, a.phif, b.phif, log = TRUE)
    //             return(list(value = logsumber + loggamma + psif))
    //     } else {
    //         Kf <- bezier2(dist.mat = dist.matf, nu = nu, phi = phif) 
    //         Kwf <- matrix(0, length(Kw), 1)
    //         Kwf[region.idx, ] <- Kf %*% wf
    //         Kw_combine <- Kw + Kwf
    //         if(sum(ag.vec) != 0) {
    //             a.vec <- sapply(1:length(Kw), function(x) ag.vec[regionlabel[x]])
    //             p_bern <- 1 / (1 + exp(-(a.vec + Kw_combine)))
    //         } else {
    //             p_bern <- 1 / (1 + exp(-Kw_combine))
    //         }
    //         logsumber <- sum(dbinom(ga, size = 1, prob = p_bern, log = TRUE))
    //             loggamma <- dgamma(phif, a.phif, b.phif, log = TRUE)
    //             return(list(value = logsumber + loggamma + psif, Kf = Kf, Kwf = Kwf))
    //     }
    // }
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_multires_cpp_main(arma::mat Yvec, arma::mat Xr, 
                           arma::vec w, arma::vec ga, double phi, double psi,
                           arma::vec wf, double phif, double psif,
                           int nmcmc, int burn, int thin, 
                           bool adapt, List tunelst, List keeplst, 
                           List tunelstf, List keeplstf, 
                           List tunepsi, List keeppsi, 
                           int tunelen, 
                           arma::vec regionlabel, 
                           arma::vec regionlabelfiner, 
                           arma::uvec regionidx, 
                           int nu, 
                           double atau, double btau,
                           double aphi, double bphi,
                           double atauf, double btauf,
                           double aphif, double bphif, 
                           arma::vec agvec, 
                           arma::mat distmat, arma::mat distmatf, 
                           bool iscplx,
                           Function get_tune, arma::mat draws) {
    
    
    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int Nvoxels = Yvec.n_rows;
    const int n = Yvec.n_cols;
    const int p = 1;
    int G = regionlabel.max();
    arma::uvec uniquelabel = find_unique(mat(regionlabelfiner));
    int Gf = uniquelabel.size() - 1;  // exclude 0
    // Rcout << "\r" << "Gf: " <<  Gf << endl;  
    // vector<int>::iterator ip;
    // arma::vec last = std::unique(regionlabelfiner.begin(), regionlabelfiner.end());
    // vector<int> vecf = regionlabelfiner;
    // ip = std::unique(vecf.begin(), vecf.end());
    
    // int Gf = last.size();
    
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_cpp(distmat, nu = nu, phi = phi);
    arma::mat Kw = Ker * w;
    arma::mat Kerf = bezier2_cpp(distmatf, nu = nu, phi = phif);
    // arma::mat Kwf = Kerf * wf;
    arma::mat Kwf(Nvoxels, 1, fill::zeros);
    Kwf.rows(regionidx - 1) = Kerf * wf;  // minus one to let idx start with zero
    
    List MM = marginal_yga_cpp(Yvec, XtY, XtX);
    arma::vec M0 = MM["M0"];
    arma::vec M1 = MM["M1"];
    
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    Rcpp::List keep = keeplst;
    Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    
    // # for w
    Rcpp::List keepf = keeplstf;
    Rcpp::List keeptmpf = keepf;  // track MH accpetance rate for adaptive tuning
    
    // # for psi
    // # keep.psi <- list(psi = 0)
    Rcpp::List keeptmppsi = keeppsi;
    
    
    int Tb = tunelen;  //# frequency of adaptive tuning
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    Rcout << "KC-S Multi-resolution Begins" << endl;
    
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
        
        // # Update tuning parameter of psi
        // #------------------------
        if(adapt == TRUE & fmod(i , Tb) == 0) {  
            // # Adaptive tuning
            // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
            keeptmppsi["psif"] = as<double>(keeptmppsi["psif"]) / Tb;
            tunepsi["psif"] = callgettune(tunepsi["psif"], keeptmppsi["psif"], i, get_tune);
            keeptmppsi["psif"] = 0;
        }
        
        // # Update parameters
        // #----------------------------
        // 
        // # ========== sample wc =============
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


            // logcond_wc_cpp(arma::vec w, double wg_new, arma::vec ga,
            //                double atau, double btau, arma::mat K, arma::mat Kw,
            //                arma::mat Kwf, int gidx, arma::vec agvec,
            //                bool isold, arma::vec regionlabel)


            List logcon_new_cpp = logcond_wc_cpp(w, wg_new, ga, atau, btau,
                                                     Ker, Kw, Kwf, g, agvec,
                                                     FALSE, regionlabel);
            List logcon_old_cpp = logcond_wc_cpp(w, wg_new, ga, atau, btau,
                                                     Ker, Kw, Kwf, g, agvec,
                                                     TRUE, regionlabel);

            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) -
                as<double>(logcon_old_cpp["value"]))) {
                w[g] = wg_new;
                Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                keep[g] = as<int>(keep[g]) + 1;
                keeptmp[g] = as<int>(keeptmp[g]) + 1;
            }
        }

        
        
        // # ========== sample wf =============
        for (int g = 0; g < Gf; ++g) {
            // Rcout << "\r" << "g: " << g + 1 << endl;
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmpf[g] = as<double>(keeptmpf[g]) / Tb;
                tunelstf[g] = callgettune(tunelstf[g], keeptmpf[g], i, get_tune);
                keeptmpf[g] = 0;
            }

            // #----------------------------
            // # Random walk Normal proposal
            // #----------------------------

            // # use random walk proposal newS = S + N(0, c) to update S
            double wfg_new = rnorm(1, wf[g], tunelstf[g])[0];



            List logcon_newf_cpp = logcond_wf_cpp(wf, wfg_new, ga, atauf, btauf,
                                                 Kerf, Kw, Kwf, g, agvec,
                                                 FALSE, regionlabel, regionidx);
            List logcon_oldf_cpp = logcond_wf_cpp(wf, wfg_new, ga, atauf, btauf,
                                                 Kerf, Kw, Kwf, g, agvec,
                                                 TRUE, regionlabel, regionidx);

            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_newf_cpp["value"]) -
                as<double>(logcon_oldf_cpp["value"]))) {
                wf[g] = wfg_new;
                Kwf = as<arma::mat>(logcon_newf_cpp["Kwf"]);
                keepf[g] = as<int>(keepf[g]) + 1;
                keeptmpf[g] = as<int>(keeptmpf[g]) + 1;
            }
        }
        
        
        // # ============ sample ga ==============
        arma::mat Kw_combine = Kw + Kwf;
        
        for (int v = 0; v < Nvoxels; ++v) {

            ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw_combine, n,
                                      regionlabel, agvec, iscplx);
        }
        
        // # ============ sample psi and hence phi ==============
        double new_psi = rnorm(1, psi, tunepsi["psi"])[0]; 
        
        // List logcond_psic_cpp(double psi, arma::vec w, arma::mat Kw, arma::mat Kwf,
        //                       double aphi, double bphi, arma::mat distmat, 
        //                       double nu, arma::vec ga, arma::vec agvec, 
        //                       bool isold, arma::vec regionlabel) 
            
        List logcon_new_psi_cpp = logcond_psic_cpp(new_psi, w, Kw, Kwf,
                                                   aphi, bphi, distmat, nu, ga,
                                                   agvec, FALSE, regionlabel);
        List logcon_old_psi_cpp = logcond_psic_cpp(psi, w, Kw, Kwf,
                                                   aphi, bphi, distmat, nu, ga,
                                                   agvec, TRUE, regionlabel);


        if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) -
            as<double>(logcon_old_psi_cpp["value"]))) {
            psi = new_psi;
            phi = exp(psi);
            Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
            Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
            keeppsi["psi"] = as<int>(keeppsi["psi"]) + 1;
            keeptmppsi["psi"] = as<int>(keeptmppsi["psi"]) + 1;
        }
        
        
        // # ============ sample psif and hence phif ==============
        double new_psif = rnorm(1, psif, tunepsi["psif"])[0]; 
        
        
        // List logcond_psif_cpp(double psif, arma::vec wf, arma::mat Kw, arma::mat Kwf,
        //                       double aphif, double bphif, arma::mat distmatf, 
        //                       double nu, arma::vec ga, arma::vec agvec, 
        //                       bool isold, arma::vec regionlabel, arma::uvec regionidx)
        List logcon_newf_psi_cpp = logcond_psif_cpp(new_psif, wf, Kw, Kwf,
                                                   aphif, bphif, distmatf, nu, ga,
                                                   agvec, FALSE, regionlabel,
                                                   regionidx);
        List logcon_oldf_psi_cpp = logcond_psif_cpp(psif, wf, Kw, Kwf,
                                                   aphif, bphif, distmatf, nu, ga,
                                                   agvec, TRUE, regionlabel,
                                                   regionidx);


        if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_newf_psi_cpp["value"]) -
            as<double>(logcon_oldf_psi_cpp["value"]))) {
            psif = new_psif;
            phif = exp(psif);
            Kerf = as<arma::mat>(logcon_newf_psi_cpp["Kf"]);
            Kwf  = as<arma::mat>(logcon_newf_psi_cpp["Kwf"]);
            keeppsi["psif"] = as<int>(keeppsi["psif"]) + 1;
            keeptmppsi["psif"] = as<int>(keeptmppsi["psif"]) + 1;
        }
        
        
        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                // Rcout << "\r" << "span(0, Nvoxels - 1): " <<  Nvoxels - 1 << endl;
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                // Rcout << "\r" << "span(Nvoxels, Nvoxels + G - 1): " <<  Nvoxels + G - 1 << endl;
                draws(d - 1, span(Nvoxels + G, Nvoxels + G + Gf - 1)) = wf.t();
                // Rcout << "\r" << "span(Nvoxels + G, Nvoxels + Gf - 1): " <<  Nvoxels + Gf - 1 << endl;
                draws(d - 1, Nvoxels + G + Gf) = phi;
                // Rcout << "\r" << "Nvoxels + Gf: " <<  Nvoxels + Gf << endl;
                draws(d - 1, draws.n_cols - 1) = phif;
                // Rcout << "\r" << "draws.n_cols - 1: " <<  draws.n_cols - 1 << endl;
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
        _["accept_keep"] = keep, 
        _["accept_keep_psi"] = keeppsi, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}













