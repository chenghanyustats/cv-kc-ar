#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double callgettune(double tune, double keep, int k, Function f) {
    double res;
    res = as<double>(f(tune, keep, k));
    return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat bezier2_phi_cpp(arma::mat distmat, double nu, arma::vec phivec) {
    int N = distmat.n_rows;
    int M = distmat.n_cols;
    arma::mat K = zeros(N, M);
    for (int i=0; i < N; i++) {
        for (int j=0; j < M; j++) {
            if (distmat(i, j) < phivec[j]) {
                K(i, j) = pow(1 - pow(distmat(i, j) / phivec[j], 2), nu);
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
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
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
List logcond_psig_cpp(arma::vec psi, arma::vec w, arma::mat Kw, 
                      double aphi, double bphi, arma::mat distmat, 
                      double nu, arma::vec ga, arma::vec agvec, 
                      bool isold, arma::vec regionlabel, int gidx) {
    
    arma::vec Kw_vec = vectorise(Kw);
    int N = Kw.n_rows;
    arma::vec phi = exp(psi);
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double logsumber;
    double loggamma;
    arma::vec logber(N);
    arma::mat K;
    
    loggamma = R::dgamma(phi[gidx], aphi, 1/bphi, 1);
    
    if (isold) { 
        if(sum(agvec) != 0) {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi[gidx]
        );
    } else {
        K = bezier2_phi_cpp(distmat, nu, phi);
        Kw = K * w;
        Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma + psi[gidx],
            Rcpp::Named("K")     = K,
            Rcpp::Named("Kw")    = Kw
        );
    }
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List logcond_psi_block_cpp(arma::vec psi, arma::vec psi_new, arma::vec w, 
                           arma::mat Kw, double aphi, double bphi, 
                           arma::mat distmat, double nu, arma::vec ga, 
                           arma::vec agvec, bool isold, arma::vec regionlabel, 
                           arma::uvec phiidx) {
    
    arma::vec Kw_vec = vectorise(Kw);
    int N = Kw.n_rows;
    arma::vec phi = exp(psi);
    arma::vec p_bern(N);
    // arma::vec a_vec(N);
    double logsumber;
    double loggamma_sum = 0;
    double psi_sum = 0;
    arma::vec logber(N);
    arma::mat K;
    

    // loggamma_sum = sum(R::dgamma(phi[phiidx], aphi, 1/bphi, 1));
    // loggamma_sum = sum(R::dgamma(phi, aphi, 1/bphi, 1));
    
    if (isold) { 
        if(sum(agvec) != 0) {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        for(int j = 0; j < phiidx.n_elem; ++j) {
            loggamma_sum = loggamma_sum + R::dgamma(phi[phiidx[j]], aphi, 1/bphi, 1);
            psi_sum = psi_sum + psi[phiidx[j]];
        }
        logsumber = sum(logber);
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma_sum + psi_sum
        );
    } else {
        psi.elem(phiidx) = psi_new;
        K = bezier2_phi_cpp(distmat, nu, phi);
        Kw = K * w;
        Kw_vec = vec(Kw);
        if(sum(agvec) != 0) {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(agvec[regionlabel[i] - 1] + Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        } else {
            for (int i=0; i < N; ++i) {
                p_bern[i] = 1 / (1 + exp(-(Kw_vec[i])));
                if (p_bern[i] == 0) p_bern[i] = 1e-8;
                if (p_bern[i] == 1) p_bern[i] = 1 - 1e-8;
                logber[i] = R::dbinom(ga[i], 1, p_bern[i], 1);
            }
        }
        logsumber = sum(logber);
        for(int j = 0; j < phiidx.n_elem; ++j) {
            loggamma_sum = loggamma_sum + R::dgamma(phi[phiidx[j]], aphi, 1/bphi, 1);
            psi_sum = psi_sum + psi[phiidx[j]];
        }
        return List::create(
            Rcpp::Named("value") = logsumber + loggamma_sum + psi_sum,
            Rcpp::Named("K")     = K,
            Rcpp::Named("Kw")    = Kw
        );
    }
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_phi_main_cpp(arma::mat Yvec, arma::mat Xr, 
                       arma::vec w, arma::vec ga, arma::vec phi, arma::vec psi,
                       int nmcmc, int burn, int thin, 
                       bool adapt, List tunelst, List keeplst, 
                       List tunepsilst, List keeppsilst, int tunelen, 
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
    const int p = 1;
    int G = regionlabel.max();
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_phi_cpp(distmat, nu, phi);
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
    Rcpp::List keeppsi = keeppsilst;
    Rcpp::List keeptmppsi = keeppsi;
    int Tb = tunelen;  //# frequency of adaptive tuning
    
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    Rcout << "KC-S-phi Begins" << endl;
    
    for (int i = 0; i < nmcmc; ++i) {
        
        // print iteration
        
        
        // # Update tuning parameter of psi
//         // #------------------------
//         if(adapt == TRUE & fmod(i , Tb) == 0) {  
//             // # Adaptive tuning
//             // as<int>(keeptmppsi["psi"]) = as<int>(keeptmppsi["psi"]) / Tb;
//             keeptmppsi["psi"] = as<double>(keeptmppsi["psi"]) / Tb;
//             tunepsi["psi"] = callgettune(tunepsi["psi"], keeptmppsi["psi"], i, get_tune);
//             keeptmppsi["psi"] = 0;
//         }
        
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
                keeptmp[g] = as<int>(keeptmp[g]) + 1;
            }
        }
        
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw, n, regionlabel, agvec, iscplx);
        }
        
        // # ============ sample psi and hence phi ==============
        for(int g = 0; g < G; ++g) {
            // # Update tuning parameter of psi
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmppsi[g] = as<double>(keeptmppsi[g]) / Tb;
                tunepsilst[g] = callgettune(tunepsilst[g], keeptmppsi[g], i, get_tune);
                keeptmppsi[g] = 0;
            }
                
            double new_psig = rnorm(1, psi[g], tunepsilst[g])[0];
            arma::vec new_psi = psi;
            new_psi[g] = new_psig;
            
            List logcon_new_psi_cpp = logcond_psig_cpp(new_psi, w, Kw, aphi, bphi, 
                                                      distmat, nu, ga, agvec, FALSE, 
                                                      regionlabel, g);
            List logcon_old_psi_cpp = logcond_psig_cpp(psi, w, Kw, aphi, bphi, 
                                                      distmat, nu, ga, agvec, TRUE, 
                                                      regionlabel, g);
            
            
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) - 
                as<double>(logcon_old_psi_cpp["value"]))) {
                psi[g] = new_psig;
                Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
                Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
                keeppsi[g] = as<int>(keeppsi[g]) + 1;
                keeptmppsi[g] = as<int>(keeptmppsi[g]) + 1;
            }
        }
        phi = exp(psi);

        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                draws(d - 1, span(Nvoxels + G, draws.n_cols - 1)) = phi.t();
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



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_phi_block_main_cpp(arma::mat Yvec, arma::mat Xr, 
                                 arma::vec w, arma::vec ga, arma::vec phi, arma::vec psi,
                                 int nmcmc, int burn, int thin, 
                                 bool adapt, List tunelst, List keeplst, 
                                 List tunepsilst, List keeppsilst, 
                                 List tunepsilst_block, List keeppsilst_block,
                                 int tunelen, 
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
    const int p = 1;
    int G = regionlabel.max();
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_phi_cpp(distmat, nu, phi);
    arma::mat Kw = Ker * w;

    
    List MM = marginal_yga_cpp(Yvec, XtY, XtX);
    arma::vec M0 = MM["M0"];
    arma::vec M1 = MM["M1"];
    arma::mat Wmat(nmcmc, G);
    arma::mat w_new;
    arma::uvec widx;
    arma::mat Psimat(nmcmc, G);
    arma::mat psi_new;
    arma::uvec phiidx;
    
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    Rcpp::List keep = keeplst;
    Rcpp::List keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    
    // # for psi
    // # keep.psi <- list(psi = 0)
    Rcpp::List keeppsi = keeppsilst;
    Rcpp::List keeptmppsi = keeppsi;
    
    Rcpp::List keeppsi_block = keeppsilst_block;
    Rcpp::List keeptmppsi_block = keeppsi_block;
    
    int Tb = tunelen;  //# frequency of adaptive tuning
    
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    Rcout << "CV-KC-phig Begins" << endl;
    
    for (int i = 0; i < nmcmc; ++i) {
        
        // #----------------------------
        // # Update parameters
        // #----------------------------
        
        // # ========== sample w =============
        int sqrtG = sqrt(G);
        for(int j = 0; j < sqrtG; ++j) {
            // Rcout << "check w" << endl;
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmp[j] = as<double>(keeptmp[j]) / Tb;
                tunelst[j] = callgettune(tunelst[j], keeptmp[j], i, get_tune);
                keeptmp[j] = 0;
            }
            // Rcout << "check w_1" << endl;
            widx = sqrtG * j + linspace<uvec>(0, sqrtG - 1, sqrtG);
            // Rcout << " widx " << widx.n_elem << endl;
            if (i > 1000) {
                // cov(W[(i-1000):(i-1), w.idx])
                arma::mat Sigma_w = cov(Wmat(span(i - 1000, i - 1), 
                                             span(widx[0], widx[sqrtG - 1])));
                Sigma_w = (Sigma_w + Sigma_w.t()) / 2;
                arma::mat sigma = as<double>(tunelst[j]) * Sigma_w;
                w_new = mvrnormArma(1, w.elem(widx), sigma);
                // Rcout << "w_new " << w_new.n_elem << endl;
            } else {
                w_new = mvrnormArma(1, w.elem(widx), eye<mat>(sqrtG,sqrtG));
                // w_new = mvrnormArma(1, w.elem(widx), 
                                    // as<double>(tunelst[j]) * eye<mat>(sqrtG,sqrtG));
                // Rcout << "w_new " << w_new.n_elem << endl;
            }
            // Rcout << "check w_2" << endl;
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
            // Rcout << "check w_3" << endl;
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
                as<double>(logcon_old_cpp["value"]))) {
                w.elem(widx) = w_new;
                Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                keep[j] = as<int>(keep[j]) + 1;
                keeptmp[j] = as<int>(keeptmp[j]) + 1;
            }
            // Rcout << "check w_4" << endl;
            Wmat(i, span(widx[0], widx[sqrtG - 1])) = mat(w.elem(widx)).t();
        }
        
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw, n, regionlabel, agvec, iscplx);
        }
        
        // // # ============ sample psi and hence phi ==============
        // // Rcout << "check psi" << endl;
        // for(int j = 0; j < sqrtG; ++j) {
        //     // # Update tuning parameter
        //     // #-------------------------
        //     if(adapt == TRUE & fmod(i , Tb) == 0) {
        //         // # Adaptive tuning
        //         // # Adaptive tuning
        //         keeptmppsi[j] = as<double>(keeptmppsi[j]) / Tb;
        //         tunepsilst[j] = callgettune(tunepsilst[j], keeptmppsi[j], i, get_tune);
        //         keeptmppsi[j] = 0;
        //     }
        //     // Rcout << "check psi_1" << endl;
        //     phiidx = sqrtG * j + linspace<uvec>(0, sqrtG - 1, sqrtG);
        //     // if (i == 100) {
        //     //     Rcout << " phiidx " << phiidx << "\r" << endl;
        //     // }
        //     if (i > 20000) {
        //         // cov(W[(i-1000):(i-1), w.idx])
        //         Rcout << cov(Psimat(span(i - 2000, i - 1),
        //                             span(phiidx[0], phiidx[sqrtG - 1]))) << endl;
        //         arma::mat Sigma_phi = cov(Psimat(span(i - 2000, i - 1),
        //                                          span(phiidx[0], phiidx[sqrtG - 1])));
        //         Sigma_phi = (Sigma_phi + Sigma_phi.t()) / 2;
        //         arma::mat sigma = as<double>(tunepsilst[j]) * Sigma_phi;
        //         psi_new = mvrnormArma(1, psi.elem(phiidx), sigma);
        //         // Rcout << "w_new " << w_new.n_elem << endl;
        //     } else {
        // 
        //         psi_new = mvrnormArma(1, psi.elem(phiidx),
        //                               as<double>(tunepsilst[j]) * eye<mat>(sqrtG,sqrtG));
        //         // Rcout << "check psi_2" << endl;
        //         // Rcout << "w_new " << w_new.n_elem << endl;
        //     }
        // 
        //     arma::vec psi_new_vec = vectorise(psi_new);
        //     // Rcout << "w_new_vec" << w_new_vec.n_elem << endl;
        //     List logcon_new_psi_cpp = logcond_psi_block_cpp(psi, psi_new_vec, w, Kw, aphi, bphi,
        //                                                     distmat, nu, ga,
        //                                                     agvec, FALSE, regionlabel,
        //                                                     phiidx);
        //     List logcon_old_psi_cpp = logcond_psi_block_cpp(psi, psi_new_vec, w, Kw, aphi, bphi,
        //                                                     distmat, nu, ga,
        //                                                     agvec, TRUE, regionlabel,
        //                                                     phiidx);
        //     // Rcout << "check psi_3" << endl;
        //     if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) -
        //         as<double>(logcon_old_psi_cpp["value"]))) {
        //         psi.elem(phiidx) = psi_new;
        //         Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
        //         Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
        //         keeppsi[j] = as<int>(keeppsi[j]) + 1;
        //         keeptmppsi[j] = as<int>(keeptmppsi[j]) + 1;
        //     }
        // 
        //     // Rcout << "phiidx" << phiidx << endl;
        //     // Rcout << "Psi" << psi.elem(phiidx) << endl;
        //     Psimat(i, span(phiidx[0], phiidx[sqrtG - 1])) =
        //         mat(psi.elem(phiidx)).t();
        // }
        // phi = exp(psi);

        
        if (i < nmcmc) {
            for(int g = 0; g < G; ++g) {
                // # Update tuning parameter of psi
                // #-------------------------
                if(adapt == TRUE & fmod(i , Tb) == 0) {
                    // # Adaptive tuning
                    keeptmppsi[g] = as<double>(keeptmppsi[g]) / Tb;
                    tunepsilst[g] = callgettune(tunepsilst[g], keeptmppsi[g], i, get_tune);
                    keeptmppsi[g] = 0;
                }
                
                double new_psig = rnorm(1, psi[g], tunepsilst[g])[0];
                arma::vec new_psi = psi;
                new_psi[g] = new_psig;
                
                List logcon_new_psi_cpp = logcond_psig_cpp(new_psi, w, Kw, aphi, bphi,
                                                           distmat, nu, ga, agvec, FALSE,
                                                           regionlabel, g);
                List logcon_old_psi_cpp = logcond_psig_cpp(psi, w, Kw, aphi, bphi,
                                                           distmat, nu, ga, agvec, TRUE,
                                                           regionlabel, g);
                
                
                if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) -
                    as<double>(logcon_old_psi_cpp["value"]))) {
                    psi[g] = new_psig;
                    Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
                    Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
                    keeppsi[g] = as<int>(keeppsi[g]) + 1;
                    keeptmppsi[g] = as<int>(keeptmppsi[g]) + 1;
                }
                Psimat.row(i) = mat(psi).t();
            }
        } else {
            int B = 1;
            for(int j = 0; j < B; ++j) {
                // # Update tuning parameter
                // #-------------------------
                if(adapt == TRUE & fmod(i , Tb) == 0) {
                    // # Adaptive tuning
                    // # Adaptive tuning
                    keeptmppsi_block[j] = as<double>(keeptmppsi_block[j]) / Tb;
                    tunepsilst_block[j] = callgettune(tunepsilst_block[j], keeptmppsi_block[j], i, get_tune);
                    keeptmppsi_block[j] = 0;
                }
                // Rcout << "check psi_1" << endl;
                // phiidx = B * j + linspace<uvec>(0, B - 1, B);
                phiidx = linspace<uvec>(0, G - 1, G);
//                 Rcout << cov(Psimat(span(i - 2000, i - 1),
//                                     span(phiidx[0], phiidx[sqrtG - 1]))) << endl;
                // arma::mat Sigma_phi = cov(Psimat(span(i - 2000, i - 1),
                //                                  span(phiidx[0], phiidx[sqrtG - 1])));
                arma::mat Sigma_phi = cov(Psimat(span(i - 2000, i - 1),
                                                 span(0, G - 1)));
                Sigma_phi = (Sigma_phi + Sigma_phi.t()) / 2;
                arma::mat sigma = as<double>(tunepsilst_block[j]) * Sigma_phi;
                // psi_new = mvrnormArma(1, psi.elem(phiidx), sigma);
                psi_new = mvrnormArma(1, psi, sigma);
                
                arma::vec psi_new_vec = vectorise(psi_new);
                // Rcout << "w_new_vec" << w_new_vec.n_elem << endl;
                List logcon_new_psi_cpp = logcond_psi_block_cpp(psi, psi_new_vec, w, Kw, aphi, bphi,
                                                                distmat, nu, ga,
                                                                agvec, FALSE, regionlabel,
                                                                phiidx);
                List logcon_old_psi_cpp = logcond_psi_block_cpp(psi, psi_new_vec, w, Kw, aphi, bphi,
                                                                distmat, nu, ga,
                                                                agvec, TRUE, regionlabel,
                                                                phiidx);
                // Rcout << "check psi_3" << endl;
                if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) -
                    as<double>(logcon_old_psi_cpp["value"]))) {
                    psi.elem(phiidx) = psi_new;
                    Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
                    Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
                    keeppsi_block[j] = as<int>(keeppsi_block[j]) + 1;
                    keeptmppsi_block[j] = as<int>(keeptmppsi_block[j]) + 1;
                }
                
                // Rcout << "phiidx" << phiidx << endl;
                // Rcout << "Psi" << psi.elem(phiidx) << endl;
                // Psimat(i, span(phiidx[0], phiidx[sqrtG - 1])) =
                //     mat(psi.elem(phiidx)).t();
                Psimat.row(i) = mat(psi).t();
                // Psimat(i, span(0, G - 1)) =
                //     mat(psi.elem(phiidx)).t();
            }
        }

        phi = exp(psi);

        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                draws(d - 1, span(Nvoxels + G, draws.n_cols - 1)) = phi.t();
            } 
        }
        if (fmod(i , 1000) == 0) {
            Rcout << "MCMC Iter: " << i + 1 << "\r" << std::flush;
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
        _["accept_keep_psi_b"] = keeppsi_block, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_pc1_phi_block_main_cpp_vec(arma::mat Yvec, arma::mat Xr, 
                                 arma::vec w, arma::vec ga, arma::vec phi, arma::vec psi,
                                 int nmcmc, int burn, int thin, 
                                 bool adapt, arma::vec tunevec, arma::vec keepvec, 
                                 arma::vec tunepsivec, arma::vec keeppsivec, 
                                 List tunepsilst_block, List keeppsilst_block,
                                 int tunelen, 
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
    const int p = 1;
    int G = regionlabel.max();
    arma::mat XtX = Xr.t() * Xr;
    arma::mat XtY = Xr.t() * Yvec.t();
    
    // Rcout << "Loading Data" << endl;
    arma::mat Ker = bezier2_phi_cpp(distmat, nu, phi);
    arma::mat Kw = Ker * w;
    
    
    List MM = marginal_yga_cpp(Yvec, XtY, XtX);
    arma::vec M0 = MM["M0"];
    arma::vec M1 = MM["M1"];
    arma::mat Wmat(nmcmc, G);
    arma::mat w_new;
    arma::uvec widx;
    arma::mat Psimat(nmcmc, G);
    arma::mat psi_new;
    arma::uvec phiidx;
    
    // #-------------------------------
    // # Adaptive tuning
    // #-------------------------------
    // # for w
    arma::vec keep = keepvec;
    arma::vec keeptmp = keep;  // track MH accpetance rate for adaptive tuning
    
    // # for psi
    // # keep.psi <- list(psi = 0)
    arma::vec keeppsi = keeppsivec;
    arma::vec keeptmppsi = keeppsi;
    
    Rcpp::List keeppsi_block = keeppsilst_block;
    Rcpp::List keeptmppsi_block = keeppsi_block;
    
    int Tb = tunelen;  //# frequency of adaptive tuning
    
    
    
    // #-------------------------------
    // # M-H algorithm
    // #------------------------------- 
    // Begin algorithm
    Rcout << "CV-KC-phig Begins" << endl;
    
    for (int i = 0; i < nmcmc; ++i) {
        
        // #----------------------------
        // # Update parameters
        // #----------------------------
        
        // # ========== sample w =============
        int sqrtG = sqrt(G);
        for(int j = 0; j < sqrtG; ++j) {
            // Rcout << "check w" << endl;
            // # Update tuning parameter
            // #-------------------------
            if(adapt == TRUE & fmod(i , Tb) == 0) {
                // # Adaptive tuning
                keeptmp[j] = keeptmp[j] / Tb;
                tunevec[j] = callgettune(tunevec[j], keeptmp[j], i, get_tune);
                keeptmp[j] = 0;
            }
            // Rcout << "check w_1" << endl;
            widx = sqrtG * j + linspace<uvec>(0, sqrtG - 1, sqrtG);
            // Rcout << " widx " << widx.n_elem << endl;
            if (i > 1000) {
                // cov(W[(i-1000):(i-1), w.idx])
                arma::mat Sigma_w = cov(Wmat(span(i - 1000, i - 1), 
                                             span(widx[0], widx[sqrtG - 1])));
                Sigma_w = (Sigma_w + Sigma_w.t()) / 2;
                arma::mat sigma = tunevec[j] * Sigma_w;
                w_new = mvrnormArma(1, w.elem(widx), sigma);
                // Rcout << "w_new " << w_new.n_elem << endl;
            } else {
                w_new = mvrnormArma(1, w.elem(widx), eye<mat>(sqrtG,sqrtG));
                // w_new = mvrnormArma(1, w.elem(widx), 
                // as<double>(tunelst[j]) * eye<mat>(sqrtG,sqrtG));
                // Rcout << "w_new " << w_new.n_elem << endl;
            }
            // Rcout << "check w_2" << endl;
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
            // Rcout << "check w_3" << endl;
            if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_cpp["value"]) - 
                as<double>(logcon_old_cpp["value"]))) {
                w.elem(widx) = w_new;
                Kw = as<arma::mat>(logcon_new_cpp["Kw"]);
                keep[j] = keep[j] + 1;
                keeptmp[j] = keeptmp[j] + 1;
            }
            // Rcout << "check w_4" << endl;
            Wmat(i, span(widx[0], widx[sqrtG - 1])) = mat(w.elem(widx)).t();
        }
        
        // # ============ sample ga ==============
        for (int v = 0; v < Nvoxels; ++v) {
            ga[v] = ga_update_pc1_cpp(v, M0[v], M1[v], Kw, n, regionlabel, agvec, iscplx);
        }
        
        // // # ============ sample psi and hence phi ==============
        // // Rcout << "check psi" << endl;
        // for(int j = 0; j < sqrtG; ++j) {
        //     // # Update tuning parameter
        //     // #-------------------------
        //     if(adapt == TRUE & fmod(i , Tb) == 0) {
        //         // # Adaptive tuning
        //         // # Adaptive tuning
        //         keeptmppsi[j] = as<double>(keeptmppsi[j]) / Tb;
        //         tunepsilst[j] = callgettune(tunepsilst[j], keeptmppsi[j], i, get_tune);
        //         keeptmppsi[j] = 0;
        //     }
        //     // Rcout << "check psi_1" << endl;
        //     phiidx = sqrtG * j + linspace<uvec>(0, sqrtG - 1, sqrtG);
        //     // if (i == 100) {
        //     //     Rcout << " phiidx " << phiidx << "\r" << endl;
        //     // }
        //     if (i > 20000) {
        //         // cov(W[(i-1000):(i-1), w.idx])
        //         Rcout << cov(Psimat(span(i - 2000, i - 1),
        //                             span(phiidx[0], phiidx[sqrtG - 1]))) << endl;
        //         arma::mat Sigma_phi = cov(Psimat(span(i - 2000, i - 1),
        //                                          span(phiidx[0], phiidx[sqrtG - 1])));
        //         Sigma_phi = (Sigma_phi + Sigma_phi.t()) / 2;
        //         arma::mat sigma = as<double>(tunepsilst[j]) * Sigma_phi;
        //         psi_new = mvrnormArma(1, psi.elem(phiidx), sigma);
        //         // Rcout << "w_new " << w_new.n_elem << endl;
        //     } else {
        // 
        //         psi_new = mvrnormArma(1, psi.elem(phiidx),
        //                               as<double>(tunepsilst[j]) * eye<mat>(sqrtG,sqrtG));
        //         // Rcout << "check psi_2" << endl;
        //         // Rcout << "w_new " << w_new.n_elem << endl;
        //     }
        // 
        //     arma::vec psi_new_vec = vectorise(psi_new);
        //     // Rcout << "w_new_vec" << w_new_vec.n_elem << endl;
        //     List logcon_new_psi_cpp = logcond_psi_block_cpp(psi, psi_new_vec, w, Kw, aphi, bphi,
        //                                                     distmat, nu, ga,
        //                                                     agvec, FALSE, regionlabel,
        //                                                     phiidx);
        //     List logcon_old_psi_cpp = logcond_psi_block_cpp(psi, psi_new_vec, w, Kw, aphi, bphi,
        //                                                     distmat, nu, ga,
        //                                                     agvec, TRUE, regionlabel,
        //                                                     phiidx);
        //     // Rcout << "check psi_3" << endl;
        //     if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) -
        //         as<double>(logcon_old_psi_cpp["value"]))) {
        //         psi.elem(phiidx) = psi_new;
        //         Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
        //         Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
        //         keeppsi[j] = as<int>(keeppsi[j]) + 1;
        //         keeptmppsi[j] = as<int>(keeptmppsi[j]) + 1;
        //     }
        // 
        //     // Rcout << "phiidx" << phiidx << endl;
        //     // Rcout << "Psi" << psi.elem(phiidx) << endl;
        //     Psimat(i, span(phiidx[0], phiidx[sqrtG - 1])) =
        //         mat(psi.elem(phiidx)).t();
        // }
        // phi = exp(psi);
        
        
        if (i < nmcmc) {
            for(int g = 0; g < G; ++g) {
                // # Update tuning parameter of psi
                // #-------------------------
                if(adapt == TRUE & fmod(i , Tb) == 0) {
                    // # Adaptive tuning
                    keeptmppsi[g] = keeptmppsi[g] / Tb;
                    tunepsivec[g] = callgettune(tunepsivec[g], keeptmppsi[g], i, get_tune);
                    keeptmppsi[g] = 0;
                }
                
                double new_psig = rnorm(1, psi[g], tunepsivec[g])[0];
                arma::vec new_psi = psi;
                new_psi[g] = new_psig;
                
                List logcon_new_psi_cpp = logcond_psig_cpp(new_psi, w, Kw, aphi, bphi,
                                                           distmat, nu, ga, agvec, FALSE,
                                                           regionlabel, g);
                List logcon_old_psi_cpp = logcond_psig_cpp(psi, w, Kw, aphi, bphi,
                                                           distmat, nu, ga, agvec, TRUE,
                                                           regionlabel, g);
                
                
                if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) -
                    as<double>(logcon_old_psi_cpp["value"]))) {
                    psi[g] = new_psig;
                    Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
                    Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
                    keeppsi[g] = keeppsi[g] + 1;
                    keeptmppsi[g] = keeptmppsi[g] + 1;
                }
                Psimat.row(i) = mat(psi).t();
            }
        } else {
            int B = 1;
            for(int j = 0; j < B; ++j) {
                // # Update tuning parameter
                // #-------------------------
                if(adapt == TRUE & fmod(i , Tb) == 0) {
                    // # Adaptive tuning
                    // # Adaptive tuning
                    keeptmppsi_block[j] = as<double>(keeptmppsi_block[j]) / Tb;
                    tunepsilst_block[j] = callgettune(tunepsilst_block[j], keeptmppsi_block[j], i, get_tune);
                    keeptmppsi_block[j] = 0;
                }
                // Rcout << "check psi_1" << endl;
                // phiidx = B * j + linspace<uvec>(0, B - 1, B);
                phiidx = linspace<uvec>(0, G - 1, G);
                //                 Rcout << cov(Psimat(span(i - 2000, i - 1),
                //                                     span(phiidx[0], phiidx[sqrtG - 1]))) << endl;
                // arma::mat Sigma_phi = cov(Psimat(span(i - 2000, i - 1),
                //                                  span(phiidx[0], phiidx[sqrtG - 1])));
                arma::mat Sigma_phi = cov(Psimat(span(i - 2000, i - 1),
                                                 span(0, G - 1)));
                Sigma_phi = (Sigma_phi + Sigma_phi.t()) / 2;
                arma::mat sigma = as<double>(tunepsilst_block[j]) * Sigma_phi;
                // psi_new = mvrnormArma(1, psi.elem(phiidx), sigma);
                psi_new = mvrnormArma(1, psi, sigma);
                
                arma::vec psi_new_vec = vectorise(psi_new);
                // Rcout << "w_new_vec" << w_new_vec.n_elem << endl;
                List logcon_new_psi_cpp = logcond_psi_block_cpp(psi, psi_new_vec, w, Kw, aphi, bphi,
                                                                distmat, nu, ga,
                                                                agvec, FALSE, regionlabel,
                                                                phiidx);
                List logcon_old_psi_cpp = logcond_psi_block_cpp(psi, psi_new_vec, w, Kw, aphi, bphi,
                                                                distmat, nu, ga,
                                                                agvec, TRUE, regionlabel,
                                                                phiidx);
                // Rcout << "check psi_3" << endl;
                if (log(runif(1, 0, 1)[0]) < (as<double>(logcon_new_psi_cpp["value"]) -
                    as<double>(logcon_old_psi_cpp["value"]))) {
                    psi.elem(phiidx) = psi_new;
                    Ker = as<arma::mat>(logcon_new_psi_cpp["K"]);
                    Kw  = as<arma::mat>(logcon_new_psi_cpp["Kw"]);
                    keeppsi_block[j] = as<int>(keeppsi_block[j]) + 1;
                    keeptmppsi_block[j] = as<int>(keeptmppsi_block[j]) + 1;
                }
                
                // Rcout << "phiidx" << phiidx << endl;
                // Rcout << "Psi" << psi.elem(phiidx) << endl;
                // Psimat(i, span(phiidx[0], phiidx[sqrtG - 1])) =
                //     mat(psi.elem(phiidx)).t();
                Psimat.row(i) = mat(psi).t();
                // Psimat(i, span(0, G - 1)) =
                //     mat(psi.elem(phiidx)).t();
            }
        }
        
        phi = exp(psi);
        
        
        // #  Save samples 
        // # -----------------------------------------------------
        if ((i + 1) > burn) {
            if (fmod(i + 1 , thin) == 0) {
                int d = (i + 1 - burn) / thin;
                draws(d - 1, span(0, Nvoxels - 1)) = ga.t();
                draws(d - 1, span(Nvoxels, Nvoxels + G - 1)) = w.t();
                draws(d - 1, span(Nvoxels + G, draws.n_cols - 1)) = phi.t();
            } 
        }
        if (fmod(i , 1000) == 0) {
            Rcout << "MCMC Iter: " << i + 1 << "\r" << std::flush;
        }
    }
    // draws = join_rows(join_rows(Gamma, W), Phi);
    
    Rcout << "Done! " << endl;
    
    // # Write output
    // #--------------
    return List::create(
        _["draws"] = draws, 
        _["accept_keep"] = keep / nmcmc, 
        _["accept_keep_psi"] = keeppsi / nmcmc, 
        _["accept_keep_psi_b"] = keeppsi_block, 
        _["burn"] = burn, 
        _["thin"] = thin, 
        _["nmcmc"] = nmcmc
    );
}







