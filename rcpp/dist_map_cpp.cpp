//dist_mat_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat dist_mat_cpp(arma::mat coord, arma::mat grid) {
    int Nvoxels = coord.n_rows;
    int G = grid.n_rows;
    arma::mat distmat = zeros(Nvoxels, G); 
    for (int i=0; i < Nvoxels; i++) {
        for (int j=0; j < G; j++) {
            distmat(i, j) = sqrt(sum(square(coord.row(i) - grid.row(j))));
        }
    }  
    return (distmat);
}

