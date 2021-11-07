


# library(clusterGeneration)
# library(microbenchmark)
# library(compiler)
# library(inline)

pkgs <- c("mvtnorm", "Matrix", "pscl", "matrixcalc", "neuRosim", "pryr",
          "pROC", "tikzDevice", "abind", "fields", "Hmisc", 
          "Rcpp", "RcppArmadillo", "emulator", "abind")

lapply(pkgs, require, character.only = TRUE)
