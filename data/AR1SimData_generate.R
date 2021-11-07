################################################################################
# Complex-valued fMRI Simulation Study: AR(1) Data Generating                  #
# "./data/AR1SimData_generate.R"                                               #
# Cheng-Han Yu, Marquette University                                           #
################################################################################
# Experimental design (same as the design without AR(1) structure)
# =====================
source("./data/SimExpDesign.R")

## function that create region labels (same as the design without AR(1) structure)
source("./data/regionlabel.R")

# AR coefficient, error variance and SNR, CNR
# =====================
no.data <- 100
# rho_vec <- c(0.1, 0.5, 0.9)
rho_vec <- c(0.5)
sig2_old <- 3
# SNR <- c(1, 300) / sqrt(sig2_old) # same as IID setting
SNR <- c(300) / sqrt(sig2_old) # same as IID setting
CNR <- c(2, 5) / sqrt(sig2_old) # same as IID setting
library(stringr)
rho_split <- str_split(rho_vec, pattern = "")
rho_char <- unlist(lapply(rho_split, function(x) x[3]))
no.rho <- length(rho_char)
len_snr_cnr <- length(SNR) * length(CNR)
which_snr_cnr <- c(1, 2) ## which (snr, cnr) are selected
no.YY <- length(which_snr_cnr) 
# Design matrix and regression coefficients with different (SNR, CNR)
# (same as the design without AR(1) structure) 

# will produce new sig2
# =====================
source("./data/SimDesignMatRegCoefAR.R")


# Generate data set
# =====================
for(r in 1:no.rho) {
    rho <- rho_vec[r]
    sig2 <- sig2_old * (1 - rho ^ 2)
    for (k in 1:no.data) {
        source("./data/AR1SimData.R")
        filename <- paste("./data/Spa_AR1_sim", k, 
                          "_rho", rho_char[r], ".RData", sep = "")
        print(filename)
        # save(YY, YYmod, file = filename)
        save(YY, YYmod, YYobs, YYvec, YYmodvec, YYobsvec, file = filename)
        # if (!exists(filename)) {
        #     # save(YY, YYmod, file = filename)
        #     save(YY, YYmod, YYobs, YYvec, YYmodvec, YYobsvec, file = filename)
        # }
    }
}