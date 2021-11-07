

################################################################################
# Load required packages and self-created functions
################################################################################
source("./load_packages.R")
source("./my_fcns.R")
source("./mcmc_fcns.R")
source("./mcmc_cpp_fcns.R")


##############################
## Simulated Data
##############################
# Experimental design (same as the design without AR(1) structure)
source("./data/SimExpDesign.R")

## function that create region labels (same as the design without AR(1) structure)
source("./data/regionlabel.R")

rho <- 0.5
sig2_old <- 3
SNR <- 300 / sqrt(sig2_old)
CNR <- 5 / sqrt(sig2_old)

# will produce new sig2
source("./data/SimDesignMatRegCoefAR.R")

## Load simulated data
load("./data/Spa_AR1_sim1_rho5.RData", verbose = TRUE)


##############################
## Compute the EB (MLE) of AR(1) coefficient for each voxel
##############################
EB_rho_arr <- array(0, dim = c(N, N))
for (j in 1:N) {
    for (i in 1:N){
        rho_re <- suppressWarnings(ar.mle(YY[[2]][i, j, 1:n], 
                                          aic = FALSE,
                                          order.max = 1, 
                                          demean = FALSE)$ar)
        rho_im <- suppressWarnings(ar.mle(YY[[2]][i, j, (n+1):(2*n)], 
                                          aic = FALSE,
                                          order.max = 1, 
                                          demean = FALSE)$ar)
        EB_rho_arr[i, j] <- (rho_re + rho_im) / 2
    }
}

EB_rho_arr_vec <- array(EB_rho_arr, dim = c(1, N*N))

##############################
# Objects that are used for all methods below
##############################
G.vec <- c(G_4, G_16, G_25, G_100)
N.voxels <- N * N
p <- 1
regionlabel.lst <- list(regionlabel_4 = regionlabel_4, 
                        regionlabel_16 = regionlabel_16,
                        regionlabel_25 = regionlabel_25,
                        regionlabel_100 = regionlabel_100)

regionmat.lst <- list(region_mat_4 = region_mat_4, 
                      region_mat_16 = region_mat_16,
                      region_mat_25 = region_mat_25,
                      region_mat_100 = region_mat_100)

active.g.lst <- list(act.4 = c(1, 2, 3, 4),
                     act.16 = c(1, 2, 4, 5, 6, 10, 11, 14, 15),
                     act.25 = c(1, 2, 5, 6, 7, 18, 19, 23, 24),
                     act.100 = c(2, 9, 10, 12, 13, 14, 19, 20,
                                 22, 23, 24, 33, 66, 67, 75, 76, 77, 78,
                                 85, 86, 87, 96))



##############################
# Algorithm
##############################

## parameters for CV-GP-AR
Keep.Lst.GP <- list()
Tune.Lst.GP <- list()
Start.Lst.GP <- list()
Name.Par.GP <- list()
for (i in 1:length(G.vec)) {
    listname <- c(paste("S.", 1:G.vec[i], sep = ""))
    keep.lst.GP <- sapply(rep(0, G.vec[i]), as.list)
    tune.lst.GP <- sapply(rep(1, G.vec[i]), as.list)
    names(keep.lst.GP) <- listname
    names(tune.lst.GP) <- listname
    Keep.Lst.GP[[i]] <- keep.lst.GP
    Tune.Lst.GP[[i]] <- tune.lst.GP
    Name.Par.GP[[i]] <- c(paste0("gamma", 1:(p * N.voxels)),
                          paste0("S", 1:(p * G.vec[i])), "r")
    Start.Lst.GP[[i]] <- list(ga = rep(1, N.voxels),
                              S = t(rmvnorm(1, rep(0, G.vec[i]), .5 * diag(G.vec[i]))),
                              r = 5)
}
keep.r.GP <- list(r = 0)
tune.r.GP <- list(r = 1)


## parameters for CV-KC-AR
Keep.Lst.KC <- list()
Tune.Lst.KC <- list()
Start.Lst.KC <- list()
Name.Par.KC <- list()
for (i in 1:length(G.vec)) {
    listname <- c(paste("S.", 1:G.vec[i], sep = ""))
    keep.lst.KC <- sapply(rep(0, G.vec[i]), as.list)
    tune.lst.KC <- sapply(rep(1, G.vec[i]), as.list)
    names(keep.lst.KC) <- listname
    names(tune.lst.KC) <- listname
    Keep.Lst.KC[[i]] <- keep.lst.KC
    Tune.Lst.KC[[i]] <- tune.lst.KC
    Name.Par.KC[[i]] <- c(paste0("gamma", 1:(p * N.voxels)),
                          paste0("w", 1:(p * G.vec[i])), "phi")
    Start.Lst.KC[[i]] <- list(ga = rep(1, N.voxels),
                              w = t(rmvnorm(1, rep(0, G.vec[i]), .5 * diag(G.vec[i]))),
                              phi = 10)
}
keep.psi.KC <- list(psi = 0)
tune.psi.KC <- list(psi = 1)


# mcmc settings
n.mcmc <- 51000
burn <- 1000
thin <- 25
# sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
## choose the number of knots
G_idx <- 3 # 1:G_4 2:G_16 3:G_25, 4:G_100


## running algorithms
gp_ar <- mcmc_gp_ar_cpp_full(YYvec[[2]], X_cplx, EB_rho_arr_vec, 
                             start.lst = Start.Lst.GP[[G_idx]], 
                             nmcmc = n.mcmc, 
                             burn = burn, thin = thin, name.par = Name.Par.GP[[G_idx]], 
                             adapt = TRUE, tunelst = Tune.Lst.GP[[G_idx]], 
                             tuner = tune.r.GP, keeplst = Keep.Lst.GP[[G_idx]],
                             keepr = keep.r.GP, tunelen = 50, 
                             regionlabel = regionlabel.lst[[G_idx]], 
                             region_mat = regionmat.lst[[G_idx]], df = 8, 
                             agvec = rep(0, G.vec[G_idx]),
                             target.accept.rate = 0.35, iscplx = TRUE)

kc_ar <- mcmc_kcs_ar_block_cpp_full(YYvec[[2]], X_cplx, EB_rho_arr_vec,
                                    start.lst = Start.Lst.KC[[G_idx]], 
                                    nmcmc = n.mcmc,
                                    burn = burn, thin = thin, 
                                    name.par = Name.Par.KC[[G_idx]],
                                    adapt = TRUE, tunelst = Tune.Lst.KC[[G_idx]],
                                    tunepsi = tune.psi.KC, 
                                    keeplst = Keep.Lst.KC[[G_idx]],
                                    keeppsi = keep.psi.KC, tunelen = 50,
                                    regionlabel = regionlabel.lst[[G_idx]],
                                    region_mat = regionmat.lst[[G_idx]], nu = 2,
                                    atau = 1/2, btau = 1/2,
                                    aphi = 1/2, bphi = 1/2,
                                    agvec = rep(0, G.vec[G_idx]),
                                    target.accept.rate = 0.35, iscplx = TRUE,
                                    kerneltype = "bezier")


###########################################################
### Evaluation Metrics
###########################################################
finaltf <- final > 0
thre_vec <- c(0.5, 0.8722)

pprob_gpar <- compute_posprob_pc(gp_ar)
Rate_GP_AR <- as.matrix(tpfpDatafcn_mcc(thre_vec, pprob_gpar, finaltf))
pprob_kcar <- compute_posprob_pc(kc_ar)
Rate_KC_AR <- as.matrix(tpfpDatafcn_mcc(thre_vec, pprob_kcar, finaltf))

###########################################################
### posterior probability map
###########################################################
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
image.plot(matrix(pprob_gpar, 20, 20), main = "CV-GP-AR", axes = FALSE)
image.plot(matrix(pprob_kcar, 20, 20), main = "CV-KC-AR", axes = FALSE)








