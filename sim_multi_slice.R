################################################################################
# Multiple-slice analysis                                                      #
# ./sim_multi_slice.R                                                          #
# Cheng-Han Yu, Marquette University                                           #
################################################################################
## Load packages
source("./load_packages.R")

## Load functions
source("./my_fcns.R")
source("./mcmc_cpp_fcns.R")

## Load data
source("./data/load_multi_slice_data.R")
# Xmat, multi_slice_data_im, multi_slice_data_re
source("./data/clean_multi_slice_data.R")


################################################################################
# Know the data: Plotting
################################################################################
par(mfrow = c(2, 2))
par(mar = c(.1, .1, 1, .1))
# image at t = 1
for (k in 1:7) {
    slice_k <- lapply(multi_slice_data_lst, function(x) x[[k]])
    invisible(lapply(names(slice_k), img_time_plot, time = 1, 
                     datalist = slice_k, slice_no = k))
}
for (k in 1:7) {
    slice_k <- lapply(multi_slice_data_mean_lst, function(x) x[[k]])
    invisible(lapply(names(slice_k), img_plot,
                     datalist = slice_k, slice_no = k))
}

# Activated v.s. Nonactivated voxel
par(mfrow = c(2, 4))
par(mar = rep(2, 4))
voxel_plot <- function(dataM, dataA, dataR, dataI, site1, site2) {
    par(mfrow = c(2, 4))
    par(mar = rep(2, 4))
    hist(dataM[, site1, site2], breaks = 100, 
         main = paste("y(", site1, ", ", site2, ") Mod", sep = ""))
    hist(dataA[, site1, site2], breaks = 100, main = "Arg")
    hist(dataR[, site1, site2], breaks = 100, main = "Re")
    hist(dataI[, site1, site2], breaks = 100, main = "Im") 
    plot(dataM[, site1, site2], type = "l",
         main = paste("time series y(", site1, ", ", site2, ") Mod", sep = ""))
    plot(dataA[, site1, site2], type = "l")
    plot(dataR[, site1, site2], type = "l")
    plot(dataI[, site1, site2], type = "l")
}
# an activated voxel at (33, 40)
for (i in 43:47) {
    for (j in 43:47) {
        voxel_plot(dataM_lst[[1]], dataA_lst[[1]], 
                   dataR_lst[[1]], dataI_lst[[1]], i, j)
        title(paste("site = (", i, j, ")"))
    }
}

################################################################################
# generate the covariate and reference function WITH time trend 
################################################################################
n_sim2 <- 510
N_sim2 <- 96
Nvoxels_sim2 <- 96 * 96
ones_sim2 <- rep(1, n_sim2 - 20) # intercept
x1_sim2 <- 1:(n_sim2-20)  # trend

totaltime_sim2 <- n_sim2
onsets_sim2 <- seq(30, n_sim2, by = 30)
dur_sim2 <- 15
s_sim2 <- neuRosim::stimfunction(totaltime = totaltime_sim2, 
                                 onsets = onsets_sim2, 
                                 durations = dur_sim2, accuracy = 1)
# HRF
canonical_sim2 <- neuRosim::specifydesign(totaltime = totaltime_sim2, 
                                          onsets = list(onsets_sim2),
                                          durations = list(dur_sim2), 
                                          effectsize = 1, TR = 1, 
                                          conv = "double-gamma")
par(mfrow = c(1, 1))
plot(canonical_sim2, type = 'l', main = "BOLD signal w/ canonical HRF", 
     xlab = "time", lwd = 2)
# design matrix
x_sim2 <- cbind(ones_sim2, x1_sim2 / (n_sim2 - 20), canonical_sim2[21:n_sim2])
X_sim2 <- as.matrix(bdiag(x_sim2, x_sim2))
x_re_sim2 <- canonical_sim2[21:n_sim2]
# first dimension is time
# combine real and imaginary part by time

Y_lst <- vector("list", length = 7)

for (k in 1:7) {
    Y_lst[[k]] <- abind::abind(dataR_lst[[k]], dataI_lst[[k]], along = 1)
}


# centering
xstar_sim2 <- cbind(ones_sim2, x1_sim2)
Xstar_sim2 <- X_sim2[, c(1, 2, 4, 5)]
HX <- emulator::quad.tform.inv(crossprod(Xstar_sim2), x = Xstar_sim2)
Hx <- emulator::quad.tform.inv(crossprod(xstar_sim2), x = xstar_sim2)

IIc_sim2 <- (diag(2*(n_sim2-20)) - HX)
II_sim2 <- (diag(n_sim2-20) - Hx)
X_cplx_sim2 <- (diag(2*(n_sim2-20)) - HX) %*% X_sim2[, c(3, 6)]
x_re_sim2 <- (diag(n_sim2-20) - Hx) %*% x_re_sim2
x_mod_sim2 <- x_re_sim2
Y_mod_lst <- dataM_lst

for (k in 1:7) {
    for (i in 1:N_sim2) {
        for (j in 1:N_sim2) {
            Y_lst[[k]][, i, j] <- IIc_sim2 %*% Y_lst[[k]][, i, j]
            Y_mod_lst[[k]][, i, j] <- II_sim2 %*% Y_mod_lst[[k]][, i, j]
        }
    }
}


################################################################################
# Activation regions and Parcellation
################################################################################
multi_slice_beta1 <- list()
multi_slice_gamma1 <- list()
for (i in 1:7) {
    beta1 <- paste0("beta1s", i, ".txt")
    path_beta1 <- paste0("./data/", beta1)
    multi_slice_beta1[[i]] <- read.table(path_beta1, sep = ",")
    
    gamma1 <- paste0("gamma1s", i, ".txt")
    path_gamma1 <- paste0("./data/", gamma1)
    multi_slice_gamma1[[i]] <- read.table(path_gamma1, sep = ",")
}

multi_slice_beta1 <- lapply(multi_slice_beta1, rotate)
multi_slice_gamma1 <- lapply(multi_slice_gamma1, rotate)

par(mfrow = c(1, 1))
par(omi = c(0, 0, 0, 0))
par(mar = c(.1, .1, 1, .1))
xx = 1.006
yy = -0.006
xlimit = c(yy, xx)

max_beta <- max(unlist(multi_slice_beta1))
max_gamma <- max(unlist(multi_slice_gamma1))
par(mfrow = c(2, 7))
for (k in 1:7) {
    image(multi_slice_data_mean_lst$Mod[[k]], 
          col = gray((0:96) / 96), axes = FALSE,
          main = paste("Beta1 --", k), xlim = xlimit,
          ylim = xlimit)
    # two activation regions
    for (i in 43:47) {
        for (j in 43:47) {
            alpha <- floor(255 * (multi_slice_beta1[[k]][i, j] / max_beta))
            points(((i)/96)-0.006, ((j)/96)-0.006, 
                   col = rgb(255, 0, 0, max = 255, alpha = alpha),
                   pch = 15, cex = 0.8)
        }
    }
    for (i in 31:35) {
        for (j in 38:42) {
            alpha <- floor(255 * (multi_slice_beta1[[k]][i, j] / max_beta))
            points(((i)/96)-0.006, ((j)/96)-0.006, 
                   col = rgb(255, 0, 0, max = 255, alpha = alpha),
                   pch = 15, cex = 0.8)
        }
    }
}
for (k in 1:7) {
    image(multi_slice_data_mean_lst$Mod[[k]], 
          col = gray((0:96) / 96), axes = FALSE,
          main = paste("Gamma1 --", k), xlim = xlimit,
          ylim = xlimit)
    # two activation regions
    for (i in 43:47) {
        for (j in 43:47) {
            alpha <- floor(255 * (multi_slice_gamma1[[k]][i, j] / max_gamma))
            points(((i)/96)-0.006, ((j)/96)-0.006, 
                   col = rgb(0, 255, 0, max = 255, alpha = alpha),
                   pch = 15, cex = 0.8)
        }
    }
    for (i in 31:35) {
        for (j in 38:42) {
            alpha <- floor(255 * (multi_slice_gamma1[[k]][i, j] / max_gamma))
            points(((i)/96)-0.006, ((j)/96)-0.006, 
                   col = rgb(0, 255, 0, max = 255, alpha = alpha),
                   pch = 15, cex = 0.8)
        }
    }
}

## grid examples
## 4 x 4 grid
for (i in 1:3) {
    xline(yy + (xx/4)*i, col = "green")
    yline(yy + (xx/4)*i, col = "green")
}

## 6 x 6 grid
for (i in 1:5) {
    xline(yy + (xx/6)*i, col = "yellow")
    yline(yy + (xx/6)*i, col = "yellow")
}

## 32 x 32 grid
for (i in 1:31) {
    xline(yy + (xx/32)*i, col = "yellow")
    yline(yy + (xx/32)*i, col = "yellow")
}


## activation binary indiators
final_slice <- matrix(0, 96, 96)
final_slice[43:47, 43:47] <- 1
final_slice[31:35, 38:42] <- 1

image(final_slice)
final_slice_tf <- final_slice > 0

########################
# Cut the region outside the brain
#########################
cut.row.idx <- c(1:17, 82:96)
cut.col.idx <- c(1:17, 82:96)

multi_slice_data_mean_lst$Mod[[k]][-cut.row.idx, -cut.col.idx]

multi_slice_beta1_cut <- lapply(multi_slice_beta1, 
                                function(x) x[-cut.row.idx, -cut.col.idx])

final_slice_tf_cut <- final_slice_tf[-cut.row.idx, -cut.col.idx]

dim(multi_slice_beta1_cut[[2]])
# dimension 64 x 64
par(mfrow = c(1, 7))
par(omi = c(0, 0, 0, 0))
par(mar = c(.1, .1, 1, .1))
for (k in 1:7) {
    image(multi_slice_data_mean_lst$Mod[[k]][-cut.row.idx, -cut.col.idx], 
          col = gray((0:96) / 96), axes = FALSE,
          main = paste(k), xlim = xlimit,
          ylim = xlimit)
    # two activation regions
    for (i in 14:18) {
        for (j in 21:25) {
            alpha <- floor(255 * (multi_slice_beta1_cut[[k]][i, j] / max_beta))
            points(((i)/64)-0.006, ((j)/64)-0.006, 
                   col = rgb(255, 0, 0, max = 255, alpha = alpha),
                   pch = 15, cex = 0.8)
        }
    }
    for (i in 26:30) {
        for (j in 26:30) {
            alpha <- floor(255 * (multi_slice_beta1_cut[[k]][i, j] / max_beta))
            points(((i)/64)-0.006, ((j)/64)-0.006,
                   col = rgb(255, 0, 0, max = 255, alpha = alpha),
                   pch = 15, cex = 0.8)
        }
    }
    for (i in 1:7) {
        yline(yy + (xx/8)*i, col = "green", lwd = 0.2)
    }
    for (i in 1:7) {
        xline(yy + (xx/8)*i, col = "green", lwd = 0.2)
    }
    
    loc.y <- seq(0.05, 0.98, by = 0.125)
    loc.x <- seq(0.05, 0.98, by = 0.125)
    for (i in 1:8) {
        for (j in 1:8) {
            points(loc.x[i], loc.y[j], pch = 16, col = "gold", cex = 0.8)
        } 
    }
}
########################
# Label region numbers
########################
N_cut_x <- 64
N_cut_y <- 64
Nvoxels_cut <- N_cut_x * N_cut_y
regionlabel_cut <- rep(0, Nvoxels_cut)
for (k in 8) {
    for (i in 1:(N_cut_x / k)) {
        for(j in 1:(N_cut_y / k)) {
            for(m in 1:k) {
                regionlabel_cut[c((j - 1) * N_cut_x * k + 
                                      N_cut_x * (m - 1) + 1:k) + (i - 1) * k] = 
                    i + (j - 1) * (N_cut_x / k)
            }
        }
    }
}

region_mat_cut <- matrix(regionlabel_cut, ncol = N_cut_y)
par(mfrow = c(1, 1))
image.plot(region_mat_cut, axes = FALSE)
G_cut <- max(region_mat_cut)  # Number of groups

##########################
# cut the margin for the data
##########################
Y_vec_lst <- vector("list", length = 7)
Y_mod_vec_lst <- vector("list", length = 7)
Y_cplx_lst <- vector("list", length = 7)
Y_cplx_vec_lst <- vector("list", length = 7)

for (k in 1:7) {
    # dim(Y_human)
    Y_lst[[k]] <- aperm::aperm(Y_lst[[k]], c(2, 3, 1))
    # dim(Y_human)
    Y_lst[[k]] <- Y_lst[[k]][-cut.row.idx, -cut.col.idx, ]
    
    Y_mod_lst[[k]] <- aperm::aperm(Y_mod_lst[[k]], c(2, 3, 1))
    Y_mod_lst[[k]] <- Y_mod_lst[[k]][-cut.row.idx, -cut.col.idx, ]
    
    # dim(Y_human)
    Y_vec_lst[[k]] <- matrix(Y_lst[[k]], N_cut_x * N_cut_y, 980)
    Y_mod_vec_lst[[k]] <- matrix(Y_mod_lst[[k]], N_cut_x * N_cut_y, 490)
    
    Y_cplx_lst[[k]] <- makecplx(Y_lst[[k]])
    Y_cplx_vec_lst[[k]] <- matrix(Y_cplx_lst[[k]] , nrow = N_cut_x * N_cut_y, 
                                  ncol = 490)
}


####################################################################
# CV-KC
####################################################################
# initial tuning parameter
p <- 1
sqrtG <- sqrt(G_cut)
listname <- c(paste("w.", 1:sqrtG, sep = ""))
keep.lst <- sapply(rep(0, sqrtG), as.list)
tune.lst <- sapply(rep(1, sqrtG), as.list)
names(keep.lst) <- listname
names(tune.lst) <- listname
name.par <- c(paste0("gamma", 1:(p * Nvoxels_cut)), 
              paste0("w", 1:(p * G_cut)), "phi")
keep.psi <- list(psi = 0)
tune.psi <- list(psi = 1)
start.lst <- list(ga = rep(1, Nvoxels_cut), 
                  w = t(mvtnorm::rmvnorm(1, rep(0, G_cut), .5 * diag(G_cut))),
                  phi = 10)


####################################################################
##### simple multi-core computing
####################################################################
library(parallel)
library(doParallel)
detectCores()
cl <- makeCluster(8, type = "FORK")
registerDoParallel(cl)
getDoParWorkers()

####################################################################
## CV-KC
####################################################################
n.mcmc <- 51000
burn <- 1000
thin <- 25
MCMC_KC_mc <- function(s, iscplx = TRUE) {
    if(iscplx) {
        Yvec <- Y_vec_lst
        Xr <- X_cplx_sim2
    } else {
        Yvec <- Y_mod_vec_lst   ## take L2 norm from the complex-valued data
        # Yvec <- YYmagvec   ## simulate as real-valued data
        Xr <- x_mod_sim2
    }
    mcmc_kcs_block_cpp_full(Yvec[[s]], Xr, start.lst = start.lst, 
                            nmcmc = n.mcmc, burn = burn, thin = thin, 
                            name.par = name.par, adapt = TRUE, 
                            tunelst = tune.lst, tunepsi = tune.psi, 
                            keeplst = keep.lst, keeppsi = keep.psi, 
                            tunelen = 50, regionlabel = regionlabel_cut,
                            region_mat = region_mat_cut, nu = 2,
                            atau = 1/2, btau = 1/2,
                            aphi = 1/2, bphi = 1/2,
                            agvec = rep(0, G_cut),
                            target.accept.rate = 0.35, iscplx = iscplx,
                            kerneltype = "bezier")
}

## 2044.744 
system.time(CV_KC_lst <- parallel::mclapply(2:6, MCMC_KC_mc, 
                                            iscplx = TRUE, mc.cores = 8))

## MO-KC
# 2058.121 
system.time(MO_KC_lst <- parallel::mclapply(2:6, MCMC_KC_mc, 
                                            iscplx = FALSE, mc.cores = 8))


####################################################################
## CV-GP
####################################################################
listname <- c(paste("S.", 1:G_cut, sep = ""))
keep.lst <- sapply(rep(0, G_cut), as.list)
tune.lst <- sapply(rep(1, G_cut), as.list)
names(keep.lst) <- listname
names(tune.lst) <- listname
name.par <- c(paste0("gamma", 1:(p * Nvoxels_cut)), 
              paste0("S", 1:(p * G_cut)), "r")
keep.r <- list(r = 0)
tune.r <- list(r = 1)
start.lst <- list(ga = rep(1, Nvoxels_cut), 
                  S = t(mvtnorm::rmvnorm(1, rep(0, G_cut), .5 * diag(G_cut))), 
                  r = 5)

MCMC_GP_mc <- function(s, iscplx = TRUE) {
    if(iscplx) {
        Yvec <- Y_vec_lst
        Xr <- X_cplx_sim2
    } else {
        Yvec <- Y_mod_vec_lst   ## take L2 norm from the complex-valued data
        # Yvec <- YYmagvec   ## simulate as real-valued data
        Xr <- x_mod_sim2
    }
    mcmc_gp_cpp_full(Yvec[[s]], Xr, start.lst = start.lst, 
                     nmcmc = n.mcmc, burn = burn, thin = thin, 
                     name.par = name.par, adapt = TRUE, 
                     tunelst = tune.lst, keeplst = keep.lst, 
                     tuner = tune.r, keepr = keep.r,
                     tunelen = 50, regionlabel = regionlabel_cut, 
                     region_mat = region_mat_cut, df = 8, 
                     agvec = rep(0, G_cut),
                     target.accept.rate = 0.35, iscplx = iscplx)
}

# 512.13
system.time(CV_GP_lst <- parallel::mclapply(2:6, MCMC_GP_mc, 
                                            iscplx = TRUE, mc.cores = 8))

## MO-GP
# 422.795 
system.time(MO_GP_lst <- parallel::mclapply(2:6, MCMC_GP_mc, 
                                            iscplx = FALSE, mc.cores = 8))
stopCluster(cl)




####################################################################
## MCMC results
####################################################################
sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
plot_parameter(MCMCsample = CV_KC_lst[[4]], draws = NULL, 
               parameter.idx = 4097:4099, 
               sampleidx = sampleidx, name.par = name.par)
par(mfrow = c(2, 2)) 
# par(mfrow = c(1, 1))
par(omi = c(0, 0, 0, 0))
par(mar = c(.5, .5, 1, .5))
for (j in 1:5) {
    posprob_kcs_sim2 <- compute_posprob_pc(CV_KC_lst[[j]], N.voxels = Nvoxels_cut)
    
    brain_grid_plot(multi_slice_data_mean_lst$Mod[[j]][-cut.row.idx, -cut.col.idx], 
                    posprob_kcs_sim2, N = 64, G = 64, 
                    mdltype = "", thre = 0.8722, is.color.bar = TRUE,
                    title = paste("CV-KC Posterior Prob Slice", j + 1))
    posprob_gp_sim2 <- compute_posprob_pc(CV_GP_lst[[j]], N.voxels = Nvoxels_cut)
    
    brain_grid_plot(multi_slice_data_mean_lst$Mod[[j]][-cut.row.idx, -cut.col.idx], 
                    posprob_gp_sim2, N = 64, G = 64, 
                    mdltype = "", thre = 0.8722, 
                    title = paste("CV-GP Posterior Prob Slice", j + 1))
    
    posprob_kcs_sim2 <- compute_posprob_pc(MO_KC_lst[[j]], N.voxels = Nvoxels_cut)
    
    brain_grid_plot(multi_slice_data_mean_lst$Mod[[j]][-cut.row.idx, -cut.col.idx], 
                    posprob_kcs_sim2, N = 64, G = 64, 
                    mdltype = "", thre = 0.8722, 
                    title = paste("MO-KC Posterior Prob Slice", j + 1))
    posprob_gp_sim2 <- compute_posprob_pc(MO_GP_lst[[j]], N.voxels = Nvoxels_cut)
    
    brain_grid_plot(multi_slice_data_mean_lst$Mod[[j]][-cut.row.idx, -cut.col.idx], 
                    posprob_gp_sim2, N = 64, G = 64, 
                    mdltype = "", thre = 0.8722, 
                    title = paste("MO-GP Posterior Prob Slice", j + 1))
    
}



par(mfcol = c(4, 5)) 
# par(mfrow = c(1, 1))
par(omi = c(0, 0, 0, 0))
par(mar = c(.5, .5, 1, .5))

for (j in 1:5) {
    if (j == 5) {
        is.color.bar = TRUE
    } else {
        is.color.bar = FALSE
    }
    posprob_kcs_sim2 <- compute_posprob_pc(CV_KC_lst[[j]], N.voxels = Nvoxels_cut)
    
    brain_grid_plot(multi_slice_data_mean_lst$Mod[[j]][-cut.row.idx, -cut.col.idx], 
                    posprob_kcs_sim2, N = 64, G = 64, is.points = TRUE,
                    mdltype = "", thre = 0.8722, is.color.bar = TRUE,
                    title = paste("CV-KC Posterior Prob Slice", j + 1))
    posprob_gp_sim2 <- compute_posprob_pc(CV_GP_lst[[j]], N.voxels = Nvoxels_cut)
    
    brain_grid_plot(multi_slice_data_mean_lst$Mod[[j]][-cut.row.idx, -cut.col.idx], 
                    posprob_gp_sim2, N = 64, G = 64, is.grid = TRUE,
                    mdltype = "", thre = 0.8722, is.color.bar = is.color.bar,
                    title = paste("CV-GP Posterior Prob Slice", j + 1))
    
    posprob_kcs_sim2 <- compute_posprob_pc(MO_KC_lst[[j]], N.voxels = Nvoxels_cut)
    
    brain_grid_plot(multi_slice_data_mean_lst$Mod[[j]][-cut.row.idx, -cut.col.idx], 
                    posprob_kcs_sim2, N = 64, G = 64, is.points = TRUE,
                    mdltype = "", thre = 0.8722, is.color.bar = is.color.bar,
                    title = paste("MO-KC Posterior Prob Slice", j + 1))
    posprob_gp_sim2 <- compute_posprob_pc(MO_GP_lst[[j]], N.voxels = Nvoxels_cut)
    
    brain_grid_plot(multi_slice_data_mean_lst$Mod[[j]][-cut.row.idx, -cut.col.idx], 
                    posprob_gp_sim2, N = 64, G = 64, is.grid = TRUE,
                    mdltype = "", thre = 0.8722, is.color.bar = is.color.bar,
                    title = paste("MO-GP Posterior Prob Slice", j + 1))
    
}


### compute rate
thre.vec <- c(0.5, 0.8722)
no.slice <- 5
Rate_CV_GP_slice <- array(NA, dim = c(no.slice, length(thre.vec), 10))
Rate_MO_GP_slice <- array(NA, dim = c(no.slice, length(thre.vec), 10))
Rate_CV_KC_slice <- array(NA, dim = c(no.slice, length(thre.vec), 10))
Rate_MO_KC_slice <- array(NA, dim = c(no.slice, length(thre.vec), 10))


for (j in 1:5) {
    posprob_kcs_sim2 <- compute_posprob_pc(CV_KC_lst[[j]], N.voxels = Nvoxels_cut)
    Rate_CV_KC_slice[j, , ] <- as.matrix(tpfpDatafcn_mcc(thre.vec, posprob_kcs_sim2, final_slice_tf_cut))
    
    posprob_gp_sim2 <- compute_posprob_pc(CV_GP_lst[[j]], N.voxels = Nvoxels_cut)
    Rate_CV_GP_slice[j, , ] <- as.matrix(tpfpDatafcn_mcc(thre.vec, posprob_gp_sim2, final_slice_tf_cut))
    
    
    posprob_kcs_sim2 <- compute_posprob_pc(MO_KC_lst[[j]], N.voxels = Nvoxels_cut)
    Rate_MO_KC_slice[j, , ] <- as.matrix(tpfpDatafcn_mcc(thre.vec, posprob_kcs_sim2, final_slice_tf_cut))
    
    posprob_gp_sim2 <- compute_posprob_pc(MO_GP_lst[[j]], N.voxels = Nvoxels_cut)
    Rate_MO_GP_slice[j, , ] <- as.matrix(tpfpDatafcn_mcc(thre.vec, posprob_gp_sim2, final_slice_tf_cut))
}

c('Sensitivity', "FNR", 'FPR', 'Specificity', "Precision", "FDR", "Accuracy", "F1", "MCC", "Thre")

(Rate_CV_KC_slice_mean <- apply(Rate_CV_KC_slice, c(2, 3), mean))
(Rate_CV_GP_slice_mean <- apply(Rate_CV_GP_slice, c(2, 3), mean))
(Rate_MO_KC_slice_mean <- apply(Rate_MO_KC_slice, c(2, 3), mean))
(Rate_MO_GP_slice_mean <- apply(Rate_MO_GP_slice, c(2, 3), mean))










































