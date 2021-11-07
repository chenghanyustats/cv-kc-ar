################################################################################
# Complex-valued fMRI data analysis                                            #
# ./sim_multi_slice.R                                                          #
# Cheng-Han Yu, Marquette University                                           #
################################################################################
################################################################################
# Load required packages and self-created functions
################################################################################
source("./load_packages.R")
source("./my_fcns.R")
source("./mcmc_cpp_fcns.R")

################################################################################
## Load data
################################################################################
load("./data/fmridata.RData", verbose = TRUE)


################################################################################
# generate the covariate and reference function 
################################################################################
n_human <- 510
N_human <- 96
Nvoxels_human <- N_human * N_human
ones_human <- rep(1, n_human - 20) # intercept
x1_human <- 1:(n_human - 20)  # trend
totaltime_human <- n_human
onsets_human <- seq(30, n_human, by = 30)
dur_human <- 15
s_human <- stimfunction(totaltime = totaltime_human, onsets = onsets_human, 
                        durations = dur_human, accuracy = 1)
# HRF
canonical_human <- specifydesign(totaltime = totaltime_human, 
                                 onsets = list(onsets_human),
                                 durations = list(dur_human), effectsize = 1, 
                                 TR = 1, conv = "double-gamma")
# design matrix
x_human <- cbind(ones_human, x1_human/(n_human-20), canonical_human[21:n_human])
X_human <- as.matrix(bdiag(x_human, x_human))

data_list_human = list('Re' = fmriR, 'Im' = fmriI, 'Mod' = fmriM, 
                       'Arg' = fmriA)
data_list_human = lapply(data_list_human, aperm, c(3, 1, 2))
dataM_f = dataA_f = dataR_f = dataI_f = array(0, c(n_human, N_human, N_human))
for (t in 1:510) {
    dataM_f[t, , ] = f(data_list_human$Mod[t, , ])
    dataA_f[t, , ] = f(data_list_human$Arg[t, , ])
    dataR_f[t, , ] = f(data_list_human$Re[t, , ])
    dataI_f[t, , ] = f(data_list_human$Im[t, , ])
}

# combine real and imaginary part by time
Y_human <- abind::abind(dataR_f[21:n_human, , ], dataI_f[21:n_human, , ], 
                        along = 1)
Xstar_human <- X_human[, c(1,2,4,5)]
H_human <- emulator::quad.tform.inv(crossprod(Xstar_human), Xstar_human)
IIc_human <- (diag(2*(n_human-20)) - H_human)
X_cplx_human <- (diag(2*(n_human-20)) - H_human) %*% X_human[, c(3, 6)]
for (i in 1:N_human) {
    for (j in 1:N_human) {
        # Complex Data
        Y_human[, i, j] <- IIc_human %*% Y_human[, i, j]
    }
}
x_re_human <- canonical_human[21:n_human]
xstar_human <- cbind(ones_human, x1_human)
h_human <- emulator::quad.tform.inv(crossprod(xstar_human), xstar_human)
x_re_human <- (diag(n_human -20) - h_human) %*% x_re_human 
x_mod_human  <- x_re_human 

################################################################################
# data visualization 
################################################################################
# image at t = 1
data_lst_human = list('Re' = dataR_f, 'Im' = dataI_f, 'Mod' = dataM_f, 
                      'Arg' = dataA_f)
img_time_plot <- function (x, time, datalist) {
    image.plot(datalist[[x]][time, , ], axes = F, 
               main = paste(x, "time = ", time))
}
par(mfrow = c(2, 2))
par(mar = c(1, 1, 2, 4))
invisible(lapply(names(data_list_human), img_time_plot, time = 1, 
                 datalist = data_lst_human))

# image of the mean over time
img_mean <- function(x) {
    apply(x, c(2, 3), mean)
}
data_mean_list_human <- lapply(data_lst_human, img_mean)
img_plot <- function (x, datalist) {
    image.plot(datalist[[x]], axes = F, 
               main = paste(x, "mean"))
}
img_plot <- function (x, datalist) {
    image(datalist[[x]], axes = F, col = tim.colors(),
          main = paste(x))
}
par(mfrow = c(2, 2))
par(mar = c(1, 1, 2, 1))
invisible(lapply(names(data_list_human), img_plot, 
                 datalist = data_mean_list_human))


################################################################################
# Parcellation
################################################################################
par(mfrow = c(1, 1))
par(omi = c(0, 0, 0, 0))
par(mar = c(.5, .5, 2, .5))
xx = 1.006
yy = -0.006
xlimit = c(yy, xx)
image(data_mean_list_human$Mod, col = gray((0:96) / 96), axes = FALSE,
      main = "True activation map", xlim = xlimit,
      ylim = xlimit)

for (i in 1:5) {
    xline(yy + (xx/6)*i, col = "green")
    yline(yy + (xx/6)*i, col = "green")
}

for (i in 1:4) {
    segments(yy + (xx/12) + (xx/6)*i, yy + (xx/6)*1, 
             yy + (xx/12) + (xx/6)*i, yy + (xx/6)*5, col= 'gold', lwd = 2)
    segments(yy + (xx/6)*1, yy + (xx/12) + (xx/6)*i, 
             yy + (xx/6)*5, yy + (xx/12) + (xx/6)*i, col= 'gold', lwd = 2)
}

image(data_mean_list_human$Mod, col = gray((0:96) / 96), axes = FALSE,
      main = "True activation map", xlim = xlimit,
      ylim = xlimit)
for (i in 1:3) {
    xline(yy + (xx/4)*i, col = "green")
    yline(yy + (xx/4)*i, col = "green")
}

for (i in 1:5) {
    segments(yy + (xx/4) + (xx/12)*i, yy + (xx/4)*1, 
             yy + (xx/4) + (xx/12)*i, yy + (xx/4)*3, col= 'gold', lwd = 2)
    segments(yy + (xx/4)*1, yy + (xx/4) + (xx/12)*i,
             yy + (xx/4)*3, yy + (xx/4) + (xx/12)*i, col= 'gold', lwd = 2)
}


# =======
# create region label
########################
regionlabel_human <- rep(0, N_human * N_human)
for (k in c(16)) {
    for (i in 1:(N_human / k)) {
        for(j in 1:(N_human / k)) {
            for(m in 1:k) {
                regionlabel_human[c((j - 1) * N_human * k + N_human * (m - 1) + 1:k) + (i - 1) * k] = 
                    i + (j - 1) * (N_human / k)
            }
        }
    }
}

########################
# Create regionlable matrix
###########################
region_mat_human <- matrix(regionlabel_human, ncol = N_human)
image.plot(region_mat_human, axes = FALSE)
G_human <- max(region_mat_human)  # Number of groups

##########################
# manipulate data so it can be used in MCMC
dim(Y_human)
Y_human_p <- aperm(Y_human, c(2, 3, 1))
dim(Y_human_p)
Y_humanvec <- matrix(Y_human_p, N_human * N_human, 980)
makecplx <- function(YY) {
    n = dim(YY)[3]
    YY[, , 1:(n / 2)] + 1i * YY[, , ((n / 2) + 1):n]
}
YY_human_cplx <- makecplx(Y_human_p)
YY_humanvec_cplx <- matrix(YY_human_cplx , nrow = N_human * N_human, ncol = 490)

################################################################################
# Cut the margin
################################################################################
cut.row.idx <- c(1:20, 77:96)
cut.col.idx <- c(1:20, 77:96)

# dimension 56 x 56
image(data_mean_list_human$Mod[-cut.row.idx, -cut.col.idx])
Mod_cut <- data_mean_list_human$Mod[-cut.row.idx, -cut.col.idx]
# str(data_mean_list_human$Mod[-cut.row.idx, -cut.col.idx])

xx <- 1.016
yy <- -0.007
xlimit = c(yy, xx)
image(data_mean_list_human$Mod[-cut.row.idx, -cut.col.idx], 
      col = gray((0:96) / 96), axes = FALSE,
      main = "True activation map", xlim = xlimit,
      ylim = xlimit)
# gridline(k = 8, 56, lwd = 2)
# 
for (i in 1:7) {
    xline(yy + (xx/8)*i, col = "green", lwd = 2)
    yline(yy + (xx/8)*i, col = "green", lwd = 2)
}


########################
N_human_cut <- 56
Nvoxels_human_cut <- N_human_cut * N_human_cut
regionlabel_human_cut <- rep(0, N_human_cut * N_human_cut)
for (k in c(7)) {
    for (i in 1:(N_human_cut / k)) {
        for(j in 1:(N_human_cut / k)) {
            for(m in 1:k) {
                regionlabel_human_cut[c((j - 1) * N_human_cut * k + 
                                            N_human_cut * (m - 1) + 1:k) + (i - 1) * k] = 
                    i + (j - 1) * (N_human_cut / k)
            }
        }
    }
}

region_mat_human_cut <- matrix(regionlabel_human_cut, ncol = N_human_cut)
image.plot(region_mat_human_cut, axes = FALSE)
G_human_cut <- max(region_mat_human_cut)  # Number of groups


# cut the margin for the data
##########################
# manipulate data so it can be used in Bezener MCMC
Y_human_cut <- Y_human_p[-cut.row.idx, -cut.col.idx, ]
dim(Y_human_cut)
Y_humanvec_cut <- matrix(Y_human_cut, dim(Y_human_cut)[1] ^ 2, 980)
dim(Y_humanvec_cut)
makecplx <- function(YY) {
    n = dim(YY)[3]
    YY[, , 1:(n / 2)] + 1i * YY[, , ((n / 2) + 1):n]
}
YY_human_cplx_cut <- makecplx(Y_human_cut)
YY_humanvec_cplx_cut <- matrix(YY_human_cplx_cut, nrow = N_human_cut * N_human_cut, 
                               ncol = 490)


####################################################################
# KC-AR
####################################################################
EB_rho_vec_cut <- matrix(0, N_human_cut * N_human_cut)
for (i in 1:N_human_cut * N_human_cut){
    rho_re <- suppressWarnings(ar.mle(Y_humanvec_cut[i, 1:490], 
                                      aic = FALSE,
                                      order.max = 1, 
                                      demean = FALSE)$ar)
    rho_im <- suppressWarnings(ar.mle(Y_humanvec_cut[i, (491):(980)], 
                                      aic = FALSE,
                                      order.max = 1, 
                                      demean = FALSE)$ar)
    EB_rho_vec_cut[i] <- (rho_re + rho_im) / 2
}

listname <- c(paste("w.", 1:G_human_cut, sep = ""))
keep.lst <- sapply(rep(0, G_human_cut), as.list)
tune.lst <- sapply(rep(1, G_human_cut), as.list)
names(keep.lst) <- listname
names(tune.lst) <- listname
name.par <- c(paste0("gamma", 1:(p * Nvoxels_human_cut)), 
              paste0("w", 1:(p * G_human_cut)), "phi")
keep.psi <- list(psi = 0)
tune.psi <- list(psi = 1)
start.lst <- list(ga = rep(1, Nvoxels_human_cut), 
                  w = t(rmvnorm(1, rep(0, G_human_cut), .5 * diag(G_human_cut))),
                  phi = 10)

n.mcmc <- 61000
burn <- 1000
thin <- 10

kcs_human_cut_ar <- mcmc_kcs_ar_block_cpp_full(Y_humanvec_cut, X_cplx_human, 
                                               EB_rho_vec_cut, 
                                               start.lst = start.lst, 
                                               nmcmc = n.mcmc,
                                               burn = burn, thin = thin, 
                                               name.par = name.par,
                                               adapt = TRUE, tunelst = tune.lst,
                                               tunepsi = tune.psi, 
                                               keeplst = keep.lst,
                                               keeppsi = keep.psi, tunelen = 50,
                                               regionlabel = regionlabel_human_cut, 
                                               region_mat = region_mat_human_cut, 
                                               nu = 2,
                                               atau = 1, btau = 1,
                                               aphi = 1, bphi = 1,
                                               agvec = rep(0, G_human_cut),
                                               target.accept.rate = 0.35, 
                                               iscplx = TRUE,
                                               kerneltype = "bezier")

# save(kcs_human_cut_ar, file = "kcs_human_cut_ar.RData")
# load("kcs_human_cut_ar.RData", verbose = TRUE)
sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
plot_parameter(MCMCsample = kcs_human_cut_ar, draws = NULL,
               parameter.idx = 3137:3139,
               sampleidx = sampleidx, name.par = name.par,
               main.title = "")

# # # Activation and Posterior Probability maps
# # # ===========================================================================
posprob_kcs_human_cut_ar <- compute_posprob_pc(kcs_human_cut_ar,
                                               N.voxels = Nvoxels_human_cut)
image_act(0.8722, posprob_kcs_human_cut_ar, kcs_human_cut_ar, is.cplx = TRUE)
image_prob(posprob_kcs_human_cut_ar, kcs_human_cut, is.cplx = TRUE)

brain_grid_plot(Mod_cut, posprob_kcs_human_cut_ar, N = 56, G = 64,
                mdltype = "CV-KC-AR", thre = 0.001, cex = 1.1, 
                is.color.bar = FALSE)

# # #############################
# # # posterior spatial effects
# # #############################
# # ===========
postmean_kcs_cut_ar <- apply(kcs_human_cut_ar$draws, 2, mean)

region_mat_human_cut <- assign_region_label(7, 56, TRUE)
centroid_cut <- compute_centroid(region_mat_human_cut)
Coord_cut <- cbind(rep(1:N_human_cut, N_human_cut),
                   rep(1:N_human_cut, each = N_human_cut))
dist.mat_cut <- dist_mat(coord = Coord_cut, grid = centroid_cut)
# ---
K_kcs_cut_ar <- bezier2(dist.mat_cut, nu = 2, postmean_kcs_cut_ar[3201])
W_kcs_cut_ar <- postmean_kcs_cut_ar[3137:3200]

S_kcs_cut_ar <- K_kcs_cut_ar %*% W_kcs_cut_ar
par(mfrow = c(1, 2))
# par(mar = c(0.1, 1, 2, 4))
par(mar = c(0.1, 0.5, 2, 0.5))  # 1200 x 560
image.plot(matrix(S_kcs_cut_ar, N_human_cut, N_human_cut), axes = FALSE,
           main = "CV-KC-AR voxel spatial effect S = Kw",
           zlim = c(-40, 10),  cex.main = 1.1, col = tim.colors()[1:50])
gridline(7, 56, col = "white", lwd = 2)
Mod_cut <- data_mean_list_human$Mod[-cut.row.idx, -cut.col.idx]
image(Mod_cut, col = gray.colors(96, alpha = 0.4), axes = FALSE, add = TRUE)
image.plot(matrix(1 / (1 + exp(-(S_kcs_cut_ar))), N_human_cut, N_human_cut),
           axes = FALSE, zlim = c(0, 1),
           main = "CV-KC-AR: 1 / (1 + exp(-S))", cex.main = 1.1)
gridline(7, 56, col = "white", lwd = 2)
image(Mod_cut, col = gray.colors(96, alpha = 0.4), axes = FALSE, add = TRUE)


image.plot(matrix(S_kcs_cut_ar, N_human_cut, N_human_cut), axes = FALSE,
           main = "KC-AR voxel spatial effect S = Kw",
           zlim = c(-50, 25),  cex.main = 1)
gridline(7, 56, col = "purple", lwd = 1)
image.plot(matrix(1 / (1 + exp(-(S_kcs_cut_ar))), N_human_cut, N_human_cut),
           axes = FALSE, zlim = c(0, 1),
           main = "KC-AR: 1 / (1 + exp(-S)) 50%", cex.main = 1)
gridline(7, 56, col = "purple", lwd = 1)



image.plot(matrix(1 / (1 + exp(-(S_kcs_cut_ar))), N_human_cut, N_human_cut),
           axes = FALSE, zlim = c(0, 1),
           main = "", cex.main = 1)
gridline(7, 56, col = "yellow", lwd = 1)
image(Mod_cut, col = gray.colors(96, alpha = 0.4), axes = FALSE, add = TRUE)

























