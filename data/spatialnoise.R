################################################################################
# # Spatial Noises from Neurosim package spatialnoise()                        #
# Research/Bezener/Analysis/Data/spatialnoise.R                                #
# Cheng-Han Yu, UC Santa Cruz                                                  #
################################################################################
# library(neuRosim)
# library(fields)
# 
# d <- c(20, 20)
# n.corr <- spatialnoise(dim = d, sigma = 15, nscan = 100,
#                        method = "corr", rho = 0.9)
# par(mfrow = c(5, 2))
# for (i in 1:10) {
#     image.plot(n.corr[, , 10*i], main = paste("corr time =", 10*i, "gaussian"))
# }
# for (i in 1:10) {
#     image.plot(n.corr[, , i], main = paste("corr time =", i))
# }
# 
# 
# n.corr <- spatialnoise(dim = d, sigma = 15, nscan = 100,
#                        method = "corr", type = "rician", rho = 0.7)
# par(mfrow = c(5, 2))
# for (i in 1:10) {
#     image.plot(n.corr[, , 10*i], main = paste("corr time =", 10*i, "rician"))
# }
# 
# 
# n.gaus <- spatialnoise(dim = d, sigma = 15, nscan = 100,
#                        method = "gaussRF", FWHM = 4)
# # the larger FWHM is, the smoother (more correlated) 
# for (i in 1:10) {
#     image.plot(n.gaus[, , 10*i], main = paste("gaussRF time =", 10*i))
# }
# 
# 
# n.gamma <- spatialnoise(dim = d, sigma = 15, nscan = 100,
#                         method = "gammaRF", FWHM = 4, 
#                         gamma.shape = 3, gamma.rate = 2)
# for (i in 1:10) {
#     image.plot(n.gamma[, , 10*i], main = paste("gammaRF time =", 10*i))
# }
################################################################################
# Observed values
#################

ones <- rep(1, n) # intercept
x_no = cbind(ones, canonical)
X_no = as.matrix(bdiag(x_no, x_no))
xstar_no = x_no[, 1]
II_no <- (diag(n) - xstar_no %*% solve(t(xstar_no) %*% xstar_no) %*% t(xstar_no))
x_re_no <-  II_no %*% x_no[, 2]
x_im_no <- x_re_no
x_mod_no <- x_re_no
x_cp_no = cbind(c(x_re_no, rep(0, n)), c(rep(0, n), x_im_no))
Xstar_no <- X_no[, c(1, 3)]
IIc_no <- (diag(2*n) - Xstar_no %*% solve(t(Xstar_no)%*%Xstar_no) %*% t(Xstar_no))
X_cplx_no <- IIc_no %*% X_no[, c(2, 4)]

YYobs_s = list()
YY_s = list()
YYmod_s = list()
phi = 0.9
noise_im_s <- spatialnoise(dim = c(N, N), sigma = sqrt(sig2), nscan = n, 
                           rho = c(phi), method = "corr")
noise_re_s <- spatialnoise(dim = c(N, N), sigma = sqrt(sig2), nscan = n, 
                           rho = c(phi), method = "corr")

## SIMULATE CIRCULAR TIME SERIES DATA...
fmritimeseries_re_no_s = array(NA, c(N, N, n))
fmritimeseries_im_no_s = array(NA, c(N, N, n))
noise_re_no_s = array(NA, c(N, N, n))
noise_im_no_s = array(NA, c(N, N, n))
Y_no_s <- array(NA, c(N, N, 2 * n))
Yobs_no_s <- array(NA, c(N, N, n))
Y_mod_no_s <- array(NA, c(N, N, n))

for (s in 1:(length(SNR) * length(CNR))) {
    for (j in 1:N) {
        for (i in 1:N) {
            if (final[i, j] != 0) {
                # Voxel is activated and so:
                fmritimeseries_re_no_s[i, j, ] = Beta_re_no_s[[s]][1] + 
                    (Beta_re_no_s[[s]][2] * final[i, j]) * canonical
                
                fmritimeseries_im_no_s[i, j, ] = Beta_im_no_s[[s]][1] + 
                    (Beta_im_no[[s]][2] * final[i, j]) * canonical
            } else {
                fmritimeseries_re_no_s[i, j, ] = Beta_re_no_s[[s]][1]
                
                fmritimeseries_im_no_s[i, j, ] = Beta_im_no_s[[s]][1]
            }
            Yobs_no_s[i, j, ] = (fmritimeseries_re_no_s[i, j, ] + 
                                     noise_re_s[i, j, ]) + 
                1i * (fmritimeseries_im_no_s[i, j, ] + noise_im_s[i, j, ])
            Y_mod_no_s[i, j, ] <- Mod(Yobs_no_s[i, j, ] )
            Y_no_s[i, j, ] <- c(fmritimeseries_re_no_s[i, j, ] + noise_re_s[i, j, ], 
                              fmritimeseries_im_no_s[i, j, ] + noise_im_s[i, j, ])
            Yobs_no_s[i, j, ] <- II_no %*% Yobs_no_s[i, j, ]
            Y_no_s[i, j, ] <- IIc_no %*% Y_no_s[i, j, ]
            Y_mod_no_s[i, j, ] <- II_no %*% Y_mod_no_s[i, j, ]
        }
    }
    YYobs_s[[s]] = Yobs_no_s
    YY_s[[s]] = Y_no_s
    YYmod_s[[s]] = Y_mod_no_s
}

Yvec_s <- lapply(YY_s, matrix, nrow = N * N, ncol = 2 * n)
Yobsvec_s <- lapply(YYobs_s, matrix, nrow = N * N, ncol = n)
Ymodvec_s <- lapply(YYmod_s, matrix, nrow = N * N, ncol = n)
no_YY_s = length(YY_s)

# par(mfrow = c(1, 3))
# par(mar = c(4, 4, 2, 1))
# image(YYmod[[1]][, , 21], axes = FALSE, col = tim.colors())
# image(YYmod[[1]][, , 22], axes = FALSE, col = tim.colors())
# image(YYmod[[1]][, , 23], axes = FALSE, col = tim.colors())
# plot(YYmod[[1]][1, 1, ], type = "l", main = "Y[1, 1] nonact")
# lines(canonical, col = "red")
# plot(YYmod[[1]][5, 5, ], type = "l", main = "Y[5, 5] act")
# lines(canonical, col = "red")




