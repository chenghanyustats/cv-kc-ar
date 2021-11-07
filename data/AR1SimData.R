################################################################################
# Complex-valued fMRI Simulation Study: AR(1) Data Generating                  #
# "./data/AR1SimData.R"                                                        #
# Cheng-Han Yu, Marquette University                                           #
################################################################################

## SIMULATE CIRCULAR TIME SERIES DATA
########################
YYobs <- YY <- YYmod <- list()
fmritimeseries.re <- array(NA, c(N, N, n))
fmritimeseries.im <- array(NA, c(N, N, n))
Y <- array(NA, c(N, N, 2 * n))
Yobs <- array(NA, c(N, N, n))
Y.mod <- array(NA, c(N, N, n))

# AR(1) noise
# =========================================================
# change seeds.re[2] from 6422 to 6423
seeds.re <- c(8710, 6423, 5403, 1074, 2114, 5065, 8086, 4862, 2898, 8941, 7328, 
              5974, 8922, 7589, 1699, 7910, 7761, 8942, 4659, 8816)
seeds.re <- c(seeds.re,
              c(1068, 1132, 1152, 1199, 1406, 1731, 1824, 1837, 1881, 1888, 
                1894, 1904, 2124, 2249, 2846, 3080, 3285, 3671, 3693, 3768, 
                3799, 3823, 3976, 4003, 4072, 4158, 4206, 4272, 4342, 4361,
                4370, 4474, 4485, 4702, 4858, 4902, 4963, 5089, 5120, 5320,
                5525, 5612, 5616, 5686, 5796, 5869, 5987, 6106, 6114, 6236,
                6438, 6533, 6618, 6745, 6758, 7285, 7477, 7488, 7605, 7655,
                7834, 7874, 8009, 8012, 8061, 8180, 8215, 8261, 8561, 8696,
                9026, 9056, 9260, 9509, 9522, 9607, 9673, 9801, 9832, 9848))


seeds.im <- c(3064, 8184, 7789, 7701, 9378, 5511, 5167, 8086, 8556, 2305, 9538,
              8106, 1345, 7284, 5238, 3839, 9708, 6960, 9466, 1511)

seeds.im <- c(seeds.im, 
              c(1019, 1074, 1201, 1213, 1228, 1355, 1420, 1423, 1569, 1664,
                1778, 1863, 2193, 2335, 2374, 2540, 2545, 2612, 2672, 2823,
                2933, 3080, 3176, 3373, 3393, 3472, 3989, 4069, 4111, 4449,
                4529, 4591, 4647, 4808, 4827, 5132, 5163, 5297, 5311, 5454,
                5479, 5697, 5738, 5862, 6007, 6088, 6163, 6168, 6172, 6185,
                6505, 6647, 6704, 6715, 6807, 6930, 7007, 7280, 7552, 7623,
                7771, 7899, 8103, 8413, 8426, 8474, 8608, 8710, 8739, 8782,
                8787, 8794, 8803, 8897, 8921, 8944, 9219, 9460, 9748, 9855))

# rho <- 0.5
set.seed(seeds.re[k])
noise.temp.re <- temporalnoise(dim = c(N, N), sigma = sqrt(sig2), nscan = n, 
                               rho = c(rho))
set.seed(seeds.im[k])
noise.temp.im <- temporalnoise(dim = c(N, N), sigma = sqrt(sig2), nscan = n, 
                               rho = c(rho))

# Create data sets
# =========================================================
for (s in 1:(length(SNR) * length(CNR))) {
    for (j in 1:N) {
        for (i in 1:N) {
            if (final[i, j] != 0) {
                # Voxel is activated and so:
                fmritimeseries.re[i, j, ] <- Beta.re[[s]][1] + 
                    (Beta.re[[s]][2] * final[i, j]) * canonical
                fmritimeseries.im[i, j, ] <- Beta.im[[s]][1] + 
                    (Beta.im[[s]][2] * final[i, j]) * canonical
            } else {
                fmritimeseries.re[i, j, ] <- Beta.re[[s]][1]
                fmritimeseries.im[i, j, ] <- Beta.im[[s]][1]
            }
            Yobs[i, j, ] <- (fmritimeseries.re[i, j, ] + noise.temp.re[i, j, ]) + 
                1i * (fmritimeseries.im[i, j, ] + noise.temp.im[i, j, ])
            Y.mod[i, j, ] <- Mod(Yobs[i, j, ])
            Y[i, j, ] <- c(fmritimeseries.re[i, j, ] + noise.temp.re[i, j, ], 
                           fmritimeseries.im[i, j, ] + noise.temp.im[i, j, ])
            Y[i, j, ] <- IIc %*% Y[i, j, ]
            Y.mod[i, j, ] <- II %*% Y.mod[i, j, ]
        }
    }
    YY[[s]] <- Y
    YYobs[[s]] <- Yobs
    YYmod[[s]] <- Y.mod
}
# no.YY <- length(YY)
YYvec <- lapply(YY, matrix, nrow = N * N, ncol = 2 * n)
YYobsvec <- lapply(YYobs, matrix, nrow = N * N, ncol = n)
YYmodvec <- lapply(YYmod, matrix, nrow = N * N, ncol = n)



# =========================================================

par(mfrow = c(2, 1))
par(mar = c(2,4,4,4))
plot(noise.temp.re[12, 12, ], type = "l", 
     main = paste("noise.temp.re[12, 12, ] rho = ", rho))
plot(noise.temp.im[12, 12, ], type = "l", 
     main = paste("noise.temp.re[12, 12, ] rho = ", rho))
plot(noise.temp.re[3, 3, ], type = "l", 
     main = paste("noise.temp.re[3, 3, ] rho = ", rho))
plot(noise.temp.im[3, 3, ], type = "l", 
     main = paste("noise.temp.re[3, 3, ] rho = ", rho))
par(mfrow = c(2, 1))
plot(scale(YYmod[[1]][1, 1, ]), type = 'l', ylab = "y(1, 1)", xlab = "time",
     main = "Nonactivated time series AR1", lwd = 2)
lines(scale(canonical), type = 'l', main = "BOLD signal w/ canonical HRF", 
      xlab = "time", lwd = 2, col = "red")
# # activation region
plot(scale(YYmod[[1]][5, 5, ]), type = 'l', ylab = "y(5, 5)", xlab = "time", 
     main = "Activated time series AR1", lwd = 2)
lines(scale(canonical), type = 'l', main = "BOLD signal w/ canonical HRF", 
      xlab = "time", lwd = 2, col = "red")
