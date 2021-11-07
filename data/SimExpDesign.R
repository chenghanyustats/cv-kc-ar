################################################################################
# CV-fMRI Simulation Study: Experimental design                                #
# Circular Case                                                                #
# ./data/SimExpDesign.R                                                             #
# Cheng-Han Yu, Marquette University                                           #
################################################################################


# load required packages and functions 
######################################
library(neuRosim)
library(pryr)
library(Matrix) 
library(mvtnorm) #rmvnorm() generate a ROW vector
library(pscl)
library(matrixcalc) # is.positive.definite(x, tol=1e-8)
library(clusterGeneration)
library(cmvnorm)
library(fields)
library(MCMCpack)
library(scales)

#############################################################################
# Generate activations using neuRosim Package
############################################# 
n <- 50
N <- 20 # Change dim to N = 24 for spatial data; N = 20 originally
totaltime <- n
onsets <- seq(1, n, by = 50)
dur <- 15
sti <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, 
                 accuracy = 1)
canonical <- specifydesign(totaltime = totaltime, onsets = list(onsets),
                          durations = list(dur), effectsize = 1, TR = 1, 
                          conv = "double-gamma")
# plot(s, type = "l")


# define activation regions
coord1 <- matrix(c(3:4, 3:8, 3:8, 3:7, 3:7, 5:6, 5:6, 
                   rep(2, 2), rep(3, 6), rep(4, 6), rep(5, 5), rep(6, 5), rep(7, 2), 
                   rep(8, 2)), 
                 ncol = 2, byrow = FALSE)
coord2 <- matrix(c(rep(c(17, 18, 19), 3), rep(c(2, 3, 4), each = 3)), 
                 ncol = 2, byrow = FALSE)
coord3 <- matrix(c(12, 11:13, 10:14, 9:15, 10:14, 11:13, 12,
                   13, rep(14, 3), rep(15, 5), rep(16, 7), rep(17, 5), rep(18, 3),
                   rep(19, 1)), ncol = 2, byrow = FALSE)
reg1 <- specifyregion(dim = c(N, N), coord = coord1, 
                      form = "manual")
reg2 <- specifyregion(dim = c(N, N), coord = coord2, 
                      form = "manual")
reg3 <- specifyregion(dim = c(N, N), coord = coord3, 
                      form = "manual")
final <- reg1 + reg2 + reg3

################################################################################
#-------------------------------------------------------------------------------
# Figure: 
# Left: BOLD signal obtained from convolving the stimulus indicator signal 
# with the canonical hemodynamic function. 
# Right: activation regions and strengths.
#-------------------------------------------------------------------------------
# png("/figures/bold_act.png", height = 300, width = 600)
par(mfrow = c(1, 2))
par(mar = c(3, 2, 1.5, .5))
plot(canonical, type = 'l', main = "BOLD signals with canonical HRF", 
     xlab = "", lwd = 3, col = "brown", cex.main = 1, cex.axis = 0.9,
     cex.lab = 0.5, ylab = "")
title(xlab = "Time", cex.lab = 0.9, line = 2)
par(mar = c(.5, .5, 1, .5))
image(final, axes = FALSE, main = "Activation map: G = 25 size 4 x 4", 
      col = c("navy", "red"))
gridline(4, N, col = "yellow")
loc.x <- 0.1 * seq(1, 9, by = 2)
loc.y <- loc.x
for (i in 1:5) {
    for (j in 1:5) {
        text(loc.x[i], loc.y[j], i+5*(j-1), col = "yellow", cex = 1.5)
    } 
}
# dev.off()

# par(mfrow = c(1, 1))
# par(mar = c(4, 4, 3, 1))
# layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = FALSE))
# plot(s, type = 'l', xlab = "time", ylab = "stimulus", lwd = 2,
#      main = "Experimental Design", yaxt = "n")
# axis(2, at = 0:1, labels = 0:1)
# plot(canonical, type = 'l', main = "BOLD signal w/ canonical HRF", 
#      xlab = "time", lwd = 2)
# par(mar = c(1,1,1,2))

# png("./figures/grid_act.png", height = 400, width = 600)
par(mfrow = c(2, 3))
par(mar = c(.5, .5, 1, .5))
image(final, axes = FALSE, main = "Activation map", col = c("navy", "red"))
image(final, axes = FALSE, main = "G = 400 size 1 x 1", xlim = c(-0.025, 1.025),
           ylim = c(-0.025, 1.025), col = tim.colors()[1:55])

gridline(1, N, col = "yellow")
# image(final, axes = TRUE, main = "Activation map")
# for (i in 1:19) {
#     xline(-0.025 + (1.05/20)*i, col = "white")
#     yline(-0.025 + (1.05/20)*i, col = "white")
# }
image(final, axes = FALSE, main = "G = 100 size 2 x 2", 
      xlim = c(-0.025, 1.025),
           ylim = c(-0.025, 1.025), col = tim.colors()[1:55])
gridline(2, N, col = "yellow")
# image(final, axes = TRUE, main = "Activation map")
# for (i in 1:9) {
#     xline(-0.025 + (1.05/10)*i, col = "white")
#     yline(-0.025 + (1.05/10)*i, col = "white")
# }
image(final, axes = FALSE, main = "G = 25 size 4 x 4", 
           xlim = c(-0.025, 1.025),
           ylim = c(-0.025, 1.025),col = tim.colors()[1:55])
gridline(4, N, col = "yellow")
# image(final, axes = TRUE, main = "Activation map")
# for (i in 1:4) {
#     xline(-0.025 + (1.05/5)*i, col = "white")
#     yline(-0.025 + (1.05/5)*i, col = "white")
# }
image(final, axes = FALSE, main = "G = 16 size 5 x 5", 
      xlim = c(-0.025, 1.025),
      ylim = c(-0.025, 1.025),col = tim.colors()[1:55])
gridline(5, N, col = "yellow")

image(final, axes = FALSE, main = "G = 4 size 10 x 10", 
      xlim = c(-0.025, 1.025),
      ylim = c(-0.025, 1.025),col = tim.colors()[1:55])
gridline(10, N, col = "yellow")
# dev.off()

# irregular grid
################
# par(mfrow = c(2, 3))
# par(mar = c(.5, .5, 1, .5))
# image(final, axes = FALSE, main = "Connectivity", xlim = c(-0.025, 1.025),
#            ylim = c(-0.025, 1.025), col = tim.colors()[1:55])
# # image(final, axes = TRUE, main = "Activation map")
# # segments(-0.025 + 0.105, -0.025+ 0.105/2, -0.025 + 0.105, -0.025 + 0.105*4, col = "white")
# # segments(-0.025 + 0.105, -0.025, -0.025 + 0.105, -0.025 + 0.105*4, col = "white")
# segments(-0.025, -0.025 + 0.105*4, -0.025 + 0.105*4, -0.025 + 0.105*4, col = "white")
# # segments(-0.025 + 0.105*4, -0.025 + 0.105*4, -0.025 + 0.105*4, -0.025, col = "white")
# segments(-0.025 + 0.105*4, -0.025 + 0.105*10, -0.025 + 0.105*4, -0.025, col = "white")
# # segments(-0.025 + 0.105, -0.025 + 0.105/2, -0.025 + 0.105*4, -0.025+ 0.105/2, col = "white")
# 
# segments(-0.025 + 0.105*8, -0.025 + 0.105*10, -0.025 + 0.105*8, -0.025, col = "white")
# # segments(-0.025 + 0.105*8, -0.025 + 0.105/2, -0.025 + 0.105*9.5, -0.025+ 0.105/2, col = "white")
# segments(-0.025 + 0.105*8, -0.025+ 0.105*2, -0.025 + 0.105*10, -0.025+ 0.105*2, col = "white")
# # segments(-0.025 + 0.105*9.5, -0.025 + 0.105*2, -0.025 + 0.105*9.5, -0.025+ 0.105/2, col = "white")
# 
# # segments(-0.025 + 0.105*4.5, -0.025+ 0.105*6.5, -0.025 + 0.105*4.5, -0.025 + 0.105*9.5, col = "white")
# # segments(-0.025 + 0.105*4.5, -0.025 + 0.105*9.5, -0.025 + 0.105*7.5, -0.025 + 0.105*9.5, col = "white")
# # segments(-0.025 + 0.105*7.5, -0.025 + 0.105*9.5, -0.025 + 0.105*7.5, -0.025+ 0.105*6.5, col = "white")
# segments(-0.025 + 0.105*8, -0.025+ 0.105*6, -0.025 + 0.105*4, -0.025+ 0.105*6, col = "white")
# 
# text(0.08, 0.08, 1, col = "yellow", cex = 5, font = 2)
# text(0.6, 0.2, 2, col = "yellow", cex = 5, font = 2)
# text(0.92, 0.08, 3, col = "yellow", cex = 5, font = 2)
# text(0.2, 0.8, 4, col = "yellow", cex = 5, font = 2)
# text(0.6, 0.8, 5, col = "yellow", cex = 5, font = 2)
# text(0.92, 0.6, 6, col = "yellow", cex = 5, font = 2)
# 
# segments(0.08, 0.08, 0.92, 0.08, col = alpha("green", 0.8), 
#          lwd = 5)
# segments(0.08, 0.08, 0.6, 0.8, col = alpha("green", 0.8), 
#          lwd = 5)
# segments(0.92, 0.08, 0.6, 0.8, col = alpha("green", 0.8), 
#          lwd = 5)
# 

# segments(0.38, 0.12, 0.88, 0.64, col = alpha("white", 0.8), 
#          lwd = abs(Prec_mat_sim[2, 12]))
# segments(0.64, 0.12, 0.88, 0.88, col = alpha("white", 0.8),
#          lwd = abs(Prec_mat_sim[3, 16]))
# segments(0.12, 0.38, 0.88, 0.64, col = alpha("white", 0.8),
#          lwd = abs(Prec_mat_sim[5, 12]))
# segments(0.38, 0.38, 0.88, 0.88, col = alpha("white", 0.8),
#          lwd = abs(Prec_mat_sim[6, 16]))
# segments(0.12, 0.64, 0.88, 0.64, col = alpha("white", 0.8),
#          lwd = abs(Prec_mat_sim[9, 12]))


# 
# 
# 
# loc.x <- c(0.08, 0.29, 0.5, 0.71, 0.92) - .06
# loc.y <- loc.x
# for (i in 1:5) {
#     for (j in 1:5) {
#         text(loc.x[i], loc.y[j], i+5*(j-1), col = "yellow", cex = 1.5,
#              font = 2)
#     } 
# }
# 
# 
# 
# 
# reg_label <- rep(0, 400)
# reg1 <- as.vector(sapply(20 * (0:7), function(x) x + (1:8)))
# reg2 <- as.vector(sapply(20 * (0:11), function(x) x + (9:16)))
# reg3 <- as.vector(sapply(20 * (0:3), function(x) x + (17:20)))
# reg4 <- as.vector(sapply(20 * (8:19), function(x) x + (1:8)))
# reg5 <- as.vector(sapply(20 * (12:19), function(x) x + (9:16)))
# reg6 <- as.vector(sapply(20 * (4:19), function(x) x + (17:20)))
# 
# reg_label[reg1] <- 1
# reg_label[reg2] <- 2
# reg_label[reg3] <- 3
# reg_label[reg4] <- 4
# reg_label[reg5] <- 5
# reg_label[reg6] <- 6
# 
# image.plot(matrix(reg_label, 20, 20))

################################################################################
# multiresolution
par(mar = c(0.1, .1, 2, .1))
lwd = 3
image(final, axes = FALSE, main = "Multi-resolution grids", 
      xlim = c(-0.025, 1.025),
      ylim = c(-0.025, 1.025),col = tim.colors()[1:55])
loc.x <- c(0.08, 0.29, 0.5, 0.71, 0.92)
loc.y <- loc.x
for (i in 1:5) {
    for (j in 1:5) {
        points(loc.x[i], loc.y[j], pch = 16, col = "white", cex = 1.2)
    } 
}


loc.x <- seq(0.025, 0.97, by = 0.105)
loc.y <- loc.x
for (i in 1:4) {
    for (j in 1:4) {
        points(loc.x[i], loc.y[j], pch = 16, col = "green", cex = 1.5)
    } 
}

for (i in 9:10) {
    for (j in 1:2) {
        points(loc.x[i], loc.y[j], pch = 16, col = "green", cex = 1.5)
    } 
}

for (i in 5:8) {
    for (j in 7:10) {
        points(loc.x[i], loc.y[j], pch = 16, col = "green", cex = 1.5)
    } 
}

gridline(4, N, col = "yellow", lwd = lwd)

segments(-0.025 + 0.105, -0.025, -0.025 + 0.105, -0.025 + 0.105*4, col = "white",lwd = lwd)
segments(-0.025 + 0.105*3, -0.025 + 0.105*4, -0.025 + 0.105*3, -0.025, col = "white",lwd = lwd)
segments(-0.025, -0.025 + 0.105, -0.025 + 0.105*4, -0.025 + 0.105, col = "white",lwd = lwd)
segments(-0.025, -0.025 + 0.105*3, -0.025 + 0.105*4, -0.025 + 0.105*3, col = "white",lwd = lwd)

segments(-0.025 + 0.105*9, -0.025 + 0.105*2, -0.025 + 0.105*9, -0.025, col = "white",lwd = lwd)
segments(-0.025 + 0.105*8, -0.025 + 0.105, -0.025 + 0.105*10, -0.025 + 0.105, col = "white",lwd = lwd)

segments(-0.025 + 0.105*5, -0.025 + 0.105*6, -0.025 + 0.105*5, -0.025 + 0.105*10, col = "white",lwd = lwd)
segments(-0.025 + 0.105*7, -0.025 + 0.105*6, -0.025 + 0.105*7, -0.025 + 0.105*10, col = "white",lwd = lwd)
segments(-0.025 + 0.105*4, -0.025 + 0.105*7, -0.025 + 0.105*8, -0.025 + 0.105*7, col = "white",lwd = lwd)
segments(-0.025 + 0.105*4, -0.025 + 0.105*9, -0.025 + 0.105*8, -0.025 + 0.105*9, col = "white",lwd = lwd)


par(mfrow = c(1, 3))
par(mar = c(3, 2, 1.5, .5))
# layout(matrix(c(1, 2, 3, 3, 4, 4), 2, 3, byrow = FALSE))
# plot(sti, type = 'l', xlab = "time", ylab = "stimulus", lwd = 2,
#      main = "Experimental Design", yaxt = "n")
plot(canonical, type = 'l', main = "Expected BOLD with canonical HRF", 
     xlab = "", lwd = 3, col = "brown", cex.main = 1, cex.axis = 0.9,
     cex.lab = 0.5, ylab = "")
title(xlab = "Time", cex.lab = 0.9, line = 2)
par(mar = c(.5, .5, 1, .5))
image(final, axes = FALSE, main = "G = 400 size 1 x 1", xlim = c(-0.025, 1.025),
      ylim = c(-0.025, 1.025), col = tim.colors()[1:55])
gridline(1, N, col = "white", lwd = 0.8)
image(final, axes = FALSE, main = "G = 25 size 4 x 4", 
      xlim = c(-0.025, 1.025),
      ylim = c(-0.025, 1.025),col = tim.colors()[1:55])
gridline(1, N, col = "white", lwd = 0.3)
gridline(4, N, col = "white", lwd = 0.3)
# loc.x <- 0.1 * seq(1, 9, by = 2)
loc.x <- c(0.08, 0.29, 0.5, 0.71, 0.92)
loc.y <- loc.x
for (i in 1:5) {
    for (j in 1:5) {
        points(loc.x[i], loc.y[j], pch = 16, col = "yellow", cex = 2)
    } 
}
# loc.x <- c(0.08, 0.29, 0.5, 0.71, 0.92) - .06
# loc.y <- loc.x
# for (i in 1:5) {
#     for (j in 1:5) {
#         text(loc.x[i], loc.y[j], i+5*(j-1), col = "yellow", cex = 1.5,
#              font = 2)
#     } 
# }



par(mfrow = c(1, 2))
par(mar = c(3, 2, 1.5, .5))
# layout(matrix(c(1, 2, 3, 3, 4, 4), 2, 3, byrow = FALSE))
# plot(sti, type = 'l', xlab = "time", ylab = "stimulus", lwd = 2,
#      main = "Experimental Design", yaxt = "n")
# plot(canonical, type = 'l', main = "Expected BOLD with canonical HRF", 
#      xlab = "", lwd = 3, col = "brown", cex.main = 1, cex.axis = 0.9,
#      cex.lab = 0.5, ylab = "")
# title(xlab = "Time", cex.lab = 0.9, line = 2)
par(mar = c(.5, .5, 1, .5))
image(final, axes = FALSE, main = "G = 400 size 1 x 1", xlim = c(-0.025, 1.025),
      ylim = c(-0.025, 1.025), col = tim.colors()[1:55])
gridline(1, N, col = "white", lwd = 0.8)
image(final, axes = FALSE, main = "G = 25 size 4 x 4", 
      xlim = c(-0.025, 1.025),
      ylim = c(-0.025, 1.025),col = tim.colors()[1:55])
gridline(1, N, col = "white", lwd = 0.3)
gridline(4, N, col = "white", lwd = 4)
# loc.x <- 0.1 * seq(1, 9, by = 2)
loc.x <- c(0.08, 0.29, 0.5, 0.71, 0.92)
loc.y <- loc.x
for (i in 1:5) {
    for (j in 1:5) {
        points(loc.x[i], loc.y[j], pch = 16, col = "yellow", cex = 1.2)
    } 
}
loc.x <- c(0.08, 0.29, 0.5, 0.71, 0.92) - .06
loc.y <- loc.x
for (i in 1:5) {
    for (j in 1:5) {
        text(loc.x[i], loc.y[j], i+5*(j-1), col = "yellow", cex = 1.5,
             font = 2)
    } 
}

par(mfrow = c(1, 2))
par(mar = c(3, 2, 1.5, .5))
# layout(matrix(c(1, 2, 3, 3, 4, 4), 2, 3, byrow = FALSE))
# plot(sti, type = 'l', xlab = "time", ylab = "stimulus", lwd = 2,
#      main = "Experimental Design", yaxt = "n")
plot(canonical, type = 'l', main = "Expected BOLD with canonical HRF", 
     xlab = "", lwd = 3, col = "brown", cex.main = .9, cex.axis = 0.9,
     cex.lab = 0.5, ylab = "")
title(xlab = "Time", cex.lab = 0.9, line = 2)
# par(mar = c(.5, .5, 1, .5))
# image(final, axes = FALSE, main = "G = 400 size 1 x 1", xlim = c(-0.025, 1.025),
#       ylim = c(-0.025, 1.025), col = tim.colors()[1:55])
# gridline(1, N, col = "white", lwd = 0.8)
par(mar = c(.5, .5, 1.2, .5))
image(final, axes = FALSE, main = "20 x 20 true activation", cex.main = 1.2,
      xlim = c(-0.025, 1.025),
      ylim = c(-0.025, 1.025),col = tim.colors()[1:55])
gridline(1, N, col = "white", lwd = 0.6)
gridline(4, N, col = "white", lwd = 0.6)
# loc.x <- 0.1 * seq(1, 9, by = 2)
loc.x <- c(0.08, 0.29, 0.5, 0.71, 0.92)
loc.y <- loc.x
for (i in 1:5) {
    for (j in 1:5) {
        points(loc.x[i], loc.y[j], pch = 16, col = "yellow", cex = 1.2)
    } 
}