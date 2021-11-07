################################################################################
# All self-created functions for the CV-KC-AR paper                                #
# ./my_fcns.R                                                                         #
# Cheng-Han Yu, Marquette University                                           #
################################################################################
# Assign label numbers
assign_region_label <- function(k, N, is.mat = TRUE) {
    # Assign region labels
    #
    # Args:
    #   k: The length of a region
    #   N: The length of an image
    #   is.mat: Return a matrix or not? 
    # Returns:
    #   A matrix with region numbers if is.mat is true, a vector otherwise
    
    regionlabel.vec <- rep(0, N * N)
    for (i in 1:(N / k)) {
        for (j in 1:(N / k)) {
            for (m in 1:k) {
                regionlabel.vec[c((j - 1) * N * k + N * (m - 1) + 1:k) + (i - 1) * k] <- 
                    i + (j - 1) * (N / k)
            }
        }
    }
    if (is.mat == TRUE) {
        return(matrix(regionlabel.vec, ncol = N))
    } else {
        return(regionlabel.vec)
    }
}

# get the number of voxels
get_voxel_number <- function(region_mat) {
    G <- max(region_mat)
    voxel_number <- rep(0, G)
    for(i in 1:G) {
        voxel_number[i] <- length(which(region_mat == i))
    }
    return(voxel_number)
}


# compute_posprob <- function(model, sampleidx) {
#     # posprob <- rep(0, N*N)
#     len = length(sampleidx)
#     # not faster than for loop
#     # test.posprob = apply(ssvsconjsnr05$psi[1, sampleidx, 1, , ], c(2, 3), 
#     #                      function(x) sum(x)/length(x))
#     GammaData <- model$Gamma[sampleidx, , 1]
#     posprob <- apply(GammaData, 2, sum) / len
#     return(posprob)
# }

compute_posprob <- function(model) {
    # posprob <- rep(0, N*N)
    gammaidx <- grep("gamma", colnames(model$draws))
    len <- nrow(model$draws)
    # not faster than for loop
    # test.posprob = apply(ssvsconjsnr05$psi[1, sampleidx, 1, , ], c(2, 3), 
    #                      function(x) sum(x)/length(x))
    GammaData <- model$draws[, gammaidx]
    posprob <- colSums(GammaData) / len
    return(posprob)
}


compute_posprob_new <- function(model, sampleidx, N, skip_idx) {
    len = length(sampleidx)
    GammaData <- model$Gamma[sampleidx, , 1]
    GammaDatafull <- rep(0, N * N)
    GammaDatafull[-skip_idx] <- GammaData
    posprob <- apply(GammaDatafull, 2, sum) / len
    return(posprob)
}


# Function to plot an activation map 
image_act <- function(aa, posprob, model, main = deparse(substitute(model)), 
                      is.cplx, cex.main = 1) {
    # aa: threshold value vector
    title <- ifelse(is.cplx, "CV", "MO")
    posprob_mat <- matrix(posprob, ncol = sqrt(length(posprob)))
    for (i in aa) {
        image(posprob_mat > i,
              # main = main,
              # main = paste0(title, " Act", " thr = ", i, "\n", main),
              main = paste0(title, main, " thr = ", i),
              axes = FALSE,
              cex.main = cex.main,
              col = c("navy", "red1") )
    }
}

# Function to plot an posterior probability map 
image_prob <- function(posprob, model, main = deparse(substitute(model)), 
                       is.cplx, cex.main = 1) {
    title = ifelse(is.cplx, "CV", "MO")
    posprob_mat <- matrix(posprob, ncol = sqrt(length(posprob)))
    image.plot(posprob_mat, 
               # main = main,
               main = paste(title, main, sep = "-"), 
               # main = main,
               cex.main = cex.main,
               axes = FALSE, col = tim.colors()[1:55])
}

# Function to compute true positives and false positives
# compute true positive rate and false positive rate
tpfpDatafcn_mcc <- function(cutoff, prob, truedata) {
    tpfp <- function(cutoff, prob, truedata) {
        # truedata: a true/false matrix
        # compute true positive rate and false positive rate
        tbl <- table(prob > cutoff, truedata)
        marg.col <- margin.table(tbl, 2)
        marg.row <- margin.table(tbl, 1)
        if ("FALSE" %in% rownames(tbl) && "TRUE" %in% rownames(tbl)) {
            # true positive (recall)
            tp <- tbl[2, 2] / marg.col[2]
            # false positive
            fp <- tbl[2, 1] / marg.col[1]
            # precision
            precision <- tbl[2, 2] / marg.row[2]
            # accuracy
            acc <- (tbl[1, 1] + tbl[2, 2]) / sum(tbl)
            # F-1
            f1 <- 2 * tbl[2, 2] / (2 * tbl[2, 2] + tbl[2, 1] + tbl[1, 2])
            # Matthews correlation coefficient 
            mcc <- (tbl[1, 1] * tbl[2, 2] - tbl[1, 2] * tbl[2, 1]) / 
                (sqrt((tbl[2, 1] + tbl[2, 2]) * (tbl[1, 2] + tbl[2, 2]) *
                          (tbl[1, 1] + tbl[2, 1]) * (tbl[1, 1] + tbl[1, 2])))
        } else if (!("FALSE" %in% rownames(tbl))) {
            # true positive
            tp <- tbl[1, 2] / marg.col[2]
            # false positive
            fp <- tbl[1, 1] / marg.col[1]
            # precision
            precision <- tbl[1, 2] / marg.row[1]
            # accuracy
            acc = tbl[1, 2] / sum(tbl)
            # F-1
            f1 <- 2 * tbl[1, 2] / (2 * tbl[1, 2] + tbl[1, 1])
            # Matthews correlation coefficient 
            mcc <- 0
        } else if  (!("TRUE" %in% rownames(tbl))) {
            # true positive
            tp <- 0
            # false positive
            fp <- 0
            # precision
            precision <- 0
            # accuracy
            acc = tbl[1, 1] / sum(tbl)
            # F-1
            f1 <- 0
            # Matthews correlation coefficient 
            mcc <- 0
        }
        rate <- c(tp, 1 - tp, fp, 1 - fp, precision, 1 - precision, acc, f1, mcc)
        names(rate) <- c('TPR/sensitivity/recall', "FNR=1-TPR",
                         'FPR', 
                         'TNR=1-FPR(specificity)',
                         "Precision(PPV)", "FDR=1-PPV", "Accuracy(ACC)",
                         "F1", "MCC")
        return(rate)
    }
    tpfpData <- as.data.frame(
        do.call(rbind, lapply(cutoff, tpfp, prob, truedata)))
    tpfpData$threshold <- cutoff
    return(tpfpData)
}



# check_parameter <- function(MCMCsample, sampleidx, variable, Sno = 1) {
#     par(mfrow = c(1, 4))
#     par(mar = c(4,4,4,1))
#     len = length(sampleidx)
#     breaks = floor(len / 20)
#     col = "navy"
#     dataname <- deparse(substitute(MCMCsample))
#     
#     if (variable == "S") {
#         # Trace plot
#         plot(MCMCsample$S[sampleidx, Sno, 1], type = "l", ylab = paste("S", Sno), 
#              main = expression(paste("Trace of ", S[Sno])))
#         
#         # ACF
#         acf(MCMCsample$S[sampleidx, Sno, 1], 
#             main = expression(paste("ACF of ", S[Sno])))
#         
#         # Distribution
#         hist(MCMCsample$S[sampleidx, Sno, 1], main = dataname, freq = FALSE,
#              breaks = breaks, xlab = expression(S[Sno]), 
#              col = col, border = FALSE)
#         
#         # Ergodic mean
#         plot(cumsum(MCMCsample$S[sampleidx, Sno, 1]) / (1:(len)), type = "l",
#              xlab = "iteration", ylab = paste("S", Sno), lwd = 2,
#              main = expression(paste("Ergodic mean ", S[Sno])))
#     } else {
#         # Trace plot
#         plot(MCMCsample$R[sampleidx], type = "l", ylab = "r", 
#              main = expression(paste("Trace of ", r)))
#         
#         # ACF
#         acf(MCMCsample$R[sampleidx], main = expression(paste("ACF of ", r)))
#         
#         # Distribution
#         hist(MCMCsample$R[sampleidx], main = dataname, freq = FALSE,
#              breaks = breaks, xlab = expression(r), col = col, border = FALSE)
#         
#         # Ergodic mean
#         plot(cumsum(MCMCsample$R[sampleidx]) / (1:(len)), type = "l",
#              xlab = "iteration", ylab = expression(r), lwd = 2,
#              main = expression(paste("Ergodic mean ", r)))
#     }
#     
# }
plot_parameter <- function(MCMCsample, draws = NULL, parameter.idx = NULL, 
                           sampleidx, name.par, hist.title = "", isquan = TRUE) {
    if (is.null(draws)) draws <- MCMCsample$draws
    ifelse(is.null(parameter.idx), par(mfrow = c(length(parameter.name), 4)),
           par(mfrow = c(length(parameter.idx), 4)))
    par(mar = c(4, 4, 4, 1))
    len <- length(sampleidx)
    breaks <- floor(len / 50)
    col <- "navy"
    if (!is.null(MCMCsample)) dataname <- deparse(substitute(MCMCsample))
    
    #     if (is.null(parameter.name)) {
    #         Names <- gsub("[[:digit:]]","", name.par)
    #     }
    Names <- gsub("[[:digit:]]","", name.par)
    parameter.name <- Names[parameter.idx]
    name.count <- table(parameter.name)
    number.idx <- gsub("[a-zA-Z]","", name.par)
    first.idx <- substr(number.idx, 1, 1)
    
    for (i in seq_along(parameter.name)) {
        if (name.count[names(name.count) == parameter.name[i]] != 1 &&
            ifelse(i == 1, TRUE, (parameter.name[i] != parameter.name[i - 1]))) {
            
            # find if subscript 0 in the names
            if (any(first.idx == 0)) {
                if (parameter.idx[i] != 1 && 
                    Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
                    # count <- parameter.idx[i] - 1
                    count <- as.numeric(number.idx[parameter.idx[i]]) - 1
                } else {
                    count <- 0
                }
            } else {
                if (parameter.idx[i] != 1 && 
                    Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
                    # count <- parameter.idx[i]
                    count <- as.numeric(number.idx[parameter.idx[i]])
                } else {
                    count <- 1
                }
            }
            # if (any(unlist(gregexpr(pattern ='0', colnames(draws))) > 0)) {
            #     if (parameter.idx[i] != 1 && 
            #         Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
            #         count <- parameter.idx[i] - 1
            #     } else {
            #         count <- 0
            #     }
            # } else {
            #     if (parameter.idx[i] != 1 && 
            #         Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
            #         count <- parameter.idx[i]
            #     } else {
            #         count <- 1
            #     }
            # }
            
            label <- as.expression(substitute(A[x], list(A = as.name(parameter.name[i]),
                                                         x = count)))
            main.trace <- substitute(paste("Trace of ", A[x]),
                                     list(A = as.name(parameter.name[i]), x = count))
            main.acf <- substitute(paste("ACF of ", A[x]),
                                   list(A = as.name(parameter.name[i]), x = count))
            main.hist <- substitute(paste("Histogram of ", A[x]),
                                    list(A = as.name(parameter.name[i]), x = count))
            main.ergo <- substitute(paste("Ergodic mean of ", A[x]),
                                    list(A = as.name(parameter.name[i]), x = count))
            count <- count + 1
            
        } else if (name.count[names(name.count) == parameter.name[i]] != 1 &&
                   (parameter.name[i] == parameter.name[i - 1])) {
            label <- as.expression(substitute(A[x], list(A = as.name(parameter.name[i]),
                                                         x = count)))
            main.trace <- substitute(paste("Trace of ", A[x]),
                                     list(A = as.name(parameter.name[i]), x = count))
            main.acf <- substitute(paste("ACF of ", A[x]),
                                   list(A = as.name(parameter.name[i]), x = count))
            main.hist <- substitute(paste("Histogram of ", A[x]),
                                    list(A = as.name(parameter.name[i]), x = count))
            main.ergo <- substitute(paste("Ergodic mean of ", A[x]),
                                    list(A = as.name(parameter.name[i]), x = count))
            count <- count + 1
        }
        
        if (name.count[names(name.count) == parameter.name[i]] == 1) {
            if (parameter.idx[i] != 1 && 
                Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
                if (any(first.idx) > 0) {
                    count <- count + parameter.idx[i] - 1
                    # count <- parameter.idx[i] - 1
                } else {
                    count <- count + parameter.idx[i]
                    # count <- parameter.idx[i]
                }
                label <- as.expression(substitute(A[x], list(A = as.name(parameter.name[i]),
                                                             x = count)))
                main.trace <- substitute(paste("Trace of ", A[x]),
                                         list(A = as.name(parameter.name[i]), x = count))
                main.acf <- substitute(paste("ACF of ", A[x]),
                                       list(A = as.name(parameter.name[i]), x = count))
                main.hist <- substitute(paste("Histogram of ", A[x]),
                                        list(A = as.name(parameter.name[i]), x = count))
                main.ergo <- substitute(paste("Ergodic mean of ", A[x]),
                                        list(A = as.name(parameter.name[i]), x = count))
                count <- count + 1
            } else {
                label <- as.expression(substitute(A, list(A = as.name(parameter.name[i]))))
                main.trace <- substitute(paste("Trace of ", A),
                                         list(A = as.name(parameter.name[i])))
                main.acf <- substitute(paste("ACF of ", A),
                                       list(A = as.name(parameter.name[i])))
                main.hist <- substitute(paste("Histogram of ", A),
                                        list(A = as.name(parameter.name[i])))
                main.ergo <- substitute(paste("Ergodic mean of ", A),
                                        list(A = as.name(parameter.name[i])))
            }
        }
        
        # ------------------------------
        plot(draws[, parameter.idx[i]], type = "l", ylab = label,
             main = main.trace)
        
        # ----------
        # ACF
        acf(draws[, parameter.idx[i]], main = main.acf)
        
        # ----------
        # Histogram
        hist(draws[, parameter.idx[i]], main = paste(main.hist, hist.title), 
             freq = FALSE, breaks = breaks, xlab = "value", col = col, border = FALSE)
        #         hist(draws[, parameter.idx[i]], main = paste(hist.title, "Histogram of", label), 
        #              freq = FALSE, breaks = breaks, xlab = label, col = col, border = FALSE)
        #         hist(draws[, parameter.idx[i]], main = paste(dataname, i), freq = FALSE,
        #              breaks = breaks, xlab = label, col = col, border = FALSE)
        #             hist(draws[, parameter.idx[i]], main = paste(dataname, i), freq = FALSE,
        #                  breaks = breaks, xlab = label, col = col, border = FALSE)
        lines(density(draws[, parameter.idx[i]]), col = "red") 
        
        # Ergodic mean
        plot(cumsum(draws[, parameter.idx[i]]) / (1:len), type = "l",
             xlab = "iteration", ylab = label, lwd = 2, main = main.ergo,
             ylim = c(min(draws[, parameter.idx[i]]), max(draws[, parameter.idx[i]])))
        if (isquan) {
            lines(x = 1:len, y = cquantile(draws[, parameter.idx[i]])[, 1],
                  col = "red", lty = 2)
            lines(x = 1:len, y = cquantile(draws[, parameter.idx[i]])[, 2],
                  col = "red", lty = 2)
        }
    }
    par(mfrow = c(1, 1))
}


plot_parameter <- function(MCMCsample, draws = NULL, parameter.idx = NULL, 
                           sampleidx, name.par, main.title = NULL, isquan = TRUE,
                           is.trace = TRUE,
                           is.acf = TRUE, is.hist = TRUE, is.ergodic = TRUE) {
    if (is.null(draws)) draws <- MCMCsample$draws
    col_num <- sum(is.trace, is.acf, is.hist, is.ergodic)
    ifelse(is.null(parameter.idx), par(mfrow = c(length(parameter.name), col_num)),
           par(mfrow = c(length(parameter.idx), col_num)))
    
    # ifelse(is.null(parameter.idx), par(mfrow = c(col_num,length(parameter.name))),
    #        par(mfrow = c(col_num, length(parameter.idx))))
    
    if(is.null(main.title)) {
        par(mar = c(4, 4, 4, 1), oma = c(0, 0, 0, 0))
    } else {
        par(mar = c(4, 4, 4, 1), oma = c(0, 0, 2, 0))
    }
    
    len <- length(sampleidx)
    breaks <- floor(len / 50)
    col <- "navy"
    if (!is.null(MCMCsample)) dataname <- deparse(substitute(MCMCsample))
    
    #     if (is.null(parameter.name)) {
    #         Names <- gsub("[[:digit:]]","", name.par)
    #     }
    Names <- gsub("[[:digit:]]","", name.par)
    parameter.name <- Names[parameter.idx]
    name.count <- table(parameter.name)
    number.idx <- gsub("[a-zA-Z]","", name.par)
    first.idx <- substr(number.idx, 1, 1)
    
    for (i in seq_along(parameter.name)) {
        if (name.count[names(name.count) == parameter.name[i]] != 1 &&
            ifelse(i == 1, TRUE, (parameter.name[i] != parameter.name[i - 1]))) {
            
            # find if subscript 0 in the names
            if (any(first.idx == 0)) {
                if (parameter.idx[i] != 1 && 
                    Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
                    # count <- parameter.idx[i] - 1
                    count <- as.numeric(number.idx[parameter.idx[i]]) - 1
                } else {
                    count <- 0
                }
            } else {
                if (parameter.idx[i] != 1 && 
                    Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
                    # count <- parameter.idx[i]
                    count <- as.numeric(number.idx[parameter.idx[i]])
                } else {
                    count <- 1
                }
            }
            # if (any(unlist(gregexpr(pattern ='0', colnames(draws))) > 0)) {
            #     if (parameter.idx[i] != 1 && 
            #         Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
            #         count <- parameter.idx[i] - 1
            #     } else {
            #         count <- 0
            #     }
            # } else {
            #     if (parameter.idx[i] != 1 && 
            #         Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
            #         count <- parameter.idx[i]
            #     } else {
            #         count <- 1
            #     }
            # }
            
            label <- as.expression(substitute(A[x], list(A = as.name(parameter.name[i]),
                                                         x = count)))
            main.trace <- substitute(paste("Trace of ", A[x]),
                                     list(A = as.name(parameter.name[i]), x = count))
            main.acf <- substitute(paste("ACF of ", A[x]),
                                   list(A = as.name(parameter.name[i]), x = count))
            main.hist <- substitute(paste("Histogram of ", A[x]),
                                    list(A = as.name(parameter.name[i]), x = count))
            main.ergo <- substitute(paste("Ergodic mean of ", A[x]),
                                    list(A = as.name(parameter.name[i]), x = count))
            count <- count + 1
            
        } else if (name.count[names(name.count) == parameter.name[i]] != 1 &&
                   (parameter.name[i] == parameter.name[i - 1])) {
            label <- as.expression(substitute(A[x], list(A = as.name(parameter.name[i]),
                                                         x = count)))
            main.trace <- substitute(paste("Trace of ", A[x]),
                                     list(A = as.name(parameter.name[i]), x = count))
            main.acf <- substitute(paste("ACF of ", A[x]),
                                   list(A = as.name(parameter.name[i]), x = count))
            main.hist <- substitute(paste("Histogram of ", A[x]),
                                    list(A = as.name(parameter.name[i]), x = count))
            main.ergo <- substitute(paste("Ergodic mean of ", A[x]),
                                    list(A = as.name(parameter.name[i]), x = count))
            count <- count + 1
        }
        
        if (name.count[names(name.count) == parameter.name[i]] == 1) {
            if (parameter.idx[i] != 1 && 
                Names[parameter.idx[i]] == Names[parameter.idx[i] - 1]) {
                if (any(first.idx) > 0) {
                    count <- count + parameter.idx[i] - 1
                    # count <- parameter.idx[i] - 1
                } else {
                    count <- count + parameter.idx[i]
                    # count <- parameter.idx[i]
                }
                label <- as.expression(substitute(A[x], list(A = as.name(parameter.name[i]),
                                                             x = count)))
                main.trace <- substitute(paste("Trace of ", A[x]),
                                         list(A = as.name(parameter.name[i]), x = count))
                main.acf <- substitute(paste("ACF of ", A[x]),
                                       list(A = as.name(parameter.name[i]), x = count))
                main.hist <- substitute(paste("Histogram of ", A[x]),
                                        list(A = as.name(parameter.name[i]), x = count))
                main.ergo <- substitute(paste("Ergodic mean of ", A[x]),
                                        list(A = as.name(parameter.name[i]), x = count))
                count <- count + 1
            } else {
                label <- as.expression(substitute(A, list(A = as.name(parameter.name[i]))))
                main.trace <- substitute(paste("Trace of ", A),
                                         list(A = as.name(parameter.name[i])))
                main.acf <- substitute(paste("ACF of ", A),
                                       list(A = as.name(parameter.name[i])))
                main.hist <- substitute(paste("Histogram of ", A),
                                        list(A = as.name(parameter.name[i])))
                main.ergo <- substitute(paste("Ergodic mean of ", A),
                                        list(A = as.name(parameter.name[i])))
            }
        }
        
        # ------------------------------
        if(is.trace) {
            plot(draws[, parameter.idx[i]], type = "l", ylab = label,
                 main = main.trace)
        }
        
        
        # ----------
        # ACF
        if(is.acf) {
            acf(draws[, parameter.idx[i]], main = main.acf)   
        }
        
        
        # ----------
        # Histogram
        # hist(draws[, parameter.idx[i]], main = paste(main.hist, hist.title), 
        #      freq = FALSE, breaks = breaks, xlab = "value", col = col, border = FALSE)
        if(is.hist) {
            hist(draws[, parameter.idx[i]], main = main.hist, 
                 freq = FALSE, breaks = breaks, xlab = "value", col = col, border = FALSE)
        }
        
        #         hist(draws[, parameter.idx[i]], main = paste(hist.title, "Histogram of", label), 
        #              freq = FALSE, breaks = breaks, xlab = label, col = col, border = FALSE)
        #         hist(draws[, parameter.idx[i]], main = paste(dataname, i), freq = FALSE,
        #              breaks = breaks, xlab = label, col = col, border = FALSE)
        #             hist(draws[, parameter.idx[i]], main = paste(dataname, i), freq = FALSE,
        #                  breaks = breaks, xlab = label, col = col, border = FALSE)
        lines(density(draws[, parameter.idx[i]]), col = "red") 
        
        # Ergodic mean
        if(is.ergodic) {
            plot(cumsum(draws[, parameter.idx[i]]) / (1:len), type = "l",
                 xlab = "iteration", ylab = label, lwd = 2, main = main.ergo,
                 ylim = c(min(draws[, parameter.idx[i]]), max(draws[, parameter.idx[i]])))
            if (isquan) {
                lines(x = 1:len, y = cquantile(draws[, parameter.idx[i]])[, 1],
                      col = "red", lty = 2)
                lines(x = 1:len, y = cquantile(draws[, parameter.idx[i]])[, 2],
                      col = "red", lty = 2)
            }
        }
    }
    mtext(main.title, outer = TRUE, cex = 1.2)
    par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
}



createlabel <- function(k, N) {
    regionlabel_vec <- rep(0, N * N)
    for (i in 1:(N / k)) {
        for(j in 1:(N / k)) {
            for(m in 1:k) {
                regionlabel_vec[c((j - 1) * N * k + N * (m - 1) + 1:k) + (i - 1) * k] = 
                    i + (j - 1) * (N / k)
            }
        }
    }
    return(regionlabel_vec)
}


gridline <- function(k, N, col = "white", lwd = 0.5) {
    total = (N / k) - 1
    u = N / k
    L = floor(k)
    U = ceiling(k)
    if (L != k) {
        idx = cumsum(rep(c(L, U), (u / 2)))
        print(idx)
        for (i in 1:u) {
            xline(-0.025 + (1.05 / N) * idx[i], col = col, lwd = lwd)
            # xline(-0.025 + (1.05 / u) * i, col = col)
            yline(-0.025 + (1.05 / N) * idx[i], col = col, lwd = lwd)
            # yline(-0.025 + (1.05 / u) * i, col = col)
        }
    } else {
        for (i in 1:total) {
            xline(-0.025 + (1.05 / u) * i, col = col, lwd = lwd)
            # xline(-0.025 + (1.05 / u) * i, col = col)
            yline(-0.025 + (1.05 / u) * i, col = col, lwd = lwd)
            # yline(-0.025 + (1.05 / u) * i, col = col)
        }
    }
}

assign_region_label <- function(k, N, is.mat = TRUE) {
    # Assign region labels
    #
    # Args:
    #   k: The length of a region
    #   N: The length of an image
    #   is.mat: Return a matrix or not? 
    # Returns:
    #   A matrix with region numbers if is.mat is true, a vector otherwise
    
    regionlabel.vec <- rep(0, N * N)
    for (i in 1:(N / k)) {
        for (j in 1:(N / k)) {
            for (m in 1:k) {
                regionlabel.vec[c((j - 1) * N * k + N * (m - 1) + 1:k) + (i - 1) * k] <- 
                    i + (j - 1) * (N / k)
            }
        }
    }
    if (is.mat == TRUE) {
        return(matrix(regionlabel.vec, ncol = N))
    } else {
        return(regionlabel.vec)
    }
}

get_voxel_number <- function(region_mat) {
    G <- max(region_mat)
    voxel_number <- rep(0, G)
    for(i in 1:G) {
        voxel_number[i] <- length(which(region_mat == i))
    }
    return(voxel_number)
}

rotate <- function(m) t(m)[, nrow(m):1] 


# image of the mean over time
img_mean <- function(x) {
    apply(x, c(2, 3), mean)
}


img_time_plot <- function (x, time, datalist, slice_no = NULL) {
    if (!is.null(slice_no)) {
        title <- paste(x, "time = ", time, "slice = ", slice_no)
    } else {
        title <- paste(x, "time = ", time)
    }
    image.plot(datalist[[x]][time, , ], axes = F, 
               main = title)
}

img_plot <- function (x, datalist, slice_no = NULL) {
    if (!is.null(slice_no)) {
        title <- paste(x, "mean", "slice = ", slice_no)
    } else {
        title <- paste(x, "mean")
    }
    image.plot(datalist[[x]], axes = F, 
               main = title)
}

makecplx <- function(YY) {
    n <- dim(YY)[3]
    YY[, , 1:(n / 2)] + 1i * YY[, , ((n / 2) + 1):n]
}

# Activation and Posterior Probability maps
color.bar <- function(lut, min, max = -min, nticks = 11, 
                      ticks = seq(min, max, len = nticks), title = '') {
    scale = (length(lut) - 1) / (max - min)
    
    # dev.new(width=1.75, height=5)
    #     par(mfrow = c(1, 1))
    # par(mar = c(10,2,3,1))
    plot(c(0, 1), c(min, max), type = 'n', bty = 'n', xaxt = 'n', 
         xlab = '', yaxt = 'n', ylab = '', main = title)
    axis(4, round(ticks, 4), las = 1, cex.axis = .8)
    for (i in 1:(length(lut) - 1)) {
        y = (i - 1) / scale + min
        rect(0, y, 1, y + 1 / scale, col = lut[i], border = NA)
    }
}

# par(mar = c(2, 0, 0, 0))
# plot(c(0, 1), c(min, max), type = 'n', bty = 'n', xaxt = 'n', 
#      xlab = '', yaxt = 'n', ylab = '', main = title)
# axis(1, round(ticks, 4), las = 1, cex.axis = .8)
# for (i in 1:(length(lut) - 1)) {
#     x = (i - 1) / scale + min
#     rect(x, 0, x + 1 / scale, 1, col = lut[i], border = NA)
# }



# Function to compute hermitian of a matrix
herm <- function(cplxmat) {
    return (t(Conj(cplxmat)))
}

# Function to compute power exponential covariance function
power_exp_cov_fcn <- function(centroid_mat, r) {
    # compute power exponential covariance function
    # C(dist) = sig2 * exp(-|dist / phi| ^ nu)
    # grid <- expand.grid(idx, idx)
    distance <- as.matrix(dist(centroid_mat))
    # powerexp <- exp(-(distance / r))
    # return(powerexp)
    return(exp(-(distance / r)))
}

# Function to compute the centroid of each region
compute_centroid <- function(region_mat) {
    G <- max(region_mat)
    centroid_mat <- matrix(0, nrow = G, ncol = 2)
    for (g in 1:G) {
        xycoord <- which(region_mat == g, arr.ind = TRUE)
        # centroid_mat[g, ] <- apply(xycoord, 2, mean)
        centroid_mat[g, ] <- colMeans(xycoord)
    }
    return(centroid_mat)
}


bezierplot <- function(x.vec, y.vec, nu, phi) {
    x.len <- length(x.vec)
    y.len <- length(y.vec)
    dist.sq <- matrix(0, nrow = x.len, ncol = y.len)
    K <- matrix(0, nrow = x.len, ncol = y.len)
    for (i in 1:x.len) {
        for (j in 1:y.len) {
            dist.sq[i, j] <- x.vec[i] ^ 2 + y.vec[j] ^ 2
        }
    }
    less.mat <- sqrt(dist.sq) < phi
    K[less.mat] <- (1 - dist.sq[less.mat] / (phi ^ 2)) ^ nu
    #     rowsum <- rowSums(K)
    #     cond <- rowsum != 0
    #     K[cond, ] <- K[cond, ] / rowsum[cond]
    return(K)
} 

Gaussianplot <- function(x.vec, y.vec, phi) {
    x.len <- length(x.vec)
    y.len <- length(y.vec)
    dist.sq <- matrix(0, nrow = x.len, ncol = y.len)
    K <- matrix(0, nrow = x.len, ncol = y.len)
    for (i in 1:x.len) {
        for (j in 1:y.len) {
            dist.sq[i, j] <- x.vec[i] ^ 2 + y.vec[j] ^ 2
        }
    }
    K <- exp(- dist.sq / (2*phi))
    return(K)
} 

# Bezier Kernel 
dist_mat <- function(coord, grid) {
    N.voxels <- dim(coord)[1]  # number of voxels
    G <- dim(grid)[1]
    dist.mat <- matrix(0, nrow = N.voxels, ncol = G)
    for (i in 1:N.voxels) {
        for (j in 1:G) {
            dist.mat[i, j] <- sqrt(sum((coord[i, ] - grid[j, ]) ^ 2))
        }
    }
    return(dist.mat)
} 

bezier2 <- function(dist.mat, nu, phi) {
    # 2-d spherical Bezier kernel
    # nu: smoothness
    # range parameter 
    # a bend matrix
    # coord: resolution of data (voxel observation site)
    # grid: resolution of process convolution (region site)
    K <- array(0, dim = dim(dist.mat))
    less.mat <- dist.mat < phi
    K[less.mat] <- (1 - (dist.mat[less.mat] / phi) ^ 2) ^ nu
    return(K)
}

normalized_bezier2 <- function(dist.mat, nu, phi) {
    K <- bezier2(dist.mat, nu, phi)
    rowsum <- rowSums(K)
    cond <- rowsum != 0
    K[cond, ] <- K[cond, ] / rowsum[cond]
    return(K)
}


compute_posprob_pc <- function(model, N.voxels = 400) {
    # posprob <- rep(0, N*N)
    if(is.null(colnames(model$draws))) {
        gammaidx <- 1:N.voxels
    } else {
        gammaidx <- grep("gamma", colnames(model$draws))
    }
    len <- nrow(model$draws)
    # not faster than for loop
    # test.posprob = apply(ssvsconjsnr05$psi[1, sampleidx, 1, , ], c(2, 3), 
    #                      function(x) sum(x)/length(x))
    GammaData <- model$draws[, gammaidx]
    posprob <- colSums(GammaData) / len
    return(posprob)
}


cquantile <- function(z, probs = c(0.025, 0.975)) {
    cquant <- matrix(0, nrow = length(z), length(probs))
    for (i in seq_along(z)) {
        cquant[i, ] <- quantile(z[1:i], probs = probs)
    }
    cquant <- as.data.frame(cquant)
    return(cquant)
}

quantile_parameter <- function(MCMCsample, draws = NULL, parameter.idx = NULL, 
                               parameter.name = NULL, prob) {
    if (is.null(draws)) draws <- MCMCsample$draws
    if (is.null(parameter.idx)) {
        draws <- MCMCsample$draws[, parameter.name]
    } else {
        draws <- MCMCsample$draws[, parameter.idx, drop = FALSE]
    }
    round(apply(draws, 2, quantile, prob = prob), 4)
}

summary_parameter <- function(MCMCsample, parameter.idx = NULL, 
                              parameter.name = NULL, prob, 
                              title = MCMCsample) {
    if (is.null(parameter.idx)) {
        draws <- MCMCsample$draws[, parameter.name]
    } else {
        draws <- MCMCsample$draws[, parameter.idx]
    }
    quan <- quantile_parameter(MCMCsample = MCMCsample, 
                               parameter.idx = parameter.idx, 
                               parameter.name = parameter.name, 
                               prob = prob)
    result <- round(rbind(apply(draws, 2, effectiveSize), 
                          apply(draws, 2, mean), 
                          quan), 4)
    rownames(result) <- c("effective_size", "post_mean", 
                          paste0(rownames(quan), "_quantile"))
    print(paste("posterior summary of", deparse(substitute(title))))
    print(result)
}


Rhat <- function(x) {
    m <- dim(x)[2]  # number of chains
    n <- dim(x)[1]  # number of iterations per chain
    phij.mean <- colMeans(x)
    B <- n / (m - 1) * sum((phij.mean - mean(phij)) ^ 2)
    W <- 1 / (m * (n - 1)) * sum(apply(t(x) - phij.mean, 1, function(a) sum(a ^ 2)))
    var.phi <- (n - 1) / n * W + 1 / n * B
    return(sqrt(var.phi / W))
}

# library(kernlab)
GaussKer <- function(dist.mat, phi = 1) {
    # 2-d Gaussian Kernel (radial basis function, RBF)
    # defined as exp(-||x - x'|| ^ 2 / (2 * sigma ^ 2) )
    # phi = 1 / sigma ^ 2, i.e., precision
    # coord: resolution of data (voxel observation site)
    # grid: resolution of process convolution (region site)
    ker <- exp(- dist.mat ^ 2 / (2 * phi))
    return(ker)
}



get_post_sample <- function(result, name.par, idx) {
    # idx <- 1:Nvoxels
    draws <- result$draws[, -idx]
    colnames(draws) <- name.par[-idx]
    return(draws)
}

color.bar <- function(lut, min, max = -min, nticks = 11, 
                      ticks = seq(min, max, len = nticks), title = '') {
    scale = (length(lut) - 1) / (max - min)
    
    # dev.new(width=1.75, height=5)
    #     par(mfrow = c(1, 1))
    # par(mar = c(10,2,3,1))
    plot(c(0, 1), c(min, max), type = 'n', bty = 'n', xaxt = 'n', 
         xlab = '', yaxt = 'n', ylab = '', main = title)
    axis(4, round(ticks, 4), las = 1, cex.axis = .8)
    for (i in 1:(length(lut) - 1)) {
        y = (i - 1) / scale + min
        rect(0, y, 1, y + 1 / scale, col = lut[i], border = NA)
    }
}


## human
brain_grid_plot <- function(brainmod, posprob, N, G, mdltype, thre,
                            pch = 15, cex = .69, isbar = TRUE, isact = FALSE,
                            istitle = TRUE) {
    xx = 1.006
    yy = -0.006
    xlimit = c(yy, xx)
    library(Hmisc)
    par(mfrow = c(1, 1))
    par(omi = c(0, 0, 0, 0))
    if (isbar) {
        par(mar = c(.5, .5, 2, 5))
    } else {
        par(mar = c(.5, .5, 2, .5))
    }
    if(istitle) {
        if(isact) {
            main = paste(mdltype, "Activation Map")
        } else {
            main = paste(mdltype, "Posterior Probability Map")
        }
    } else {
        main = ""
        par(mar = c(.5, .5, .5, .5))
    }
    
    image(brainmod, col = gray((0:N) / N), axes = FALSE, 
          main = main)
    k <- sqrt(G)
    for (i in 1:(k-1)) {
        xline(yy + (xx/k)*i, col = "yellow")
        yline(yy + (xx/k)*i, col = "yellow")
    }
    for (i in 1:N) {
        ll <- length(which(matrix(posprob, nrow = N)[i, ] > thre))
        idx <- which(matrix(posprob, nrow = N)[i, ] > thre)
        brk <- 150
        if (isact) {
            if (ll > 0) {
                points(rep((i - .5)/N, ll), 
                       (which(matrix(posprob, nrow = N)[i, ] > thre) - .5)/N, 
                       col = "red", pch = pch, cex = cex)
            }
        } else {
            rbPal <- colorRampPalette(tim.colors()[seq(1, 55)])
            Col1 <- rbPal(brk)[as.numeric(cut(posprob, breaks = brk))]
            Color <- matrix(Col1, nc = N, nr = N)
            if (ll > 0) {
                points(rep((i - .5)/N, ll), 
                       (which(matrix(posprob, nrow = N)[i, ] > thre) - .5)/N, 
                       col = Color[i, idx], pch = pch, cex = cex)
            }
        }
    }
    mini = min(posprob)
    maxi = max(posprob)
    if(isbar) {
        subplot(color.bar(rbPal(brk), min = mini, max = maxi), 
                x = .995, y = .5, size = c(.5, 5))
    }
}

## simulated
brain_grid_plot <- function(brainmod, posprob, N, G, mdltype, thre,
                            pch = 15, cex = .69, title = NULL,
                            is.color.bar = FALSE, is.grid = FALSE,
                            is.points = FALSE) {
    library(Hmisc)
    # par(mfrow = c(1, 1))
    # par(omi = c(0, 0, 0, 0))
    # par(mar = c(.85, .5, 2, 5))
    if (is.null(title)) {
        title <- paste(mdltype, "Posterior Probability")
    } 
    image(brainmod, col = gray((0:N) / N), axes = FALSE, main = title)
    # image(brainmod, col = gray((0:N) / N), axes = FALSE, 
    #       main = paste(mdltype, "Posterior Probability"))
    k <- sqrt(G)
    if (is.grid) {
        for (i in 1:(k-1)) {
            xline(yy + (xx/k)*i, col = "yellow")
            yline(yy + (xx/k)*i, col = "yellow")
        }
    }
    
    if (is.points) {
        loc.y <- seq(0.05, 0.98, by = 0.125)
        loc.x <- seq(0.05, 0.98, by = 0.125)
        for (i in 1:k) {
            for (j in 1:k) {
                points(loc.x[i], loc.y[j], pch = 16, col = "gold", cex = 0.8)
            } 
        }
    }
    
    for (i in 1:N) {
        ll <- length(which(matrix(posprob, nrow = N)[i, ] > thre))
        idx <- which(matrix(posprob, nrow = N)[i, ] > thre)
        brk <- 150
        rbPal <- colorRampPalette(tim.colors()[seq(1, 55)])
        Col1 <- rbPal(brk)[as.numeric(cut(posprob, breaks = brk))]
        Color <- matrix(Col1, nc = N, nr = N)
        if (ll > 0) {
            points(rep((i - .5)/N, ll), 
                   (which(matrix(posprob, nrow = N)[i, ] > thre) - .5)/N, 
                   col = Color[i, idx], pch = pch, cex = cex)
        }
    }
    mini = min(posprob)
    maxi = max(posprob)
    # subplot(color.bar(rbPal(brk), min = mini, max = maxi), 
    #         x = .995, y = .5, size = c(.5, 8))
    if(is.color.bar) {
        subplot(color.bar(rbPal(brk), min = mini, max = maxi), 
                x = 1.01, y = .5, size = c(.5, 5.5))
    }
}





