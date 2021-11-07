################################################################################
# Load data for Multiple-slice analysis                                        #
# ./clean_multi_slice_data.R                                                   #
# Cheng-Han Yu, Marquette University                                           #
################################################################################
## Clean the data loaded from ./load_multi_slice_data.R and create data objects 
## for running the proposed algorithms

if (!exists(multi_slice_data_re)) {
    source("./data/load_multi_slice_data.R")
}

multi_slice_data_im_arr <- lapply(multi_slice_data_im, function(x){
    x <- as.matrix(x)
    return(array(x, dim = c(96, 96, 490)))
})

multi_slice_data_re_arr <- lapply(multi_slice_data_re, function(x){
    x <- as.matrix(x)
    return(array(x, dim = c(96, 96, 490)))
})

## 490 x 96 x 96
multi_slice_data_im_perm <- lapply(multi_slice_data_im_arr, function(x) {
    return(aperm::aperm(x, c(3, 1, 2)))
})

multi_slice_data_re_perm <- lapply(multi_slice_data_re_arr, function(x) {
    return(aperm::aperm(x, c(3, 1, 2)))
})

par(mfcol = c(7, 2))
lapply(1:7, function(k) {
    par(mar = c(.1, .1, 1, 1))
    image.plot(rotate(multi_slice_data_re_perm[[k]][1, , ]), 
               axes = F, main = paste("Re", k))
    image.plot(rotate(multi_slice_data_im_perm[[k]][1, , ]), 
               axes = F, main = paste("Im", k))
})

multi_slice_data_cplx_perm <- vector("list", length = 7)
multi_slice_data_cplx_perm <- lapply(multi_slice_data_cplx_perm, 
                                     function(x) array(0, dim = c(490, 96, 96)))

for (k in 1:7) {
    for (i in 1:96) {
        for (j in 1:96) {
            multi_slice_data_cplx_perm[[k]][, i, j] <- 
                multi_slice_data_re_perm[[k]][, i, j] + 
                1i * multi_slice_data_im_perm[[k]][, i, j]
        }
    }
}

multi_slice_data_mod_perm <- lapply(multi_slice_data_cplx_perm, Mod)
multi_slice_data_arg_perm <- lapply(multi_slice_data_cplx_perm, Arg)


dataM_lst <- dataA_lst <- dataR_lst <- dataI_lst <- 
    lapply(vector("list", length = 7), function(x) array(0, c(490, 96, 96)))

for (k in 1:7) {
    for (t in 1:490) {
        dataM_lst[[k]][t, , ] = rotate(multi_slice_data_mod_perm[[k]][t, , ])
        dataA_lst[[k]][t, , ] = rotate(multi_slice_data_arg_perm[[k]][t, , ])
        dataR_lst[[k]][t, , ] = rotate(multi_slice_data_re_perm[[k]][t, , ])
        dataI_lst[[k]][t, , ] = rotate(multi_slice_data_im_perm[[k]][t, , ])
    }
}

multi_slice_data_lst <- list('Re' = dataR_lst, 'Im' = dataI_lst, 
                             'Mod' = dataM_lst, 'Arg' = dataA_lst)


multi_slice_data_mean_lst <- lapply(multi_slice_data_lst, 
                                    function(x) lapply(x, img_mean))

