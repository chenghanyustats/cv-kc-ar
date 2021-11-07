################################################################################
# Load data for Multiple-slice analysis                                        #
# ./load_multi_slice_data.R                                                    #
# Cheng-Han Yu, Marquette University                                           #
################################################################################
# Create data Xmat, multi_slice_data_im and multi_slice_data_re
Xmat <- read.table("./data/Xmat.txt", sep = ",")
multi_slice_data_im <- list()
multi_slice_data_re <- list()

## 7 slices
for (i in 1:7) {
    data_im <- paste0("data", i, "I.txt")
    path_im <- paste0("./data/", data_im)
    multi_slice_data_im[[i]] <- read.table(path_im, sep = ",")
    
    data_re <- paste0("data", i, "R.txt")
    path_re <- paste0("./data/", data_re)
    multi_slice_data_re[[i]] <- read.table(path_re, sep = ",")
}




