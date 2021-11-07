################################################################################
# Create region labels                                                         #
# ./data/regionlabel.R                                                              #
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

regionlabel_100 <- assign_region_label(k = 2, N = N, is.mat = FALSE)
regionlabel_25 <- assign_region_label(k = 4, N = N, is.mat = FALSE)
regionlabel_16 <- assign_region_label(k = 5, N = N, is.mat = FALSE)
regionlabel_4 <- assign_region_label(k = 10, N = N, is.mat = FALSE)

region_mat_100 <- assign_region_label(k = 2, N = N, is.mat = TRUE)
region_mat_25 <- assign_region_label(k = 4, N = N, is.mat = TRUE)
region_mat_16 <- assign_region_label(k = 5, N = N, is.mat = TRUE)
region_mat_4 <- assign_region_label(k = 10, N = N, is.mat = TRUE)

G_100 <- max(region_mat_100)  # Number of groups
G_25 <- max(region_mat_25)  # Number of groups
G_16 <- max(region_mat_16)  # Number of groups
G_4 <- max(region_mat_4)  # Number of groups

get_voxel_number <- function(region_mat) {
    G <- max(region_mat)
    voxel_number <- rep(0, G)
    for(i in 1:G) {
        voxel_number[i] <- length(which(region_mat == i))
    }
    return(voxel_number)
}
