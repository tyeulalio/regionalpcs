#' Get significant principal components
#'
#' @param x A data frame or matrix of methylation values;
#' rows = features, columns = samples
#' @param method String indicating the method for estimating dimension;
#' "gd" = Gavish-Donoho (default), "mp" = Marchenko-Pastur
#' @param verbose Boolean; print output statements
#'
#' @return List containing four elements;
#' sig_pcs = significant PCs,
#' percent_var = percent variance explained,
#' loadings = PC loadings,
#' est_dim = estimated dimension
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' x <- diag(4)
#' get_sig_pcs(x, "gd")
get_sig_pcs <- function(x, method = "gd", verbose = FALSE) {
    # Check input validity
    stopifnot(method %in% c("gd", "mp"))
    stopifnot(is.data.frame(x) | is.matrix(x))

    # Compute PCs using PCAtools package
    pca_res <- PCAtools::pca(x)

    # Extract eigenvalues and square them to get variance explained
    eig_sq <- pca_res$sdev^2

    # Placeholder for noise selection
    noise_select <- 1

    # Compute significant dimensions
    dims_res <- compute_dimension(x, eig_sq, noise_select, method)

    # Extract rotated PCs and loadings
    pcs <- pca_res$rotated
    loads <- pca_res$loadings

    # Extract estimated dimension
    est_dim <- dims_res[[1]]

    # Filter significant PCs
    sig_pcs <- pcs[, seq_len(est_dim), drop = FALSE] %>%
        as.data.frame()

    # Filter significant loadings
    sig_loads <- loads[, seq_len(est_dim), drop = FALSE] %>%
        as.data.frame()

    # Compute percent variance explained
    pct_var <- eig_sq / sum(eig_sq)
    subset_pct_var <- data.frame(
        percent_variance_explained = pct_var[seq_len(est_dim)]) %>%
        t() %>%
        as.data.frame()

    names(subset_pct_var) <- names(sig_pcs)

    # Package results into a list and return
    final_res <- list(
        sig_pcs = sig_pcs, percent_var = subset_pct_var,
        loadings = sig_loads, est_dim = est_dim
        )
    return(final_res)
}


#' Compute significant dimensions of a matrix
#' using the Marchenko-Pastur or Gavish-Donoho methods
#'
#' @param x A data frame or matrix of methylation values;
#' rows = features, columns = samples
#' @param var_explained A numeric vector containing the
#' variance explained by successive PCs,
#' sorted in decreasing order. (Used for PCAtools)
#' @param noise_select Numeric scalar specifying the
#' variance of the random noise (Used for PCAtools)
#' @param method String indicating the method for estimating dimension;
#' "gd" = Gavish-Donoho (default), "mp" = Marchenko-Pastur
#' @param verbose Boolean indicating whether to print
#' statements while running, default = FALSE
#'
#' @return Numeric scalar representing the optimal
#' number of PCs to retain using the specified method
#' @export
#'
#' @examples
#' x <- diag(4)
#' pca_res <- PCAtools::pca(x) # Run PCA
#' eig_sq <- pca_res$sdev^2 # Compute variance explained
#' compute_dimension(x, eig_sq, 1, "gd")
compute_dimension <- function(x, var_explained, noise_select,
                                method, verbose = FALSE) {
    # Initialize result variable
    dims_res <- NA

    # Select method for dimension estimation
    if (method == "gd") {
    if (verbose) message("Using Gavish-Donoho method")
        dims_res <- PCAtools::chooseGavishDonoho(x,
            var.explained = var_explained,
            noise = noise_select
            )
    }

    if (method == "mp") {
    if (verbose) message("Using Marchenko-Pastur method")
        dims_res <- PCAtools::chooseMarchenkoPastur(x,
            var.explained = var_explained,
            noise = noise_select
            )
    }

    # Return computed dimensions
    return(dims_res)
}
