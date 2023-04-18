#' Get significant principal components
#'
#' @param x A data frame or matrix of methylation value, rows=features, columns=samples
#' @param method String indicating the method for estimating dimension; "gd"=Gavish-Donoho (default), "mp"=Marchenko-Pastur
#' @param verbose Boolean; print output statements
#'
#' @return List containing four elements; sig_pcs = significant PCs, percent_var = percent variance explained, loadings = PC loadings, est_dim = estimated dimension
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' x <- diag(4)
#' get_sig_pcs(x, "gv")
get_sig_pcs <- function(x, method='gd', verbose=FALSE){
  # check input
  if (!method %in% c('gd', 'mp')){
    print("Error: method must be either 'gd' for Gavish-Donoho or 'mp' for Marchenko-Pastur" )
    return(-1)
  }
  if (!(is.data.frame(x) | is.matrix(x))){
    print("Error: x must be either a dataframe or matrix")
    return(-1)
  }

  # compute PCs of x
  # rows = features
  # cols = samples
  pca_res <- PCAtools::pca(x)

  # -- get dimension using method
  # compute variance explained
  eig_sq <- pca_res$sdev^2
  # select noise
  noise_select <- 1

  # compute significant dimension here
  dims_res <- compute_dimension(x, eig_sq, noise_select, method)


  # get the PCs
  pcs <- pca_res$rotated

  # get the loadings
  loads <- pca_res$loadings

  # get the estimated dimension
  est_dim <- dims_res[[1]]

  # filter to sig PCs
  sig_pcs <- pcs[,1:est_dim,drop=FALSE] %>%
    as.data.frame()
  sig_loads <- loads[,1:est_dim,drop=FALSE] %>%
    as.data.frame()

  # get percent variance explained
  pct_var <- eig_sq / sum(eig_sq)
  subset_pct_var <- data.frame(percent_variance_explained=pct_var[1:est_dim, drop=FALSE]) %>%
    t() %>%
    as.data.frame()
  names(subset_pct_var) <- names(sig_pcs)

  # return all of the results
  final_res <- list(sig_pcs=sig_pcs,
                    percent_var=subset_pct_var,
                    loadings=sig_loads,
                    est_dim=est_dim)
  final_res
}

#' Compute significant dimensions of a matrix using the Marchenko-Pastur or Gavish-Donoho methods
#'
#' @param x A data frame or matrix of methylation value, rows=features, columns=samples
#' @param var_explained A numeric vector containing the variance explained by successive PCs. This should be sorted in decreasing order. (Used for PCAtools)
#' @param noise_select Numeric scalar specifying the variance of the random noise (Used for PCAtools)
#' @param method String indicating the method for estimating dimension; "gd"=Gavish-Donoho (default), "mp"=Marchenko-Pastur
#' @param verbose Boolean indicating whether to print statements while running, default=FALSE
#'
#' @return Numeric scalar representing the optimal number of PCs to retain using the specified method
#' @export
#'
#' @examples
#' x <- diag(4)
#' pca_res <- PCAtools::pca(x) # run PCA
#' eig_sq <- pca_res$sdev^2 # compute variance explained
#' compute_dimension(x, eig_sq, 1, "gv")
compute_dimension <-function(x, var_explained, noise_select, method, verbose=FALSE){
  dims_res <- NA
  if (method == "gd"){
    if (verbose) print("Using Gavish-Donoho method")
    dims_res <- PCAtools::chooseGavishDonoho(x,
                                             var.explained=var_explained,
                                             noise=noise_select)
  }
  if (method == 'mp'){
    if (verbose) print("Using Marchenko-Pastur method")
    dims_res <- PCAtools::chooseMarchenkoPastur(x,
                                                var.explained=var_explained,
                                                noise=noise_select)
  }
  dims_res
}
