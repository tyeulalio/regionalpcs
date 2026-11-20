#' Compute regional principal components for methylation data
#'
#' @param meth Numeric matrix or data frame of methylation values,
#' with uniquely named CpGs in rows and samples in columns. Values must be
#' finite. The input scale affects dimension selection; see the vignette.
#' @param region_map Data frame mapping CpGs to regions. The first two columns
#' contain region IDs and CpG IDs, respectively. All mapped CpGs must occur in
#' the row names of `meth`. Region and within-region CpG order are preserved,
#' including repeated mappings.
#' @param pc_method Method to use for PC computation,
#' either 'gd' (Gavish-Donoho) or 'mp' (Marchenko-Pastur)
#' @param verbose Logical, should progress messages be displayed?
#'
#' @return A list with `regional_pcs` (retained PCs by samples),
#' `percent_variance` (fractions of within-region variance, not percentages),
#' `loadings` (a named list of CpG-by-PC data frames), and `info`
#' (a metadata matrix with region, estimated dimension, and CpG count).
#' @export
#'
#' @examples
#' # Create synthetic methylation data
#' meth_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(meth_data) <- paste0("CpG", 1:100)
#' colnames(meth_data) <- paste0("Sample", 1:10)
#'
#' # Create a synthetic region map
#' region_map_data <- data.frame(
#'     region_id = rep(c("Gene1", "Gene2"), each = 50),
#'     cpg_id = rownames(meth_data)
#' )
#'
#' # Run the function
#' compute_regional_pcs(meth_data, region_map_data, pc_method = 'gd')
#'
compute_regional_pcs <- function(meth, region_map,
    pc_method = c('gd', 'mp'), verbose = FALSE) {
    region_map <- .prepare_region_map(meth, region_map)
    pc_method <- match.arg(pc_method)
    msg <- switch(
        pc_method,
        gd = "Using Gavish-Donoho method",
        mp = "Using Marchenko-Pastur method"
    )
    message(msg)

    regions <- unique(region_map[['region_id']])
    num_regions <- length(regions)
    if (verbose) message("Number of regions:", num_regions)

    # Group indices once, retaining region order and repeated mappings.
    row_indices <- match(region_map$cpg_id, rownames(meth))
    region_rows <- split(row_indices,
        factor(region_map$region_id, levels = regions))
    region_res <- lapply(seq_along(regions), function(i) {
        .summarize_region_rows(regions[[i]], region_rows[[i]],
            meth, pc_method, verbose)
    })

    region_pcs <- do.call(rbind, lapply(region_res, combine_results, 'sig_pcs'))
    region_pct_var <- do.call(rbind,
        lapply(region_res, combine_results, 'percent_var'))

    region_loadings <- lapply(region_res, function(x) x$loadings)
    names(region_loadings) <- regions

    region_info <- do.call(
        rbind,
        lapply(region_res, function(x) x[c('region', 'est_dim', 'num_cpgs')])
        )
    return(list(
        regional_pcs = region_pcs,
        percent_variance = region_pct_var,
        loadings = region_loadings,
        info = region_info
        ))
}


#' Combine results dataframes across regions
#'
#' @param res List of lists; contains summarized region results
#' @param df_name String;
#' name of result being combined (sig_pcs or percent_var)
#'
#' @return Data Frame containing results
#' @export
#'
#' @examples
#' # Create example data for 'sig_pcs' and 'percent_var'
#'     sig_pcs_example <- data.frame(pcs = c("PC1", "PC2"),
#' value = c(0.2, 0.4))
#' percent_var_example <- data.frame(pcs = c("PC1", "PC2"),
#' value = c(0.7, 0.3))
#'
#' # Create 'res' list containing both 'sig_pcs' and 'percent_var'
#' res <- list(region = "Region1", sig_pcs = sig_pcs_example,
#' percent_var = percent_var_example)
#'
#' # Example function use: Combine 'sig_pcs' across regions
#' combined_sig_pcs <- combine_results(res, df_name = "sig_pcs")
#' print(combined_sig_pcs)
combine_results <- function(res, df_name) {
    formatted_df <- as.data.frame(t(res[[df_name]]))
    if (nrow(formatted_df) > 0L) {
        rownames(formatted_df) <- paste(res$region,
            rownames(formatted_df), sep = '-')
    }
    return(formatted_df)
}


#' Summarize a region using regional principal components
#'
#' @param region String; name of region being processed
#' @param meth Data frame or matrix;
#' Methylation values to summarize; rows=CpGs, columns=samples
#' @param region_map Data frame;
#' Mapping of CpGs to regions,
#' column 1 should be regions,
#' column 2 should be CpGs with the same names as the rows of meth.
#' Column names are ignored; all mapped CpGs must occur in meth.
#' @param pc_method String;
#' indicating the method for estimating dimension;
#' "gd"=Gavish-Donoho, "mp"=Marchenko-Pastur
#' @param verbose Boolean; print output statements
#'
#' @return list containing PC results
#' @export
#'
#' @examples
#' # Create the region map with just one region containing 10 CpGs
#' region_map <- data.frame(region_id = rep(1, 10), cpg_id = seq(1, 10))
#'
#' # Create methylation data frame
#' set.seed(123)
#' meth <- as.data.frame(matrix(runif(10 * 20, min = 0, max = 1), nrow = 10))
#' rownames(meth) <- seq(1, 10)
#'
#' # Call the function
#' summarize_region(1, region_map, meth, 'gd')
summarize_region <- function(region, region_map,
                                meth, pc_method, verbose = FALSE) {
    region_map <- .prepare_region_map(meth, region_map)
    pc_method <- match.arg(pc_method, c('gd', 'mp'))
    rows <- match(region_map$cpg_id[region_map$region_id == region],
        rownames(meth))
    if (length(rows) == 0L) {
        stop("No CpGs are mapped to the requested region.")
    }
    .summarize_region_rows(region, rows, meth, pc_method, verbose)
}

.summarize_region_rows <- function(region, rows, meth, pc_method, verbose) {
    if (verbose) {
        message("Processing ", region)
    }
    # Subset before conversion to avoid copying the full matrix.
    region_meth <- as.data.frame(meth[rows, , drop = FALSE])

    num_cpgs <- nrow(region_meth)

    pc_res <- get_sig_pcs(region_meth, pc_method = pc_method, verbose)

    pc_res$num_cpgs <- num_cpgs
    pc_res$region <- region

    return(pc_res)
}

.prepare_region_map <- function(meth, region_map) {
    if (!(is.matrix(meth) || is.data.frame(meth)) ||
        (is.matrix(meth) && !is.numeric(meth)) ||
        (is.data.frame(meth) && !all(vapply(meth, is.numeric, logical(1))))) {
        stop("'meth' must be a numeric matrix or data frame.")
    }
    if (nrow(meth) == 0L || ncol(meth) < 2L) {
        stop("'meth' must contain at least one CpG and two samples.")
    }
    ids <- rownames(meth)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        stop("'meth' must have unique, non-missing CpG row names.")
    }
    if (!is.data.frame(region_map) || ncol(region_map) < 2L ||
        nrow(region_map) == 0L) {
        stop("'region_map' must be a non-empty data frame with two columns.")
    }
    region_map <- as.data.frame(region_map[, seq_len(2), drop = FALSE])
    names(region_map) <- c('region_id', 'cpg_id')
    region_map[] <- lapply(region_map, as.character)
    if (anyNA(region_map) || any(!nzchar(region_map$region_id)) ||
        any(!nzchar(region_map$cpg_id))) {
        stop("'region_map' IDs must be non-missing and non-empty.")
    }
    if (anyNA(match(region_map$cpg_id, ids))) {
        stop("All CpGs in 'region_map' must occur in rownames(meth). ",
            "Filter the map to the measured CpGs before computing PCs.")
    }
    region_map
}
