#' Compute regional principal components for methylation data
#'
#' @param meth Data frame of methylation beta values,
#' with CpGs in rows and samples in columns
#' @param region_map Data frame mapping CpGs to gene regions
#' @param pc_method Method to use for PC computation,
#' either 'gd' (Gavish-Donoho) or 'mp' (Marchenko-Pastur)
#' @param verbose Logical, should progress messages be displayed?
#'
#' @return A list containing several elements,
#' including the regional PCs, percent variance, and other information
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
    # Perform input checks
    stopifnot(length(intersect(rownames(meth), region_map[, 2])) > 0)
    pc_method <- match.arg(pc_method) # matching argument
    # output method being used for user
    msg <- switch(
        pc_method,
        gd = "Using Gavish-Donoho method",
        mp = "Using Marchenko-Pastur method"
    )
    message(msg)

    # Set column names for region_map
    colnames(region_map)[seq_len(2)] <- c('region_id', 'cpg_id')
    region_map$region_id <- as.character(region_map$region_id)

    # Get the unique regions
    regions <- unique(region_map[['region_id']])
    num_regions <- length(regions)
    if (verbose) message("Number of regions:", num_regions)

    # Summarize each of the region types
    region_res <- lapply(
        regions, summarize_region, region_map,
        meth, pc_method, verbose
    )

    # Combine results into data frames
    region_pcs <- do.call(rbind, lapply(region_res, combine_results, 'sig_pcs'))
    region_pct_var <- do.call(rbind,
        lapply(region_res, combine_results, 'percent_var'))

    # Combine loadings into one list
    region_loadings <- lapply(region_res, function(x) x$loadings)
    names(region_loadings) <- regions

    # Combine the other information
    region_info <- do.call(
        rbind,
        lapply(region_res, function(x) x[c('region', 'est_dim', 'num_cpgs')])
        )
    # Return all results in a list
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
    .data <- dplyr::`.data`
    # Transform and format data frame according to df_name
    formatted_df <- res[[df_name]] |>
        t() |>
        as.data.frame() |>
        tibble::rownames_to_column('pc') |>
        dplyr::mutate(pc = paste(res$region, .data$pc, sep = '-')) |>
        tibble::column_to_rownames('pc')

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
#' column 2 should be CpGs with the same names as the rows of meth
#' @param pc_method String;
#' indicating the method for estimating dimension;
#' "gd"=Gavish-Donoho (default), "mp"=Marchenko-Pastur
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
    # Output processing message if verbose is TRUE
    if (verbose) {
        message("Processing ", region)
    }
    .data <- dplyr::`.data`
    # Filter CpGs for the given region
    region_cpgs <- region_map |>
        dplyr::filter(.data$region_id == region)

    # Subset methylation data to CpGs in the region
    region_meth <- as.data.frame(meth)[region_cpgs$cpg_id, , drop = FALSE]

    # Get number of CpGs
    num_cpgs <- nrow(region_meth)

    # Get significant PCs for this region
    pc_res <- get_sig_pcs(region_meth, pc_method = pc_method, verbose)

    # Add metadata to the result
    pc_res$num_cpgs <- num_cpgs
    pc_res$region <- region

    return(pc_res)
}
