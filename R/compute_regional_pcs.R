#' Compute regional PCs using methylation data for all regions included in the region map
#'
#' @param meth Data frame or matrix of methylation values to summarize; rows=CpGs, columns=samples
#' @param region_map Mapping of CpGs to regions, column 1 should be regions, column 2 should be CpGs with the same names as the rows of meth
#' @param pc_method String indicating the method for estimating dimension; "gd"=Gavish-Donoho (default), "mp"=Marchenko-Pastur
#' @param verbose Boolean; print output statements
#'
#' @return list containing PC results
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
compute_regional_pcs <- function(meth, region_map, pc_method='gd', verbose=FALSE){
  # perform input checks
  # check that rownames of meth == column 2 of region_map
  if (length(intersect(rownames(meth), region_map[,2])) == 0){
    print("Error: no matching values in rownames(meth) and column 2 of region map.")
    return(-1)
  }
  if (!pc_method %in% c('gd', 'mp')){
    print("Error: pc_method must be either 'gd' or 'mp'")
    return(-1)
  }
  if (pc_method == "gd"){
    print("Using Gavish-Donoho method")
  }
  if (pc_method == 'mp'){
    print("Using Marchenko-Pastur method")
  }

  # set column names here
  colnames(region_map)[1:2] <- c('region_id', 'cpg_id')
  region_map$region_id <- as.character(region_map$region_id)

  # get the regions
  regions <- unique(region_map[['region_id']])
  num_regions <- length(regions)
  if (verbose) print(paste("Number of regions:", num_regions))


  # summarize each of the region types
  region_res <- lapply(regions, summarize_region, region_map, meth, pc_method, verbose)


  # combine results across regions
  combine_dfs <- function(res, df_name){
    formatted_df <- res[[df_name]] %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column('pc') %>%
      dplyr::mutate(pc = paste(res$region, pc, sep='-')) %>%
      tibble::column_to_rownames('pc')
  }

  # combine into data frames
  region_pcs <- do.call(rbind, lapply(region_res, combine_dfs, 'sig_pcs'))
  region_pct_var <- do.call(rbind, lapply(region_res, combine_dfs, 'percent_var'))

  # combine loadings into one list
  region_loadings <- lapply(region_res, function(x) x$loadings)
  names(region_loadings) <- regions

  # combine the rest of the info
  region_info <- do.call(rbind, lapply(region_res,
                                       function(x) x[c('region', 'est_dim', 'num_cpgs')]))

  # return everythign in a list
  return_res <- list(regional_pcs=region_pcs,
                     percent_variance=region_pct_var,
                     loadings=region_loadings,
                     info=region_info)


  return(return_res)
}

#' Summarize a region using regional principal components
#'
#' @param region String; name of region being processed
#' @param meth Data frame or matrix; Methylation values to summarize; rows=CpGs, columns=samples
#' @param region_map Data frame; Mapping of CpGs to regions, column 1 should be regions, column 2 should be CpGs with the same names as the rows of meth
#' @param pc_method String; indicating the method for estimating dimension; "gd"=Gavish-Donoho (default), "mp"=Marchenko-Pastur
#' @param verbose Boolean; print output statements
#'
#' @return list containing PC results
#' @export
#'
#' @examples
summarize_region <- function(region, region_map, meth, pc_method, verbose){
  if (verbose){
    print(paste("Processing", region))
  }

  # summarize each region
  region_cpgs <- region_map %>%
    dplyr::filter(region_id == region)

  # subset methylation to the cpgs in the region
  region_meth <- as.data.frame(meth)[region_cpgs$cpg_id,,drop=FALSE]

  num_cpgs <- nrow(region_meth)

  # get the significant pcs for this region
  pc_res <- get_sig_pcs(region_meth, method=pc_method, verbose)

  # add number of cpgs and region name to results
  pc_res$num_cpgs <- num_cpgs
  pc_res$region <- region
  pc_res
}
