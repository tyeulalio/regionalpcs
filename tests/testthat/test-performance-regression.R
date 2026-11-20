# Reference implementation for testing region indexing and result assembly.
# Both implementations use the same PCA routine.
legacy_regional_pcs <- function(meth, region_map, method) {
    names(region_map)[1:2] <- c("region_id", "cpg_id")
    region_map$region_id <- as.character(region_map$region_id)
    regions <- unique(region_map$region_id)
    res <- lapply(regions, function(region) {
        selected <- dplyr::filter(region_map, .data$region_id == region)
        x <- as.data.frame(meth)[selected$cpg_id, , drop = FALSE]
        ans <- get_sig_pcs(x, method)
        ans$num_cpgs <- nrow(x)
        ans$region <- region
        ans
    })
    combine <- function(x, field) {
        out <- as.data.frame(t(x[[field]]))
        out <- tibble::rownames_to_column(out, "pc")
        out$pc <- paste(x$region, out$pc, sep = "-")
        tibble::column_to_rownames(out, "pc")
    }
    list(
        regional_pcs = do.call(rbind, lapply(res, combine, "sig_pcs")),
        percent_variance = do.call(rbind, lapply(res, combine, "percent_var")),
        loadings = stats::setNames(lapply(res, `[[`, "loadings"), regions),
        info = do.call(rbind, lapply(res,
            function(x) x[c("region", "est_dim", "num_cpgs")]))
    )
}

regional_fixture <- function() {
    set.seed(904)
    x <- matrix(rnorm(12 * 30, sd = 0.1), 12, 30)
    x <- x + outer(seq_len(12), rnorm(30))
    dimnames(x) <- list(paste0("cg", seq_len(12)),
        paste0("sample-", seq_len(30)))
    map <- data.frame(region = c(rep("z-region", 5), rep("a", 5)),
        cpg = c("cg4", "cg2", "cg1", "cg4", "cg8",
            "cg3", "cg7", "cg2", "cg9", "cg11"))
    list(meth = x, map = map)
}

test_that("indexed computation preserves all results and ordering", {
    f <- regional_fixture()
    for (method in c("gd", "mp")) {
        for (x in list(f$meth, as.data.frame(f$meth))) {
            expected <- legacy_regional_pcs(x, f$map, method)
            actual <- compute_regional_pcs(x, f$map, method)
            expect_equal(actual, expected, tolerance = 1e-12)
            expect_identical(names(actual$loadings), c("z-region", "a"))
            expect_identical(colnames(actual$regional_pcs), colnames(x))
        }
    }
})

test_that("factor IDs and extra annotation columns work", {
    f <- regional_fixture()
    expected <- compute_regional_pcs(f$meth, f$map)
    f$map[] <- lapply(f$map, factor)
    f$map$annotation <- "promoter"
    expect_equal(compute_regional_pcs(f$meth, f$map), expected)
})

test_that("standalone region summaries match batched output", {
    f <- regional_fixture()
    single <- summarize_region("z-region", f$map, f$meth, "gd")
    all <- compute_regional_pcs(f$meth, f$map)
    expect_equal(single$loadings, all$loadings[["z-region"]])
    expect_equal(single$num_cpgs, 5L)
    expect_error(summarize_region("missing", f$map, f$meth, "gd"),
        "No CpGs")
})

test_that("invalid mapping and methylation inputs fail clearly", {
    f <- regional_fixture()
    bad <- f$map
    bad$cpg[1] <- "unmeasured"
    expect_error(compute_regional_pcs(f$meth, bad), "All CpGs")
    bad$cpg[1] <- NA_character_
    expect_error(compute_regional_pcs(f$meth, bad), "non-missing")
    expect_error(compute_regional_pcs(f$meth, f$map[FALSE, ]), "non-empty")
    expect_error(compute_regional_pcs(unname(f$meth), f$map), "row names")
    bad <- f$meth
    rownames(bad)[1] <- rownames(bad)[2]
    expect_error(compute_regional_pcs(bad, f$map), "unique")
    expect_error(compute_regional_pcs(f$meth[, 1, drop = FALSE], f$map),
        "two samples")
    expect_error(compute_regional_pcs(matrix("x", 12, 30), f$map), "numeric")
})

test_that("combining no retained PCs preserves the empty output shape", {
    empty <- data.frame(row.names = c("sample1", "sample2"))
    out <- combine_results(list(region = "region", sig_pcs = empty), "sig_pcs")
    expect_equal(dim(out), c(0L, 2L))
    expect_identical(names(out), c("sample1", "sample2"))
})

test_that("a single-CpG region keeps sample and loading dimensions", {
    f <- regional_fixture()
    map <- data.frame(region = "only-region", cpg = "cg1")
    for (method in c("gd", "mp")) {
        actual <- compute_regional_pcs(f$meth, map, method)
        expect_equal(actual, legacy_regional_pcs(f$meth, map, method))
        expect_equal(dim(actual$regional_pcs), c(1L, 30L))
        expect_equal(dim(actual$loadings[[1]]), c(1L, 1L))
    }
})
