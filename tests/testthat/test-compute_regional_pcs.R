# Create synthetic methylation data
meth_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(meth_data) <- paste0("CpG", 1:100)
colnames(meth_data) <- paste0("Sample", 1:10)

# Create a synthetic region map
region_map_data <- data.frame(
  region_id = rep(c("Gene1", "Gene2"), each = 50),
  cpg_id = rownames(meth_data)
)

# Now add some tests
test_that("compute_regional_pcs works as expected", {

  # Test that the function returns a list
  res <- compute_regional_pcs(meth_data, region_map_data, pc_method = 'gd')
  expect_type(res, "list")

  # Test that the list has the expected named elements
  expect_named(res, c("regional_pcs", "percent_variance", "loadings", "info"))

  # Test that the function produces an error for unsupported pc_method
  expect_error(compute_regional_pcs(meth_data,
                                    region_map_data,
                                    pc_method = 'xyz'))

  # Test that the function works for Marchenko-Pastur method as well
  res_mp <- compute_regional_pcs(meth_data,
                                 region_map_data,
                                 pc_method = 'mp')
  expect_type(res_mp, "list")
})
