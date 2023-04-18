# test_that("function runs", {
#   car_map <- data.frame(gear=paste0('gear', mtcars$gear),
#                         car=rownames(mtcars))
#   expect_that(compute_regional_pcs(mtcars, car_map, 'gd', verbose=TRUE)$regional_pcs,
#               is.data.frame()
#               )
# })

# test_that("count_regions", {
#   car_map <- data.frame(gear=mtcars$gear,
#                         car=rownames(mtcars))
#   expect_equal(compute_regional_pcs(mtcars, car_map), 1)
# })
