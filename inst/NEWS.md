# Changes in version 1.11.1

* Reduced repeated mapping scans and matrix conversions in
  `compute_regional_pcs()`.
* Added input checks and consistent handling of region-map columns in
  `summarize_region()` and `compute_regional_pcs()`.
* Moved RMTstat to Imports and dplyr and tibble to Suggests.
* Updated the vignette and function documentation, including input scales,
  variance fractions, and component-selection assumptions.
* Added regression tests for both selection methods and regional output order.

# Changes in version 0.99.0 (2023-08-25)

## NEW FEATURES

* Added `compute_regional_pcs()` function to compute principal components for given gene regions.

## IMPROVEMENTS

* Improved performance of `compute_regional_pcs()` by using a more efficient algorithm.

* Updated vignette with more in-depth explanations and examples.

## BUG FIXES

* Fixed bug in `compute_regional_pcs()` that produced incorrect results for certain edge cases.

* Corrected typos in man pages.

## DOCUMENTATION

* Added a detailed vignette for working with `compute_regional_pcs()`.

* Improved man pages with runnable examples for all exported functions.

# Changes in version 0.0.1

* Initial release.
