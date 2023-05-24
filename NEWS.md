# SPIAT 1.3.2

NOTES

* Fixed typo in citation and SPIAT overview diagram.
* Moved the following packages from Imports to Suggests: alphahull, plotly.

# SPIAT 1.3.1

* Added citation to the newly published paper.

# SPIAT 1.3.0

Development version on Bioconductor 3.18.

# SPIAT 1.1.6

BUG FIXES

* Fixed the mixing score and normalised mixing score calculation. Each reference-reference interaction is now counted once (was treated directional and counted twice) and the fraction of normalised mixing score is fixed.

# SPIAT 1.1.5

SIGNIFICANT USER-VISIBLE CHANGES

* Return message instead of error when there are no cells of interest present in the image (identify_neighborhoods()).
* Removed the option of manually defining tumour regions in `identify_bordering_cells()`. Removed parameters `n_of_polygons` and `draw`.

NOTES

* Moved the following packages from Imports to Suggestions: graphics, umap, Rtsne, rlang, ComplexHeatmap and elsa. SpatialExperiment requires version >= 1.8.0. Removed xROI.

# SPIAT 1.1.4

BUG FIXES
* Fixed error when there are only one cell in the clusters. (`identify_neighborhoods()`).

# SPIAT 1.1.3

BUG FIXES
* The calculation of cell types of interest to `All_cells_in_the_structure` in `calculate_proportions_of_cells_in_structure()` was incorrect. Now fixed.

# SPIAT 1.1.2

SIGNIFICANT USER-VISIBLE CHANGES

* Re-organised the vignettes.  

# SPIAT 1.1.1

BUG FIXES
* Fix bug when Cell.ID column is missing from the spe_object in `identify_neighborhood()`.

# SPIAT 1.1.0

Development version on Bioconductor 3.17.

# SPIAT 0.99.14

BUG FIXES

* Fixed bugs in `identify_neighborhoods()`:
1) Assign "Free_cell" to the cells of interest when the number of clustered cells are smaller than `min_cluster_size` in each cluster; 
2) Fixed spe_object output (Adding "Neighborhood" column had a bug previously).

# SPIAT 0.99.13

BUG FIXES

* Minor bug in `average_nearest_neighbor_index()` - the variable `output` was 
not defined when there is no reference cell found in the image.

# SPIAT 0.99.12

SIGNIFICANT USER-VISIBLE CHANGES

* `average_nearest_neighbor_index()` now returns the index along with the 
pattern type and the p value.

# SPIAT 0.99.11

BUG FIXES

* Renamed the file of `R_BC()`.
* Fixed a minor bug in `identify_bordering_cells()` that causes error.

# SPIAT 0.99.10

SIGNIFICANT USER-VISIBLE CHANGES

* `mixing_score_summary()` ensures returning data in any situation. There will 
be difference between returning `NA` and `0`. See updated documentation.

# SPIAT 0.99.9

BUG FIXES

* Fixed `format_image_to_spe()` "general" format. Assigned rownames (markers) 
and colnames (Cell IDs) to the assay of the formatted spe object.

# SPIAT 0.99.8

BUG FIXES

* Fixed `format_image_to_spe()` "general" format. The NA intensity value removal 
was for markers. Instead, it should remove cells that have NA marker intensities.
* Fixed `predict_phenotypes()` plot error.

# SPIAT 0.99.7

BUG FIXES

* Fixed formatting image error when `intensity_matrix` is NULL under "general" 
option.

# SPIAT 0.99.6

SIGNIFICANT USER-VISIBLE CHANGES

* `identify_bordering_cells()` emit Warning when no bordering cells are detected.

BUG FIXES

* Fixed the message displayed for NA removal when formatting image using 
"general" format.
* Fixed plot functions to display only one plot. `plot_marker_level_heatmap()`, 
and `plot_distance_heatmap()`.


# SPIAT 0.99.5

SIGNIFICANT USER-VISIBLE CHANGES

* Generalised functions in tumour structure families to suit other tissue and cell
types, not restricted to tumour and immune cells. The affected functions include:
    1. `identify_bordering_cells()` - deleted "tumour" from the plot title;
    2. `plot_cell_categories()` - cell categories changed to not contain "immune" when `feature_colname == "Structure"`;
    3. `calculate_distance_to_tumour_margin()` renamed to `calculate_distance_to_margin()`;
    4. `R_BT()` renamed to `R_BC()`;
    5. `calculate_summary_distances_of_cells_to_borders()` - one column in the returned data frame has name change;
    6. `defined_strcture()` - parameter `name_of_immune_cells` renamed to `cell_types_of_interest`;
* `identify_bordering_cells()` can return the number of clusters of the 
specified cell type.

BUG FIXES

* Fixed the calculation of the normalised cross-K AUC in `AUC_of_cross_funcion()`.
* `identify_neighborhood()` returns an SPE object instead of sending `ERROR` 
when the cells of interest do no form any clusters.

# SPIAT 0.99.4

BUG FIXES

* Fixed internal function `bind_info()` to avoid duplicated columns.

# SPIAT 0.99.3

BUG FIXES

* Added checks for `dimensionality_reduction_plot()` to return error when the 
cell count in an image is too low.

# SPIAT 0.99.2

SIGNIFICANT USER-VISIBLE CHANGES

* Added `perplexity` parameter to `dimensionality_reduction_plot()`.
* Added `plot_final_border` parameter to `identify_bordering_cells()`.

BUG FIXES

* `image_splitter()` returns `NULL` for the sub-images that do not contain any 
cells. 
* `calculate_spatial_autocorrelation()`
* `identify_bordering_cells()` only plots the bordering cells when bordering 
cells are detected and the user chooses to plot them.
* The calculation of one cell type to `All_cells_in_the_structure` in `calculate_proportions_of_cells_in_structure()` was incorrect. Now fixed.

# SPIAT 0.99.1

SIGNIFICANT USER-VISIBLE CHANGES

* Updated the main object class from `SingleCellExperiment` to `SpatialExperiment`.
* Replace "sce" with "spe" in function names.
* Reading data by `format_image_to_spe()` using `general` format is recommended.
For other data platforms, `format_inform_to_spe()`, `format_halo_to_spe()`, `format_codex_to_spe()`, and `format_cellprofiler_to_spe()` are also available.
* Deleted "Visium" format option in `format_image_to_spe()`.
* `image_splitter()` now returns a list of spe objects.

# SPIAT 0.99.0

Version submitted to Nature Communications. Access the release v0.99.0 [here](https://github.com/TrigosTeam/SPIAT/releases/tag/v0.99.0).

SIGNIFICANT USER-VISIBLE CHANGES

* Systematically renamed parameter names. For example, the parameter for selecting a column of interest is changed from "column" to "feature_colname".
* Added a function `entropy_gradient_aggregated()` for calculating cell colocalisation. 
* Added a function `dimensionality_reduction_plot()` for visualisation and quality control.
* Added examples to documentation. The example datasets were simulated by [**spaSim**](https://github.com/TrigosTeam/spaSim/releases/tag/v0.99.1).

# SPIAT 0.4

SIGNIFICANT USER-VISIBLE CHANGES

* Added a new parameter "column" to most functions to select a column of interest.
* Users can define cell types based on certain marker combinations or cell properties.
* Improved tumour border detection (`identify_bordering_cells()`).
* Added functions to define tumour structure and characterise immune populations in different tissue regions.
* Added functions to characterise spatial heterogeneity.
* Added functions for cell colocalisation, including mixing score, normalised mixing score, cross K function.

# SPIAT 0.2

* Version for the [paper](https://www.biorxiv.org/content/10.1101/2020.05.28.122614v1) submitted to BioRxiv.
